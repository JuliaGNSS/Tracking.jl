"""
$(SIGNATURES)

Main tracking function that processes one or more `BandMeasurement`s and updates
the tracking state. Performs downconversion, correlation, and Doppler
estimation for all satellites in the track state. Returns an updated
`TrackState` with new phase/Doppler estimates and decoded bits.

Three input shapes for the first positional argument:

| Argument                                | Meaning                                                                        |
|:--------------------------------------- |:------------------------------------------------------------------------------ |
| `AbstractVecOrMat`                      | Bare sample buffer. Single-band TrackState only.                               |
| `BandMeasurement`                       | One band's bundled buffer + sample rate. Single-band TrackState.               |
| `NamedTuple{...}` of `BandMeasurement`s | Multi-band: one `BandMeasurement` per band id (see `GNSSSignals.get_band_id`). |

The bare-buffer form `track(buf, state, fs; intermediate_frequency = ...)`
is preserved as a thin wrapper that builds a single-entry
`NamedTuple{(get_band_id(band),)}` internally. The two-phase inner loop
(downconvert+correlate across all groups, then estimate across the whole
TrackState) is the same shape regardless of how many measurements are
passed.

The returned `TrackState` is *structurally* detached from the input:
each group's key set and slot vector are copied, so
[`add_satellite!`](@ref) / [`remove_satellite!`](@ref) and tracking
itself on either state never affect the other's satellites. The copy
is shallow, however — per-satellite scratch vectors (each signal's
`filtered_prompts` and `correlator_outputs`, the soft-bit buffer, and
the CN0 estimator's prompt buffer) are shared with the input and are
overwritten by the next `track` call on either state. Treat the input
as a stale handle after the call; `deepcopy` it first if you need to
snapshot those buffers. The same applies to a bare
`downconvert_and_correlate`: the returned state's
`correlator_outputs` alias the input's, so reuse of one input state
across several calls appends to the same buffers.

For real-time loops processing many chunks of signal in sequence, **construct
the correlator once outside the loop** and pass it via the
`downconvert_and_correlator` keyword argument:

```julia
dc = CPUThreadedDownconvertAndCorrelator()
while got_chunk(rx)
    chunk = read_chunk!(rx)
    track_state =
        track(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
end
```

The default kwarg value builds a fresh correlator (with fresh per-thread
scratch buffers) on every call, which is fine for one-shot use but
defeats the allocation-free design in tight loops. See also [`track!`](@ref)
for the in-place variant that avoids rebuilding `track_state` per call.

The coherent-integration length is a **per-signal** setting that lives on each
[`TrackedSignal`](@ref) (its `preferred_num_code_blocks_to_integrate` field,
addressed by `(group, prn, signal)`), not a `track!` argument. Set it with
[`set_preferred_num_code_blocks_to_integrate!`](@ref); the actual length is
capped per integration by the signal's bit/secondary-code period and held at 1
until bit/secondary sync. Defaults to 1 (1 ms for GPS L5I / L1 C/A). Different
satellites — and different signals on one satellite — can therefore integrate
for different lengths:

```julia
set_preferred_num_code_blocks_to_integrate!(track_state, :gps_l5, 1, GPSL5I, 10)  # PRN 1 L5I: 10 ms
```

The conventional estimator auto-scales each signal's loop bandwidth by `1/N`
for its integration length `N`, so longer integration stays stable without
re-tuning (see [`ConventionalPLLAndDLL`](@ref)).
"""
function track(
    measurements::BandMeasurements,
    track_state::TS;
    kwargs...,
) where {TS<:TrackState}
    # Detach the slot storage from the input once — keys (`Indices`) *and*
    # values (`_detach_groups_slot_vectors`, #123) — so a later
    # `add_satellite!`/`remove_satellite!` on the returned state cannot
    # corrupt the input's key set. Then run the fully in-place pipeline on
    # the detached copy, which avoids re-copying the slot vectors twice per
    # chunk iteration (issue #133). The copy is otherwise shallow:
    # per-satellite scratch vectors are shared with the input — see the
    # docstring above.
    detached =
        TrackState(track_state; groups = _detach_groups_slot_vectors(track_state.groups))
    track!(measurements, detached; kwargs...)::TS
end

# Wrap a bare buffer / single `BandMeasurement` into the one-entry
# `BandMeasurements` NamedTuple keyed by the TrackState's only band. Shared
# by the `track` and `track!` convenience wrappers; errors (via
# `_single_band`) on multi-band TrackStates.
@inline function _single_band_measurements(
    measurement::BandMeasurement,
    track_state::TrackState,
)
    key = get_band_id(_single_band(track_state))
    NamedTuple{(key,)}((measurement,))
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    m = BandMeasurement(signal, sampling_frequency, intermediate_frequency)
    track(_single_band_measurements(m, track_state), track_state; kwargs...)
end

# Single-BandMeasurement convenience: still a single-band path.
function track(measurement::BandMeasurement, track_state::TrackState; kwargs...)
    track(_single_band_measurements(measurement, track_state), track_state; kwargs...)
end

"""
$(SIGNATURES)

In-place version of [`track`](@ref). Mutates `track_state` by overwriting the
`Vector{TrackedSat}` slots inside each group instead of rebuilding new
immutable wrappers. Returns the same `track_state` object.

After one warmup call (which seats each satellite's `filtered_prompts`
buffer capacity), the single-threaded path is fully allocation-free.

The threaded path (`CPUThreadedDownconvertAndCorrelator`) keeps a small
residual — about 64 B per **completed code-block integration** per group —
but only when the process runs with more than one thread *and* the group
holds more than one satellite (so Polyester's `@batch` actually distributes
work). `track!` launches one `@batch` per code block, so this residual
scales with the chunk length rather than staying flat per call; for a
real-time loop with fixed-size chunks it is bounded per call and, in
practice, dwarfed by the input sample buffer the caller allocates each
chunk. The single-threaded backend has none of it.

The root cause is not the per-satellite state per se — it is that the
GNSS **signal** object is not an `isbits` type: `GPSL1CA`, for example, holds
a `Matrix{Int16}` code table and a (also non-`isbits`) `SignalLUT`. Polyester
roots bare `Array`s and `isbits` values into its worker tasks for free (that
is why a `Vector{Float64}` kernel is allocation-free), but it pays a small
per-launch allocation to root any *other* non-`isbits` object touched inside
the `@batch` region. Each satellite's `downconvert_and_correlate` reaches its
signal (for code-replica generation), so the parallel loop touches that
non-`isbits` object once per launch. It could be removed by generating the code
from the LUT's *bare arrays* (a prototype doing so inside the loop measures 0 B
at full parallel throughput) — which needs a GNSSSignals-side bare-array
`gen_code!` — or, in-tree, by a serial code-gen pre-pass (allocation-free but
slower, since it serializes ~30 % of the work). Neither is currently worth
~64 B/block. See [`CPUThreadedDownconvertAndCorrelator`](@ref) for the full
analysis.

For real-time loops, **construct the correlator once outside the loop**
and pass it via the `downconvert_and_correlator` keyword argument:

```julia
track_state = TrackState(signal, initial_sats)
dc = CPUThreadedDownconvertAndCorrelator()      # hoist!

while got_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
end
```

The correlator holds long-lived per-thread scratch buffers that grow on
first use; rebuilding it via the default kwarg value would re-grow them
every call.
"""
function track!(
    measurements::BandMeasurements,
    track_state::TrackState;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    doppler_update_interval = nothing,
)
    _validate_measurements(track_state, measurements)
    reset_start_sample_and_bit_buffer!(track_state)
    # Resolve the Doppler-update / chunk interval. `nothing` => auto: the
    # smallest primary-code period across all tracked signals, so a default
    # chunk holds one code period of the shortest signal. The measurement is
    # walked chunk by chunk: each chunk correlates (collecting every completed
    # correlator output per signal into its `correlator_outputs` buffer) with
    # the NCO Doppler held fixed, then the estimator folds over those outputs
    # and updates every sat's NCO once — a common epoch across all sats.
    chunk_duration = _resolve_doppler_update_interval(doppler_update_interval, track_state)
    _validate_doppler_update_interval(chunk_duration, measurements)
    # One correlate pass + one estimate per chunk. The pass runs each satellite
    # from wherever it stands to its last completed code-block boundary inside
    # the chunk (`stop_before_partial` — the chunk-clamped trailing partial is
    # NOT integrated); the estimator then folds the collected outputs and
    # writes the new NCO Doppler. The residue is picked up by the NEXT chunk's
    # pass, which therefore covers boundary → boundary in a single kernel
    # window, entirely at the just-updated Doppler. So every completed
    # integration is produced by a single NCO Doppler and each correction takes
    # effect right at its completing boundary (the classic per-completion loop
    # timing), while the estimator still runs once per chunk at a common epoch
    # — without splitting each code period into two kernel invocations.
    #
    # `samples_unchanged`: the measurement buffers are fixed for the whole
    # call, so sample-derived backend caches (the bit backends' shared band
    # pack) are built on the very first pass and reused ever after.
    chunk_index = 0
    while _chunks_left(chunk_duration, chunk_index, measurements)
        downconvert_and_correlate!(
            downconvert_and_correlator,
            measurements,
            track_state;
            chunk_index,
            chunk_duration,
            stop_before_partial = true,
            samples_unchanged = chunk_index > 0,
        )
        estimate_dopplers_and_filter_prompt!(track_state, measurements)
        chunk_index += 1
    end
    # Drain the buffer: consume every satellite's trailing partial — from its
    # last completed boundary to the buffer end — into its live accumulator (at
    # the final chunk's Doppler), so the integration carries into the next
    # `track!` call. A boundary landing exactly on the buffer end completes
    # here, so fold once more; a no-op (per-sat early return) otherwise.
    downconvert_and_correlate!(
        downconvert_and_correlator,
        measurements,
        track_state;
        samples_unchanged = chunk_index > 0,
    )
    estimate_dopplers_and_filter_prompt!(track_state, measurements)
    return track_state
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track!(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    m = BandMeasurement(signal, sampling_frequency, intermediate_frequency)
    track!(_single_band_measurements(m, track_state), track_state; kwargs...)
end

# Single-BandMeasurement convenience: still a single-band path.
function track!(measurement::BandMeasurement, track_state::TrackState; kwargs...)
    track!(_single_band_measurements(measurement, track_state), track_state; kwargs...)
end

# Loop termination: the chunk grid is walked until the previous chunk's end
# already reached the buffer end on every band (`_chunk_last_sample` clamps at
# `num_samples`, so the count is finite and independent of per-sat progress —
# satellites lagging behind a chunk boundary are caught up by later passes and
# by `track!`'s final buffer-draining pass). For `chunk_index == 0` the
# convention `_chunk_last_sample(…, -1, …) == 0` makes this "is the buffer
# non-empty on any band".
@inline function _chunks_left(chunk_duration, chunk_index::Int, measurements)
    for m in measurements
        n = get_num_samples(m)
        _chunk_last_sample(chunk_duration, chunk_index - 1, m.sampling_frequency, n) < n &&
            return true
    end
    false
end

# Resolve the per-chunk update interval to a concrete time. `nothing` => auto:
# the smallest primary-code period across every signal in every group, so the
# default chunk holds exactly one code period of the shortest signal (e.g. 1 ms
# for a GPS L1 C/A + Galileo E1B track state).
@inline _resolve_doppler_update_interval(doppler_update_interval, ::TrackState) =
    doppler_update_interval
@inline function _resolve_doppler_update_interval(::Nothing, track_state::TrackState)
    _smallest_code_period(track_state)
end

# Primary-code period (a time) of one signal.
@inline _code_period(signal::AbstractGNSSSignal) =
    get_code_length(signal) / get_code_frequency(signal)

@inline _min_signal_code_period(::Tuple{}, m) = m
@inline _min_signal_code_period(signals::Tuple, m) =
    _min_signal_code_period(Base.tail(signals), min(m, _code_period(first(signals))))

@inline _min_group_code_period(::Tuple{}, m) = m
@inline _min_group_code_period(groups::Tuple, m) = _min_group_code_period(
    Base.tail(groups),
    _min_signal_code_period(first(groups).signals, m),
)

function _smallest_code_period(track_state::TrackState)
    groups = Tuple(track_state.groups)
    # Every TrackState has at least one group and every group at least one
    # signal, so seeding from the first signal's period is safe and keeps the
    # reduction type-stable across heterogeneous signal tuples.
    init = _code_period(first(first(groups).signals))
    _min_group_code_period(groups, init)
end

# A chunk must cover at least one sample on every band, otherwise the chunk
# grid could fail to advance and `track!` would not terminate. The default
# (smallest code period) is always many samples; this only guards against a
# user-supplied `doppler_update_interval` shorter than a sample period. The dimension
# check turns a plain number (e.g. `doppler_update_interval = 1e-3`) into a clear
# ArgumentError instead of a cryptic Unitful conversion error.
function _validate_doppler_update_interval(chunk_duration, measurements::BandMeasurements)
    dimension(chunk_duration) == dimension(1.0s) || throw(
        ArgumentError(
            "doppler_update_interval must be a time quantity, e.g. `1u\"ms\"` or `1e-3u\"s\"` " *
            "(with `using Unitful`); got $chunk_duration.",
        ),
    )
    for m in measurements
        samples_per_chunk = uconvert(NoUnits, chunk_duration * m.sampling_frequency)
        samples_per_chunk >= 1 || throw(
            ArgumentError(
                "doppler_update_interval $chunk_duration is shorter than one sample period " *
                "at sampling frequency $(m.sampling_frequency); pick a longer interval.",
            ),
        )
    end
    nothing
end
