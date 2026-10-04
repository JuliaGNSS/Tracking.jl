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

The bare-buffer form `track(buf, state, fs; intermediate_frequency = ...)` is a
thin wrapper that builds the single-entry `NamedTuple{(get_band_id(band),)}`.

The returned `TrackState` is *structurally* detached from the input: each
group's key set and slot vector are copied, so [`add_satellite!`](@ref) /
[`remove_satellite!`](@ref) and tracking on either state never affect the
other's satellites. The copy is shallow, however: per-satellite scratch vectors
(each signal's `filtered_prompts` and `correlator_outputs`, the soft-bit buffer,
the CN0 estimator's prompt buffer) are shared and overwritten by the next `track`
call on either state. Treat the input as a stale handle; `deepcopy` it first to
snapshot those buffers. The same holds for a bare `downconvert_and_correlate`.

Each signal's noise estimator is shared the same way and advanced in place (a
[`CorrelatorNoiseEstimator`](@ref)'s window and RNG stream), so branched states
share one noise reference and race on it when advanced concurrently. Build a
separate `TrackState` per thread — see `downconvert_and_correlate`.

For real-time loops, **construct the correlator once outside the loop** and pass
it via the `downconvert_and_correlator` keyword argument:

```julia
dc = CPUThreadedDownconvertAndCorrelator()
while got_chunk(rx)
    chunk = read_chunk!(rx)
    track_state =
        track(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
end
```

The default kwarg value builds a fresh correlator (and scratch buffers) on every
call, which defeats the allocation-free design in tight loops. See also
[`track!`](@ref), the in-place variant.

The coherent-integration length is a **per-signal** setting on each
[`TrackedSignal`](@ref) (`preferred_num_code_blocks_to_integrate`), set with
[`set_preferred_num_code_blocks_to_integrate!`](@ref). It is capped by the
signal's bit/secondary-code period, held at 1 until bit/secondary sync, and
defaults to [`default_num_code_blocks_to_integrate`](@ref):

```julia
set_preferred_num_code_blocks_to_integrate!(track_state, :gps_l5, 1, GPSL5I, 10)  # PRN 1 L5I: 10 ms
```

Longer integration needs no loop re-tuning: the conventional estimator caps
the bandwidths for stability (see [`ConventionalPLLAndDLL`](@ref)).
"""
function track(
    measurements::BandMeasurements,
    track_state::TS;
    kwargs...,
) where {TS<:TrackState}
    # Detach keys and slot values once (#123), then run the in-place pipeline
    # on the copy rather than re-copying per chunk (#133). Shallow otherwise —
    # see the docstring.
    detached =
        TrackState(track_state; groups = _detach_groups_slot_vectors(track_state.groups))
    track!(measurements, detached; kwargs...)::TS
end

# Wrap a single `BandMeasurement` into the one-entry `BandMeasurements` keyed by
# the TrackState's only band; errors (via `_single_band`) on multi-band states.
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
residual of about 96 B per **completed code-block integration** per group, only
with more than one thread *and* more than one work item in the group's loop (a
measured noise density adds a work item, so a single-satellite group pays it
too). It scales with the chunk length but is dwarfed by the caller's sample
buffer. The single-threaded backend has none. The cause is noted at the
definition of `CPUThreadedDownconvertAndCorrelator`.

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
    # Walk the measurement chunk by chunk: one correlate pass up to each sat's
    # last completed boundary (`stop_before_partial`), then one estimate at a
    # common epoch. See "Chunked Doppler updates" in docs/src/track.md.
    chunk_duration = _resolve_doppler_update_interval(doppler_update_interval, track_state)
    _validate_doppler_update_interval(chunk_duration, measurements)
    # `samples_unchanged`: the buffers are fixed for the whole call, so
    # sample-derived backend caches are built on the first pass only.
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
    # Drain each satellite's trailing partial into its live accumulator so it
    # carries into the next `track!` call; a boundary exactly at the buffer end
    # completes here, hence the final fold. Noise is measured here only when no
    # chunk ran, otherwise every sample would enter the noise window twice.
    downconvert_and_correlate!(
        downconvert_and_correlator,
        measurements,
        track_state;
        samples_unchanged = chunk_index > 0,
        measure_noise = chunk_index == 0,
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

# Loop termination: walk the chunk grid until the previous chunk's end reached
# the buffer end on every band, independent of per-sat progress (lagging sats
# catch up in later passes and the final drain). For `chunk_index == 0`,
# `_chunk_last_sample(…, -1, …) == 0` makes this "is any band non-empty".
@inline function _chunks_left(chunk_duration, chunk_index::Int, measurements)
    for m in measurements
        n = get_num_samples(m)
        _chunk_last_sample(chunk_duration, chunk_index - 1, m.sampling_frequency, n) < n &&
            return true
    end
    false
end

# Resolve the per-chunk update interval. `nothing` => auto: the smallest
# primary-code period across all signals (1 ms for GPS L1 C/A + Galileo E1B).
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

# A chunk must cover at least one sample on every band, or the chunk grid could
# fail to advance and `track!` would not terminate. The dimension check turns a
# unitless interval into a clear ArgumentError instead of a Unitful error.
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
