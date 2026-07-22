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
`filtered_prompts`, the soft-bit buffer, and the CN0 estimator's
prompt buffer) are shared with the input and are overwritten by the
next `track` call on either state. Treat the input as a stale handle
after the call; `deepcopy` it first if you need to snapshot those
buffers.

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
)
    _validate_measurements(track_state, measurements)
    reset_start_sample_and_bit_buffer!(track_state)
    # The measurement buffers are fixed for the whole call, so sample-derived
    # backend caches (the bit backends' shared band pack) are built on the
    # first iteration and reused ever after (`samples_unchanged`).
    first_iteration = true
    while true
        _all_groups_reached_end(track_state, measurements) && break

        downconvert_and_correlate!(
            downconvert_and_correlator,
            measurements,
            track_state;
            samples_unchanged = !first_iteration,
        )
        first_iteration = false
        estimate_dopplers_and_filter_prompt!(track_state, measurements)
    end
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

# Loop termination: every group has consumed its band's measurement to the
# end. Each group's `signal_start_sample` advances to `num_samples + 1`
# when that group's measurement is fully integrated; the outer `while`
# exits once every group has reached its own band's measurement end.
@inline function _all_groups_reached_end(
    track_state::TrackState,
    measurements::BandMeasurements,
)
    _check_all_groups_at_end(Tuple(track_state.groups), measurements)
end

# Recursive tuple-walk: each step has fully concrete types.
@inline _check_all_groups_at_end(::Tuple{}, ::BandMeasurements) = true
@inline function _check_all_groups_at_end(t::Tuple, measurements::BandMeasurements)
    g = first(t)
    m = measurements[get_band_id(g.band)]
    target = get_num_samples(m) + 1
    @inbounds for sat in g.satellites.values
        sat.signal_start_sample == target || return false
    end
    _check_all_groups_at_end(Base.tail(t), measurements)
end
