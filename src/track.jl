"""
$(SIGNATURES)

Main tracking function that processes one or more `Measurement`s and updates
the tracking state. Performs downconversion, correlation, and Doppler
estimation for all satellites in the track state. Returns an updated
`TrackState` with new phase/Doppler estimates and decoded bits.

Three input shapes for the first positional argument:

| Argument                              | Meaning                                                            |
|---------------------------------------|--------------------------------------------------------------------|
| `AbstractVecOrMat`                    | Today's bare buffer. Single-band TrackState only.                  |
| `Measurement`                         | One band's bundled buffer + sample rate. Single-band TrackState.   |
| `NamedTuple{...}` of `Measurement`s   | Multi-band: one `Measurement` per band key (see [`band_key`](@ref)). |

The bare-buffer form `track(buf, state, fs; intermediate_frequency = ...)`
is preserved as a thin wrapper that builds a single-entry
`NamedTuple{(band_key(state),)}` internally. The two-phase inner loop
(downconvert+correlate across all groups, then estimate across the whole
TrackState) is the same shape regardless of how many measurements are
passed.

For real-time loops processing many chunks of signal in sequence, **construct
the correlator once outside the loop** and pass it via the
`downconvert_and_correlator` keyword argument:

```julia
dc = CPUThreadedDownconvertAndCorrelator()
while got_chunk(rx)
    chunk = read_chunk!(rx)
    track_state = track(chunk, track_state, sampling_frequency;
                        downconvert_and_correlator = dc)
end
```

The default kwarg value builds a fresh correlator (with fresh per-thread
scratch buffers) on every call, which is fine for one-shot use but
defeats the allocation-free design in tight loops. See also [`track!`](@ref)
for the in-place variant that avoids rebuilding `track_state` per call.
"""
function track(
    measurements::Measurements,
    track_state::TS;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    preferred_num_code_blocks_to_integrate = 1,
) where {TS<:TrackState}
    _validate_measurements(track_state, measurements)
    track_state = reset_start_sample_and_bit_buffer(track_state)::TS
    while true
        _all_groups_reached_end(track_state, measurements) && break

        track_state = downconvert_and_correlate(
            downconvert_and_correlator,
            measurements,
            track_state,
            preferred_num_code_blocks_to_integrate,
        )::TS
        track_state = estimate_dopplers_and_filter_prompt(
            track_state,
            measurements,
            preferred_num_code_blocks_to_integrate,
        )::TS
    end
    return track_state
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    band = _single_band(track_state)
    key = band_key(band)
    m = Measurement(signal, sampling_frequency, intermediate_frequency)
    measurements = NamedTuple{(key,)}((m,))
    track(measurements, track_state; kwargs...)
end

# Single-Measurement convenience: still a single-band path.
function track(
    measurement::Measurement,
    track_state::TrackState;
    kwargs...,
)
    band = _single_band(track_state)
    key = band_key(band)
    measurements = NamedTuple{(key,)}((measurement,))
    track(measurements, track_state; kwargs...)
end

"""
$(SIGNATURES)

In-place version of [`track`](@ref). Mutates `track_state` by overwriting the
`Vector{TrackedSat}` slots inside each group instead of rebuilding new
immutable wrappers. Returns the same `track_state` object.

After one warmup call (which seats each satellite's `filtered_prompts`
buffer capacity), the single-threaded path is fully allocation-free; the
threaded path keeps an irreducible Polyester `@batch` closure-capture
allocation of about 160 B per GNSS system per call.

For real-time loops, **construct the correlator once outside the loop**
and pass it via the `downconvert_and_correlator` keyword argument:

```julia
track_state = TrackState(signal, initial_sats)
dc = CPUThreadedDownconvertAndCorrelator()      # hoist!

while got_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_frequency;
           downconvert_and_correlator = dc)
end
```

The correlator holds long-lived per-thread scratch buffers that grow on
first use; rebuilding it via the default kwarg value would re-grow them
every call.
"""
function track!(
    measurements::Measurements,
    track_state::TrackState;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    preferred_num_code_blocks_to_integrate = 1,
)
    _validate_measurements(track_state, measurements)
    reset_start_sample_and_bit_buffer!(track_state)
    while true
        _all_groups_reached_end(track_state, measurements) && break

        downconvert_and_correlate!(
            downconvert_and_correlator,
            measurements,
            track_state,
            preferred_num_code_blocks_to_integrate,
        )
        estimate_dopplers_and_filter_prompt!(
            track_state,
            measurements,
            preferred_num_code_blocks_to_integrate,
        )
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
    band = _single_band(track_state)
    key = band_key(band)
    m = Measurement(signal, sampling_frequency, intermediate_frequency)
    measurements = NamedTuple{(key,)}((m,))
    track!(measurements, track_state; kwargs...)
end

# Single-Measurement convenience: still a single-band path.
function track!(
    measurement::Measurement,
    track_state::TrackState;
    kwargs...,
)
    band = _single_band(track_state)
    key = band_key(band)
    measurements = NamedTuple{(key,)}((measurement,))
    track!(measurements, track_state; kwargs...)
end

# Loop termination: every group has consumed its band's measurement to the
# end. Each group's `signal_start_sample` advances to `num_samples + 1`
# when that group's measurement is fully integrated; the outer `while`
# exits once every group has reached its own band's measurement end.
@inline function _all_groups_reached_end(
    track_state::TrackState, measurements::Measurements,
)
    _check_all_groups_at_end(Tuple(track_state.groups), measurements)
end

# Recursive tuple-walk: each step has fully concrete types.
@inline _check_all_groups_at_end(::Tuple{}, ::Measurements) = true
@inline function _check_all_groups_at_end(t::Tuple, measurements::Measurements)
    g = first(t)
    m = measurements[band_key(g.band)]
    target = get_num_samples(m) + 1
    @inbounds for sat in g.satellites.values
        sat.signal_start_sample == target || return false
    end
    _check_all_groups_at_end(Base.tail(t), measurements)
end
