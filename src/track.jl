"""
$(SIGNATURES)

Main tracking function that processes a signal and updates the tracking state.
Performs downconversion, correlation, and Doppler estimation for all satellites
in the track state. Returns an updated TrackState with new phase/Doppler estimates
and decoded bits.

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
    signal::AbstractVecOrMat,
    track_state::TS,
    sampling_frequency;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    intermediate_frequency = 0.0Hz,
    preferred_num_code_blocks_to_integrate = 1,
) where {TS<:TrackState}
    track_state = reset_start_sample_and_bit_buffer(track_state)::TS
    num_samples_signal = get_num_samples(signal)
    while true
        has_integration_reached_signal_end_for_all_satellites(
            track_state,
            num_samples_signal,
        ) && break

        track_state = downconvert_and_correlate(
            downconvert_and_correlator,
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )::TS
        track_state = estimate_dopplers_and_filter_prompt(
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )::TS
    end
    return track_state
end

"""
$(SIGNATURES)

In-place version of [`track`](@ref). Mutates `track_state` by overwriting the
`Vector{TrackedSat}` slots inside each system instead of rebuilding new
immutable wrappers. Returns the same `track_state` object.

After one warmup call (which seats each satellite's `filtered_prompts`
buffer capacity), the single-threaded path is fully allocation-free; the
threaded path keeps an irreducible Polyester `@batch` closure-capture
allocation of about 160 B per GNSS system per call.

For real-time loops, **construct the correlator once outside the loop**
and pass it via the `downconvert_and_correlator` keyword argument:

```julia
track_state = TrackState(system, initial_sats)
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
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    intermediate_frequency = 0.0Hz,
    preferred_num_code_blocks_to_integrate = 1,
)
    reset_start_sample_and_bit_buffer!(track_state)
    num_samples_signal = get_num_samples(signal)
    while true
        has_integration_reached_signal_end_for_all_satellites(
            track_state,
            num_samples_signal,
        ) && break

        downconvert_and_correlate!(
            downconvert_and_correlator,
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )
        estimate_dopplers_and_filter_prompt!(
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )
    end
    return track_state
end
