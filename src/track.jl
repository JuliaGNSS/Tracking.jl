"""
$(SIGNATURES)

Main tracking function that processes a signal and updates the tracking state.
Performs downconversion, correlation, and Doppler estimation for all satellites
in the track state. Returns an updated TrackState with new phase/Doppler estimates
and decoded bits.
"""
function track(
    signal::AbstractVecOrMat,
    track_state::TS,
    sampling_frequency;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(
        Val(sampling_frequency),
    ),
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

In steady state (no satellites added or removed, `filtered_prompts` capacity
already grown — see [`prewarm!`](@ref)), each individual stage
(`reset_start_sample_and_bit_buffer!`, `downconvert_and_correlate!`,
`estimate_dopplers_and_filter_prompt!`) is allocation-free. The full call
currently shows a small per-call residual (≤ a few hundred bytes) due to
internal closures used in `bit_buffer`'s find-bit logic; those allocations
are short-lived and stay on the young-generation pool, so they do not
trigger GC pauses in the hot SDR loop.
"""
function track!(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(
        Val(sampling_frequency),
    ),
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

"""
$(SIGNATURES)

Pre-grow each tracked satellite's `filtered_prompts` buffer so that the first
calls to [`track!`](@ref) do not pay the geometric-doubling reallocations
that `push!` would otherwise incur. After `prewarm!`, steady-state `track!`
calls are allocation-free.

`max_prompts_per_track_call` is the upper bound on how many integrations can
complete inside a single `track!` call. For one-millisecond code periods and
a few-millisecond chunk size this is small (≤ a handful); for longer chunks
size accordingly.
"""
function prewarm!(track_state::TrackState, max_prompts_per_track_call::Integer)
    for system_sats_state in track_state.multiple_system_sats_state
        for tracked_sat in system_sats_state.states.values
            sizehint!(tracked_sat.sat_state.filtered_prompts, max_prompts_per_track_call)
        end
    end
    return track_state
end
