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
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(
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
