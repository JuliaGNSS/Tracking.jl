function track(
    signal::AbstractVecOrMat,
    track_state::TS,
    sampling_frequency;
    intermediate_frequency = 0.0Hz,
    preferred_num_code_blocks_to_integrate = 1,
    maximum_expected_sampling_frequency = Val(sampling_frequency),
) where {TS<:TrackState}
    track_state = reset_start_sample(track_state)::TS
    num_samples_signal = get_num_samples(signal)
    while true
        has_integration_reached_signal_end_for_all_satellites(
            track_state,
            num_samples_signal,
        ) && break

        track_state = downconvert_and_correlate(
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
            num_samples_signal,
            maximum_expected_sampling_frequency,
        )::TS
        track_state = estimate_dopplers_and_filter_prompt(
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )::TS
        track_state = post_process(track_state)::TS
    end
    return track_state
end
