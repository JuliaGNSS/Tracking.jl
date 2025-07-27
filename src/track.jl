function track(
    signal::AbstractVecOrMat,
    track_state::TS,
    sampling_frequency;
    intermediate_frequency = 0.0Hz,
    preferred_num_code_blocks_to_integrate = 1,
    min_integration_time = 0.75ms,
    maximum_expected_sampling_frequency = Val(sampling_frequency),
) where {TS<:TrackState}
    multiple_system_sats_state = track_state.multiple_system_sats_state
    sat_sample_params = init_sample_params(
        multiple_system_sats_state,
        preferred_num_code_blocks_to_integrate,
    )
    num_samples_signal = get_num_samples(signal)
    while true
        sat_sample_params = calc_sample_params(
            multiple_system_sats_state,
            sat_sample_params,
            num_samples_signal,
            sampling_frequency,
            preferred_num_code_blocks_to_integrate,
        )
        all(x -> all(y -> y.signal_samples_to_integrate == 0, x), sat_sample_params) &&
            break

        track_state = downconvert_and_correlate(
            signal,
            track_state,
            sat_sample_params,
            sampling_frequency,
            intermediate_frequency,
            num_samples_signal,
            maximum_expected_sampling_frequency,
        )::TS
        track_state = estimate_dopplers_and_filter_prompt(
            track_state,
            sat_sample_params,
            sampling_frequency,
            min_integration_time,
        )::TS
        track_state = post_process(track_state)::TS
    end
    return track_state
end
