function track(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = 0.0Hz,
    preferred_num_code_blocks_to_integrate = 1,
    min_integration_time = 0.75ms,
)
    system_sats_states = track_state.system_sats_states
    doppler_estimator = track_state.doppler_estimator
    sat_sample_params =
        init_sample_params(system_sats_states, preferred_num_code_blocks_to_integrate)
    num_samples_signal = get_num_samples(signal)
    while true
        sat_sample_params = calc_sample_params(
            system_sats_states,
            sat_sample_params,
            num_samples_signal,
            sampling_frequency,
            preferred_num_code_blocks_to_integrate,
        )
        all(x -> all(y -> y.signal_samples_to_integrate == 0, x), sat_sample_params) &&
            break
        correlators = downconvert_and_correlate(
            track_state.downconvert_and_correlator,
            signal,
            sampling_frequency,
            intermediate_frequency,
            system_sats_states,
            sat_sample_params,
        )
        system_sats_states = update(
            system_sats_states,
            sampling_frequency,
            intermediate_frequency,
            correlators,
            sat_sample_params,
            num_samples_signal,
        )
        doppler_estimator, dopplers_and_prompts = estimate_dopplers_and_filter_prompt(
            doppler_estimator,
            system_sats_states,
            sat_sample_params,
            sampling_frequency,
            min_integration_time,
        )
        system_sats_states =
            update(system_sats_states, dopplers_and_prompts, sat_sample_params)
        track_state = TrackState(
            track_state.post_process,
            system_sats_states,
            doppler_estimator,
            track_state.downconvert_and_correlator,
            track_state.num_samples,
        )
        track_state = post_process(track_state)
    end
    return track_state
end
