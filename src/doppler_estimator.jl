function update(
    system::AbstractGNSS,
    sat_state::SatState,
    sample_params::SampleParams,
    carrier_doppler,
    code_doppler,
    filtered_prompt,
    doppler_estimator,
)
    cn0_estimator = update(get_cn0_estimator(sat_state), filtered_prompt)
    bit_buffer = buffer(
        system,
        sat_state.bit_buffer,
        sample_params.num_code_blocks_to_integrate,
        found(sat_state.sc_bit_detector),
        filtered_prompt,
    )
    sc_bit_detector = find(system, sat_state.sc_bit_detector, filtered_prompt)
    correlator = zero(sat_state.correlator)
    integrated_samples = 0

    return SatState(
        sat_state,
        code_doppler,
        carrier_doppler,
        integrated_samples,
        correlator,
        sat_state.correlator,
        filtered_prompt,
        sample_params.signal_start_sample + sample_params.signal_samples_to_integrate - 1,
        sc_bit_detector,
        cn0_estimator,
        bit_buffer,
        doppler_estimator,
    )
end