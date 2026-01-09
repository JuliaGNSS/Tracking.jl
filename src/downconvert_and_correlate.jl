function update(
    system::AbstractGNSS,
    sat_state::SatState,
    integrated_samples::Int,
    intermediate_frequency,
    sampling_frequency,
    correlator,
    is_integration_completed::Bool,
)
    carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
    code_frequency = sat_state.code_doppler + get_code_frequency(system)
    carrier_phase = update_carrier_phase(
        integrated_samples,
        carrier_frequency,
        sampling_frequency,
        sat_state.carrier_phase,
    )
    code_phase = update_code_phase(
        system,
        integrated_samples,
        code_frequency,
        sampling_frequency,
        sat_state.code_phase,
        has_bit_or_secondary_code_been_found(sat_state),
    )
    total_integrated_samples = sat_state.integrated_samples + integrated_samples

    SatState(
        sat_state;
        code_phase,
        carrier_phase,
        integrated_samples = total_integrated_samples,
        signal_start_sample = sat_state.signal_start_sample + integrated_samples,
        is_integration_completed,
        correlator,
    )
end
