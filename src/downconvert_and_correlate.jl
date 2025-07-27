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
        found(sat_state.sc_bit_detector),
    )
    total_integrated_samples = sat_state.integrated_samples + integrated_samples

    SatState(
        sat_state;
        code_phase,
        carrier_phase,
        integrated_samples = total_integrated_samples,
        signal_start_sample = sat_state.signal_start_sample + integrated_samples,
        correlator = is_integration_completed ? zero(correlator) : correlator,
        last_fully_integrated_correlator = is_integration_completed ? correlator :
                                           sat_state.last_fully_integrated_correlator,
    )
end

"""
$(SIGNATURES)

Downconvert and correlate a single satellite on the CPU.
"""
function downconvert_and_correlate!(
    system,
    signal,
    correlator,
    code_replica,
    code_phase,
    carrier_replica,
    carrier_phase,
    downconverted_signal,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn,
    maximum_expected_sampling_frequency,
)
    gen_code_replica!(
        code_replica,
        system,
        code_frequency,
        sampling_frequency,
        code_phase,
        signal_start_sample,
        num_samples_left,
        correlator.shifts,
        prn,
        maximum_expected_sampling_frequency,
    )
    gen_carrier_replica!(
        carrier_replica,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        signal_start_sample,
        num_samples_left,
    )
    downconvert!(
        downconverted_signal,
        signal,
        carrier_replica,
        signal_start_sample,
        num_samples_left,
    )
    correlate(
        correlator,
        downconverted_signal,
        code_replica,
        signal_start_sample,
        num_samples_left,
    )
end
