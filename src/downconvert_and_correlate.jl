function update(
    system::AbstractGNSS,
    sat_state::SatState,
    sat_params::SampleParams,
    intermediate_frequency,
    sampling_frequency,
    correlator,
    num_samples_signal::Int,
)
    carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
    code_frequency = sat_state.code_doppler + get_code_frequency(system)
    carrier_phase = update_carrier_phase(
        sat_params.signal_samples_to_integrate,
        carrier_frequency,
        sampling_frequency,
        sat_state.carrier_phase,
    )
    code_phase = update_code_phase(
        system,
        sat_params.signal_samples_to_integrate,
        code_frequency,
        sampling_frequency,
        sat_state.code_phase,
        found(sat_state.sc_bit_detector),
    )
    sample_of_last_fully_integrated_correlator =
        sat_state.sample_of_last_fully_integrated_correlator -
        (sat_params.signal_start_sample == 1 ? num_samples_signal : 0)
    SatState(
        sat_state,
        code_phase,
        carrier_phase,
        sat_state.integrated_samples + sat_params.signal_samples_to_integrate,
        correlator,
        sample_of_last_fully_integrated_correlator,
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
