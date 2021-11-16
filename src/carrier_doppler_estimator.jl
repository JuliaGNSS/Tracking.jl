abstract type AbstractCarrierDopplerEstimatorState end
abstract type AbstractCarrierDopplerEstimatorVariable end

struct CostasLoopState{LF <: AbstractLoopFilter} <: AbstractCarrierDopplerEstimatorState
    loop_filter::LF
end

struct CostasLoopBandwidth <: AbstractCarrierDopplerEstimatorVariable
    bandwidth::typeof(1.0Hz)
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(correlator, correlator_sample_shifts)
    p = get_prompt(correlator, correlator_sample_shifts)
    atan(imag(p) / real(p))
end


function est_carrier_doppler(
    costas_loop_state::CostasLoopState,
    costas_loop_bandwidth::CostasLoopBandwidth,
    correlator,
    filtered_correlator,
    carrier_phase,
    carrier_frequency,
    correlator_sample_shifts,
    sampling_frequency,
    integration_time,
    integrated_samples_wrt_signal,
    init_carrier_doppler,
)
    pll_discriminator = pll_disc(
        filtered_correlator,
        correlator_sample_shifts
    )
    carrier_freq_update, loop_filter = filter_loop(
        costas_loop_state.loop_filter,
        pll_discriminator,
        integration_time,
        costas_loop_bandwidth.bandwidth
    )
    carrier_doppler = carrier_freq_update + init_carrier_doppler
    carrier_doppler, carrier_phase, CostasLoopState(loop_filter)
end