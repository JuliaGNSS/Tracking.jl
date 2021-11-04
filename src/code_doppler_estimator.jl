abstract type AbstractCodeDopplerEstimatorState end
abstract type AbstractCodeDopplerEstimatorVariable end

struct EarlyPromptLateLoopState{LF <: AbstractLoopFilter} <: AbstractCodeDopplerEstimatorState
    loop_filter::LF
end

struct EarlyPromptLateLoopBandwidth <: AbstractCodeDopplerEstimatorVariable
    bandwidth::typeof(1.0Hz)
end

"""
$(SIGNATURES)

Calculates the code phase error in chips.
"""
function dll_disc(
    correlator,
    correlator_sample_shifts,
    early_late_index_shift,
    code_phase_delta
)
    E = abs(get_early(correlator, correlator_sample_shifts, early_late_index_shift))
    L = abs(get_late(correlator, correlator_sample_shifts, early_late_index_shift))
    distance_between_early_and_late =
        get_early_late_sample_spacing(correlator_sample_shifts, early_late_index_shift) *
        code_phase_delta
    (E - L) / (E + L) / (2 * (2 - distance_between_early_and_late))
end


function est_code_doppler(
    epl_loop_state::EarlyPromptLateLoopState,
    epl_loop_bandwidth::EarlyPromptLateLoopBandwidth,
    correlator,
    filtered_correlator,
    correlator_sample_shifts,
    early_late_index_shift,
    code_frequency,
    sampling_frequency,
    integration_time,
    init_code_doppler,
)
    dll_discriminator = dll_disc(
        filtered_correlator,
        correlator_sample_shifts,
        early_late_index_shift,
        code_frequency / sampling_frequency
    )
    code_freq_update, loop_filter = filter_loop(
        epl_loop_state.loop_filter,
        dll_discriminator,
        integration_time,
        epl_loop_bandwidth.bandwidth
    )
    code_doppler = code_freq_update + init_code_doppler
    code_doppler, EarlyPromptLateLoopState(loop_filter)
end