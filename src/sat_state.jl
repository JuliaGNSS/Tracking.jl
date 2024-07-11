struct SatState{
    C<:AbstractCorrelator,
    DE<:Maybe{<:AbstractSatDopplerEstimator},
    DC<:Maybe{<:AbstractSatDownconvertAndCorrelator},
    P<:Maybe{AbstractSatPostProcess},
}
    prn::Int
    code_phase::Float64
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    carrier_doppler::typeof(1.0Hz)
    integrated_samples::Int
    correlator::C
    last_fully_integrated_correlator::C
    last_fully_integrated_filtered_prompt::ComplexF64
    sample_of_last_fully_integrated_correlator::Int
    sc_bit_detector::SecondaryCodeOrBitDetector
    cn0_estimator::MomentsCN0Estimator
    bit_buffer::BitBuffer
    doppler_estimator::DE
    downconvert_and_correlator::DC
    post_process::P
end

get_prn(s::SatState) = s.prn
get_num_ants(s::SatState{<:AbstractCorrelator{M}}) where {M} = M
get_code_phase(s::SatState) = s.code_phase
get_code_doppler(s::SatState) = s.code_doppler
get_carrier_phase(s::SatState) = s.carrier_phase * 2π
get_carrier_doppler(s::SatState) = s.carrier_doppler
get_integrated_samples(s::SatState) = s.integrated_samples
get_correlator(s::SatState) = s.correlator
get_last_fully_integrated_correlator(s::SatState) = s.last_fully_integrated_correlator
get_last_fully_integrated_filtered_prompt(s::SatState) =
    s.last_fully_integrated_filtered_prompt
get_sample_of_last_fully_integrated_correlator(s::SatState) =
    s.sample_of_last_fully_integrated_correlator
get_secondary_code_or_bit_detector(s::SatState) = s.sc_bit_detector
get_cn0_estimator(s::SatState) = s.cn0_estimator
get_bit_buffer(s::SatState) = s.bit_buffer
get_doppler_estimator(s::SatState) = s.bit_buffer
get_downconvert_and_correlator(s::SatState) =
    s.downconvert_and_correlator, get_post_process(s::SatState) = s.post_process

function SatState(
    system::AbstractGNSS,
    prn::Int,
    sampling_frequency,
    code_phase,
    carrier_doppler;
    num_ants::NumAnts = NumAnts(1),
    correlator = get_default_correlator(system, sampling_frequency, num_ants),
    carrier_phase = 0.0,
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(system),
    num_prompts_for_cn0_estimation::Int = 100,
    doppler_estimator::Maybe{AbstractSatDopplerEstimator} = ConventionalPLLAndDLL(
        carrier_doppler,
        code_doppler,
    ),
    downconvert_and_correlator::Maybe{AbstractSatDownconvertAndCorrelator} = nothing,
    post_process::Maybe{AbstractSatPostProcess} = NoSatPostProcess(),
)
    SatState(
        prn,
        float(code_phase),
        float(code_doppler),
        float(carrier_phase) / 2π,
        float(carrier_doppler),
        0,
        correlator,
        correlator,
        complex(0.0, 0.0),
        -1,
        SecondaryCodeOrBitDetector(),
        MomentsCN0Estimator(num_prompts_for_cn0_estimation),
        BitBuffer(),
        doppler_estimator,
        downconvert_and_correlator,
        post_process,
    )
end

function SatState(acq::AcquisitionResults; args...)
    SatState(
        acq.system,
        acq.prn,
        acq.sampling_frequency,
        acq.code_phase,
        acq.carrier_doppler;
        args...,
    )
end

function SatState(
    sat_state::SS,
    code_phase,
    carrier_phase,
    integrated_samples::Int,
    correlator::AbstractCorrelator,
    sample_of_last_fully_integrated_correlator::Int,
) where {SS<:SatState}
    SS(
        sat_state.prn,
        code_phase,
        sat_state.code_doppler,
        carrier_phase,
        sat_state.carrier_doppler,
        integrated_samples,
        correlator,
        sat_state.last_fully_integrated_correlator,
        sat_state.last_fully_integrated_filtered_prompt,
        sample_of_last_fully_integrated_correlator,
        sat_state.sc_bit_detector,
        sat_state.cn0_estimator,
        sat_state.bit_buffer,
        sat_state.doppler_estimator,
        sat_state.downconvert_and_correlator,
        sat_state.post_process,
    )
end

function SatState(
    sat_state::SS,
    code_doppler,
    carrier_doppler,
    integrated_samples::Int,
    correlator::AbstractCorrelator,
    last_fully_integrated_correlator,
    last_fully_integrated_filtered_prompt,
    sample_of_last_fully_integrated_correlator,
    sc_bit_detector,
    cn0_estimator,
    bit_buffer,
    doppler_estimator::AbstractSatDopplerEstimator,
) where {SS<:SatState}
    SS(
        sat_state.prn,
        sat_state.code_phase,
        code_doppler,
        sat_state.carrier_phase,
        carrier_doppler,
        integrated_samples,
        correlator,
        last_fully_integrated_correlator,
        last_fully_integrated_filtered_prompt,
        sample_of_last_fully_integrated_correlator,
        sc_bit_detector,
        cn0_estimator,
        bit_buffer,
        doppler_estimator,
        sat_state.downconvert_and_correlator,
        sat_state.post_process,
    )
end

struct SystemSatsState{
    S<:AbstractGNSS,
    SS<:SatState,
    DE<:Maybe{AbstractSystemDopplerEstimator},
    DC<:Maybe{AbstractSystemDownconvertAndCorrelator},
    P<:Maybe{AbstractSystemPostProcess},
    I,
}
    system::S
    states::Dictionary{I,SS}
    doppler_estimator::DE
    downconvert_and_correlator::DC
    post_process::P
end

const MultipleSystemSatsState{N,S,SS,DE,DC,P,I} =
    TupleLike{<:NTuple{N,SystemSatsState{<:S,<:SS,<:DE,<:DC,<:P,<:I}}}

function merge_sats(
    multiple_system_sats_state::MultipleSystemSatsState{N},
    system_idx,
    new_sat_states::Dictionary{I,<:SatState},
    num_samples::Int,
) where {N,I}
    system_sats_state = get_system_sats_state(multiple_system_sats_state, system_idx)
    initiated_new_sat_states = map(
        sat_state -> initiate_downconvert_and_correlator(
            system_sats_state.system,
            sat_state,
            num_samples,
        ),
        new_sat_states,
    )
    @set multiple_system_sats_state[system_idx].states =
        merge(system_sats_state.states, initiated_new_sat_states)
end

function to_dictionary(sat_states::Dictionary{I,<:SatState}) where {I}
    sat_states
end

function to_dictionary(sat_states::Vector{<:SatState})
    Dictionary(map(get_prn, sat_states), sat_states)
end

function to_dictionary(sat_state::SatState)
    dictionary((get_prn(sat_state) => sat_state,))
end

function SystemSatsState(
    system::AbstractGNSS,
    states;
    doppler_estimator::Maybe{AbstractSystemDopplerEstimator} = SystemConventionalPLLAndDLL(),
    downconvert_and_correlator::Maybe{AbstractSystemDownconvertAndCorrelator} = nothing,
    post_process::Maybe{AbstractSystemPostProcess} = NoSystemPostProcess(),
)
    SystemSatsState(
        system,
        to_dictionary(states),
        doppler_estimator,
        downconvert_and_correlator,
        post_process,
    )
end

function SystemSatsState(
    system_sats_state::SystemSatsState,
    states::Dictionary{I,<:SatState};
    doppler_estimator = system_sats_state.doppler_estimator,
    downconvert_and_correlator = system_sats_state.downconvert_and_correlator,
    post_process = system_sats_state.post_process,
) where {I}
    SystemSatsState(
        system_sats_state.system,
        states,
        doppler_estimator,
        downconvert_and_correlator,
        post_process,
    )
end

get_system(sss::SystemSatsState) = sss.system
get_states(sss::SystemSatsState) = sss.states
get_sat_state(sss::SystemSatsState, identifier) = sss.states[identifier]
get_downconvert_and_correlator(sss::SystemSatsState) = sss.downconvert_and_correlator
get_post_process(sss::SystemSatsState) = sss.post_process

function estimate_cn0(sss::SystemSatsState, sat_identifier)
    system = sss.system
    estimate_cn0(
        get_cn0_estimator(get_sat_state(sss, sat_identifier)),
        get_code_length(system) / get_code_frequency(system),
    )
end
