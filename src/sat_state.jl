struct SatState{C<:AbstractCorrelator}
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

function SatState(
    system::AbstractGNSS,
    prn::Int,
    sampling_frequency,
    code_phase,
    carrier_doppler;
    num_ants::NumAnts = NumAnts(1),
    carrier_phase = 0.0,
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(system),
    num_prompts_for_cn0_estimation::Int = 100,
)
    SatState(
        prn,
        float(code_phase),
        float(code_doppler),
        float(carrier_phase) / 2π,
        float(carrier_doppler),
        0,
        get_default_correlator(system, sampling_frequency, num_ants),
        get_default_correlator(system, sampling_frequency, num_ants),
        complex(0.0, 0.0),
        -1,
        SecondaryCodeOrBitDetector(),
        MomentsCN0Estimator(num_prompts_for_cn0_estimation),
        BitBuffer(),
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

struct SystemSatsState{S<:AbstractGNSS,SS<:SatState}
    system::S
    states::Vector{SS}
end

get_system(sss::SystemSatsState) = sss.system
get_states(sss::SystemSatsState) = sss.states
get_sat_state(sss::SystemSatsState, sat_idx::Integer) = sss.states[sat_idx]

function estimate_cn0(sss::SystemSatsState, sat_idx::Integer)
    system = sss.system
    estimate_cn0(
        get_cn0_estimator(get_sat_state(sss, sat_idx)),
        get_code_length(system) / get_code_frequency(system),
    )
end
