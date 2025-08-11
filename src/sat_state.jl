struct SatState{C<:AbstractCorrelator,PCF<:AbstractPostCorrFilter}
    prn::Int
    code_phase::Float64
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    carrier_doppler::typeof(1.0Hz)
    integrated_samples::Int
    signal_start_sample::Int
    correlator::C
    last_fully_integrated_correlator::C
    last_fully_integrated_filtered_prompt::ComplexF64
    cn0_estimator::MomentsCN0Estimator
    bit_buffer::BitBuffer
    post_corr_filter::PCF
end

get_prn(s::SatState) = s.prn
get_num_ants(s::SatState{<:AbstractCorrelator{M}}) where {M} = M
get_code_phase(s::SatState) = s.code_phase
get_code_doppler(s::SatState) = s.code_doppler
get_carrier_phase(s::SatState) = s.carrier_phase * 2π
get_carrier_doppler(s::SatState) = s.carrier_doppler
get_integrated_samples(s::SatState) = s.integrated_samples
get_signal_start_sample(s::SatState) = s.signal_start_sample
get_correlator(s::SatState) = s.correlator
get_last_fully_integrated_correlator(s::SatState) = s.last_fully_integrated_correlator
get_last_fully_integrated_filtered_prompt(s::SatState) =
    s.last_fully_integrated_filtered_prompt
get_post_corr_filter(s::SatState) = s.post_corr_filter
get_cn0_estimator(s::SatState) = s.cn0_estimator
get_bit_buffer(s::SatState) = s.bit_buffer
@inline has_bit_or_secondary_code_been_found(s::SatState) =
    has_bit_or_secondary_code_been_found(get_bit_buffer(s))

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
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
)
    SatState(
        prn,
        float(code_phase),
        float(code_doppler),
        float(carrier_phase) / 2π,
        float(carrier_doppler),
        0,
        1,
        correlator,
        correlator,
        complex(0.0, 0.0),
        MomentsCN0Estimator(num_prompts_for_cn0_estimation),
        BitBuffer(),
        post_corr_filter,
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
    sat_state::SatState{C,PCF};
    prn = nothing,
    code_phase = nothing,
    code_doppler = nothing,
    carrier_phase = nothing,
    carrier_doppler = nothing,
    integrated_samples = nothing,
    signal_start_sample = nothing,
    correlator::Maybe{C} = nothing,
    last_fully_integrated_correlator = nothing,
    last_fully_integrated_filtered_prompt = nothing,
    cn0_estimator = nothing,
    bit_buffer = nothing,
    post_corr_filter::Maybe{PCF} = nothing,
) where {C<:AbstractCorrelator,PCF<:AbstractPostCorrFilter}
    SatState{C,PCF}(
        isnothing(prn) ? sat_state.prn : prn,
        isnothing(code_phase) ? sat_state.code_phase : code_phase,
        isnothing(code_doppler) ? sat_state.code_doppler : code_doppler,
        isnothing(carrier_phase) ? sat_state.carrier_phase : carrier_phase,
        isnothing(carrier_doppler) ? sat_state.carrier_doppler : carrier_doppler,
        isnothing(integrated_samples) ? sat_state.integrated_samples : integrated_samples,
        isnothing(signal_start_sample) ? sat_state.signal_start_sample :
        signal_start_sample,
        isnothing(correlator) ? sat_state.correlator : correlator,
        isnothing(last_fully_integrated_correlator) ?
        sat_state.last_fully_integrated_correlator : last_fully_integrated_correlator,
        isnothing(last_fully_integrated_filtered_prompt) ?
        sat_state.last_fully_integrated_filtered_prompt :
        last_fully_integrated_filtered_prompt,
        isnothing(cn0_estimator) ? sat_state.cn0_estimator : cn0_estimator,
        isnothing(bit_buffer) ? sat_state.bit_buffer : bit_buffer,
        isnothing(post_corr_filter) ? sat_state.post_corr_filter : post_corr_filter,
    )
end

function reset_start_sample_and_bit_buffer(sat_state)
    SatState(sat_state; signal_start_sample = 1, bit_buffer = reset(sat_state.bit_buffer))
end

struct SystemSatsState{S<:AbstractGNSS,SS<:SatState,I}
    system::S
    states::Dictionary{I,SS}
end

const MultipleSystemSatsState{N,I,S,SS} =
    TupleLike{<:NTuple{N,SystemSatsState{<:S,<:SS,<:I}}}

function merge_sats(
    multiple_system_sats_state::MultipleSystemSatsState{N},
    system_idx,
    new_sat_states::Dictionary{I,<:SatState},
) where {N,I}
    system_sats_state = get_system_sats_state(multiple_system_sats_state, system_idx)
    @set multiple_system_sats_state[system_idx].states =
        merge(system_sats_state.states, new_sat_states)
end

function filter_out_sats(
    multiple_system_sats_state::MultipleSystemSatsState,
    system_idx::Union{Symbol,Integer},
    identifiers,
)
    filtered_sat_states = map(
        last,
        filter(
            ((id,),) -> !in(id, identifiers),
            pairs(multiple_system_sats_state[system_idx].states),
        ),
    )
    @set multiple_system_sats_state[system_idx].states = filtered_sat_states
end

function reset_start_sample_and_bit_buffer(
    multiple_system_sats_state::MultipleSystemSatsState,
)
    map(multiple_system_sats_state) do system_sats_state
        new_sat_states = map(reset_start_sample_and_bit_buffer, system_sats_state.states)
        SystemSatsState(system_sats_state, new_sat_states)
    end
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

function SystemSatsState(system::AbstractGNSS, states;)
    SystemSatsState(system, to_dictionary(states))
end

function SystemSatsState(
    system_sats_state::SystemSatsState,
    states::Dictionary{I,<:SatState};
) where {I}
    SystemSatsState(system_sats_state.system, states)
end

get_system(sss::SystemSatsState) = sss.system
get_states(sss::SystemSatsState) = sss.states
get_sat_state(sss::SystemSatsState, identifier) = sss.states[identifier]

function estimate_cn0(sss::SystemSatsState, sat_identifier)
    system = sss.system
    estimate_cn0(
        get_cn0_estimator(get_sat_state(sss, sat_identifier)),
        get_code_length(system) / get_code_frequency(system),
    )
end
