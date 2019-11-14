abstract type AbstractCorrelator{T} end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator for the three correlators: Early, Prompt and Late
"""
struct EarlyPromptLateCorrelator{T} <: AbstractCorrelator{T}
    early::T
    prompt::T
    late::T
end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator constructor without parameters assumes a single antenna.
"""
function EarlyPromptLateCorrelator()
    EarlyPromptLateCorrelator(NumAnts(1))
end

function EarlyPromptLateCorrelator(num_ants::NumAnts{1})
    EarlyPromptLateCorrelator(
        complex(0.0, 0.0),
        complex(0.0, 0.0),
        complex(0.0, 0.0)
    )
end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator constructor that considers multiple antennas. The number of
antennas has to be specified by `num_ants::NumAnts{N}` where N is the number of antenna
elements.
"""
function EarlyPromptLateCorrelator(num_ants::NumAnts{N}) where N
    EarlyPromptLateCorrelator(
        zero(SVector{N, ComplexF64}),
        zero(SVector{N, ComplexF64}),
        zero(SVector{N, ComplexF64})
    )
end

"""
$(SIGNATURES)

Get the early correlator
"""
@inline get_early(correlator::EarlyPromptLateCorrelator) = correlator.early

"""
$(SIGNATURES)

Get the prompt correlator
"""
@inline get_prompt(correlator::EarlyPromptLateCorrelator) = correlator.prompt

"""
$(SIGNATURES)

Get the late correlator
"""
@inline get_late(correlator::EarlyPromptLateCorrelator) = correlator.late

"""
$(SIGNATURES)

Reset the correlator
"""
function zero(correlator::EarlyPromptLateCorrelator{T}) where T
    EarlyPromptLateCorrelator(zero(T), zero(T), zero(T))
end

"""
$(SIGNATURES)

Filter the correlator by the function `post_corr_filter`.
"""
function filter(post_corr_filter, correlator::EarlyPromptLateCorrelator)
    EarlyPromptLateCorrelator(
        post_corr_filter(get_early(correlator)),
        post_corr_filter(get_prompt(correlator)),
        post_corr_filter(get_late(correlator))
    )
end

"""
$(SIGNATURES)

Calculate the shift between the early and late in samples.
"""
function get_early_late_sample_shift(
    ::Type{S},
    correlator::EarlyPromptLateCorrelator,
    sample_frequency,
    preferred_code_shift
) where S <: AbstractGNSSSystem
    round(Int, preferred_code_shift * sample_frequency / get_code_frequency(S))
end

"""
$(SIGNATURES)

Normalize the correlator.
"""
function normalize(correlator::EarlyPromptLateCorrelator, integrated_samples)
    EarlyPromptLateCorrelator(
        get_early(correlator) / integrated_samples,
        get_prompt(correlator) / integrated_samples,
        get_late(correlator) / integrated_samples
    )
end

"""
$(SIGNATURES)

Perform a correlation iteration.
"""
@inline function correlate_iteration(
    ::Type{S},
    correlator::EarlyPromptLateCorrelator,
    current_signal,
    carrier,
    code_register::Integer,
    early_late_sample_shift::Integer,
    early_code
) where S <: AbstractGNSSSystem
    prompt_code = ((code_register & 1 << early_late_sample_shift) > 0) * 2 - 1
    late_code = ((code_register & 1 << (2 * early_late_sample_shift)) > 0) * 2 - 1
    early = get_early(correlator) + current_signal * carrier * early_code
    prompt = get_prompt(correlator) + current_signal * carrier * prompt_code
    late = get_late(correlator) + current_signal * carrier * late_code
    EarlyPromptLateCorrelator(early, prompt, late)
end

struct VeryEarlyPromptLateCorrelator{T} <: AbstractCorrelator{T}
    very_early::T
    early::T
    prompt::T
    late::T
    very_late::T
end

function VeryEarlyPromptLateCorrelator(num_ants::NumAnts{1})
    VeryEarlyPromptLateCorrelator(
        complex(0.0, 0.0),
        complex(0.0, 0.0),
        complex(0.0, 0.0),
        complex(0.0, 0.0),
        complex(0.0, 0.0)
    )
end

function VeryEarlyPromptLateCorrelator(num_ants::NumAnts{N}) where N
    VeryEarlyPromptLateCorrelator(
        zero(SVector{N, ComplexF64}),
        zero(SVector{N, ComplexF64}),
        zero(SVector{N, ComplexF64}),
        zero(SVector{N, ComplexF64}),
        zero(SVector{N, ComplexF64})
    )
end

@inline get_very_early(correlator::VeryEarlyPromptLateCorrelator) = correlator.early
@inline get_early(correlator::VeryEarlyPromptLateCorrelator) = correlator.early
@inline get_prompt(correlator::VeryEarlyPromptLateCorrelator) = correlator.prompt
@inline get_late(correlator::VeryEarlyPromptLateCorrelator) = correlator.late
@inline get_very_late(correlator::VeryEarlyPromptLateCorrelator) = correlator.late

function zero(correlator::VeryEarlyPromptLateCorrelator{T}) where T
    EarlyPromptLateCorrelator(zero(T), zero(T), zero(T), zero(T), zero(T))
end

function filter(post_corr_filter, correlator::VeryEarlyPromptLateCorrelator)
    EarlyPromptLateCorrelator(
        post_corr_filter(get_very_early(correlator)),
        post_corr_filter(get_early(correlator)),
        post_corr_filter(get_prompt(correlator)),
        post_corr_filter(get_late(correlator)),
        post_corr_filter(get_very_late(correlator))
    )
end

function normalize(correlator::VeryEarlyPromptLateCorrelator, integrated_samples)
    EarlyPromptLateCorrelator(
        get_very_early(correlator) / integrated_samples,
        get_early(correlator) / integrated_samples,
        get_prompt(correlator) / integrated_samples,
        get_late(correlator) / integrated_samples,
        get_very_late(correlator) / integrated_samples
    )
end

# TODO: correlate and dump for very early, very late
#=
Base.@propagate_inbounds @inline function correlate_iteration(
    ::Type{S},
    correlator::VeryEarlyPromptLateCorrelator,
    current_signal,
    early_late_sample_shift,
    carrier,
    prn,
    total_code_length,
    prompt_code_phase
) where S <: AbstractGNSSSystem
    early_code_phase = prompt_code_phase + code_phase_delta * early_late_sample_shift
    early_code_phase += (early_code_phase < 0) * total_code_length
    late_code_phase = prompt_code_phase - code_phase_delta * early_late_sample_shift
    late_code_phase -= (late_code_phase >= total_code_length) * total_code_length
    early_code = get_code_unsafe(S, early_code_phase, prn)
    prompt_code = get_code_unsafe(S, code_phase, prn)
    late_code = get_code_unsafe(S, late_code_phase, prn)
    early = get_early(correlator) + current_signal * carrier * early_code
    prompt = get_prompt(correlator) + current_signal * carrier * prompt_code
    late = get_late(correlator) + current_signal * carrier * late_code
    VeryEarlyPromptLateCorrelator(early, prompt, late)
end
=#
