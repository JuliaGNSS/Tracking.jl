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
        zero(Complex{Float64}),
        zero(Complex{Float64}),
        zero(Complex{Float64})
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
        zero(SVector{N, Complex{Float64}}),
        zero(SVector{N, Complex{Float64}}),
        zero(SVector{N, Complex{Float64}})
    )
end

"""
$(SIGNATURES)

Get number of antennas from correlator
"""
get_num_ants(correlator::EarlyPromptLateCorrelator{Complex{T}}) where T = 1
get_num_ants(correlator::EarlyPromptLateCorrelator{SVector{N, Complex{T}}}) where {N, T} = N

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

Calculate the replica phase offset required for the correlator taps with
respect to the prompt correlator, expressed in samples. The shifts are 
ordered from latest to earliest replica.
"""
function get_correlator_sample_shifts(
    system::AbstractGNSS,
    correlator::EarlyPromptLateCorrelator,
    sampling_frequency,
    preferred_code_shift
)
    @SVector [i*round(Int, preferred_code_shift * sampling_frequency / get_code_frequency(system)) for i in -1:1]
end

function get_early_late_sample_spacing(correlator::EarlyPromptLateCorrelator, correlator_sample_shifts::SVector{3})
    correlator_sample_shifts[3] - correlator_sample_shifts[1]
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

Perform a correlation.
"""
function correlate(
    correlator::EarlyPromptLateCorrelator,
    downconverted_signal,
    code,
    correlator_sample_shifts::SVector{3},
    start_sample,
    num_samples_left,
    agc_attenuation,
    agc_bits,
    carrier_bits::Val{NC}
) where NC
    late = zero(Complex{Int32})
    prompt = zero(Complex{Int32})
    early = zero(Complex{Int32})
    @inbounds for i = start_sample:num_samples_left + start_sample - 1
        late = late + downconverted_signal[i] * code[i]
    end
    @inbounds for i = start_sample:num_samples_left + start_sample - 1
        prompt = prompt + downconverted_signal[i] * code[i + correlator_sample_shifts[2]-correlator_sample_shifts[1]]
    end
    @inbounds for i = start_sample:num_samples_left + start_sample - 1
        early = early + downconverted_signal[i] * code[i + correlator_sample_shifts[3]-correlator_sample_shifts[1]]
    end
    EarlyPromptLateCorrelator(
        get_early(correlator) + early * agc_attenuation / 1 << (agc_bits + NC),
        get_prompt(correlator) + prompt * agc_attenuation / 1 << (agc_bits + NC),
        get_late(correlator) + late * agc_attenuation / 1 << (agc_bits + NC)
    )
end

function correlate(
    correlator::EarlyPromptLateCorrelator{<: SVector{N}},
    downconverted_signal::AbstractMatrix,
    code,
    correlator_sample_shifts::SVector{3},
    start_sample,
    num_samples_left,
    agc_attenuation,
    agc_bits,
    carrier_bits::Val{NC}
) where {N, NC}
    late = zero(MVector{N, Complex{Int32}})
    prompt = zero(MVector{N, Complex{Int32}})
    early = zero(MVector{N, Complex{Int32}})
    @inbounds for j = 1:length(late), i = start_sample:num_samples_left + start_sample - 1
        late[j] = late[j] + downconverted_signal[i,j] * code[i]
    end
    @inbounds for j = 1:length(late), i = start_sample:num_samples_left + start_sample - 1
        prompt[j] = prompt[j] + downconverted_signal[i,j] * code[i + correlator_sample_shifts[2] - correlator_sample_shifts[1]]
    end
    @inbounds for j = 1:length(late), i = start_sample:num_samples_left + start_sample - 1
        early[j] = early[j] + downconverted_signal[i,j] * code[i + correlator_sample_shifts[3] - correlator_sample_shifts[1]]
    end
    EarlyPromptLateCorrelator(
        get_early(correlator) + early .* agc_attenuation / 1 << (agc_bits + NC),
        get_prompt(correlator) + prompt .* agc_attenuation / 1 << (agc_bits + NC),
        get_late(correlator) + late .* agc_attenuation / 1 << (agc_bits + NC)
    )
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
        zero(Complex{Float64}),
        zero(Complex{Float64}),
        zero(Complex{Float64}),
        zero(Complex{Float64}),
        zero(Complex{Float64})
    )
end

function VeryEarlyPromptLateCorrelator(num_ants::NumAnts{N}) where N
    VeryEarlyPromptLateCorrelator(
        zero(SVector{N, Complex{Float64}}),
        zero(SVector{N, Complex{Float64}}),
        zero(SVector{N, Complex{Float64}}),
        zero(SVector{N, Complex{Float64}}),
        zero(SVector{N, Complex{Float64}})
    )
end

"""
$(SIGNATURES)

Get number of antennas from correlator
"""
get_num_ants(correlator::VeryEarlyPromptLateCorrelator{Complex{T}}) where T = 1
function get_num_ants(
    correlator::VeryEarlyPromptLateCorrelator{SVector{N, Complex{T}}}
) where {N, T}
    N
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
) where S <: AbstractGNSS
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

include("correlators/generic_correlator.jl")