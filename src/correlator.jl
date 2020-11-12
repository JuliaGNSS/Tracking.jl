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

Calculate the shift between the early and late in samples.
"""
function get_early_late_sample_shift(
    ::Type{S},
    correlator::EarlyPromptLateCorrelator,
    sampling_frequency,
    preferred_code_shift
) where S <: AbstractGNSSSystem
    round(Int, preferred_code_shift * sampling_frequency / get_code_frequency(S))
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
    early_late_sample_shift,
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
        prompt = prompt + downconverted_signal[i] * code[i + early_late_sample_shift]
    end
    @inbounds for i = start_sample:num_samples_left + start_sample - 1
        early = early + downconverted_signal[i] * code[i + 2 * early_late_sample_shift]
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
    early_late_sample_shift,
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
        prompt[j] = prompt[j] + downconverted_signal[i,j] * code[i + early_late_sample_shift]
    end
    @inbounds for j = 1:length(late), i = start_sample:num_samples_left + start_sample - 1
        early[j] = early[j] + downconverted_signal[i,j] * code[i + 2 * early_late_sample_shift]
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


"""
$(SIGNATURES)

Generic correlator holding a user defined number of correlation values.
"""
struct GenericCorrelator{T} <: AbstractCorrelator{Vector{T}}
    taps::Vector{T}
    early_index::Int64
    prompt_index::Int64
    late_index::Int64
end

"""
$(SIGNATURES)

GenericCorrelator constructor without parameters assumes a single antenna 
and tree correlator elements.
"""
function GenericCorrelator()
    GenericCorrelator(NumAnts(1))
end

"""
$(SIGNATURES)

GenericCorrelator constructor for single antenna. The number of correlator
taps is fixed to three.
"""
function GenericCorrelator(num_ants::NumAnts{1})
    GenericCorrelator(
        [zero(ComplexF64) for i = 1:3],
        1,
        2,
        3
    )
end

"""
$(SIGNATURES)

GenericCorrelator constructor that allows for the configuration of the
number of antenna elements using `num_ants::NumAnts{N}` The number of 
correlator taps is fixed to three.
"""
function GenericCorrelator(num_ants::NumAnts{N}) where N
    GenericCorrelator(
        num_ants,
        NumTaps(3),
        2
    )
end

"""
$(SIGNATURES)

GenericCorrelator constructor for single antenna correlation, that 
allows for the configuration of the number of correlator taps using
`num_taps::NumTaps{M}`. The early-late spacing is fixed to 2.
"""
function GenericCorrelator(num_ants::NumAnts{1}, num_taps::NumTaps{M}) where M
    prompt_idx = ceil(M/2)
    GenericCorrelator(
        num_ants,
        num_taps,
        2
    )
end

"""
$(SIGNATURES)

GenericCorrelator constructor that allows for the configuration of the
number of antenna elements using `num_ants::NumAnts{N}`, the number of 
correlator taps using `num_taps::NumTaps{M}`. The early-late spacing 
is fixed to 2.
"""
function GenericCorrelator(num_ants::NumAnts{N}, num_taps::NumTaps{M}) where {M, N}
    prompt_idx = ceil(M/2)
    GenericCorrelator(
        num_ants,
        num_taps,
        2
    )
end

"""
$(SIGNATURES)

GenericCorrelator constructor for single antenna correlation, that 
allows for the configuration of the number of correlator taps using 
`num_taps::NumTaps{M}` and the spacing of the early and late correlator 
using `el_spacing`.
"""
function GenericCorrelator(num_ants::NumAnts{1}, num_taps::NumTaps{M}, el_spacing) where M
    prompt_index = ceil(Int64, M/2)
    early_index  = prompt_index - ceil(Int64, el_spacing/2)
    late_index   = prompt_index + ceil(Int64, el_spacing/2)
    @assert late_index <= M
    GenericCorrelator(
        [zero(ComplexF64) for i = 1:M],
        early_index,
        prompt_index,
        late_index
    )
end

"""
$(SIGNATURES)

GenericCorrelator constructor that allows for the configuration of the
number of antenna elements using `num_ants::NumAnts{N}`, the number of 
correlator taps using `num_taps::NumTaps{M}` and the spacing of the 
early and late correlator using `el_spacing`.
"""
function GenericCorrelator(num_ants::NumAnts{N}, num_taps::NumTaps{M}, el_spacing) where {M, N}
    prompt_index = ceil(Int64, M/2)
    early_index  = prompt_index - ceil(Int64, el_spacing/2)
    late_index   = prompt_index + ceil(Int64, el_spacing/2)
    @assert late_index <= M
    GenericCorrelator(
        [zero(SVector{N, ComplexF64}) for i = 1:M],
        early_index,
        prompt_index,
        late_index
    )
end


"""
$(SIGNATURES)

Get number of antennas from correlator
"""
get_num_ants(correlator::GenericCorrelator{Vector{Complex{T}}}) where T = 1
get_num_ants(correlator::GenericCorrelator{Vector{SVector{N,T}}}) where {N,T} = N

"""
$(SIGNATURES)

Get number of correlator taps
"""
get_num_taps(correlator::GenericCorrelator) = length(correlator.taps)

"""
$(SIGNATURES)

Get early correlator index
"""
get_early_index(correlator::GenericCorrelator) = correlator.early_index

"""
$(SIGNATURES)

Get prompt correlator index
"""
get_prompt_index(correlator::GenericCorrelator) = correlator.prompt_index

"""
$(SIGNATURES)

Get late correlator index
"""
get_late_index(correlator::GenericCorrelator) = correlator.late_index

"""
$(SIGNATURES)

Get all correlator taps
"""
get_taps(correlator::GenericCorrelator) = correlator.taps

"""
$(SIGNATURES)

Get a specific correlator tap with `index` counted positive for late and 
negative for early correlators.
"""
function get_tap(correlator::GenericCorrelator, index::Integer) 
    correlator.taps[index+get_prompt_index(correlator)]
end

"""
$(SIGNATURES)

Get the early correlator
"""
function get_early(correlator::GenericCorrelator)
    correlator.taps[get_early_index(correlator)]
end

"""
$(SIGNATURES)

Get the prompt correlator
"""
function get_prompt(correlator::GenericCorrelator)
    correlator.taps[get_prompt_index(correlator)]
end

"""
$(SIGNATURES)

Get the late correlator
"""
function get_late(correlator::GenericCorrelator)
    correlator.taps[get_late_index(correlator)]
end

"""
$(SIGNATURES)

Reset the Correlator
"""
function zero(correlator::GenericCorrelator{Vector{T}}) where T
    GenericCorrelator(
        [zero(T) for i = 1:get_num_taps(correlator)],
        correlator.early_index,
        correlator.prompt_index,
        correlator.late_index
    )
end

"""
$(SIGNATURES)

Filter the correlator by the function `post_corr_filter`
"""
function filter(post_corr_filter, correlator::GenericCorrelator)
    GenericCorrelator(
        map(x->post_corr_filter(x), get_taps(correlator)),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end

"""
$(SIGNATURES)

Calculate the shift between the early and late in samples
"""
function get_early_late_sample_shift(
    ::Type{S},
    correlator::GenericCorrelator,
    sampling_frequency,
    preferred_code_shift
) where S <: AbstractGNSSSystem
    round(Int, preferred_code_shift * sampling_frequency / get_code_frequency(S))
end

"""
$(SIGNATURES)

Normalize the correlator
"""
function normalize(correlator::GenericCorrelator, integrated_samples)
    GenericCorrelator(
        map(x->x/integrated_samples, get_taps(correlator)),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end
"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::GenericCorrelator,
    downconverted_signal,
    code,
    early_late_sample_shift,
    start_sample, 
    num_samples_left,
    agc_attenuation,
    agc_bits,
    carrier_bits::Val{NC}
) where {NC}
    nt = get_num_taps(correlator)
    taps = Vector{Complex{Int32}}(undef,nt)
    for i = 1:nt
        c = zero(Complex{Int32})
        offset = (i-1)*early_late_sample_shift
        @inbounds for j = start_sample:num_samples_left + start_sample - 1
            c += downconverted_signal[j] .* code[j + offset]
        end
        taps[i] = c
    end
    return GenericCorrelator(
        get_taps(correlator) + taps .* agc_attenuation / 1 << (agc_bits + NC),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end

"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::GenericCorrelator{<: SVector{N}},
    downconverted_signal::AbstractMatrix,
    code,
    early_late_sample_shift,
    start_sample, 
    num_samples_left,
    agc_attenuation,
    agc_bits,
    carrier_bits::Val{NC}
) where {N,NC}
    nt = get_num_taps(correlator)
    taps = Vector{SVector{N,Complex{Int32}}}(undef,nt)
    for i = 1:nt
        c = zero(MVector{N, Complex{Int32}})
        offset = (i-1)*early_late_sample_shift
        @inbounds for j = start_sample:num_samples_left + start_sample - 1
            c += downconverted_signal[j,:] .* code[j + offset]
        end
        taps[i] = c
    end

    return GenericCorrelator(
        get_taps(correlator) + taps .* agc_attenuation / 1 << (agc_bits + NC),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end


