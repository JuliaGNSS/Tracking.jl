abstract type AbstractCorrelator{M} end
abstract type AbstractEarlyPromptLateCorrelator{M} <: AbstractCorrelator{M} end
"""
$(SIGNATURES)

EarlyPromptLateCorrelator holding a user defined number of correlation values.
"""
struct EarlyPromptLateCorrelator{
    M,
    T,
    TA <: AbstractVector{T},
    S <: AbstractVector{Int}
} <: AbstractEarlyPromptLateCorrelator{M}
    accumulators::TA
    shifts::S
    early_shift::Int
    prompt_shift::Int
    late_shift::Int
    function EarlyPromptLateCorrelator(
        accumulators::AbstractVector{SVector{M, T}},
        shifts::S,
        early_shift::Int,
        prompt_shift::Int,
        late_shift::Int,
    ) where {M, T <: Complex, S}
        idxs = SVector(early_shift, prompt_shift, late_shift)
        length(accumulators) < 3 ?
            error("Early-Prompt-Late-Correlator needs at least 3 distinct accumulators") :
            length(shifts) != length(accumulators) ?
                error("Number of accumulators and shifts must be the same.") :
                !allunique(shifts[idxs]) ?
                    error("Early, Prompt and Late must not have the same sample shift. Either increase sample rate or adjust the preferred code shift.") :
                    new{M, SVector{M, T}, typeof(accumulators), S}(accumulators, shifts, early_shift, prompt_shift, late_shift)
    end
    function EarlyPromptLateCorrelator(
        accumulators::AbstractVector{T},
        shifts::S,
        early_shift::Int,
        prompt_shift::Int,
        late_shift::Int,
    ) where {T <: Complex, S}
        idxs = SVector(early_shift, prompt_shift, late_shift)
        length(accumulators) < 3 ?
            error("Early-Prompt-Late-Correlator needs at least 3 distinct accumulators") :
            length(shifts) != length(accumulators) ?
                error("Number of accumulators and shifts must be the same.") :
                !allunique(shifts[idxs]) ?
                    error("Early, Prompt and Late must not have the same sample shift. Either increase sample rate or adjust the preferred code shift.") :
                    new{1, T, typeof(accumulators), S}(accumulators, shifts, early_shift, prompt_shift, late_shift)
    end
end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator constructor without parameters and some default parameters.
"""
function EarlyPromptLateCorrelator(
    system::AbstractGNSS,
    sampling_frequency;
    preferred_code_shift = 0.5,
    preferred_early_late_to_prompt_code_shift = 0.5, # TODO: At the momemnt this must be similar to preferred_code_shift
    num_ants::NumAnts = NumAnts(1),
    num_accumulators = NumAccumulators(3),
)
    shifts = get_correlator_sample_shifts(
        system,
        get_num_accumulators(num_accumulators),
        sampling_frequency,
        preferred_code_shift
    )
    prompt_shift = length(shifts) >> 1 + 1

    sample_shift = preferred_code_shift_to_sample_shift(
        preferred_early_late_to_prompt_code_shift,
        sampling_frequency,
        system
    )
    _, late_shift = findmin(abs, shifts .+ sample_shift)
    _, early_shift = findmin(abs, shifts .- sample_shift)

    EarlyPromptLateCorrelator(
        get_initial_accumulator(num_ants, num_accumulators),
        shifts,
        early_shift,
        prompt_shift,
        late_shift
    )
end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator holding a user defined number of correlation values.
"""
struct VeryEarlyPromptLateCorrelator{
    M,
    T,
    TA <: AbstractVector{T},
    S <: AbstractVector{Int}
} <: AbstractEarlyPromptLateCorrelator{M}
    accumulators::TA
    shifts::S
    very_early_shift::Int
    early_shift::Int
    prompt_shift::Int
    late_shift::Int
    very_late_shift::Int
    function VeryEarlyPromptLateCorrelator(
        accumulators::AbstractVector{SVector{M, T}},
        shifts::S,
        very_early_shift::Int,
        early_shift::Int,
        prompt_shift::Int,
        late_shift::Int,
        very_late_shift::Int
    ) where {M, T <: Complex, S}
        idxs = SVector(very_early_shift, early_shift, prompt_shift, late_shift, very_late_shift)
        length(accumulators) < 5 ?
            error("Very-Early-Prompt-Late-Correlator needs at least 5 distinct accumulators") :
            length(shifts) != length(accumulators) ?
                error("Number of accumulators and shifts must be the same.") :
                !allunique(shifts[idxs]) ?
                    error("Early, Prompt and Late must not have the same sample shift. Either increase sample rate or adjust the preferred code shift.") :
                    new{M, SVector{M, T}, typeof(accumulators), S}(accumulators, shifts, very_early_shift, early_shift, prompt_shift, late_shift, very_late_shift)
    end
    function VeryEarlyPromptLateCorrelator(
        accumulators::AbstractVector{T},
        shifts::S,
        very_early_shift::Int,
        early_shift::Int,
        prompt_shift::Int,
        late_shift::Int,
        very_late_shift::Int
    ) where {T <: Complex, S}
        idxs = SVector(very_early_shift, early_shift, prompt_shift, late_shift, very_late_shift)
        length(accumulators) < 5 ?
            error("Very-Early-Prompt-Late-Correlator needs at least 5 distinct accumulators") :
            length(shifts) != length(accumulators) ?
                error("Number of accumulators and shifts must be the same.") :
                !allunique(shifts[idxs]) ?
                    error("Early, Prompt and Late must not have the same sample shift. Either increase sample rate or adjust the preferred code shift.") :
                    new{1, T, typeof(accumulators), S}(accumulators, shifts, very_early_shift, early_shift, prompt_shift, late_shift, very_late_shift)
    end
end

"""
$(SIGNATURES)

VeryEarlyPromptLateCorrelator constructor without parameters and some default parameters.
Default parameters take from https://gnss-sdr.org/docs/sp-blocks/tracking/#implementation-galileo_e1_dll_pll_veml_tracking
"""
function VeryEarlyPromptLateCorrelator(
    system::AbstractGNSS,
    sampling_frequency;
    preferred_early_late_to_prompt_code_shift = 0.15,
    preferred_very_early_late_to_prompt_code_shift = 0.6,
    num_ants::NumAnts = NumAnts(1),
)
    early_late_sample_shift = preferred_code_shift_to_sample_shift(
        preferred_early_late_to_prompt_code_shift,
        sampling_frequency,
        system
    )
    very_early_late_sample_shift = preferred_code_shift_to_sample_shift(
        preferred_very_early_late_to_prompt_code_shift,
        sampling_frequency,
        system
    )
    shifts = SVector(
        -very_early_late_sample_shift,
        -early_late_sample_shift,
        0,
        early_late_sample_shift,
        very_early_late_sample_shift
    )

    VeryEarlyPromptLateCorrelator(
        get_initial_accumulator(num_ants, NumAccumulators(5)),
        shifts,
        5,
        4,
        3,
        2,
        1
    )
end

type_for_num_ants(num_ants::NumAnts{1}) = ComplexF64
type_for_num_ants(num_ants::NumAnts{N}) where N = SVector{N, ComplexF64}

function get_initial_accumulator(
    num_ants::NumAnts,
    num_accumulators::NumAccumulators{M}
) where M
    zero(SVector{M, type_for_num_ants(num_ants)})
end

function get_initial_accumulator(
    num_ants::NumAnts,
    num_accumulators::Integer
)
    [zero(type_for_num_ants(num_ants)) for i = 1:num_accumulators]
end


"""
$(SIGNATURES)

Get number of antennas from correlator
"""
get_num_ants(correlator::AbstractCorrelator{M}) where {M} = M

"""
$(SIGNATURES)

Get number of accumulators
"""
get_num_accumulators(correlator::AbstractCorrelator) = size(correlator.accumulators, 1)
get_num_accumulators(num_accumulators::Int) = num_accumulators
get_num_accumulators(num_accumulators::NumAccumulators{M}) where M = M

"""
$(SIGNATURES)

Get early correlator index
"""
@inline get_very_early_index(corr::VeryEarlyPromptLateCorrelator) = corr.very_early_shift

"""
$(SIGNATURES)

Get early correlator index
"""
@inline get_early_index(corr::AbstractEarlyPromptLateCorrelator) = corr.early_shift

"""
$(SIGNATURES)

Get prompt correlator index
"""
@inline get_prompt_index(corr::AbstractEarlyPromptLateCorrelator) = corr.prompt_shift

"""
$(SIGNATURES)

Get late correlator index
"""
@inline get_late_index(corr::AbstractEarlyPromptLateCorrelator) = corr.late_shift

"""
$(SIGNATURES)

Get late correlator index
"""
@inline get_very_late_index(corr::VeryEarlyPromptLateCorrelator) = corr.very_late_shift

"""
$(SIGNATURES)

Get all correlator accumulators
"""
get_accumulators(correlator::AbstractCorrelator) = correlator.accumulators

"""
$(SIGNATURES)

Get a specific accumulator with `index` counted negative for late and
positive for early accumulators.
"""
function get_accumulator(correlator::AbstractCorrelator, index::Integer)
    correlator.accumulators[index + get_prompt_index(correlator)]
end

"""
$(SIGNATURES)

Get very early correlator
"""
function get_very_early(correlator::VeryEarlyPromptLateCorrelator)
    correlator.accumulators[get_very_early_index(correlator)]
end

"""
$(SIGNATURES)

Get early correlator
"""
function get_early(correlator::AbstractEarlyPromptLateCorrelator)
    correlator.accumulators[get_early_index(correlator)]
end

"""
$(SIGNATURES)

Get prompt correlator
"""
function get_prompt(correlator::AbstractCorrelator)
    correlator.accumulators[get_prompt_index(correlator)]
end

"""
$(SIGNATURES)

Get late correlator
"""
function get_late(correlator::AbstractEarlyPromptLateCorrelator)
    correlator.accumulators[get_late_index(correlator)]
end

"""
$(SIGNATURES)

Get very late correlator
"""
function get_very_late(correlator::VeryEarlyPromptLateCorrelator)
    correlator.accumulators[get_very_late_index(correlator)]
end

"""
$(SIGNATURES)

Update the Correlator
"""
function update_accumulator(correlator::EarlyPromptLateCorrelator, accumulators) 
    EarlyPromptLateCorrelator(
        accumulators,
        correlator.shifts,
        correlator.early_shift,
        correlator.prompt_shift,
        correlator.late_shift
    )
end
function update_accumulator(correlator::VeryEarlyPromptLateCorrelator, accumulators) 
    VeryEarlyPromptLateCorrelator(
        accumulators,
        correlator.shifts,
        correlator.very_early_shift,
        correlator.early_shift,
        correlator.prompt_shift,
        correlator.late_shift,
        correlator.very_late_shift
    )
end

"""
$(SIGNATURES)

Zero the Correlator
"""
function zero(correlator::AbstractCorrelator)
    update_accumulator(correlator, zero(correlator.accumulators))
end

"""
$(SIGNATURES)

Filter the Correlator by the function `post_corr_filter`
"""
function apply(post_corr_filter, correlator::AbstractCorrelator)
    update_accumulator(correlator, map(post_corr_filter, get_accumulators(correlator)))
end

"""
$(SIGNATURES)

Calculate the replica phase offset required for the correlator with
respect to the prompt correlator, expressed in samples. The shifts are
ordered from latest to earliest replica.
"""
function get_correlator_sample_shifts(
    system::AbstractGNSS,
    num_correlators::Int,
    sampling_frequency,
    preferred_code_shift
)
    half_num_correlators = num_correlators >> 1
    sample_shift = preferred_code_shift_to_sample_shift(
        preferred_code_shift,
        sampling_frequency,
        system
    )
    (-half_num_correlators:half_num_correlators) .* sample_shift
end

function preferred_code_shift_to_sample_shift(
    preferred_code_shift,
    sampling_frequency,
    system
)
    sample_shift = round(Int, preferred_code_shift * sampling_frequency / get_code_frequency(system))
    max(1, sample_shift)
end

"""
$(SIGNATURES)

Calculates the index for the early and late samples
"""
function get_early_late_index_shift(
    system,
    correlator_sample_shifts,
    sampling_frequency,
    preferred_code_shift
)
    sample_shift = preferred_code_shift_to_sample_shift(
        preferred_code_shift,
        sampling_frequency,
        system
    )
    min_val, min_idx = findmin(abs, correlator_sample_shifts .- sample_shift)
    min_idx
end

"""
$(SIGNATURES)

Calculate the total spacing between early and late correlator in samples.
"""
function get_early_late_sample_spacing(
    correlator::AbstractEarlyPromptLateCorrelator
)
    correlator.shifts[get_early_index(correlator)] -
    correlator.shifts[get_late_index(correlator)]
end

"""
$(SIGNATURES)

Normalize the correlator
"""
function normalize(correlator::AbstractCorrelator, integrated_samples)
    apply(x -> x / integrated_samples, correlator)
end
"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::AbstractCorrelator{1},
    downconverted_signal::AbstractVector,
    code,
    start_sample,
    num_samples
)
    a_re = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    a_im = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    d_re = downconverted_signal.re
    d_im = downconverted_signal.im
    @avx for i = start_sample:num_samples + start_sample - 1
        for j = 1:length(a_re)
            sample_shift = correlator.shifts[j] - correlator.shifts[1]
            a_re[j] += d_re[i] * code[i + sample_shift]
            a_im[j] += d_im[i] * code[i + sample_shift]
        end
    end
    accumulators_result = complex.(a_re, a_im)
    update_accumulator(correlator, map(+, get_accumulators(correlator), accumulators_result))
end

function zero_accumulators(accumulators::SVector, signal)
    zeros(MVector{length(accumulators), real(eltype(signal))})
end
function zero_accumulators(accumulators::Vector, signal)
    zeros(real(eltype(signal)), length(accumulators))
end

"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::AbstractCorrelator{M},
    downconverted_signal::AbstractMatrix,
    code,
    start_sample,
    num_samples,
) where {M}
    a_re = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    a_im = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    d_re = downconverted_signal.re
    d_im = downconverted_signal.im
    @avx for i = start_sample:num_samples + start_sample - 1
        for k = 1:size(a_re, 2)
            for j = 1:size(a_re, 1)
                shift = correlator.shifts[k] - correlator.shifts[1]
                a_re[j,k] += d_re[i,j] * code[shift + i]
                a_im[j,k] += d_im[i,j] * code[shift + i]
            end
        end
    end
    update_accumulator(correlator, add_to_previous(get_accumulators(correlator), a_re, a_im))
end

function add_to_previous(accumulators::SVector{NC, <:SVector}, a_re, a_im) where NC
    SVector{NC, eltype(accumulators)}(map(+, accumulators, eachcol(complex.(a_re, a_im))))
end
function add_to_previous(accumulators::Vector{<:SVector}, a_re, a_im)
    map(+, accumulators, eachcol(complex.(a_re, a_im)))
end

function zero_accumulators(accumulators::SVector{NC, <:SVector{NA}}, signal) where {NC, NA}
    zeros(MMatrix{NA, NC, real(eltype(signal))})
end
function zero_accumulators(accumulators::Vector{<:SVector{NA}}, signal) where NA
    zeros(real(eltype(signal)), NA, length(accumulators))
end
