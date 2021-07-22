abstract type AbstractCorrelator{T, TA} end
"""
$(SIGNATURES)

EarlyPromptLateCorrelator holding a user defined number of correlation values.
"""
struct EarlyPromptLateCorrelator{T, TA <: AbstractVector{T}} <: AbstractCorrelator{T, TA}
    accumulators::TA
    EarlyPromptLateCorrelator{T, TA}(accumulators) where {T, TA <: AbstractVector{T}} = length(accumulators) < 3 ?
        error("Early-Prompt-Late-Correlator needs at least 3 accumulators") : new{T, TA}(accumulators)
end

function EarlyPromptLateCorrelator(accumulators::AbstractVector)
    EarlyPromptLateCorrelator{typeof(first(accumulators)), typeof(accumulators)}(accumulators)
end

type_for_num_ants(num_ants::NumAnts{1}) = ComplexF64
type_for_num_ants(num_ants::NumAnts{N}) where N = SVector{N, ComplexF64}

"""
$(SIGNATURES)

EarlyPromptLateCorrelator constructor without parameters and some default parameters.
"""
function EarlyPromptLateCorrelator(
    num_ants::NumAnts = NumAnts(1),
    num_accumulators::NumAccumulators{M} = NumAccumulators(3),
) where M
    EarlyPromptLateCorrelator(
        zero(SVector{M, type_for_num_ants(num_ants)})
    )
end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator constructor for large number of accumulators.
"""
function EarlyPromptLateCorrelator(
    num_ants::NumAnts,
    num_accumulators::Integer
)
    EarlyPromptLateCorrelator(
        [zero(type_for_num_ants(num_ants)) for i = 1:num_accumulators]
    )
end


"""
$(SIGNATURES)

Get number of antennas from correlator
"""
get_num_ants(correlator::AbstractCorrelator{Complex{T}}) where {T} = 1
get_num_ants(correlator::AbstractCorrelator{SVector{N,T}}) where {N,T} = N

"""
$(SIGNATURES)

Get number of accumulators
"""
get_num_accumulators(correlator::AbstractCorrelator) = size(correlator.accumulators, 1)

"""
$(SIGNATURES)

Get early correlator index
"""
@inline get_early_index(correlator_sample_shifts, early_late_index_shift) =
    get_prompt_index(correlator_sample_shifts) + early_late_index_shift

"""
$(SIGNATURES)

Get prompt correlator index
"""
@inline get_prompt_index(correlator_sample_shifts) = findfirst(iszero, correlator_sample_shifts)

"""
$(SIGNATURES)

Get late correlator index
"""
@inline get_late_index(correlator_sample_shifts, early_late_index_shift) =
    get_prompt_index(correlator_sample_shifts) - early_late_index_shift

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
function get_accumulator(correlator::AbstractCorrelator, correlator_sample_shifts, index::Integer)
    correlator.accumulators[index + get_prompt_index(correlator_sample_shifts)]
end

"""
$(SIGNATURES)

Get early correlator
"""
function get_early(
    correlator::AbstractCorrelator,
    correlator_sample_shifts,
    early_late_index_shift
)
    correlator.accumulators[
        get_early_index(correlator_sample_shifts, early_late_index_shift)
    ]
end

"""
$(SIGNATURES)

Get prompt correlator
"""
function get_prompt(correlator::AbstractCorrelator, correlator_sample_shifts)
    correlator.accumulators[get_prompt_index(correlator_sample_shifts)]
end
CUDA.dot
"""
$(SIGNATURES)

Get late correlator
"""
function get_late(
    correlator::AbstractCorrelator,
    correlator_sample_shifts,
    early_late_index_shift
)
    correlator.accumulators[
        get_late_index(correlator_sample_shifts, early_late_index_shift)
    ]
end

"""
$(SIGNATURES)

Zero the Correlator
"""
function zero(correlator::T) where T <: AbstractCorrelator
    T(zero(correlator.accumulators))
end

"""
$(SIGNATURES)

Filter the Correlator by the function `post_corr_filter`
"""
function filter(post_corr_filter, correlator::T) where T <: AbstractCorrelator
    (T.name.wrapper)(map(x -> post_corr_filter(x), get_accumulators(correlator)))
end

"""
$(SIGNATURES)

Calculate the replica phase offset required for the correlator with
respect to the prompt correlator, expressed in samples. The shifts are
ordered from latest to earliest replica.
"""
function get_correlator_sample_shifts(
    system::AbstractGNSS,
    correlator::AbstractCorrelator{T, <:SVector{M}},
    sampling_frequency,
    preferred_code_shift
) where {T,M}
    num_corrs = floor(Int, M / 2)
    sample_shift = preferred_code_shift_to_sample_shift(
        preferred_code_shift,
        sampling_frequency,
        system
    )
    SVector{M}(-num_corrs:num_corrs) .* sample_shift
end
function get_correlator_sample_shifts(
    system::AbstractGNSS,
    correlator::AbstractCorrelator,
    sampling_frequency,
    preferred_code_shift
)
    num_corrs = floor(Int, get_num_accumulators(correlator) / 2)
    sample_shift = preferred_code_shift_to_sample_shift(
        preferred_code_shift,
        sampling_frequency,
        system
    )
    (-num_corrs:num_corrs) .* sample_shift
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
    correlator,
    sampling_frequency,
    preferred_code_shift
)
    preferred_code_phase_distance = Inf
    early_late_index_shift = 1
    for i = get_prompt_index(correlator_sample_shifts) + 1:length(correlator_sample_shifts)
        code_phase_shift = correlator_sample_shifts[i] * get_code_frequency(system) / sampling_frequency
        preferred_code_phase_distance_temp = abs(code_phase_shift - preferred_code_shift)
        if preferred_code_phase_distance_temp < preferred_code_phase_distance
            early_late_index_shift = i - get_prompt_index(correlator_sample_shifts)
            preferred_code_phase_distance = preferred_code_phase_distance_temp
        end
    end
    early_late_index_shift
end

"""
$(SIGNATURES)

Calculate the total spacing between early and late correlator in samples.
"""
function get_early_late_sample_spacing(
    correlator_sample_shifts,
    early_late_index_shift
)
    correlator_sample_shifts[get_early_index(correlator_sample_shifts, early_late_index_shift)] -
    correlator_sample_shifts[get_late_index(correlator_sample_shifts, early_late_index_shift)]
end

"""
$(SIGNATURES)

Normalize the correlator
"""
function normalize(correlator::AbstractCorrelator, integrated_samples)
    filter(x -> x / integrated_samples, correlator)
end
"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::T,
    downconverted_signal::AbstractVector,
    code,
    correlator_sample_shifts,
    start_sample,
    num_samples
) where {T <: AbstractCorrelator}
    accumulators = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    d_re = downconverted_signal.re
    d_im = downconverted_signal.im
    a_re = real.(accumulators)
    a_im = imag.(accumulators)
    @avx for i = start_sample:num_samples + start_sample - 1
        for j = 1:length(accumulators)
            sample_shift = correlator_sample_shifts[j] - correlator_sample_shifts[1]
            a_re[j] += d_re[i] * code[i + sample_shift]
            a_im[j] += d_im[i] * code[i + sample_shift]
        end
    end
    accumulators_result = complex.(a_re, a_im)
    T(map(+, get_accumulators(correlator), accumulators_result))
end

function zero_accumulators(accumulators::SVector, signal)
    zeros(MVector{length(accumulators), eltype(signal)})
end
function zero_accumulators(accumulators::Vector, signal)
    zeros(eltype(signal), length(accumulators))
end

"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::AbstractCorrelator{<: SVector{N}},
    downconverted_signal::AbstractMatrix,
    code,
    correlator_sample_shifts,
    start_sample,
    num_samples,
) where {N}
    accumulators = zero(MMatrix{N, length(correlator_sample_shifts), eltype(downconverted_signal)})
    a_re = real.(accumulators)
    a_im = imag.(accumulators)
    d_re = downconverted_signal.re
    d_im = downconverted_signal.im
    @avx for i = start_sample:num_samples + start_sample - 1
        for k = 1:size(accumulators, 2)
            for j = 1:size(accumulators, 1)
                shift = correlator_sample_shifts[k] - correlator_sample_shifts[1]
                a_re[j,k] += d_re[i,j] * code[shift + i]
                a_im[j,k] += d_im[i,j] * code[shift + i]
            end
        end
    end

    accumulators_new = SVector(eachcol(complex.(SMatrix(a_re), SMatrix(a_im)))...)
    typeof(correlator)(map(+, get_accumulators(correlator), accumulators_new))
end
