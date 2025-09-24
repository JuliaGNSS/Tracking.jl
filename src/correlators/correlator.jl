abstract type AbstractCorrelator{M} end
abstract type AbstractEarlyPromptLateCorrelator{M} <: AbstractCorrelator{M} end

type_for_num_ants(num_ants::NumAnts{1}) = ComplexF64
type_for_num_ants(num_ants::NumAnts{N}) where {N} = SVector{N,ComplexF64}

function get_initial_accumulator(
    num_ants::NumAnts,
    num_accumulators::NumAccumulators{M},
) where {M}
    zero(SVector{M,type_for_num_ants(num_ants)})
end

function get_initial_accumulator(num_ants::NumAnts, num_accumulators::Integer)
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

"""
$(SIGNATURES)

Get all correlator accumulators
"""
get_accumulators(correlator::AbstractCorrelator) = correlator.accumulators

"""
$(SIGNATURES)

Get prompt correlator index
"""
function get_prompt_index(correlator::AbstractCorrelator)
    accumulators = get_accumulators(correlator)
    div(length(accumulators) - 1, 2) + 1
end

"""
$(SIGNATURES)

Get prompt correlator
"""
function get_prompt(correlator::AbstractCorrelator)
    get_accumulators(correlator)[get_prompt_index(correlator)]
end

function get_late_accumulator_index(correlator::AbstractEarlyPromptLateCorrelator)
    max(1, get_prompt_index(correlator) - 1)
end

function get_early_accumulator_index(correlator::AbstractEarlyPromptLateCorrelator)
    min(length(get_accumulators(correlator)), get_prompt_index(correlator) + 1)
end

"""
$(SIGNATURES)

Get early correlator
"""
function get_early(correlator::AbstractEarlyPromptLateCorrelator)
    get_accumulators(correlator)[get_early_accumulator_index(correlator)]
end

"""
$(SIGNATURES)

Get late correlator
"""
function get_late(correlator::AbstractEarlyPromptLateCorrelator)
    get_accumulators(correlator)[get_late_accumulator_index(correlator)]
end

"""
$(SIGNATURES)

Calculate the total spacing between early and late correlator in samples.
"""
function get_early_late_sample_spacing(
    correlator::AbstractEarlyPromptLateCorrelator,
    sampling_frequency,
    code_frequency,
)
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    sample_shifts[get_early_accumulator_index(correlator)] -
    sample_shifts[get_late_accumulator_index(correlator)]
end

"""
$(SIGNATURES)

Zero the correlator
"""
function zero(correlator::AbstractCorrelator)
    update_accumulator(correlator, zero(correlator.accumulators))
end

"""
$(SIGNATURES)

Is zero correlator
"""
function is_zero(correlator::AbstractCorrelator)
    get_prompt(correlator)[1] == 0
end

"""
$(SIGNATURES)

Filter the correlator by the function `post_corr_filter`
"""
function apply(post_corr_filter, correlator::AbstractCorrelator)
    update_accumulator(correlator, map(post_corr_filter, get_accumulators(correlator)))
end

"""
$(SIGNATURES)

Normalize the correlator
"""
function normalize(correlator::AbstractCorrelator, integrated_samples)
    apply(x -> x / integrated_samples, correlator)
end

function calc_preferred_code_shift_to_sample_shift(
    preferred_code_shift,
    sampling_frequency,
    code_frequency,
)
    sample_shift = round(Int, preferred_code_shift * sampling_frequency / code_frequency)
    max(1, sample_shift)
end

"""
$(SIGNATURES)

Perform a correlation for single antenna systems
"""
function correlate(
    correlator::AbstractCorrelator{1},
    downconverted_signal::AbstractVector,
    sample_shifts,
    code,
    start_sample,
    num_samples,
)
    a_re = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    a_im = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    d_re = downconverted_signal.re
    d_im = downconverted_signal.im
    latest_sample_shift = minimum(sample_shifts)
    @avx for i = start_sample:num_samples+start_sample-1
        for j = 1:length(a_re)
            sample_shift = sample_shifts[j] - latest_sample_shift
            a_re[j] += d_re[i] * code[i+sample_shift]
            a_im[j] += d_im[i] * code[i+sample_shift]
        end
    end
    accumulators_result = complex.(a_re, a_im)
    update_accumulator(
        correlator,
        map(+, get_accumulators(correlator), accumulators_result),
    )
end

function zero_accumulators(accumulators::SVector, signal)
    zeros(MVector{length(accumulators),real(eltype(signal))})
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
    sample_shifts,
    code,
    start_sample,
    num_samples,
) where {M}
    a_re = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    a_im = zero_accumulators(get_accumulators(correlator), downconverted_signal)
    d_re = downconverted_signal.re
    d_im = downconverted_signal.im
    latest_sample_shift = minimum(sample_shifts)
    @avx for i = start_sample:num_samples+start_sample-1
        for k = 1:size(a_re, 2)
            for j = 1:size(a_re, 1)
                shift = sample_shifts[k] - latest_sample_shift
                a_re[j, k] += d_re[i, j] * code[shift+i]
                a_im[j, k] += d_im[i, j] * code[shift+i]
            end
        end
    end
    update_accumulator(
        correlator,
        add_to_previous(get_accumulators(correlator), a_re, a_im),
    )
end

function add_to_previous(accumulators::SVector{NC,<:SVector}, a_re, a_im) where {NC}
    SVector{NC,eltype(accumulators)}(map(+, accumulators, eachcol(complex.(a_re, a_im))))
end
function add_to_previous(accumulators::Vector{<:SVector}, a_re, a_im)
    map(+, accumulators, eachcol(complex.(a_re, a_im)))
end

function zero_accumulators(accumulators::SVector{NC,<:SVector{NA}}, signal) where {NC,NA}
    zeros(MMatrix{NA,NC,real(eltype(signal))})
end
function zero_accumulators(accumulators::Vector{<:SVector{NA}}, signal) where {NA}
    zeros(real(eltype(signal)), NA, length(accumulators))
end
