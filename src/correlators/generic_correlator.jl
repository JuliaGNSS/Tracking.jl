"""
$(SIGNATURES)

Generic correlator holding a user defined number of correlation values.
"""
struct GenericCorrelator{T, TA <: AbstractVector{T}} <: AbstractCorrelator{T}
    correlators::TA
    early_index::Int
    prompt_index::Int
    late_index::Int
end

"""
$(SIGNATURES)

GenericCorrelator constructor without parameters and some default parameters.
"""
function GenericCorrelator(;
    num_correlators::NumCorrelators{M} = NumCorrelators(3),
    num_ants = NumAnts(1),
    early_late_index_offset::Integer = 2
) where M
    prompt_index = ceil(Int, M / 2)
    early_index = prompt_index + early_late_index_offset >> 1
    late_index = prompt_index - early_late_index_offset >> 1
    GenericCorrelator(
        zero(SVector{M, type_for_num_ants(num_ants)}),
        early_index,
        prompt_index,
        late_index
    )
end

type_for_num_ants(num_ants::NumAnts{1}) = ComplexF64
type_for_num_ants(num_ants::NumAnts{N}) where N = SVector{N, ComplexF64}

"""
$(SIGNATURES)

GenericCorrelator constructor for large number of correlators.
"""
function GenericCorrelator(
    num_correlators::Integer;
    num_ants = NumAnts(1),
    early_late_index_offset::Integer = 2
) where M
    prompt_index = ceil(Int, num_correlators / 2)
    early_index = prompt_index + early_late_index_offset >> 1
    late_index = prompt_index - early_late_index_offset >> 1
    GenericCorrelator(
        [zero(type_for_num_ants(num_ants)) for i = 1:num_correlators],
        early_index,
        prompt_index,
        late_index
    )
end


"""
$(SIGNATURES)

Get number of antennas from correlator
"""
get_num_ants(correlator::GenericCorrelator{Complex{T}}) where {T} = 1
get_num_ants(correlator::GenericCorrelator{SVector{N,T}}) where {N,T} = N

"""
$(SIGNATURES)

Get number of correlators
"""
get_num_correlators(correlator::GenericCorrelator) = size(correlator.correlators, 1)

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

Get all correlator correlators
"""
get_correlators(correlator::GenericCorrelator) = correlator.correlators

"""
$(SIGNATURES)

Get a specific correlator with `index` counted negative for late and
positive for early correlators.
"""
function get_correlator(correlator::GenericCorrelator, index::Integer)
    correlator.correlators[index + get_prompt_index(correlator)]
end

"""
$(SIGNATURES)

Get the early correlator
"""
function get_early(correlator::GenericCorrelator)
    correlator.correlators[get_early_index(correlator)]
end

"""
$(SIGNATURES)

Get the prompt correlator
"""
function get_prompt(correlator::GenericCorrelator)
    correlator.correlators[get_prompt_index(correlator)]
end

"""
$(SIGNATURES)

Get the late correlator
"""
function get_late(correlator::GenericCorrelator)
    correlator.correlators[get_late_index(correlator)]
end

"""
$(SIGNATURES)

Reset the Correlator
"""
function zero(correlator::GenericCorrelator)
    GenericCorrelator(
        zero(correlator.correlators),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end

"""
$(SIGNATURES)

Filter the correlator by the function `post_corr_filter`
"""
function filter(post_corr_filter, correlator::GenericCorrelator)
    GenericCorrelator(
        map(x -> post_corr_filter(x), get_correlators(correlator)),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
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
    correlator::GenericCorrelator{T, <:SVector{M}},
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
    correlator::GenericCorrelator,
    sampling_frequency,
    preferred_code_shift
)
    num_corrs = floor(Int, get_num_correlators(correlator) / 2)
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
    sample_shift < 1 &&
        throw(ArgumentError("Sample shift between early and prompt / late and prompt is less than 1."))
    sample_shift
end

"""
$(SIGNATURES)

Calculate the total spacing between early and late correlator in samples.
"""
function get_early_late_sample_spacing(
    correlator,
    correlator_sample_shifts
)
    correlator_sample_shifts[get_early_index(correlator)] -
    correlator_sample_shifts[get_late_index(correlator)]
end

"""
$(SIGNATURES)

Normalize the correlator
"""
function normalize(correlator::GenericCorrelator, integrated_samples)
    filter(x -> x / integrated_samples, correlator)
end
"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::GenericCorrelator,
    downconverted_signal::AbstractVector,
    code,
    correlator_sample_shifts,
    start_sample,
    num_samples
)
    correlators = zeros_correlator(get_correlators(correlator), downconverted_signal)
    @inbounds @fastmath for i = start_sample:num_samples + start_sample - 1
        for j = 1:length(correlators)
            sample_shift = correlator_sample_shifts[j] - correlator_sample_shifts[1]
            correlators[j] += downconverted_signal[i] * code[i + sample_shift]
        end
    end

    return GenericCorrelator(
        map(+, get_correlators(correlator), correlators),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end

function zeros_correlator(correlators::SVector, signal)
    zeros(MVector{length(correlators), eltype(signal)})
end
function zeros_correlator(correlators::Vector, signal)
    zeros(eltype(signal), length(correlators))
end

"""
$(SIGNATURES)

Perform a correlation for multi antenna systems
"""
function correlate(
    correlator::GenericCorrelator{<: SVector{N}},
    downconverted_signal::AbstractMatrix,
    code,
    correlator_sample_shifts,
    start_sample,
    num_samples,
) where N

    correlators = map(correlator_sample_shifts) do correlator_sample_shift
        correlate_single_tap(
            NumAnts(N),
            correlator_sample_shift - correlator_sample_shifts[1],
            start_sample,
            num_samples,
            downconverted_signal,
            code
        )
    end
    
    GenericCorrelator(
        map(+, get_correlators(correlator), correlators),
        get_early_index(correlator),
        get_prompt_index(correlator),
        get_late_index(correlator)
    )
end

function correlate_single_tap(
    ::NumAnts{N},
    offset,
    start_sample,
    num_samples,
    downconverted_signal,
    code
) where N
    correlator = zero(MVector{N, eltype(downconverted_signal)})
    @inbounds @fastmath for i = start_sample:num_samples + start_sample - 1, j = 1:length(correlator)
        correlator[j] += downconverted_signal[i,j] * code[i + offset]
    end
    SVector(correlator)
end
