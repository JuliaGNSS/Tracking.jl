"""
$(SIGNATURES)

EarlyPromptLateCorrelator holding a user defined number of correlation values.
The code is shifted in samples. Hence, the specified code shift is actually a
preferred code shift, because depending on sampling frequency and 
code frequency the specified code shift might not be the actual code shift. It is as
close as possible, though. The algorithm makes sure that at least one sample is shifted.
"""
struct EarlyPromptLateCorrelator{M,T} <: AbstractEarlyPromptLateCorrelator{M}
    accumulators::SVector{3,T}
    preferred_early_late_to_prompt_code_shift::Float64
    actual_early_late_code_spacing::Float64
    function EarlyPromptLateCorrelator(
        accumulators::AbstractVector{SVector{M,T}},
        preferred_early_late_to_prompt_code_shift,
        actual_early_late_code_spacing = 2 * preferred_early_late_to_prompt_code_shift,
    ) where {M,T<:Complex}
        new{M,SVector{M,T}}(accumulators, preferred_early_late_to_prompt_code_shift, actual_early_late_code_spacing)
    end
    function EarlyPromptLateCorrelator(
        accumulators::AbstractVector{T},
        preferred_early_late_to_prompt_code_shift,
        actual_early_late_code_spacing = 2 * preferred_early_late_to_prompt_code_shift,
    ) where {T<:Complex}
        new{1,T}(accumulators, preferred_early_late_to_prompt_code_shift, actual_early_late_code_spacing)
    end
end

"""
$(SIGNATURES)

EarlyPromptLateCorrelator constructor.
"""
function EarlyPromptLateCorrelator(;
    num_ants::NumAnts = NumAnts(1),
    preferred_early_late_to_prompt_code_shift = 0.5,
)
    EarlyPromptLateCorrelator(
        get_initial_accumulator(num_ants, NumAccumulators(3)),
        preferred_early_late_to_prompt_code_shift,
    )
end

"""
$(SIGNATURES)

Update the correlator with new accumulators
"""
function update_accumulator(correlator::EarlyPromptLateCorrelator, accumulators)
    EarlyPromptLateCorrelator(
        accumulators,
        correlator.preferred_early_late_to_prompt_code_shift,
        correlator.actual_early_late_code_spacing,
    )
end

function set_actual_early_late_code_spacing(
    correlator::EarlyPromptLateCorrelator,
    spacing::Float64,
)
    EarlyPromptLateCorrelator(
        get_accumulators(correlator),
        correlator.preferred_early_late_to_prompt_code_shift,
        spacing,
    )
end

"""
$(SIGNATURES)

Get the exact code phase shifts in chips for each correlator tap.
"""
function get_correlator_code_shifts(correlator::EarlyPromptLateCorrelator)
    shift = correlator.preferred_early_late_to_prompt_code_shift
    SVector(-shift, 0.0, shift)
end

"""
$(SIGNATURES)

Calculate the replica phase offset required for the correlator with
respect to the prompt correlator, expressed in samples. The shifts are
ordered from latest to earliest replica.
"""
function get_correlator_sample_shifts(
    correlator::EarlyPromptLateCorrelator,
    sampling_frequency,
    code_frequency,
)
    sample_shift = calc_preferred_code_shift_to_sample_shift(
        correlator.preferred_early_late_to_prompt_code_shift,
        sampling_frequency,
        code_frequency,
    )
    SVector(-sample_shift, 0, sample_shift)
end