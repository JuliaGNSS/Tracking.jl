abstract type AbstractCN0Estimator end

"""
$(SIGNATURES)

MomentsCN0Estimator to estimate the CN0
"""
struct MomentsCN0Estimator{N} <: AbstractCN0Estimator
    prompt_buffer::SVector{N, ComplexF64}
    current_index::Int
    length::Int
end

function MomentsCN0Estimator(N)
    MomentsCN0Estimator{N}(zero(SVector{N, ComplexF64}), 0, 0)
end

length(cn0_state::MomentsCN0Estimator) = cn0_state.length
get_prompt_buffer(cn0_state::MomentsCN0Estimator) = cn0_state.prompt_buffer
get_current_index(cn0_state::MomentsCN0Estimator) = cn0_state.current_index

"""
$(SIGNATURES)

Updates the `cn0_estimator` to include the information of the current prompt value.
"""
function update(cn0_estimator::MomentsCN0Estimator{N}, prompt) where N
    next_index = mod(get_current_index(cn0_estimator), N) + 1
    next_prompt_buffer = setindex(get_prompt_buffer(cn0_estimator), prompt, next_index)
    next_length = min(length(cn0_estimator) + 1, N)
    MomentsCN0Estimator{N}(next_prompt_buffer, next_index, next_length)
end

"""
$(SIGNATURES)

Estimates the CN0 based on the struct `cn0_estimator`.
"""
function estimate_cn0(cn0_estimator::MomentsCN0Estimator, integration_time)
    length(cn0_estimator) == 0 && return 0.0dBHz
    abs2_prompt_buffer = abs2.(get_prompt_buffer(cn0_estimator))
    M₂ = 1 / length(cn0_estimator) * sum(abs2_prompt_buffer)
    M₄ = 1 / length(cn0_estimator) * sum(abs2_prompt_buffer .^ 2)
    Pd = sqrt(abs(2 * M₂^2 - M₄))
    SNR = Pd / (M₂ - Pd)
    dBHz(SNR / integration_time)
end
