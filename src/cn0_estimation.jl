abstract type AbstractCN0Estimator end

"""
$(SIGNATURES)

MomentsCN0Estimator to estimate the CN0
"""
struct MomentsCN0Estimator <: AbstractCN0Estimator
    prompt_buffer::Vector{ComplexF64}
    buffer_current_index::Int
    filled_buffer_length::Int
end

function MomentsCN0Estimator(N)
    MomentsCN0Estimator(zeros(ComplexF64, N), 0, 0)
end

length(estimator::MomentsCN0Estimator) = estimator.filled_buffer_length
get_prompt_buffer(estimator::MomentsCN0Estimator) = estimator.prompt_buffer
get_current_index(estimator::MomentsCN0Estimator) = estimator.buffer_current_index

"""
$(SIGNATURES)

Buffers the prompts such that they can be used to estimate the CN0
"""
function update(estimator::MomentsCN0Estimator, prompt)
    buffer_length = length(estimator.prompt_buffer)
    next_index = mod(get_current_index(estimator), buffer_length) + 1
    estimator.prompt_buffer[next_index] = prompt
    filled_buffer_length = min(length(estimator) + 1, buffer_length)
    MomentsCN0Estimator(estimator.prompt_buffer, next_index, filled_buffer_length)
end

"""
$(SIGNATURES)

Estimates the CN0 based on the struct `MomentsCN0Estimator`.
"""
function estimate_cn0(estimator::MomentsCN0Estimator, integration_time)
    length(estimator) == 0 && return 0.0dBHz
    prompt_buffer = get_prompt_buffer(estimator)
    M₂ = 1 / length(estimator) * sum(abs2, prompt_buffer)
    M₄ = 1 / length(estimator) * sum(x -> abs2(x)^2, prompt_buffer)
    Pd = sqrt(abs(2 * M₂^2 - M₄))
    noise_power = abs(M₂ - Pd)
    SNR = Pd / noise_power
    dBHz(SNR / integration_time)
end
