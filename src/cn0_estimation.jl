abstract type AbstractCN0Estimator end

"""
$(SIGNATURES)

MomentsCN0Estimator to estimate the CN0
"""
struct MomentsCN0Estimator <: AbstractCN0Estimator end

"""
$(SIGNATURES)

MomentsCN0Estimator to estimate the CN0
"""
struct PromptsBuffer
    prompt_buffer::Vector{ComplexF64}
    current_index::Int
    length::Int
end

function PromptsBuffer(N)
    PromptsBuffer(zeros(ComplexF64, N), 0, 0)
end

length(buffer::PromptsBuffer) = buffer.length
get_prompt_buffer(buffer::PromptsBuffer) = buffer.prompt_buffer
get_current_index(buffer::PromptsBuffer) = buffer.current_index

"""
$(SIGNATURES)

Updates the `cn0_estimator` to include the information of the current prompt value.
"""
function update(buffer::PromptsBuffer, prompt)
    buffer_length = length(buffer.prompt_buffer)
    next_index = mod(get_current_index(buffer), buffer_length) + 1
    buffer.prompt_buffer[next_index] = prompt
    next_length = min(length(buffer) + 1, buffer_length)
    PromptsBuffer(buffer.prompt_buffer, next_index, next_length)
end

"""
$(SIGNATURES)

Estimates the CN0 based on the struct `cn0_estimator`.
"""
function estimate_cn0(
    buffer::PromptsBuffer,
    cn0_estimator::MomentsCN0Estimator,
    integration_time,
)
    length(buffer) == 0 && return 0.0dBHz
    abs2_prompt_buffer = abs2.(get_prompt_buffer(buffer))
    M₂ = 1 / length(buffer) * sum(abs2_prompt_buffer)
    M₄ = 1 / length(buffer) * sum(abs2_prompt_buffer .^ 2)
    Pd = sqrt(abs(2 * M₂^2 - M₄))
    noise_power = (M₂ - Pd)
    noise_power_non_neg = noise_power - 2 * (noise_power < 0) * noise_power
    SNR = Pd / noise_power_non_neg
    dBHz(SNR / integration_time)
end
