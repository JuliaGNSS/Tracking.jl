"""
$(SIGNATURES)

BitBuffer to buffer bits
"""
struct BitBuffer
    buffer::UInt64
    length::Int
end

function BitBuffer()
    BitBuffer(0, 0)
end

@inline get_bits(bit_buffer::BitBuffer) = bit_buffer.buffer
@inline length(bit_buffer::BitBuffer) = bit_buffer.length

"""
$(SIGNATURES)

Buffer data bits based on the prompt accumulation and the current prompt value.
"""
function buffer(
    system,
    bit_buffer,
    prompt_accumulator,
    secondary_code_or_bit_found,
    prev_code_phase,
    code_phase,
    integration_time,
    prompt_correlator
)
    prompt_accumulator = prompt_accumulator + secondary_code_or_bit_found *
        prompt_correlator

    if secondary_code_or_bit_found &&
        (code_phase - prev_code_phase < 0 || integration_time == 1 / get_data_frequency(system))
        bit = real(prompt_accumulator) > 0
        bit_buffer = BitBuffer(
            get_bits(bit_buffer) << 1 + UInt64(bit),
            length(bit_buffer) + 1
        )
        prompt_accumulator = zero(prompt_accumulator)
    end
    bit_buffer, prompt_accumulator
end
