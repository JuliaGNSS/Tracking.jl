"""
$(SIGNATURES)

BitBuffer to buffer bits
"""
struct BitBuffer
    buffer::UInt64
    length::Int
    prompt_accumulator::ComplexF64
    prompt_accumulator_integrated_code_blocks::Int
end

function BitBuffer()
    BitBuffer(0, 0, complex(0.0, 0.0), 0)
end

@inline get_bits(bit_buffer::BitBuffer) = bit_buffer.buffer
@inline length(bit_buffer::BitBuffer) = bit_buffer.length

"""
$(SIGNATURES)

Buffer data bits based on the prompt accumulation and the current prompt value.
"""
function buffer(
    system::AbstractGNSS,
    bit_buffer,
    integrated_code_blocks,
    secondary_code_or_bit_found,
    prompt,
)
    prompt_accumulator =
        bit_buffer.prompt_accumulator + secondary_code_or_bit_found * prompt
    prompt_accumulator_integrated_code_blocks =
        bit_buffer.prompt_accumulator_integrated_code_blocks +
        secondary_code_or_bit_found * integrated_code_blocks

    num_code_blocks_that_form_a_bit = Int(
        get_code_frequency(system) / (get_code_length(system) * get_data_frequency(system)),
    )

    return if secondary_code_or_bit_found &&
              prompt_accumulator_integrated_code_blocks == num_code_blocks_that_form_a_bit
        bit = real(prompt_accumulator) > 0
        BitBuffer(
            get_bits(bit_buffer) << 1 + UInt64(bit),
            length(bit_buffer) + 1,
            zero(prompt_accumulator),
            0,
        )
    else
        BitBuffer(
            bit_buffer.buffer,
            bit_buffer.length,
            prompt_accumulator,
            prompt_accumulator_integrated_code_blocks,
        )
    end
end
