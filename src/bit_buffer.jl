"""
$(SIGNATURES)

BitBuffer to buffer bits
"""
struct BitBuffer
    code_block_buffer::UInt128
    code_block_buffer_lengh::Int
    found::Bool
    buffer::UInt128
    length::Int
    prompt_accumulator::ComplexF64
    prompt_accumulator_integrated_code_blocks::Int
end

function BitBuffer()
    BitBuffer(0, 0, false, 0, 0, complex(0.0, 0.0), 0)
end

@inline get_bits(bit_buffer::BitBuffer) = bit_buffer.buffer
@inline length(bit_buffer::BitBuffer) = bit_buffer.length
@inline has_bit_or_secondary_code_been_found(bit_buffer::BitBuffer) = bit_buffer.found

# Number of primary-code blocks that form one navigation bit.
# Returns 0 for pilot signals (`data_frequency = 0`), where the concept is
# undefined; callers must guard for that case before using the result.
@inline function _calc_num_code_blocks_that_form_a_bit(signal::AbstractGNSSSignal)
    data_freq = get_data_frequency(signal)
    iszero(data_freq) && return 0
    Int(get_code_frequency(signal) / (get_code_length(signal) * data_freq))
end

"""
$(SIGNATURES)

Buffer data bits based on the prompt accumulation and the current prompt value.
"""
function buffer(signal::AbstractGNSSSignal, bit_buffer, integrated_code_blocks, prompt)
    # The divide is deferred to the helper — pilot signals
    # (`get_data_frequency = 0`) would otherwise blow up here with `Int(Inf)`.
    # Pilots take the `_buffer_find_bit` branch with `bit_buffer.found = false`
    # forever (their per-signal `is_upcoming_integration_new_bit` returns
    # false), so the value is unused for them.
    num_code_blocks_that_form_a_bit = _calc_num_code_blocks_that_form_a_bit(signal)

    if (bit_buffer.found == false)
        return _buffer_find_bit(
            signal, bit_buffer, num_code_blocks_that_form_a_bit,
            integrated_code_blocks, prompt,
        )
    end

    prompt_accumulator = bit_buffer.prompt_accumulator + prompt
    prompt_accumulator_integrated_code_blocks =
        bit_buffer.prompt_accumulator_integrated_code_blocks + integrated_code_blocks

    if prompt_accumulator_integrated_code_blocks == num_code_blocks_that_form_a_bit
        bit = real(prompt_accumulator) > 0
        return BitBuffer(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            get_bits(bit_buffer) << 1 + UInt64(bit),
            length(bit_buffer) + 1,
            zero(prompt_accumulator),
            0,
        )
    else
        return BitBuffer(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            bit_buffer.buffer,
            bit_buffer.length,
            prompt_accumulator,
            prompt_accumulator_integrated_code_blocks,
        )
    end
end

function _buffer_find_bit(signal, bit_buffer, num_code_blocks_that_form_a_bit,
                          integrated_code_blocks, prompt)
    if (integrated_code_blocks != 1)
        error(
            "The number code blocks must be equal to 1 if bit or secondary code hasn't been found yet.",
        )
    end
    code_block_buffer = bit_buffer.code_block_buffer << 1 + UInt64(real(prompt) > 0)
    code_block_buffer_lengh = bit_buffer.code_block_buffer_lengh + 1
    bit_found = is_upcoming_integration_new_bit(
        signal,
        code_block_buffer,
        code_block_buffer_lengh,
    )
    if (bit_found == false)
        return BitBuffer(
            code_block_buffer,
            code_block_buffer_lengh,
            false,
            0,
            0,
            complex(0.0, 0.0),
            0,
        )
    end
    num_bits = min(
        div(code_block_buffer_lengh, num_code_blocks_that_form_a_bit),
        div(sizeof(code_block_buffer) * 8, num_code_blocks_that_form_a_bit),
    )
    bits = reduce(num_bits:-1:1; init = UInt128(0)) do bits, bit_index
        # Don't shadow the outer `integrated_code_blocks` argument here:
        # because this closure reassigns its own local of the same name,
        # Julia conservatively boxes the outer parameter even though this
        # closure path does not always run, leading to a per-call
        # Core.Box allocation.
        bit_sum = sum(0:num_code_blocks_that_form_a_bit-1) do code_block_index
            buffer_code_block_index =
                (bit_index - 1) * num_code_blocks_that_form_a_bit + code_block_index
            ((code_block_buffer & (1 << buffer_code_block_index)) > 0) * 2 - 1
        end
        bits << 1 + (bit_sum > 0)
    end
    return BitBuffer(
        code_block_buffer,
        code_block_buffer_lengh,
        true,
        bits,
        num_bits,
        complex(0, 0),
        0,
    )
end

function reset(bit_buffer::BitBuffer)
    BitBuffer(
        bit_buffer.code_block_buffer,
        bit_buffer.code_block_buffer_lengh,
        bit_buffer.found,
        0,
        0,
        bit_buffer.prompt_accumulator,
        bit_buffer.prompt_accumulator_integrated_code_blocks,
    )
end