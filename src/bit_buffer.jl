"""
$(SIGNATURES)

BitBuffer to buffer bits.

The `code_block_buffer` field is the sync-search sliding window — its
width `B` is chosen per signal by [`get_code_block_buffer_type`](@ref) so
that a single integer can hold the entire pre-sync search horizon (one
NH10 period for GPS L5I, 40 primary blocks for GPS L1 C/A, 1800 chips
for the GPS L1C-P overlay, etc.). After sync the field is dead state and
the decoded navigation bits accumulate in the fixed-width
`buffer::UInt128` instead.
"""
struct BitBuffer{B<:Unsigned}
    code_block_buffer::B
    code_block_buffer_lengh::Int
    found::Bool
    buffer::UInt128
    length::Int
    prompt_accumulator::ComplexF64
    prompt_accumulator_integrated_code_blocks::Int
end

# Default constructor preserves the pre-refactor `UInt128`-backed search
# buffer. Once `get_code_block_buffer_type` lands (Step 2) the per-signal
# `TrackedSignal` constructor picks the right width instead.
function BitBuffer()
    BitBuffer{UInt128}(zero(UInt128), 0, false, zero(UInt128), 0, complex(0.0, 0.0), 0)
end

# Typed empty constructor used by the per-signal `TrackedSignal` path.
function BitBuffer{B}() where {B<:Unsigned}
    BitBuffer{B}(zero(B), 0, false, zero(UInt128), 0, complex(0.0, 0.0), 0)
end

# Convenience outer constructor for the 7-arg form: pins `B` from the
# first argument and converts the rest to the field types. Used by test
# code that builds a `BitBuffer` from raw integer / Complex{Int} literals.
function BitBuffer(
    code_block_buffer::B,
    code_block_buffer_lengh::Integer,
    found::Bool,
    buffer::Integer,
    length::Integer,
    prompt_accumulator::Complex,
    prompt_accumulator_integrated_code_blocks::Integer,
) where {B<:Unsigned}
    BitBuffer{B}(
        code_block_buffer,
        Int(code_block_buffer_lengh),
        found,
        UInt128(buffer),
        Int(length),
        ComplexF64(prompt_accumulator),
        Int(prompt_accumulator_integrated_code_blocks),
    )
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
function buffer(signal::AbstractGNSSSignal, bit_buffer::BitBuffer{B}, integrated_code_blocks, prompt) where {B<:Unsigned}
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
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            get_bits(bit_buffer) << 1 + UInt64(bit),
            length(bit_buffer) + 1,
            zero(prompt_accumulator),
            0,
        )
    else
        return BitBuffer{B}(
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

function _buffer_find_bit(signal, bit_buffer::BitBuffer{B}, num_code_blocks_that_form_a_bit,
                          integrated_code_blocks, prompt) where {B<:Unsigned}
    if (integrated_code_blocks != 1)
        error(
            "The number code blocks must be equal to 1 if bit or secondary code hasn't been found yet.",
        )
    end
    code_block_buffer = bit_buffer.code_block_buffer << 1 + B(real(prompt) > 0)
    code_block_buffer_lengh = bit_buffer.code_block_buffer_lengh + 1
    bit_found = is_upcoming_integration_new_bit(
        signal,
        code_block_buffer,
        code_block_buffer_lengh,
    )
    if (bit_found == false)
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_lengh,
            false,
            zero(UInt128),
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
            ((code_block_buffer & (one(B) << buffer_code_block_index)) > 0) * 2 - 1
        end
        bits << 1 + (bit_sum > 0)
    end
    return BitBuffer{B}(
        code_block_buffer,
        code_block_buffer_lengh,
        true,
        bits,
        num_bits,
        complex(0, 0),
        0,
    )
end

function reset(bit_buffer::BitBuffer{B}) where {B<:Unsigned}
    BitBuffer{B}(
        bit_buffer.code_block_buffer,
        bit_buffer.code_block_buffer_lengh,
        bit_buffer.found,
        zero(UInt128),
        0,
        bit_buffer.prompt_accumulator,
        bit_buffer.prompt_accumulator_integrated_code_blocks,
    )
end
