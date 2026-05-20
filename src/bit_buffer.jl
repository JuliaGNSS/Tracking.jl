"""
    SyncResult

Outcome of a per-signal bit-sync / secondary-code-sync detector call.

Fields:

- `found::Bool` — whether the detector locked on this update.
- `phase::Int` — when `found = true`, the secondary-code chip offset of
  the buffer's *first* sample, in `0:secondary_code_length-1`. Zero for
  signals without a secondary code (the L1 C/A bit-edge case picks the
  bit boundary inside the detector, not a chip offset).
- `polarity::Int8` — `+1` or `-1`; which match orientation the detector
  locked. Carries through to the post-sync prompt accumulator so that a
  negative-polarity lock doesn't trip the downstream bit decoder.
"""
struct SyncResult
    found::Bool
    phase::Int
    polarity::Int8
end

"""
$(SIGNATURES)

Hamming-tolerance template matcher used by every per-signal
`is_upcoming_integration_new_bit` implementation.

Compares `code_block_bits & mask` against `template` for the "positive
polarity" hit and against `template ⊻ mask` for the "negative polarity"
hit. The first orientation whose Hamming distance does not exceed
`max_errors` wins; if neither fits within tolerance the result reports
`found = false`. `phase` is hard-coded to 0 — callers that need a
non-zero phase (L1C-P overlay search) compute it themselves and build
their own `SyncResult`.

Inlined so the per-signal template / mask / tolerance constants fold at
the call site.
"""
@inline function _try_match(
    code_block_bits::B,
    template::B,
    mask::B,
    max_errors::Int,
) where {B<:Unsigned}
    masked = code_block_bits & mask
    dist_pos = count_ones(masked ⊻ template)
    dist_pos <= max_errors && return SyncResult(true, 0, Int8(+1))
    dist_neg = count_ones(masked ⊻ (template ⊻ mask))
    dist_neg <= max_errors && return SyncResult(true, 0, Int8(-1))
    return SyncResult(false, 0, Int8(0))
end

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

"""
$(SIGNATURES)

Width of the sync-search sliding window for `signal`, returned as a
concrete `Unsigned` subtype.

Each signal's bit/secondary-code sync detector matches the buffered
primary-code-block signs against a per-signal template. This trait
picks the smallest integer width that holds one full search horizon for
that signal:

| Signal            | Returns   | Window |
|-------------------|-----------|--------|
| GPS L1 C/A        | `UInt64`  | 40 blocks (2 × 20 blocks/symbol) |
| Galileo E1B       | `UInt8`   | 8 blocks (2 × 4 blocks/symbol)   |
| GPS L5I           | `UInt32`  | 20 blocks (2 × NH10 length)      |
| GPS L1C-D         | `UInt8`   | n/a — symbol = primary period, buffer unused |
| GPS L1C-P         | `UInt1800`| 1800-chip overlay (added in step 4) |

The default for any signal not specialized below is `UInt64`. The width
flows through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF}` so the
parameter chain stays type-stable at construction.
"""
@inline get_code_block_buffer_type(::AbstractGNSSSignal) = UInt64

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
    sync = is_upcoming_integration_new_bit(
        signal,
        code_block_buffer,
        code_block_buffer_lengh,
    )
    if !sync.found
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
