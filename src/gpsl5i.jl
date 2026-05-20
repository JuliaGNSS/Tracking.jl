"""
$(SIGNATURES)

Secondary-code sync detector for GPS L5I.

Matches the low 10 bits of `code_block_bits` against the NH10 template
`0x035` (= `0000110101`). The negated-polarity template `0x3ca` is
searched in the same call via [`_try_match`](@ref). Returns
[`SyncResult`](@ref).
"""
@inline function is_upcoming_integration_new_bit(
    ::GPSL5I,
    ::Integer,           # PRN — ignored; NH10 is shared across PRNs
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    num_code_blocks < 10 && return SyncResult(false, 0, Int8(0))
    # Tolerance 2 ≈ 10 % per-block error — matches L1 C/A's per-chip
    # confidence (3/40 ≈ 7.5 %) within the shorter 10-block window.
    _try_match(code_block_bits, B(0x035), B(0x3ff), 2)
end

function get_default_correlator(gpsl5::GPSL5I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 20-block sync-search window (2 × NH10) needs at least 20 bits.
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
