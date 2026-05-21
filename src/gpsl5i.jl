"""
$(SIGNATURES)

Secondary-code sync detector for GPS L5I.

Matches the low 10 bits of `code_block_bits` against the NH10 template
`0x035` (= `0000110101`). The negated-polarity template `0x3ca` is
searched in the same call via [`_try_match`](@ref). Returns
[`SyncResult`](@ref).
"""
@inline function is_upcoming_integration_new_bit(
    signal::GPSL5I,
    ::Integer,           # PRN — ignored; NH10 is shared across PRNs
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    num_code_blocks < 10 && return SyncResult(false, 0, Int8(0))
    # Tolerance is a percentage of the 10-block window. The default
    # 2.5 % discretizes to floor(0.025 × 10) = 0 (exact match).
    # Adjustable via `get_bit_edge_or_secondary_code_tolerance(::GPSL5I)`.
    max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * 10)
    _try_match(code_block_bits, B(0x035), B(0x3ff), max_errors)
end

function get_default_correlator(gpsl5::GPSL5I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 20-block sync-search window (2 × NH10) needs at least 20 bits.
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
