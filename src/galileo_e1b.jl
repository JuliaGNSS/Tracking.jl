"""
$(SIGNATURES)

Bit-sync detector for Galileo E1B.

Matches the 8-bit `code_block_bits` window against the bit-edge template
`0x0f` (4 ones followed by 4 zeros). The negated-polarity template
`0xf0` is searched in the same call via [`_try_match`](@ref). Returns
[`SyncResult`](@ref).
"""
@inline function is_upcoming_integration_new_bit(
    ::GalileoE1B,
    ::Integer,           # PRN — ignored; same template for every PRN
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    num_code_blocks < 8 && return SyncResult(false, 0, Int8(0))
    _try_match(code_block_bits, B(0x0f), B(0xff), 0)
end

# TODO: Very early very late correlator?
function get_default_correlator(galileo_e1b::GalileoE1B, num_ants::NumAnts = NumAnts(1))
    VeryEarlyPromptLateCorrelator(; num_ants)
end

# 8-block sync-search window (2 × 4 blocks/symbol) needs at least 8 bits.
@inline get_code_block_buffer_type(::GalileoE1B) = UInt8
