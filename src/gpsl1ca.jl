"""
$(SIGNATURES)

Bit-sync detector for GPS L1 C/A.

Matches the low 40 bits of `code_block_bits` against the bit-edge
template `0xfffff` (20 ones followed by 20 zeros). The negated-polarity
template `0xfffff00000` is searched in the same call via
[`_try_match`](@ref). Returns [`SyncResult`](@ref).
"""
@inline function is_upcoming_integration_new_bit(
    ::GPSL1CA,
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    num_code_blocks < 40 && return SyncResult(false, 0, Int8(0))
    _try_match(code_block_bits, B(0xfffff), B(0xffffffffff), 0)
end

"""
$(SIGNATURES)

Get the default correlator for the given GNSS system. Returns an
EarlyPromptLateCorrelator for GPS L1 or a VeryEarlyPromptLateCorrelator
for systems like Galileo E1B that use BOC modulation.
"""
function get_default_correlator(gpsl1::GPSL1CA, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 40-block sync-search window (2 × 20 blocks/symbol) needs at least 40 bits.
@inline get_code_block_buffer_type(::GPSL1CA) = UInt64
