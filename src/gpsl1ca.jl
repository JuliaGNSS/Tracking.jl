"""
$(SIGNATURES)

Bit-sync detector for GPS L1 C/A.

Matches the low 40 bits of `code_block_bits` against the bit-edge
template `0xfffff` (20 ones followed by 20 zeros). The negated-polarity
template `0xfffff00000` is searched in the same call via
[`_try_match`](@ref). Returns [`SyncResult`](@ref).

The two blocks straddling the claimed bit edge (window bits 19 and 20)
must match exactly; the Hamming tolerance applies only to the other 38
blocks. Otherwise a single tolerated bit-flip adjacent to the edge would
let the window one block before the true edge fire first whenever the
first transition is preceded by a repeated bit, permanently misaligning
the bit grid by 1 ms (issue #124).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GPSL1CA,
    ::Integer,           # PRN — ignored; same template for every PRN
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    num_code_blocks < 40 && return SyncResult(false, 0, Int8(0))
    # Tolerance is a percentage of the 40-block window. The default
    # 2.5 % discretizes to floor(0.025 × 40) = 1 bit-flip allowed.
    # Adjustable via `get_bit_edge_or_secondary_code_tolerance(::GPSL1CA)`.
    max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * 40)
    _try_match(code_block_bits, B(0xfffff), B(0xffffffffff), max_errors, B(0x180000))
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
