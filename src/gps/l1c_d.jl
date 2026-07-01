"""
$(SIGNATURES)

Symbol-sync detector for GPS L1C-D.

L1C-D broadcasts one CNAV-2 channel symbol per 10 ms primary code
period (100 sps; IS-GPS-800G §3.2.3), so the buffer of primary-block
signs is itself the symbol stream — there is no sub-symbol boundary to
find. The detector therefore reports `found = true` from the very first
integration, leaving downstream consumers (GNSSDecoder.jl) to resolve
the residual ±1 polarity ambiguity via the CNAV-2 preamble's CRC.

Same shape as the Galileo E1B "1-block-per-symbol" case.
"""
@inline detect_bit_or_secondary_code_sync(
    signal::GPSL1C_D,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

# GPS L1C (both components) is BOC(1,1)-class — L1C-D is BOC(1,1), L1C-P is
# TMBOC(6,1,4/33), i.e. 29/33 BOC(1,1) symbols with 4/33 BOC(6,1) sharpening
# the peak slightly. The split-spectrum BOC autocorrelation has side peaks at
# ±0.5 chip that a plain early-late discriminator can false-lock onto. As for
# the equivalent Galileo E1 signals (E1B/E1C), the VeryEarlyPromptLate
# correlator's very-early/very-late taps feed the VEML discriminator that
# mitigates those side peaks, so L1C uses the same default.
function get_default_correlator(::Union{GPSL1C_D,GPSL1C_P}, num_ants::NumAnts = NumAnts(1))
    VeryEarlyPromptLateCorrelator(; num_ants)
end

# 1 channel symbol = 1 primary code period (100 sps, 10 ms primary
# period). No sub-symbol boundary to search for — the sync-search buffer
# is dead state at runtime; we still pick a concrete type to keep the
# `BitBuffer{B}` parameter chain stable. UInt8 is the smallest legal
# Unsigned.
@inline get_code_block_buffer_type(::GPSL1C_D) = UInt8
