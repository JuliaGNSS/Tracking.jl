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
@inline function detect_bit_or_secondary_code_sync(
    ::GPSL1C_D,
    ::Integer,           # PRN — ignored
    ::Unsigned,
    ::Integer,
)
    SyncResult(true, 0, Int8(+1))
end

# GPS L1C-D is BOC(1,1): its autocorrelation has nulls at ±0.25-0.3 chip and
# side-lobes at ±0.4-0.6 chip. The C/A-style 0.5-chip early-late spacing would
# place the early/late taps on those side-lobes, biasing the DLL discriminator
# (it reports a non-zero error even when perfectly aligned) and walking the code
# replica off the peak. A narrow 0.1-chip spacing keeps both taps on the BOC main
# peak where the discriminator is unbiased.
function get_default_correlator(gpsl1c_d::GPSL1C_D, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants, preferred_early_late_to_prompt_code_shift = 0.1)
end

# 1 channel symbol = 1 primary code period (100 sps, 10 ms primary
# period). No sub-symbol boundary to search for — the sync-search buffer
# is dead state at runtime; we still pick a concrete type to keep the
# `BitBuffer{B}` parameter chain stable. UInt8 is the smallest legal
# Unsigned.
@inline get_code_block_buffer_type(::GPSL1C_D) = UInt8
