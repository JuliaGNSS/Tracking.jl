# BeiDou B2b_I is the open-service data signal on the 1207.14 MHz carrier — the
# `E5b` band, which Galileo E5b shares exactly. A single BPSK(10) 10230-chip
# primary code at 10.23 Mcps (1 ms period), no secondary code, and the B-CNAV3
# message at 1000 sym/s, so the primary period already *is* the symbol period
# (BDS-SIS-ICD-B2b-1.0 §5). No open-service quadrature counterpart, hence
# `get_relative_power` 1.0.
#
# Note the ICD defines ranging codes for PRN 6-58 only; GNSSSignals returns an
# all-zero code for any other PRN index, which correlates to nothing. That is an
# upstream property of the code table, not something the tracker screens for.

"""
$(SIGNATURES)

Symbol-sync detector for BeiDou B2b_I.

B2b_I broadcasts one B-CNAV3 symbol per 1 ms primary code period
(1000 sym/s; BDS-SIS-ICD-B2b-1.0 §5 — 10230 chips at 10.23 Mcps with no
secondary code, so the symbol period *is* the primary code period), so the
buffer of primary-block signs is itself the symbol stream — there is no
sub-symbol boundary to find. The detector therefore reports `found = true`
from the very first integration, leaving downstream consumers
(GNSSDecoder.jl) to resolve the residual ±1 polarity ambiguity via the
B-CNAV3 preamble.

Same shape as the Galileo E1B / E6-B and GPS L1C-D "1-block-per-symbol"
cases.
"""
@inline detect_bit_or_secondary_code_sync(
    signal::BeiDouB2bI,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

# B2b_I is plain BPSK(10) (`LOC`), so the C/A-style EarlyPromptLate default
# applies.
function get_default_correlator(::BeiDouB2bI, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 1 symbol = 1 primary code period (1000 sym/s, 1 ms primary period). No
# sub-symbol boundary to search for — the sync-search buffer is dead state at
# runtime; we still pick a concrete type to keep the `BitBuffer{B}` parameter
# chain stable. UInt8 is the smallest legal Unsigned.
@inline get_code_block_buffer_type(::BeiDouB2bI) = UInt8
