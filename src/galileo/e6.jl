# Galileo E6 is the civil pair on the 1278.75 MHz carrier (`E6`): E6-B carries
# the C/NAV message — and with it the High Accuracy Service — at 1000 sym/s,
# E6-C is its dataless pilot under a per-SVID 100-chip CS100 overlay. Both share
# the 1 ms / 5115-chip primary code (5.115 Mcps) and are tracked as two
# independent BPSK(5) signals, the same data/pilot split as Galileo E1 and E5.
#
# The two differ in *carrier* phase, not only in code: OS SIS ICD Eq. 10 forms
# the composite as `(e_E6-B − e_E6-C)/√2`, so E6-C sits at π against the E6-B
# reference (GNSSSignals reports that via `get_carrier_phase_offset`). Tracking
# reads that offset generically when several signals on one carrier share a
# driver loop (`conventional_pll_and_dll.jl`), so nothing per-signal is needed
# here for it.

"""
$(SIGNATURES)

Symbol-sync detector for Galileo E6-B.

E6-B broadcasts one C/NAV channel symbol per 1 ms primary code period
(1000 sym/s; Galileo OS SIS ICD v2.2 Table 5 — 5115 chips at 5.115 Mcps
with no secondary code, so the symbol period *is* the primary code
period), so the buffer of primary-block signs is itself the symbol stream
— there is no sub-symbol boundary to find. The detector therefore reports
`found = true` from the very first integration, leaving downstream
consumers (GNSSDecoder.jl) to resolve the residual ±1 polarity ambiguity
via the C/NAV preamble.

Same shape as the Galileo E1B and GPS L1C-D "1-block-per-symbol" cases,
at E6-B's four-times-faster symbol rate.
"""
@inline detect_bit_or_secondary_code_sync(
    signal::GalileoE6B,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

"""
$(SIGNATURES)

Secondary-code sync detector for Galileo E6-C — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the per-SVID
100-chip CS100 secondary code (Galileo E6-B/C Codes Technical Note §2.4,
the OS SIS ICD's CS100₁₋₅₀ assigned CS100ₙ to SVID `n`) overlaid on the
1 ms primary code period, giving a 100 ms cycle. E6-C is a dataless pilot;
the CS100 overlay is its only sync feature, so the detector locks after a
single CS100 period in the worst case and reports the upcoming
integration's CS100 chip in `SyncResult.phase`. The per-PRN packed
reference comes from the generic [`_packed_secondary_code`](@ref), which
reads the signal's [`PerPRNSecondaryCode`](@ref) — E6-C draws the same
CS100₁₋₅₀ half of the table as [`GalileoE5aQ`](@ref). With `N = 100` the
trait default routes E6-C to the soft
[`_detect_secondary_code_cfar`](@ref)
([`uses_soft_secondary_code_detection`](@ref)), so this method is reached
only if a caller forces the hard path. Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE6C,
    prn::Integer,        # selects the SVID's CS100 column in the reference
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both E6 components are plain BPSK(5) (`LOC`) — no subcarrier, so no
# split-spectrum side peaks — and the C/A-style EarlyPromptLate default applies.
function get_default_correlator(
    ::Union{GalileoE6B,GalileoE6C},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# 1 channel symbol = 1 primary code period (1000 sym/s, 1 ms primary
# period). No sub-symbol boundary to search for — the sync-search buffer
# is dead state at runtime; we still pick a concrete type to keep the
# `BitBuffer{B}` parameter chain stable. UInt8 is the smallest legal
# Unsigned.
@inline get_code_block_buffer_type(::GalileoE6B) = UInt8
# E6-C: the CS100 overlay search needs a 100-block window; UInt128 is the
# smallest built-in unsigned that holds it (the rotation search masks down
# to the low 100 bits) — as for Galileo E5a-Q / E5b-Q.
@inline get_code_block_buffer_type(::GalileoE6C) = UInt128
