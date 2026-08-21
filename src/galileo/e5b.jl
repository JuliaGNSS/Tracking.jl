# Galileo E5b is the upper sideband of the wideband E5 AltBOC(15,10) signal, on
# the 1207.14 MHz carrier (`E5b`, shared exactly with BeiDou B2b). The OS SIS
# ICD blesses processing E5a and E5b "as though they were two separate QPSK
# signals" (§2.3.1.2), which is what GNSSSignals models and what Tracking.jl
# tracks: the in-phase E5b-I carries I/NAV at 250 sym/s under a 4-chip CS4
# overlay, the quadrature E5b-Q is a dataless pilot under a per-SVID 100-chip
# CS100 overlay. Both share the 1 ms / 10230-chip primary code (10.23 Mcps) and
# are handled as two independent signals, the same split as E5a and GPS L5.

"""
$(SIGNATURES)

Secondary-code sync detector for Galileo E5b-I — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the 4-chip CS4
secondary code (`1110`, shared across all SVIDs; Galileo OS SIS ICD v2.2
§3.5.1) overlaid on the 1 ms primary code period. E5b-I carries the I/NAV
data stream at 250 sym/s, so one CS4 period (4 primary blocks) is exactly
one channel symbol: the detector locks the secondary phase, and data-bit
decoding then integrates one CS4 period per symbol. This is the same
"overlay period = one symbol" shape as Galileo E5a-I's CS20 at 50 sym/s,
scaled to E5b-I's five-times-faster symbol rate.

At the default 2.5 % tolerance the hard-path error budget discretizes to 0
(exact match over the 4-chip window), but with `N = 4` the trait default
routes E5b-I to the soft, CFAR detector
([`uses_soft_secondary_code_detection`](@ref)) instead — a 4-chip hard
template match would be badly false-lock-prone. The packed reference comes
from the generic [`_packed_secondary_code`](@ref). Returns
[`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE5bI,
    prn::Integer,        # ignored by the reference; CS4 is shared across SVIDs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

"""
$(SIGNATURES)

Secondary-code sync detector for Galileo E5b-Q — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the per-SVID
100-chip CS100 secondary code (Galileo OS SIS ICD v2.2 §3.5.2, the CS100
codes 51-100 assigned CS100₍ₙ₊₅₀₎ to SVID `n`) overlaid on the 1 ms
primary code period, giving a 100 ms cycle. E5b-Q is a dataless pilot; the
CS100 overlay is its only sync feature, so the detector locks after a
single CS100 period in the worst case and reports the upcoming
integration's CS100 chip in `SyncResult.phase`. The per-PRN packed
reference comes from the generic [`_packed_secondary_code`](@ref), which
reads the signal's [`PerPRNSecondaryCode`](@ref) — the same shape as
Galileo E5a-Q, which draws the other half (CS100₁₋₅₀) of the same table.
With `N = 100` the trait default routes E5b-Q to the soft
[`_detect_secondary_code_cfar`](@ref)
([`uses_soft_secondary_code_detection`](@ref)), so this method is reached
only if a caller forces the hard path. Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE5bQ,
    prn::Integer,        # selects the SVID's CS100 column in the reference
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Modelled as two independent BPSK(10) (`LOC`) sidebands rather than one
# AltBOC(15,10), so neither component sees a split-spectrum autocorrelation and
# the C/A-style EarlyPromptLate default applies to each — same as E5a.
function get_default_correlator(
    ::Union{GalileoE5bI,GalileoE5bQ},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# E5b-I: the CS4 search window is 4 blocks. UInt32 matches the other
# short-secondary-code signals (GPS L5I/L5Q, Galileo E1C/E5a-I) — far wider
# than the horizon needs, but it keeps the hard rotation sweep available and
# the widths uniform across that family.
@inline get_code_block_buffer_type(::GalileoE5bI) = UInt32
# E5b-Q: the CS100 overlay search needs a 100-block window; UInt128 is the
# smallest built-in unsigned that holds it (the rotation search masks down
# to the low 100 bits) — as for Galileo E5a-Q.
@inline get_code_block_buffer_type(::GalileoE5bQ) = UInt128
