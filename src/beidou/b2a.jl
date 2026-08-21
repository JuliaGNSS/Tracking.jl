# BeiDou B2a is a data/pilot pair on the 1176.45 MHz carrier — the `L5` band it
# shares with GPS L5 and Galileo E5a. Both components are BPSK(10) 10230-chip
# primary codes at 10.23 Mcps (1 ms period), split 50/50 in power: the in-phase
# B2a data channel carries B-CNAV2 at 200 sym/s under a 5-chip shared secondary
# code, the quadrature B2a pilot is dataless under a per-PRN 100-chip secondary
# code (BDS-SIS-ICD-B2a-1.0 §5.2). Tracking.jl treats them as two independent
# signals, the same split as GPS L5 and Galileo E5a/E5b.

"""
$(SIGNATURES)

Secondary-code sync detector for the BeiDou B2a data component — the
generic [`_detect_secondary_code_sync`](@ref) rotation search over the
5-chip secondary code (`00010`, shared across all PRNs;
BDS-SIS-ICD-B2a-1.0 §5.2.1) overlaid on the 1 ms primary code period. B2a
carries B-CNAV2 at 200 sym/s, so one secondary period (5 primary blocks)
is exactly one channel symbol: the detector locks the secondary phase, and
data-bit decoding then integrates one secondary period per symbol — the
same "overlay period = one symbol" shape as Galileo E5a-I's CS20 and
E5b-I's CS4.

With `N = 5` the trait default routes B2a data to the soft, CFAR detector
([`uses_soft_secondary_code_detection`](@ref)); a 5-chip hard template
match would be badly false-lock-prone. The packed reference comes from the
generic [`_packed_secondary_code`](@ref). Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB2aI,
    prn::Integer,        # ignored by the reference; the 5-chip code is shared
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

"""
$(SIGNATURES)

Secondary-code sync detector for the BeiDou B2a pilot component — the
generic [`_detect_secondary_code_sync`](@ref) rotation search over the
per-PRN 100-chip secondary code (truncated length-1021 Weil codes,
BDS-SIS-ICD-B2a-1.0 §5.2.1 Table 5-4) overlaid on the 1 ms primary code
period, giving a 100 ms cycle. The B2a pilot is dataless; the overlay is
its only sync feature, so the detector locks after a single overlay period
in the worst case and reports the upcoming integration's secondary chip in
`SyncResult.phase`. The per-PRN packed reference comes from the generic
[`_packed_secondary_code`](@ref), which reads the signal's
[`PerPRNSecondaryCode`](@ref) — the same shape as Galileo E5a-Q / E5b-Q /
E6-C. With `N = 100` the trait default routes the pilot to the soft
[`_detect_secondary_code_cfar`](@ref)
([`uses_soft_secondary_code_detection`](@ref)), so this method is reached
only if a caller forces the hard path. Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB2aQ,
    prn::Integer,        # selects the PRN's secondary-code column
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both B2a components are plain BPSK(10) (`LOC`), so the C/A-style
# EarlyPromptLate default applies to each — as for GPS L5 on the same carrier.
function get_default_correlator(
    ::Union{BeiDouB2aI,BeiDouB2aQ},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# B2a data: the 5-chip secondary search window is 5 blocks. UInt32 matches the
# other short-secondary-code signals (GPS L5I/L5Q, Galileo E1C/E5a-I/E5b-I) —
# wider than the horizon needs, but it keeps the hard rotation sweep available
# and the widths uniform across that family.
@inline get_code_block_buffer_type(::BeiDouB2aI) = UInt32
# B2a pilot: the 100-chip overlay search needs a 100-block window; UInt128 is
# the smallest built-in unsigned that holds it (the rotation search masks down
# to the low 100 bits).
@inline get_code_block_buffer_type(::BeiDouB2aQ) = UInt128
