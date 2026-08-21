# BeiDou B1C is the modern data/pilot pair on the L1 carrier (1575.42 MHz),
# alongside GPS L1 C/A, GPS L1C and Galileo E1. Both components ride a
# 10230-chip primary code at 1.023 Mcps (10 ms period): B1C data carries B-CNAV1
# at 100 sym/s with no secondary code, B1C pilot is dataless under a per-PRN
# 1800-chip overlay (18 s cycle). The power split is 1:3 data:pilot
# (BDS-SIS-ICD-B1C-1.0 §4.2.2), which is what `get_relative_power` reports.
#
# Modulation: GNSSSignals generates both components as sine-phased BOC(1,1).
# For the data component that is the ICD's own modulation. The pilot is
# specified as QMBOC(6,1,4/33) — a BOC(1,1) arm and a BOC(6,1) arm in phase
# *quadrature* — whose principal, 29/33-power axis simply *is* BOC(1,1), so a
# real correlator on the pilot's own carrier phase uses a BOC(1,1) replica and
# cannot capture the orthogonal BOC(6,1) arm at all. Either way both components
# are BOC(1,1)-class here, which is what sets the correlator default below.

"""
$(SIGNATURES)

Symbol-sync detector for the BeiDou B1C data component.

B1C data broadcasts one B-CNAV1 symbol per 10 ms primary code period
(100 sym/s; BDS-SIS-ICD-B1C-1.0 Table 5-1 — 10230 chips at 1.023 Mcps with
no secondary code, so the symbol period *is* the primary code period), so
the buffer of primary-block signs is itself the symbol stream — there is
no sub-symbol boundary to find. The detector therefore reports
`found = true` from the very first integration, leaving downstream
consumers (GNSSDecoder.jl) to resolve the residual ±1 polarity ambiguity
via the B-CNAV1 preamble.

Exactly the GPS L1C-D case, at the same 10 ms period and 100 sym/s rate.
"""
@inline detect_bit_or_secondary_code_sync(
    signal::BeiDouB1C_D,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

"""
$(SIGNATURES)

Secondary-code sync detector for the BeiDou B1C pilot component.

B1C pilot broadcasts a per-PRN 1800-chip overlay code
(BDS-SIS-ICD-B1C-1.0 §5.2.2, truncated Weil codes) on top of the 10 ms
primary code, giving an 18 second cycle — dimensionally the same overlay as
GPS L1C-P's, on the same carrier. The generic
[`_detect_secondary_code_sync`](@ref) waits for the sliding
`code_block_bits` window to fill to 1800 primary periods, then runs a
single 1800-phase shifted Hamming-distance sweep
([`_secondary_code_search`](@ref)) against the PRN's overlay pattern,
accepting the best alignment within the 2.5 % tolerance
([`get_bit_edge_or_secondary_code_tolerance`](@ref), 45 errors here) and
reporting the secondary-chip offset of the *upcoming* integration in
`SyncResult.phase`.

Returns `SyncResult(false, 0, 0)` until 1800 blocks have been buffered.
Like GPS L1C-P's, the overlay reaches the sweep through the generic
[`_packed_secondary_code`](@ref): GNSSSignals exposes it as a
[`PerPRNSecondaryCode`](@ref), and the per-chip packing cost is negligible
next to the 1800-phase sweep that follows.
"""
function detect_bit_or_secondary_code_sync(
    signal::BeiDouB1C_P,
    prn::Integer,        # selects the PRN's 1800-chip overlay column
    code_block_bits::UInt1800,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both B1C components are BOC(1,1)-class (see the file header on the pilot's
# QMBOC). The split-spectrum BOC autocorrelation has side peaks at ±0.5 chip
# that a plain early-late discriminator can false-lock onto, so B1C takes the
# same VeryEarlyPromptLate default as the other BOC(1,1)-class L1 signals —
# GPS L1C-D/P and Galileo E1B/E1C — whose very-early/very-late taps feed the
# VEML discriminator that mitigates those side peaks.
function get_default_correlator(
    ::Union{BeiDouB1C_D,BeiDouB1C_P},
    num_ants::NumAnts = NumAnts(1),
)
    VeryEarlyPromptLateCorrelator(; num_ants)
end

# 1 symbol = 1 primary code period (100 sym/s, 10 ms primary period). No
# sub-symbol boundary to search for — the sync-search buffer is dead state at
# runtime; we still pick a concrete type to keep the `BitBuffer{B}` parameter
# chain stable. UInt8 is the smallest legal Unsigned.
@inline get_code_block_buffer_type(::BeiDouB1C_D) = UInt8

# The 1800-chip overlay search needs an exact-width 1800-bit container. The
# `UInt1800` alias is defined in the top-level Tracking module via
# `BitIntegers.@define_integers 1800` for GPS L1C-P's identically sized overlay,
# and is what `BitBuffer{B}` carries for B1C pilot throughout the tracker.
@inline get_code_block_buffer_type(::BeiDouB1C_P) = UInt1800

# B1C pilot keeps the hard-decision rotation sweep, for the same reason GPS
# L1C-P does: its 1800-chip / 18 s overlay is far too long to integrate
# coherently per bin (the soft CFAR detector's model), and a 1800-chip code at
# the 45-error budget is not false-lock-prone anyway. The
# `uses_soft_secondary_code_detection` default already excludes it (N = 1800
# exceeds the 100-chip cap); this makes the intent explicit and keeps B1C pilot
# on the hard path even if that default cap were ever widened.
@inline uses_soft_secondary_code_detection(::BeiDouB1C_P) = false
