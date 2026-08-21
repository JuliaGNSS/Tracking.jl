# BeiDou B3I is the legacy open-service signal on the 1268.52 MHz carrier
# (`B3I`): a single BPSK 10230-chip primary code at 10.23 Mcps (1 ms period),
# again with no open-service quadrature counterpart (`get_relative_power` 1.0).
#
# Structurally B3I is B1I at five times the chipping rate: the same 20-chip
# Neuman-Hoffman overlay on the MEO/IGSO satellites (PRN 6-58, D1 at 50 sym/s),
# the same absent overlay on the GEO ones (PRN 1-5, 59-63, D2)
# — BDS-SIS-ICD-B3I-1.0 §5.2.1. Everything the `b1i.jl` header says applies here
# verbatim: the all-ones GEO column plus the D2 symbol rate leave the sync
# search nothing to lock, so a GEO satellite ranges but stays pre-sync, and
# `get_data_frequency` reports the D1 rate for every PRN.

"""
$(SIGNATURES)

Secondary-code sync detector for BeiDou B3I — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the per-PRN
20-chip Neuman-Hoffman overlay (NH20; BDS-SIS-ICD-B3I-1.0 §5.2.1)
overlaid on the 1 ms primary code period, giving a 20 ms tiered code. The
detector behaves exactly as [`BeiDouB1I`](@ref)'s — one NH20 period is one
D1 symbol on the MEO/IGSO satellites (PRN 6-58), while the GEO satellites'
all-ones column, at their 2-block D2 symbol rate, leaves nothing for any
rotation to lock so they never sync. As for B1I, the live detector is the
soft [`_detect_secondary_code_cfar`](@ref) and this method is reached only
if a caller forces B3I onto the hard path; the per-PRN packed reference
comes from the generic [`_packed_secondary_code`](@ref) either way. Returns
[`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB3I,
    prn::Integer,        # selects the PRN's NH20 column (all-ones on GEO)
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# B3I is plain BPSK (`LOC`), so the C/A-style EarlyPromptLate default applies.
function get_default_correlator(::BeiDouB3I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Sync-search window is one NH20 period (20 blocks); UInt32 holds it.
@inline get_code_block_buffer_type(::BeiDouB3I) = UInt32
