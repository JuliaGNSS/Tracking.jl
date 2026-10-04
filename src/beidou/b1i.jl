# BeiDou B1I is the legacy open-service signal on the 1561.098 MHz carrier
# (`B1I`): a single BPSK 2046-chip primary code at 2.046 Mcps (1 ms period), no
# quadrature counterpart in the open-service ICD — hence `get_relative_power`
# of 1.0, it is its own composite.
#
# The NH20 overlay (BDS-SIS-ICD-B1I-3.0 §5.2.1) is applied only on the MEO/IGSO
# satellites (PRN 6-58, D1 at 50 sym/s: one NH20 period = one symbol). The GEO
# satellites (PRN 1-5, 59-63) carry none and broadcast D2 at 500 sym/s;
# GNSSSignals gives them all-ones columns. An all-ones reference is
# rotation-invariant and the 2-block D2 symbols leave no rotation bin standing
# out, so a GEO satellite never syncs: it tracks and ranges but decodes no bits
# (`test/beidou_b1i.jl`). Decoding it would need a per-PRN data rate, but
# `get_data_frequency` is per signal type. See docs/src/bit_sync.md.

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for BeiDou B1I over the per-PRN NH20
overlay (20 ms tiered code) — see [`_detect_secondary_code_sync`](@ref). GEO
PRNs never sync (see the file header). The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB1I,
    prn::Integer,        # selects the PRN's NH20 column (all-ones on GEO)
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# B1I is plain BPSK (`LOC`) — no subcarrier, so no split-spectrum side peaks —
# and the C/A-style EarlyPromptLate default applies.
function get_default_correlator(::BeiDouB1I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Holds one NH20 period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::BeiDouB1I) = UInt32
