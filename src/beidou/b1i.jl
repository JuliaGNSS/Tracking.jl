# BeiDou B1I is the legacy open-service signal on the 1561.098 MHz carrier
# (`B1I`): a single BPSK 2046-chip primary code at 2.046 Mcps (1 ms period), no
# quadrature counterpart in the open-service ICD — hence `get_relative_power`
# of 1.0, it is its own composite.
#
# ## The NH20 overlay is per-satellite-orbit, and that shapes the detector
#
# BDS-SIS-ICD-B1I-3.0 §5.2.1 applies the 20-bit Neuman-Hoffman overlay only on
# the MEO/IGSO satellites (PRN 6-58), which broadcast the D1 message at 50 sym/s
# — so one NH20 period is exactly one data symbol, the same "overlay period =
# one symbol" shape as Galileo E5a-I's CS20. The GEO satellites (PRN 1-5 and
# 59-63) broadcast D2 at 500 sym/s and carry *no* overlay; GNSSSignals models
# that uniformly by returning a `PerPRNSecondaryCode` whose GEO columns are
# all-ones, so a GEO PRN's tiered code equals its primary code.
#
# Tracking.jl takes that table as-is and stays on the generic secondary-code
# path for every PRN. On the D1 PRNs that is exactly right. On a GEO PRN the
# sync search simply never completes: an all-ones reference is
# rotation-invariant, so the 20 rotation hypotheses differ only in how their
# 20-block bins straddle the *data* — and a GEO satellite's D2 symbols are 2
# blocks long, so every bin averages ~10 random symbols and no rotation stands
# out. The soft CFAR detector needs a peak that beats its runner-up, so it
# correctly declines to lock rather than picking one at random. A GEO satellite
# therefore tracks and can be ranged on, but stays in the pre-sync one-block
# integration regime and decodes no bits. (The two GEO facts are what rule the
# lock out together: no overlay *and* D2. Were the overlay missing at the D1
# rate, the same bins would reduce to the bit-edge search GPS L1 C/A uses and
# lock on the symbol boundary — `test/beidou_b1i.jl` pins both halves.)
#
# Decoding a GEO satellite therefore needs a per-PRN data rate upstream:
# `get_data_frequency` is a per-signal-type accessor and reports the D1 value
# (50 Hz) for every PRN, so there is nothing that could drive a 2-block bit-edge
# search at the D2 rate. Ranging on one is unaffected.

"""
$(SIGNATURES)

Secondary-code sync detector for BeiDou B1I — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the per-PRN
20-chip Neuman-Hoffman overlay (NH20; BDS-SIS-ICD-B1I-3.0 §5.2.1)
overlaid on the 1 ms primary code period, giving a 20 ms tiered code. On
the MEO/IGSO satellites (PRN 6-58) that carry NH20, one overlay period is
exactly one D1 data symbol at 50 sym/s: the detector locks the secondary
phase, and data-bit decoding then integrates one NH20 period per symbol.

On the GEO satellites (PRN 1-5, 59-63) the overlay column is all-ones —
they carry no NH20, and they are also the D2 satellites — so there is
nothing for the rotation search to lock and the soft CFAR detector does not
sync at all (see the file header). Those satellites track and range, but
stay pre-sync.

At runtime this method is reached only if a caller forces B1I onto the hard
path; the live detector is the soft
[`_detect_secondary_code_cfar`](@ref), which
[`uses_soft_secondary_code_detection`](@ref) selects for a 20-chip overlay.
Both read the same per-PRN reference — this one via the generic
[`_packed_secondary_code`](@ref), which reads the signal's
[`PerPRNSecondaryCode`](@ref). Returns [`SyncResult`](@ref).
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

# Sync-search window is one NH20 period (20 blocks); UInt32 holds it (the
# rotation search masks down to the low 20 bits) — as for Galileo E5a-I's CS20.
@inline get_code_block_buffer_type(::BeiDouB1I) = UInt32
