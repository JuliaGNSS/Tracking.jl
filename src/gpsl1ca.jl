# GPS L1 C/A locks its bit edge with the soft-decision, maximum-energy
# CFAR detector `_detect_bit_edge_cfar` (driven by the `soft_prompts`
# history) rather than a hard-decision template match. The trait below
# routes `_buffer_find_bit` to that path; there is no
# `detect_bit_or_secondary_code_sync(::GPSL1CA, …)` method.
#
# The soft detector is edge-locked by construction: it can only ever fire
# at the energy-maximizing phase's own bit boundary, so it cannot lock one
# block early the way the old fixed-tolerance template matcher did when the
# first data transition was preceded by a repeated bit (issue #124). The
# lock latency self-paces with C/N₀ — ~40 ms for a clean signal, longer in
# noise — via `get_bit_edge_detection_confidence(::GPSL1CA)`.
@inline uses_soft_bit_edge_detection(::GPSL1CA) = true

"""
$(SIGNATURES)

Get the default correlator for the given GNSS system. Returns an
EarlyPromptLateCorrelator for GPS L1 or a VeryEarlyPromptLateCorrelator
for systems like Galileo E1B that use BOC modulation.
"""
function get_default_correlator(gpsl1::GPSL1CA, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 40-block sync-search window (2 × 20 blocks/symbol) needs at least 40 bits.
@inline get_code_block_buffer_type(::GPSL1CA) = UInt64
