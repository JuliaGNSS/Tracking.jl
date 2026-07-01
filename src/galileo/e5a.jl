"""
$(SIGNATURES)

Secondary-code sync detector for Galileo E5a-I — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the 20-chip
CS20 secondary code (Galileo OS SIS ICD Table 19) overlaid on the 1 ms
primary code period. E5a-I carries the F/NAV data stream at 50 sps, so one
CS20 period (20 primary blocks) is exactly one channel symbol: the
detector locks the secondary phase, and data-bit decoding then integrates
one CS20 period per symbol. Default 2.5 % tolerance discretizes to 0
(exact match over the 20-chip window). The packed reference comes from the
generic [`_packed_secondary_code`](@ref). Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE5aI,
    prn::Integer,        # ignored by the reference; CS20 is shared across PRNs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

"""
$(SIGNATURES)

Secondary-code sync detector for Galileo E5a-Q — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the per-PRN
100-chip CS100 secondary code (Galileo OS SIS ICD Table 20) overlaid on
the 1 ms primary code period, giving a 100 ms cycle. E5a-Q is a dataless
pilot; the CS100 overlay is its only sync feature, so the detector locks
after a single CS100 period in the worst case and reports the upcoming
integration's CS100 chip in `SyncResult.phase`. The per-PRN packed
reference comes from the generic [`_packed_secondary_code`](@ref), which
reads the signal's [`PerPRNSecondaryCode`](@ref). Returns
[`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE5aQ,
    prn::Integer,        # selects the PRN's CS100 column in the reference
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both E5a components are BPSK (`LOC`) on the L5 carrier (1176.45 MHz), so
# the C/A-style EarlyPromptLate default applies.
function get_default_correlator(
    ::Union{GalileoE5aI,GalileoE5aQ},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# E5a-I: sync-search window is one CS20 period (20 blocks); UInt32 holds it.
@inline get_code_block_buffer_type(::GalileoE5aI) = UInt32
# E5a-Q: the CS100 overlay search needs a 100-block window; UInt128 is the
# smallest built-in unsigned that holds it (the rotation search masks down
# to the low 100 bits).
@inline get_code_block_buffer_type(::GalileoE5aQ) = UInt128
