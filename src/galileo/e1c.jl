const GalileoE1CAny = Union{GalileoE1C,GalileoE1C_BOC11}

"""
$(SIGNATURES)

Secondary-code sync detector for Galileo E1C (both [`GalileoE1C`](@ref)
and the BOC(1,1) approximation [`GalileoE1C_BOC11`](@ref)).

E1C is the E1 pilot channel: no navigation data, but a 25-chip CS25
secondary code (Galileo OS SIS ICD Table 4) overlaid on the 4 ms primary
code period, giving a 100 ms cycle. The generic
[`_detect_secondary_code_sync`](@ref) rotation search locks after a single
CS25 period in the worst case (default 2.5 % tolerance discretizes to 0,
i.e. exact match over the 25-chip window) and reports the upcoming
integration's CS25 chip in `SyncResult.phase`. The packed reference comes
from the generic [`_packed_secondary_code`](@ref). Returns
[`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE1CAny,
    prn::Integer,        # ignored by the reference; CS25 is shared across PRNs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# E1C shares the E1 modulation family with E1B (CBOC / BOC(1,1)), so it
# uses the same VeryEarlyPromptLate default.
function get_default_correlator(::GalileoE1CAny, num_ants::NumAnts = NumAnts(1))
    VeryEarlyPromptLateCorrelator(; num_ants)
end

# Sync-search window is one CS25 period (25 blocks); UInt32 holds it (the
# rotation search masks down to the low 25 bits).
@inline get_code_block_buffer_type(::GalileoE1CAny) = UInt32
