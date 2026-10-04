const GalileoE1CAny = Union{GalileoE1C,GalileoE1C_BOC11}

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for the Galileo E1C pilot (both
[`GalileoE1C`](@ref) and the BOC(1,1) approximation [`GalileoE1C_BOC11`](@ref))
over the shared CS25 code (Galileo OS SIS ICD Table 4; 100 ms cycle) — see
[`_detect_secondary_code_sync`](@ref). The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
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

# Holds one CS25 period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::GalileoE1CAny) = UInt32
