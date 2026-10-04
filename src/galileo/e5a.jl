"""
$(SIGNATURES)

Hard-path secondary-code sync detector for Galileo E5a-I over the shared CS20
code (Galileo OS SIS ICD Table 19) — see [`_detect_secondary_code_sync`](@ref).
One CS20 period is exactly one 50 sps F/NAV symbol. The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
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

Hard-path secondary-code sync detector for the Galileo E5a-Q pilot over the
per-PRN CS100 code (Galileo OS SIS ICD Table 20; 100 ms cycle) — see
[`_detect_secondary_code_sync`](@ref). The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
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

# Hold one CS20 / CS100 period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::GalileoE5aI) = UInt32
@inline get_code_block_buffer_type(::GalileoE5aQ) = UInt128
