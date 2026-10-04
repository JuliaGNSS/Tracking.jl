# Galileo E5b is the upper sideband of the E5 AltBOC(15,10) signal on the
# 1207.14 MHz carrier (`E5b`, shared with BeiDou B2b), processed "as though"
# E5a and E5b "were two separate QPSK signals" (OS SIS ICD §2.3.1.2): E5b-I
# carries I/NAV at 250 sym/s under a CS4 overlay, E5b-Q is a dataless pilot
# under a per-SVID CS100 overlay. Tracked as two independent signals.

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for Galileo E5b-I over the shared CS4
code (`1110`; OS SIS ICD v2.2 §3.5.1) — see [`_detect_secondary_code_sync`](@ref).
One CS4 period is exactly one 250 sym/s I/NAV symbol. The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE5bI,
    prn::Integer,        # ignored by the reference; CS4 is shared across SVIDs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for the Galileo E5b-Q pilot over the
per-SVID CS100 code (OS SIS ICD v2.2 §3.5.2: CS100₍ₙ₊₅₀₎ for SVID `n`; 100 ms
cycle) — see [`_detect_secondary_code_sync`](@ref). The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE5bQ,
    prn::Integer,        # selects the SVID's CS100 column in the reference
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Modelled as two independent BPSK(10) (`LOC`) sidebands rather than one
# AltBOC(15,10), so neither component sees a split-spectrum autocorrelation and
# the C/A-style EarlyPromptLate default applies to each — same as E5a.
function get_default_correlator(
    ::Union{GalileoE5bI,GalileoE5bQ},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# Hold one CS4 / CS100 period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::GalileoE5bI) = UInt32
@inline get_code_block_buffer_type(::GalileoE5bQ) = UInt128
