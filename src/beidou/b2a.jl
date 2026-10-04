# BeiDou B2a is a data/pilot pair on the 1176.45 MHz carrier — the `L5` band it
# shares with GPS L5 and Galileo E5a. Both components are BPSK(10) 10230-chip
# primary codes at 10.23 Mcps (1 ms period), split 50/50 in power: the in-phase
# B2a data channel carries B-CNAV2 at 200 sym/s under a 5-chip shared secondary
# code, the quadrature B2a pilot is dataless under a per-PRN 100-chip secondary
# code (BDS-SIS-ICD-B2a-1.0 §5.2). Tracking.jl treats them as two independent
# signals, the same split as GPS L5 and Galileo E5a/E5b.

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for the BeiDou B2a data component over
the shared 5-chip code (`00010`; BDS-SIS-ICD-B2a-1.0 §5.2.1) — see
[`_detect_secondary_code_sync`](@ref). One period is exactly one 200 sym/s
B-CNAV2 symbol. The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB2aI,
    prn::Integer,        # ignored by the reference; the 5-chip code is shared
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for the BeiDou B2a pilot over the
per-PRN 100-chip code (BDS-SIS-ICD-B2a-1.0 §5.2.1 Table 5-4; 100 ms cycle) —
see [`_detect_secondary_code_sync`](@ref). The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB2aQ,
    prn::Integer,        # selects the PRN's secondary-code column
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both B2a components are plain BPSK(10) (`LOC`), so the C/A-style
# EarlyPromptLate default applies to each — as for GPS L5 on the same carrier.
function get_default_correlator(
    ::Union{BeiDouB2aI,BeiDouB2aQ},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# Hold one 5-chip / 100-chip period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::BeiDouB2aI) = UInt32
@inline get_code_block_buffer_type(::BeiDouB2aQ) = UInt128
