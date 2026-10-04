# GPS L5 is a quadrature pair on the L5 carrier (1176.45 MHz): the in-phase
# L5I channel carries the CNAV data stream under a 10-chip Neuman-Hoffman
# (NH10) overlay, the quadrature L5Q channel is a dataless pilot under a
# 20-chip Neuman-Hoffman (NH20) overlay. Both share the 1 ms / 10230-chip
# primary code (10.23 Mcps); Tracking.jl treats them as two independent
# signals (`GPSL5I`, `GPSL5Q`), same as the L1C data/pilot split.

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for GPS L5I over the shared NH10 code —
see [`_detect_secondary_code_sync`](@ref). The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GPSL5I,
    prn::Integer,        # ignored by the reference; NH10 is shared across PRNs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for the GPS L5Q pilot over the shared
NH20 code — see [`_detect_secondary_code_sync`](@ref). The default detector is
the soft one ([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GPSL5Q,
    prn::Integer,        # ignored by the reference; NH20 is shared across PRNs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both L5 components are BPSK (`LOC`) on the L5 carrier, so the C/A-style
# EarlyPromptLate default applies to each.
function get_default_correlator(::Union{GPSL5I,GPSL5Q}, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Holds one NH10 / NH20 period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
@inline get_code_block_buffer_type(::GPSL5Q) = UInt32
