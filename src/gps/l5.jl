# GPS L5 is a quadrature pair on the L5 carrier (1176.45 MHz): the in-phase
# L5I channel carries the CNAV data stream under a 10-chip Neuman-Hoffman
# (NH10) overlay, the quadrature L5Q channel is a dataless pilot under a
# 20-chip Neuman-Hoffman (NH20) overlay. Both share the 1 ms / 10230-chip
# primary code (10.23 Mcps); Tracking.jl treats them as two independent
# signals (`GPSL5I`, `GPSL5Q`), same as the L1C data/pilot split.

"""
$(SIGNATURES)

Secondary-code sync detector for GPS L5I — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the NH10
secondary code (window 10; the default 2.5 % tolerance discretizes to 0,
i.e. exact match). The negated-polarity (data-bit-1) case is handled
inside the search, so the detector locks after a single NH10 period in the
worst case and reports the upcoming integration's NH10 chip in
`SyncResult.phase`. Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GPSL5I,
    prn::Integer,        # ignored by the reference; NH10 is shared across PRNs
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# NH10 packed newest-first as `0x035` (= `0000110101`) — bit `i` holds
# NH10 chip `9 - i`, matching the prompt buffer's newest-in-bit-0 fill
# order. Shared across PRNs.
@inline _packed_secondary_code(::Type{B}, ::GPSL5I, ::Integer) where {B<:Unsigned} =
    B(0x035)

"""
$(SIGNATURES)

Secondary-code sync detector for GPS L5Q — the generic
[`_detect_secondary_code_sync`](@ref) rotation search over the NH20
secondary code (window 20; the default 2.5 % tolerance discretizes to 0,
i.e. exact match). L5Q is a pilot (no navigation data); the 20-chip
Neuman-Hoffman overlay is the only sync feature, so the detector locks
after a single NH20 period in the worst case and reports the upcoming
integration's NH20 chip in `SyncResult.phase`. Returns [`SyncResult`](@ref).

The packed reference is derived generically from `get_secondary_code`
(see [`_packed_secondary_code`](@ref)); no bespoke packing is needed.
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

# L5I sync-search window is one NH10 period (10 blocks); L5Q's is one NH20
# period (20 blocks). UInt32 holds either with room to spare (the rotation
# search masks down to the low 10 / 20 bits).
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
@inline get_code_block_buffer_type(::GPSL5Q) = UInt32
