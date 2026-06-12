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

function get_default_correlator(gpsl5::GPSL5I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Sync-search window is one NH10 period (10 blocks); UInt32 holds it with room
# to spare (the rotation search masks down to the low 10 bits).
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
