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

function get_default_correlator(gpsl5q::GPSL5Q, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Sync-search window is one NH20 period (20 blocks); UInt32 holds it (the
# rotation search masks down to the low 20 bits).
@inline get_code_block_buffer_type(::GPSL5Q) = UInt32
