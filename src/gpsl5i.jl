"""
$(SIGNATURES)

Secondary-code sync detector for GPS L5I.

Runs the generic [`_secondary_code_search`](@ref) rotation search over the
NH10 secondary code, packed newest-first as `0x035` (= `0000110101`) — bit
`i` holds NH10 chip `9 - i`, matching the prompt buffer's newest-in-bit-0
fill order. The negated-polarity (data-bit-1) case is handled inside the
search, so the detector locks after a single NH10 period in the worst case
and reports the upcoming integration's NH10 chip in `SyncResult.phase`.
Returns [`SyncResult`](@ref).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GPSL5I,
    ::Integer,           # PRN — ignored; NH10 is shared across PRNs
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    secondary_code_length = get_secondary_code_length(signal)   # 10 (NH10)
    num_code_blocks < secondary_code_length && return SyncResult(false, 0, Int8(0))
    # Tolerance is a percentage of the secondary-code window. The default
    # 2.5 % discretizes to floor(0.025 × 10) = 0 (exact match).
    # Adjustable via `get_bit_edge_or_secondary_code_tolerance(::GPSL5I)`.
    max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * secondary_code_length)
    _secondary_code_search(code_block_bits, B(0x035), secondary_code_length, max_errors)
end

function get_default_correlator(gpsl5::GPSL5I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Sync-search window is one NH10 period (10 blocks); UInt32 holds it with room
# to spare (the rotation search masks down to the low 10 bits).
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
