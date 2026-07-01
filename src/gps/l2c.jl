# GPS L2C is a time-multiplexed pair of codes sharing the L2 carrier
# (1227.60 MHz): the moderate-length L2CM code carries the CNAV data
# stream, the very long L2CL code is a dataless pilot. Tracking.jl treats
# them as two independent signals (`GPSL2CM`, `GPSL2CL`), same as the L1C
# data/pilot split.

"""
$(SIGNATURES)

Symbol-sync detector for GPS L2CM.

L2CM broadcasts one CNAV symbol per 20 ms L2CM code period (50 sps;
IS-GPS-200N §3.2.2 — the 10230-chip code at 511.5 kcps is exactly one
symbol long), so the buffer of primary-block signs is itself the symbol
stream — there is no sub-symbol boundary to find. The detector therefore
reports `found = true` from the very first integration, leaving downstream
consumers (GNSSDecoder.jl) to resolve the residual ±1 polarity ambiguity
via the CNAV preamble's CRC.

Same shape as the GPS L1C-D / Galileo E1B "1-block-per-symbol" case.
"""
@inline detect_bit_or_secondary_code_sync(
    signal::GPSL2CM,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

"""
$(SIGNATURES)

"Sync" detector for GPS L2CL — a no-op. L2CL is a dataless pilot with no
secondary/overlay code and a single 767250-chip code (1.5 s period; a
whole code period is the coherent-integration unit). There is no bit and
no secondary code to lock, so the detector never reports `found`: the
tracker keeps one code block per integration (the `calc_num_code_blocks`
cap for a signal whose `get_secondary_code_length` is 1) and simply
tracks. Returns `SyncResult(false, 0, 0)`.
"""
@inline detect_bit_or_secondary_code_sync(::GPSL2CL, ::Integer, ::Unsigned, ::Integer) =
    SyncResult(false, 0, Int8(0))

# Both L2C components are BPSK (`LOC`), so the C/A-style EarlyPromptLate
# default applies.
function get_default_correlator(::Union{GPSL2CM,GPSL2CL}, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 1 CNAV symbol = 1 L2CM code period (50 sps, 20 ms period). No sub-symbol
# boundary to search for — the sync-search buffer is dead state at runtime;
# we still pick a concrete type to keep the `BitBuffer{B}` parameter chain
# stable. UInt8 is the smallest legal Unsigned.
@inline get_code_block_buffer_type(::GPSL2CM) = UInt8
# L2CL has no bit/secondary sync feature either; its search buffer is
# likewise dead state.
@inline get_code_block_buffer_type(::GPSL2CL) = UInt8
