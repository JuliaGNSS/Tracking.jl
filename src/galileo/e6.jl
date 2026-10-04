# Galileo E6 is the civil pair on the 1278.75 MHz carrier (`E6`): E6-B carries
# the C/NAV message — and with it the High Accuracy Service — at 1000 sym/s,
# E6-C is its dataless pilot under a per-SVID 100-chip CS100 overlay. Both share
# the 1 ms / 5115-chip primary code (5.115 Mcps) and are tracked as two
# independent BPSK(5) signals, the same data/pilot split as Galileo E1 and E5.
#
# E6-C sits at π carrier phase against E6-B (OS SIS ICD Eq. 10); GNSSSignals
# reports it via `get_carrier_phase_offset`, which the shared driver loop reads
# generically.

"""
$(SIGNATURES)

Symbol-sync detector for Galileo E6-B: one C/NAV symbol per 1 ms primary code
period (1000 sym/s; OS SIS ICD v2.2 Table 5), so it locks immediately — see
[`_detect_symbol_is_code_block_sync`](@ref).
"""
@inline detect_bit_or_secondary_code_sync(
    signal::GalileoE6B,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for the Galileo E6-C pilot over the
per-SVID CS100 code (E6-B/C Codes Technical Note §2.4: CS100ₙ for SVID `n`, as
[`GalileoE5aQ`](@ref); 100 ms cycle) — see [`_detect_secondary_code_sync`](@ref).
The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::GalileoE6C,
    prn::Integer,        # selects the SVID's CS100 column in the reference
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Both E6 components are plain BPSK(5) (`LOC`) — no subcarrier, so no
# split-spectrum side peaks — and the C/A-style EarlyPromptLate default applies.
function get_default_correlator(
    ::Union{GalileoE6B,GalileoE6C},
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(; num_ants)
end

# One symbol per primary code period: nothing to search, so the buffer is dead
# state (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::GalileoE6B) = UInt8
# E6-C: holds one CS100 period.
@inline get_code_block_buffer_type(::GalileoE6C) = UInt128
