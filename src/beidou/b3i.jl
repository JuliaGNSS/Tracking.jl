# BeiDou B3I is the legacy open-service signal on the 1268.52 MHz carrier
# (`B3I`): a single BPSK 10230-chip primary code at 10.23 Mcps (1 ms period),
# again with no open-service quadrature counterpart (`get_relative_power` 1.0).
#
# Structurally B1I at five times the chipping rate, with the same NH20 overlay
# on MEO/IGSO and none on GEO (BDS-SIS-ICD-B3I-1.0 §5.2.1); the `b1i.jl` header
# applies verbatim.

"""
$(SIGNATURES)

Hard-path secondary-code sync detector for BeiDou B3I over the per-PRN NH20
overlay — see [`_detect_secondary_code_sync`](@ref); behaves as
[`BeiDouB1I`](@ref)'s. The default detector is the soft one
([`uses_soft_secondary_code_detection`](@ref)).
"""
@inline function detect_bit_or_secondary_code_sync(
    signal::BeiDouB3I,
    prn::Integer,        # selects the PRN's NH20 column (all-ones on GEO)
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# B3I is plain BPSK (`LOC`), so the C/A-style EarlyPromptLate default applies.
function get_default_correlator(::BeiDouB3I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Holds one NH20 period (see `get_code_block_buffer_type`).
@inline get_code_block_buffer_type(::BeiDouB3I) = UInt32
