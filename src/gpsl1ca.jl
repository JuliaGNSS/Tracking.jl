"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL1CA.
"""
function is_upcoming_integration_new_bit(gpsl1::GPSL1CA, code_block_bits, num_code_blocks)
    num_code_blocks < 40 && return false
    masked_bit_synchronizer = code_block_bits & 0xffffffffff # First 40 bits
    # Upcoming integration will be a new bit if masked_bit_synchronizer contains
    # 20 zeros and 20 ones or 20 ones and 20 zeros
    masked_bit_synchronizer == 0xfffff || masked_bit_synchronizer == 0xfffff00000
end

"""
$(SIGNATURES)

Get the default correlator for the given GNSS system. Returns an
EarlyPromptLateCorrelator for GPS L1 or a VeryEarlyPromptLateCorrelator
for systems like Galileo E1B that use BOC modulation.
"""
function get_default_correlator(gpsl1::GPSL1CA, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 40-block sync-search window (2 × 20 blocks/symbol) needs at least 40 bits.
@inline get_code_block_buffer_type(::GPSL1CA) = UInt64
