"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GalileoE1B.
"""
function is_upcoming_integration_new_bit(
    galileo_e1b::GalileoE1B,
    code_block_bits,
    num_code_blocks,
)
    num_code_blocks < 8 && return false
    masked_bit_synchronizer = code_block_bits & 0xff # First 8 bits
    # Upcoming integration will be a new bit if masked_bit_synchronizer contains
    # 20 zeros and 20 ones or 20 ones and 20 zeros
    masked_bit_synchronizer == 0xf || masked_bit_synchronizer == 0xf0
end

# TODO: Very early very late correlator?
function get_default_correlator(galileo_e1b::GalileoE1B, num_ants::NumAnts = NumAnts(1))
    VeryEarlyPromptLateCorrelator(; num_ants)
end
