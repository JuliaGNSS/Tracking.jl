"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL5I.
"""
function is_upcoming_integration_new_bit(gpsl5::GPSL5I, code_block_bits, num_code_blocks)
    num_code_blocks < 10 && return false
    masked_bit_synchronizer = code_block_bits & 0x3ff # First 10
    xored_bit_synchronizer = masked_bit_synchronizer ⊻ 0x35 # 0x35 == 0000110101
    # If xored_bit_synchronizer == 0 -> bit -1 and if xored_bit_synchronizer == 0x3ff -> bit 1
    xored_bit_synchronizer == 0 || xored_bit_synchronizer == 0x3ff
end

function get_default_correlator(gpsl5::GPSL5I, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 20-block sync-search window (2 × NH10) needs at least 20 bits.
@inline get_code_block_buffer_type(::GPSL5I) = UInt32
