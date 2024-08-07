"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL5.
"""
function is_upcoming_integration_new_bit(gpsl5::GPSL5, prns, num_prns)
    num_prns < 10 && return false
    masked_bit_synchronizer = prns & 0x3ff # First 10
    xored_bit_synchronizer = masked_bit_synchronizer ⊻ 0x35 # 0x35 == 0000110101
    # If xored_bit_synchronizer == 0 -> bit -1 and if xored_bit_synchronizer == 0x3ff -> bit 1
    xored_bit_synchronizer == 0 || xored_bit_synchronizer == 0x3ff
end

function get_default_correlator(
    gpsl5::GPSL5,
    sampling_frequency,
    num_ants::NumAnts = NumAnts(1),
)
    EarlyPromptLateCorrelator(gpsl5, sampling_frequency; num_ants)
end
