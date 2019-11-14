"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL1.
"""
function is_upcoming_integration_new_bit(::Type{GPSL1}, prns, num_prns)
    num_prns < 40 && return false
    masked_bit_synchronizer = prns & 0xffffffffff # First 40 bits
    # Upcoming integration will be a new bit if masked_bit_synchronizer contains
    # 20 zeros and 20 ones or 20 ones and 20 zeros
    masked_bit_synchronizer == 0xfffff || masked_bit_synchronizer == 0xfffff00000
end

function get_default_correlator(::Type{GPSL1}, num_ants::NumAnts{N}) where N
    EarlyPromptLateCorrelator(num_ants)
end
