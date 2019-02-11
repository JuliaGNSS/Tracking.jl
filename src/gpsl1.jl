"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL1.
"""
function is_upcoming_integration_new_bit(::Type{GPSL1}, synchronisation_buffer, num_bits_in_buffer)
    if num_bits_in_buffer < 40
        return false
    end
    masked_bit_synchronizer = synchronisation_buffer & 0xffffffffff # First 40 bits
    # Upcoming integration will be a new bit if masked_bit_synchronizer contains
    # 20 zeros and 20 ones or 20 ones and 20 zeros
    masked_bit_synchronizer == 0xfffff || masked_bit_synchronizer == 0xfffff00000
end
