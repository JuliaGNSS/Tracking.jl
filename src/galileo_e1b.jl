"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GalileoE1B.
"""
function is_upcoming_integration_new_bit(::Type{GalileoE1B}, synchronisation_buffer, num_bits_in_buffer)
    if num_bits_in_buffer < 8
        return false
    end
    masked_bit_synchronizer = synchronisation_buffer & 0xff # First 8 bits
    # Upcoming integration will be a new bit if masked_bit_synchronizer contains
    # 20 zeros and 20 ones or 20 ones and 20 zeros
    masked_bit_synchronizer == 0xf || masked_bit_synchronizer == 0xf0
end
