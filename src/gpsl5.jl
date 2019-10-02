"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL5.
"""
function is_upcoming_integration_new_bit(::Type{GPSL5}, synchronisation_buffer, num_bits_in_buffer)
    if num_bits_in_buffer < 10
        return false
    end
    masked_bit_synchronizer = synchronisation_buffer & 0x3ff # First 10
    xored_bit_synchronizer = masked_bit_synchronizer ⊻ 0x35 # 0x35 == 0000110101
    # If xored_bit_synchronizer == 0 -> bit -1 and if xored_bit_synchronizer == 0x3ff -> bit 1
    xored_bit_synchronizer == 0 || xored_bit_synchronizer == 0x3ff
end

"""
$(SIGNATURES)

Adjusts the code phase if the data bit has not been found yet.
"""
function adjust_code_phase(data_bits::DataBits{GPSL5}, phase)
    ifelse(found(data_bits), phase, mod(phase, get_shortest_code_length(GPSL5)))
end

"""
$(SIGNATURES)

Calculates the integrated PRNs specifically for GPSL5.
"""
function calc_integrated_prns(::Type{GPSL5}, integrated_samples, sample_freq)
    ceil(Int, integrated_samples / (sample_freq / get_code_frequency(GPSL5) * get_shortest_code_length(GPSL5)))
end
