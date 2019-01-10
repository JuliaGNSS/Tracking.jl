function is_upcoming_integration_new_bit(::Type{GPSL5}, synchronisation_buffer, num_bits_in_buffer)
    if num_bits_in_buffer < 10
        return false
    end
    masked_bit_synchronizer = synchronisation_buffer & 0x3ff # First 10
    xored_bit_synchronizer = masked_bit_synchronizer ⊻ 0x35 # 0x35 == 0000110101
    # If xored_bit_synchronizer == 0 -> bit -1 and if xored_bit_synchronizer == 0x3ff -> bit 1
    xored_bit_synchronizer == 0 || xored_bit_synchronizer == 0x3ff
end

function adjust_code_phase(system::GPSL5, data_bits, phase)
    ifelse(found(data_bits), phase, mod(phase, system.code_length_wo_neuman_hofman_code))
end

function calc_integrated_prns(system::GPSL5, integrated_samples, sample_freq)
    ceil(Int, integrated_samples / (sample_freq / system.code_freq * system.code_length_wo_neuman_hofman_code))
end
