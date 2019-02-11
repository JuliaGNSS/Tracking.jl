function DataBits(system::T) where T <: AbstractGNSSSystem
    DataBits{T}(0, 0, -1, 0, 0, 0)
end

function DataBits(data_bits::DataBits{T}, buffer, num_bits_in_buffer) where T <: AbstractGNSSSystem
    DataBits{T}(data_bits.synchronisation_buffer, data_bits.num_bits_in_synchronisation_buffer, data_bits.first_found_after_num_prns, data_bits.prompt_accumulator, buffer, num_bits_in_buffer)
end

"""
$(SIGNATURES)

Buffers the data bits.
"""
function buffer(data_bits, system::T, prompt_real, num_integrated_prns) where T<:AbstractGNSSSystem
    if found(data_bits)
        prompt_accumulator = data_bits.prompt_accumulator + prompt_real
        bitbuffer = ifelse(data_bits.first_found_after_num_prns == num_integrated_prns, data_bits.buffer << 1 + UInt(prompt_accumulator > 0), data_bits.buffer)
        num_bits_in_buffer = ifelse(data_bits.first_found_after_num_prns == num_integrated_prns, data_bits.num_bits_in_buffer + 1, data_bits.num_bits_in_buffer)
        new_prompt_accumulator = ifelse(data_bits.first_found_after_num_prns == num_integrated_prns, zero(prompt_accumulator), prompt_accumulator)
        return DataBits{T}(data_bits.synchronisation_buffer, data_bits.num_bits_in_synchronisation_buffer, data_bits.first_found_after_num_prns, new_prompt_accumulator, bitbuffer, num_bits_in_buffer)
    else
        synchronisation_buffer = data_bits.synchronisation_buffer << 1 + UInt(prompt_real > 0)
        num_bits_in_synchronisation_buffer = data_bits.num_bits_in_synchronisation_buffer + 1
        first_found_after_num_prns = is_upcoming_integration_new_bit(T, synchronisation_buffer, num_bits_in_synchronisation_buffer) ? num_integrated_prns + 1 : -1
        return DataBits{T}(synchronisation_buffer, num_bits_in_synchronisation_buffer, first_found_after_num_prns, data_bits.prompt_accumulator, data_bits.buffer, data_bits.num_bits_in_buffer)
    end
end

"""
$(SIGNATURES)

Checks if a possible data bit transition has been found yet.
"""
function found(data_bits::DataBits)
    data_bits.first_found_after_num_prns != -1
end
