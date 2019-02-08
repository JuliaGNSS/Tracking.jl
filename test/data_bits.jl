@testset "Buffer data bits" begin
    gpsl1 = GPSL1()
    gpsl5 = GPSL5()

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 0, 0)) == false
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 12 << 41 + 2^20-1, 40)) == true
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 12 << 41 + (2^20-1) << 20, 40)) == true
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 12 << 41 + 2^20-1, 20)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL5, 12 << 11 + 0x35, 10)) == true # 0x35 == 0000110101
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL5, 12 << 11 + 0x35, 7)) == false
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL5, 12 << 11 + 0x3ca, 10)) == true # 0x3ca == 1111001010

    data_bits = Tracking.DataBits(gpsl1)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl1, 5.0, 1)
    @test next_data_bits.synchronisation_buffer == 1
    @test next_data_bits.num_bits_in_synchronisation_buffer == 1
    @test next_data_bits.first_found_after_num_prns == -1
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL1}(12 << 40 + 2^19-1, 39, -1, 0.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl1, 5.0, 3)
    @test next_data_bits.synchronisation_buffer == 12 << 41 + 2^20-1
    @test next_data_bits.num_bits_in_synchronisation_buffer == 40
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL5}(12 << 10 + 0x1a, 9, -1, 0.0, 0, 0) # 0x1a == 000011010
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 3)
    @test next_data_bits.synchronisation_buffer == 12 << 11 + 0x35 # 0x35 == 0000110101
    @test next_data_bits.num_bits_in_synchronisation_buffer == 10
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 0.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 3)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 5.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 4)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 1
    @test next_data_bits.num_bits_in_buffer == 1

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, -10.0, 4)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 1

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 1, 1)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 4)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 3
    @test next_data_bits.num_bits_in_buffer == 2
end
