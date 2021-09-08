@testset "Bit buffer" begin
    bit_buffer = Tracking.BitBuffer()
    @test @inferred(get_bits(bit_buffer)) == 0
    @test @inferred(Tracking.length(bit_buffer)) == 0
    @test @allocated(Tracking.BitBuffer()) == 0

    prompt_accumulator = 0.0 + 0.0im
    secondary_code_or_bit_found = false
    prev_code_phase = 0.0
    code_phase = 1023
    prompt_correlator = 1.0 + 0.0im
    integration_time = 5ms
    gpsl1 = GPSL1()
    next_bit_buffer, next_prompt_accumulator = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        prompt_accumulator,
        secondary_code_or_bit_found,
        prev_code_phase,
        code_phase,
        integration_time,
        prompt_correlator
    )

    @test next_prompt_accumulator == 0.0 + 0.0im
    @test @inferred(get_bits(next_bit_buffer)) == 0
    @test @inferred(Tracking.length(next_bit_buffer)) == 0

    prompt_accumulator = 20.0 + 0.0im
    secondary_code_or_bit_found = true
    prev_code_phase = 0.0
    code_phase = 1023
    prompt_correlator = 1.0 + 0.0im
    integration_time = 5ms

    next_bit_buffer, next_prompt_accumulator = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        prompt_accumulator,
        secondary_code_or_bit_found,
        prev_code_phase,
        code_phase,
        integration_time,
        prompt_correlator
    )
    @test next_prompt_accumulator == 21.0 + 0.0im
    @test @inferred(get_bits(next_bit_buffer)) == 0
    @test @inferred(Tracking.length(next_bit_buffer)) == 0

    prompt_accumulator = 20.0 + 0.0im
    secondary_code_or_bit_found = true
    prev_code_phase = 9000
    code_phase = 1023
    prompt_correlator = 1.0 + 0.0im
    integration_time = 20ms

    next_bit_buffer, next_prompt_accumulator = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        prompt_accumulator,
        secondary_code_or_bit_found,
        prev_code_phase,
        code_phase,
        integration_time,
        prompt_correlator
    )
    @test next_prompt_accumulator == 0.0 + 0.0im
    @test @inferred(get_bits(next_bit_buffer)) == 1
    @test @inferred(Tracking.length(next_bit_buffer)) == 1

    prompt_accumulator = 20.0 + 0.0im
    secondary_code_or_bit_found = true
    prev_code_phase = 0
    code_phase = 0
    prompt_correlator = 1.0 + 0.0im
    integration_time = 20ms

    next_bit_buffer, next_prompt_accumulator = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        prompt_accumulator,
        secondary_code_or_bit_found,
        prev_code_phase,
        code_phase,
        integration_time,
        prompt_correlator
    )
    @test next_prompt_accumulator == 0.0 + 0.0im
    @test @inferred(get_bits(next_bit_buffer)) == 1
    @test @inferred(Tracking.length(next_bit_buffer)) == 1


end
