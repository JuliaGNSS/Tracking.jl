@testset "Bit buffer" begin
    bit_buffer = Tracking.BitBuffer()
    @test @inferred(get_bits(bit_buffer)) == 0
    @test @inferred(Tracking.length(bit_buffer)) == 0
    @test @allocated(Tracking.BitBuffer()) == 0
    @test bit_buffer.prompt_accumulator == 0.0 + 0.0im
    @test bit_buffer.prompt_accumulator_integrated_code_blocks == 0

    secondary_code_or_bit_found = false
    prompt_correlator = 1.0 + 0.0im
    integrated_code_blocks = 1
    gpsl1 = GPSL1()
    next_bit_buffer = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        integrated_code_blocks,
        secondary_code_or_bit_found,
        prompt_correlator,
    )

    @test next_bit_buffer.prompt_accumulator == 0.0 + 0.0im
    @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
    @test @inferred(get_bits(next_bit_buffer)) == 0
    @test @inferred(Tracking.length(next_bit_buffer)) == 0

    bit_buffer = Tracking.BitBuffer(0, 0, 20.0 + 0.0im, 1)
    secondary_code_or_bit_found = true

    next_bit_buffer = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        integrated_code_blocks,
        secondary_code_or_bit_found,
        prompt_correlator,
    )
    @test next_bit_buffer.prompt_accumulator == 21.0 + 0.0im
    @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 2
    @test @inferred(get_bits(next_bit_buffer)) == 0
    @test @inferred(Tracking.length(next_bit_buffer)) == 0

    integrated_code_blocks = 19

    next_bit_buffer = @inferred Tracking.buffer(
        gpsl1,
        bit_buffer,
        integrated_code_blocks,
        secondary_code_or_bit_found,
        prompt_correlator,
    )
    @test next_bit_buffer.prompt_accumulator == 0.0 + 0.0im
    @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
    @test @inferred(get_bits(next_bit_buffer)) == 1
    @test @inferred(Tracking.length(next_bit_buffer)) == 1
end
