@testset "Code replica" begin

    code = zeros(Int16, 2502)

    Tracking.gen_code_replica!(
        code,
        GPSL1,
        1023e3,
        2.5e6,
        2.0,
        11,
        2480,
        1,
        1
    )

    @test code[11:2492] == get_code.(GPSL1, (-1:2480) * 1023e3 / 2.5e6 .+ 2.0, 1)
end

@testset "Update code phase" begin
    code_phase = 10
    code_frequency = 10Hz
    sample_frequency = 100Hz
    num_samples = 2000
    bit_found = true
    phase = @inferred Tracking.update_code_phase(
        GPSL1,
        num_samples,
        code_frequency,
        sample_frequency,
        code_phase,
        bit_found
    )
    @test phase ≈ mod(10 + 0.1 * 2000, 1023)
end
