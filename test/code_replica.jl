@testset "Code replica" begin

    code = zeros(Int16, 2502)
    gpsl1 = GPSL1()
    Tracking.gen_code_replica!(
        code,
        gpsl1,
        1023e3,
        2.5e6,
        2.0,
        11,
        2480,
        SVector(-1, 0, 1),
        1
    )

    @test code[11:2492] == get_code.(gpsl1, (-1:2480) * 1023e3 / 2.5e6 .+ 2.0, 1)
end

@testset "Update code phase" begin
    gpsl1 = GPSL1()
    code_phase = 10
    code_frequency = 10Hz
    sampling_frequency = 100Hz
    num_samples = 2000
    bit_found = true
    phase = @inferred Tracking.update_code_phase(
        gpsl1,
        num_samples,
        code_frequency,
        sampling_frequency,
        code_phase,
        bit_found
    )
    @test phase ≈ mod(10 + 0.1 * 2000, 1023)
end
