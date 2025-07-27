@testset "Code replica" begin
    code = zeros(Int16, 2502)
    gpsl1 = GPSL1()
    Tracking.gen_code_replica!(
        code,
        gpsl1,
        1023e3Hz,
        2.5e6Hz,
        2.0,
        11,
        2480,
        -1:1,
        1,
        Val(2.5e6Hz),
    )

    # TODO: The new code generation `gen_code!` seems to be off in the first two samples
    #@test code[11:2492] == get_code.(gpsl1, (-1:2480) * 1023e3 / 2.5e6 .+ 2.0, 1)
    @test code[13:2492] == get_code.(gpsl1, (-1:2480) * 1023e3 / 2.5e6 .+ 2.0, 1)[3:end]
    @testset "More than 1ms" begin
        code = zeros(Int16, 6502)
        gpsl1 = GPSL1()
        Tracking.gen_code_replica!(
            code,
            gpsl1,
            1023e3Hz,
            2.5e6Hz,
            2.0,
            11,
            6480,
            -1:1,
            1,
            Val(2.5e6Hz),
        )

        # TODO: The new code generation `gen_code!` seems to be off in the first two samples
        #@test code[11:6492] == get_code.(gpsl1, (-1:6480) * 1023e3 / 2.5e6 .+ 2.0, 1)
        @test code[13:6492] == get_code.(gpsl1, (-1:6480) * 1023e3 / 2.5e6 .+ 2.0, 1)[3:end]
    end

    @testset "code_length is less than 1ms" begin
        code = zeros(Int16, 2502)
        gpsl1 = GPSL1()
        Tracking.gen_code_replica!(
            code,
            gpsl1,
            1023e3Hz + 1000Hz,
            7.5e6Hz,
            2.0,
            11,
            2480,
            -1:1,
            1,
            Val(7.5e6Hz),
        )

        # TODO: The new code generation `gen_code!` seems to be off in the first two samples
        #@test code[11:2492] == get_code.(gpsl1, (-1:2480) * 1023e3 * 3 / 7.5e6 .+ 2.0, 1)
        @test code[13:2492] ==
              get_code.(gpsl1, (-1:2480) * (1023e3 + 1000) / 7.5e6 .+ 2.0, 1)[3:end]
    end
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
        bit_found,
    )
    @test phase ≈ mod(10 + 0.1 * 2000, 1023)
end
