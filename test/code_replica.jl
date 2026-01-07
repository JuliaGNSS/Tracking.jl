module CodeReplicaTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, get_code
using Tracking: gen_code_replica!, update_code_phase

@testset "Code replica" begin
    code = zeros(Int16, 2502)
    gpsl1 = GPSL1()
    gen_code_replica!(code, gpsl1, 1023e3Hz, 2.5e6Hz, 2.0, 11, 2480, -1:1, 1, Val(2.5e6Hz))

    @test code[11:2492] == get_code.(gpsl1, (-1:2480) * 1023e3 / 2.5e6 .+ 2.0, 1)
    @testset "More than 1ms" begin
        code = zeros(Int16, 6502)
        gpsl1 = GPSL1()
        gen_code_replica!(
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

        @test code[11:6492] == get_code.(gpsl1, (-1:6480) * 1023e3 / 2.5e6 .+ 2.0, 1)
    end

    @testset "code_length is less than 1ms" begin
        code = zeros(Int16, 2502)
        gpsl1 = GPSL1()
        gen_code_replica!(
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

        @test code[11:2492] == get_code.(gpsl1, (-1:2480) * (1023e3 + 1000) / 7.5e6 .+ 2.0, 1)
    end
end

@testset "Update code phase" begin
    gpsl1 = GPSL1()
    code_phase = 10
    code_frequency = 10Hz
    sampling_frequency = 100Hz
    num_samples = 2000
    bit_found = true
    phase = @inferred update_code_phase(
        gpsl1,
        num_samples,
        code_frequency,
        sampling_frequency,
        code_phase,
        bit_found,
    )
    @test phase ≈ mod(10 + 0.1 * 2000, 1023)
end

end
