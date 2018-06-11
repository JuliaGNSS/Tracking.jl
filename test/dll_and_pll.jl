@testset "initialize PLL" begin
    correlator_output = [0.5 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im]
    PLL, sampled_carrier, phase = @inferred Tracking.init_PLL(1 / 3 * π, 50, 4e6, 18.0, 1e-3)
    @test sampled_carrier ≈ cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π)
    @test phase ≈ mod2pi((2 * π * 50 / 4e6) * 4000 + 1 / 3 * π)
    next_PLL, next_sampled_carrier, next_phase, frequency_update = @inferred PLL(correlator_output)
    @test next_phase ≈ mod2pi((2 * π * 50 / 4e6) * 8000 + 1 / 3 * π)
    @test next_sampled_carrier ≈ cis.(2 * π * 50 / 4e6 * (4001:8000) + 1 / 3 * π)
    @test frequency_update == 0.0
end

@testset "PLL aiding" begin
    correlator_output = [0.5 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im]
    PLL, sampled_carrier, phase = @inferred Tracking.init_PLL(1 / 3 * π, 50, 4e6, 18.0, 1e-3)
    next_PLL, next_sampled_carrier, next_phase, frequency_update = @inferred PLL(correlator_output, 1)
    @test next_phase ≈ mod2pi((2 * π * (50 ) / 4e6) * 4000 + 1 / 3 * π + (2 * π * (50 + 1) / 4e6) * 4000) # phase after first loop + init phase + phase after second loop with velocity aiding
    @test next_sampled_carrier ≈ cis.(2 * π * (50 + 1 + frequency_update) / 4e6 * (1:4000) + phase)
    next_PLL2, next_sampled_carrier2, next_phase2, frequency_update2 = @inferred next_PLL(correlator_output, 1)
    @test next_phase2 ≈ mod2pi((2 * π * (50 ) / 4e6) * 4000 + 1 / 3 * π + (2 * π * (50 + 1) / 4e6) * 4000 * 2) # phase after first loop + init phase + phase after second loop with velocity aiding + 3rd loop phase
    @test next_sampled_carrier2 ≈ cis.(2 * π * (50 + 1 + frequency_update) / 4e6 * (1:4000) + next_phase)
end

@testset "inizialize DLL" begin
    correlator_output = [0.5 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im]
    DLL, sampled_code, phase = @inferred Tracking.init_DLL(2, 1023e3, 4e6, 1, 1e-3, 1)
    gen_sampled_code, get_code_phase = @inferred GNSSSignals.init_gpsl1_codes()
    @test sampled_code[2] == gen_sampled_code(1:4000, 1023e3, 2, 4e6, 1)
    @test phase == get_code_phase(4000, 1023e3, 2, 4e6)
    @test sampled_code[1] ≈ gen_sampled_code(-1:3998, 1023e3, 2.0, 4e6, 1)
    next_DLL, next_sampled_code, next_phase, frequency_update = @inferred DLL(correlator_output);
    @test next_sampled_code[3] == gen_sampled_code(3:4002, 1023e3 + frequency_update, next_phase, 4e6, 1)
    @test next_sampled_code[2] == gen_sampled_code(1:4000, 1023e3 + frequency_update, next_phase, 4e6, 1)
    @test frequency_update == 0.0
end

@testset "DLL aiding" begin
    correlator_output = [0.5 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im]
    gen_sampled_code, get_code_phase = @inferred GNSSSignals.init_gpsl1_codes()
    DLL, sampled_code, phase = @inferred Tracking.init_DLL(2, 1023e3, 4e6, 1, 1e-3, 1)
    next_DLL, next_sampled_code, next_phase, frequency_update = @inferred DLL(correlator_output, 1);
    @test next_phase ≈ mod((1023e3 + frequency_update + 1) / 4e6 * 4000 + 2.0 + 1023 / 2, 1023) - 1023 / 2
    @test next_sampled_code[2] == gen_sampled_code(1:4000, 1023e3 + 1 + frequency_update, phase, 4e6, 1)
    @test next_sampled_code[3] == gen_sampled_code(3:4002, 1023e3 + 1 + frequency_update, phase, 4e6, 1)
    next_DLL2, next_sampled_code2, next_phase2, frequency_update2 = @inferred next_DLL(correlator_output, 1);
    @test next_phase2 ≈ mod((1023e3 + frequency_update2 + 1) / 4e6 * 8000 + 2.0 + 1023 / 2, 1023) - 1023 / 2
    @test next_sampled_code2[2] == gen_sampled_code(1:4000, 1023e3 + 1 + frequency_update, next_phase, 4e6, 1)
    @test next_sampled_code2[3] == gen_sampled_code(3:4002, 1023e3 + 1 + frequency_update, next_phase , 4e6, 1)
end