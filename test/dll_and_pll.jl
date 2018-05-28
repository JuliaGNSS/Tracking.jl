@testset "initialize PLL" begin
    correlator_output = [0.5 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im]
    PLL, sampled_carrier, phase = Tracking.init_PLL(1 / 3 * π, 50, 4000, 4e6, 18.0, 1e-3)
    @test sampled_carrier ≈ cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π)
    @test phase ≈ mod2π((2 * π * 50 / 4e6) * 4000 + 1 / 3 * π)
    next_PLL, next_sampled_carrier, next_phase = PLL(correlator_output)
    @test next_phase ≈ mod2π((2 * π * 50 / 4e6) * 4000 + 1 / 3 * π)
    @test next_sampled_carrier ≈ cis.(2 * π * 50 / 4e6 * (1:4000) + next_phase)
end


@testset "inizialize DLL" begin
    correlator_output = [0.5 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im]
    DLL, sampled_code, phase = Tracking.init_DLL(2, 1023e3, 4000, 4e6, 1, 1e-3, 1)
    gen_sampled_code, get_code_phase = GNSSSignals.init_gpsl1_codes()
    @test sampled_code[2] == gen_sampled_code(1:4000, 1023e3, 2, 4e6, 1)
    @test phase == get_code_phase(4000, 1023e3, 2, 4e6)
    @test sampled_code[1] == gen_sampled_code(1:4000, 1023e3, 1.5, 4e6, 1)
    next_DLL, next_sampled_code, next_phase = DLL(correlator_output);
    @test sampled_code[3] == gen_sampled_code(1:4000, 1023e3, 2.5, 4e6, 1)
    @test next_sampled_code[2] == gen_sampled_code(1:4000, 1023e3, phase, 4e6, 1)
end

