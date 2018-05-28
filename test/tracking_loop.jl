function beamform(x)
    sum(x)
end

@testset "Tracking_loop" begin
    tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = Tracking.init_tracking(init_PLL, init_DLL, 20.0, 50, 20.0, 1023e3, 1e-3, 4e6, beamform, 4000, 18.0, 1.0, 1)
    @test code_phase == 20
    @test prompt_correlated_signal = 1
    @test prompt_beamformed_signal = 2
    @test tracking_loop([1 3 45 6 72  34 7])
end