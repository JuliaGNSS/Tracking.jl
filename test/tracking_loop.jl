#For now
function beamform(x)
    beamformed_x = map(x_per_phaseshift -> [0.5 0.5 0.5 0.5] * (x_per_phaseshift), x)
    hcat(beamformed_x...)
end

@testset "Tracking_loop" begin
    incoming_signals = [cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π), cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π), cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π), cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π)]
    correlator_output = [0.5 + 0.0im 1.0 + 0.0im 0.5 + 0.0im]
    tracking_loop = Tracking.init_tracking(Tracking.init_PLL, Tracking.init_DLL, 1/3 * π, 50, 2.0, 1023e3, 1e-3, 4e6, beamform, 4000, 18.0, 1.0, 1)
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = tracking_loop(incoming_signals)
    @test code_phase == 20
    @test prompt_correlated_signal == 1
    @test prompt_beamformed_signal == 2
    """
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = next_tracking_loop(incoming_signal)
    println("Next evaluation step -------------------------------------------------------------")
    @test code_phase == 20
    @test prompt_correlated_signal == 1
    @test prompt_beamformed_signal == 2
    """
end