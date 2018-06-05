#For now
function beamform(x)
[0.5 0.5 0.5 0.5] * x
end


@testset "correlate" begin
    x = [1im 2im 3im; 4im 5im 6im]
    replica_signal = [1im 2im 3im]
    downconverted_signals = @inferred Tracking.downconvert(x, replica_signal)
    @test downconverted_signals == [1 + 0im 4 + 0im 9 + 0im; 4 + 0im 10 + 0im 18 + 0im]
end
@testset "downconvert" begin
    x = [1im 2im 3im; 4im 5im 6im]
    replica_signal = [1im 2im 3im]
    downconverted_signals = @inferred Tracking.downconvert(x, replica_signal)
    replica_codes = [1, 1, -1]
    correlated_signals = @inferred Tracking.correlate(downconverted_signals, replica_codes)
    @test correlated_signals[1] == -4 + 0im
end
@testset "Tracking_loop" begin
    test_signal = cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π)
    samples_code = GNSSSignals.gen_sat_code(1:4000,1.023e6,2.0,4e6,SATELLITE_1_CODE)
    incoming_signals = [1,1,1,1] .* (test_signal .* samples_code)'
    tracking_loop = Tracking.init_tracking(1/3 * π, 50, 2.0, 1023e3, 1e-3, 4e6, beamform, 18.0, 1.0, 1)
    next_tracking_loop, code_phase, prompts_correlated_signals, prompt_beamformed_signal = tracking_loop(incoming_signals)
    @test code_phase ≈ 2.0 atol = 0.001
    @test @inferred(beamform(prompts_correlated_signals))[1] == prompt_beamformed_signal
    next_tracking_loop, code_phase, prompts_correlated_signals, prompt_beamformed_signal = next_tracking_loop(incoming_signals)
    @test code_phase ≈ 2.0 atol=0.003
    @test @inferred(beamform(prompts_correlated_signals))[1] == prompt_beamformed_signal
end