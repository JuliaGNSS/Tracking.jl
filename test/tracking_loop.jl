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
    @test code_phase ≈ 2.0
    @test @inferred(beamform(prompts_correlated_signals))[1] == prompt_beamformed_signal

    compare_beamformed_signal = [-16.0-2.94209e-15im -24.0+5.55112e-17im -8.0+6.16174e-15im]
    compareDLL, compareReplica, comparePhase = @inferred Tracking.init_DLL(2.0, 1023e3,  4e6, 1.0, 1e-3, 1)
    next_compareDLL, next_compareReplica, next_comparePhase = @inferred compareDLL(compare_beamformed_signal)
    next_compare_beamformed_signal = [-11.4127-3.7082im -19.0211-6.18034im -11.4127-3.7082im]

    next_tracking_loop, code_phase, prompts_correlated_signals, prompt_beamformed_signal = next_tracking_loop(incoming_signals)
    @test @inferred(next_compareDLL(next_compare_beamformed_signal))[3] ≈ code_phase
    println(code_phase)
    @test @inferred(beamform(prompts_correlated_signals))[1] == prompt_beamformed_signal
end