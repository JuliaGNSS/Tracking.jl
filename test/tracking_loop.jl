#For now
function beamform(x)
[0.5 0.5 0.5 0.5] * x
end
scale_factor = 1.023e6/1575.43e6

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
    test_signal = cis.(2 * π * 50 / 4e6 * (1:32000) + 1 / 3 * π)
    samples_code = GNSSSignals.gen_sat_code(1:32000, 1.023e6, 2.0, 4e6, SATELLITE_1_CODE) * 10^(-20/20)
    incoming_signals1 = [1,1,1,1] .* (test_signal[1:4000] .* samples_code[1:4000])'
    incoming_signals2 = [1,1,1,1] .* (test_signal[4001:8000] .* samples_code[4001:8000])'
    incoming_signals3 = [1,1,1,1] .* (test_signal[8001:12000] .* samples_code[8001:12000])'
    tracking_loop = Tracking.init_tracking(1/3 * π, 50, 2.0, 1023e3, 1e-3, 4e6, beamform, 18.0, 1.0, 1, scale_factor)
    next_tracking_loop, code_phase, prompts_correlated_signals = tracking_loop(incoming_signals1)
    @test code_phase ≈ 2.0 atol = 0.001
    next_tracking_loop, code_phase, prompts_correlated_signals = next_tracking_loop(incoming_signals2)
    @test code_phase ≈ 2.0 atol = 0.001
    next_tracking_loop, code_phase, prompts_correlated_signals = next_tracking_loop(incoming_signals3)
    @test code_phase ≈ 2.0 atol = 0.001
    for i = 0:4
        incoming_signals = [1,1,1,1] .* (test_signal[(12001 +4000*i):(16000+4000*i)] .* samples_code[(12001+4000*i):(16000+4000*i)])'
        next_tracking_loop, code_phase, prompts_correlated_signals = next_tracking_loop(incoming_signals)
        @test code_phase ≈ 2.0 atol = 0.001
    end

end