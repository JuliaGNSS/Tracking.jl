#For now
# This takes an 3x4 Matrix as an input x
function beamform(x)
    beamformed_x = map(x_per_phaseshift -> [0.5 0.5 0.5 0.5] * (x_per_phaseshift), x)
    hcat(beamformed_x...)
end


@testset "correlate_and_downconvert" begin
    x = [[1im 2im 3im], [4im 5im 6im]]
    replica_signal = [1im 1im 1im]
    downconverted_signals = Tracking._downconvert(x, replica_signal)
    @test @inferred downconverted_signals == [[1 + 0im 2 + 0im 3 + 0im], [4 + 0im 5 + 0im 6 + 0im]]

    replica_codes = [[1 1 -1], [-1 -1 1]]
    correlated_signals = Tracking._correlate(downconverted_signals, replica_codes)
    @test @inferred correlated_signals[1][1] == [1+0im 2+0im 3+0im; 1+0im 2+0im 3+0im; -1+0im -2+0im -3+0im]
end

@testset "Tracking_loop" begin
    test_signal = cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π)
    incoming_signals = [test_signal, test_signal, test_signal, test_signal]
    tracking_loop = Tracking.init_tracking(1/3 * π, 50, 2.0, 1023e3, 1e-3, 4e6, beamform, 18.0, 1.0, 1)
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = tracking_loop(incoming_signals)
    @test @inferred code_phase ≈ 2.0
    @test @inferred prompt_beamformed_signal == beamform([prompt_correlated_signal,prompt_correlated_signal,prompt_correlated_signal])[2]

    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = next_tracking_loop(incoming_signals)
    @test @inferred code_phase ≈ 2.00096
    @test @inferred prompt_beamformed_signal == beamform([prompt_correlated_signal,prompt_correlated_signal,prompt_correlated_signal])[2]
end