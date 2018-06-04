#For now
function beamform(x)
[0.5 0.5 0.5 0.5] * (x)
end


@testset "correlate_and_downconvert" begin[]
    x = [1im 2im 3im; 4im 5im 6im]
    replica_signal = [1im 2im 3im]
    downconverted_signals = Tracking._downconvert(x, replica_signal)
    @test @inferred downconverted_signals == [1 + 0im 4 + 0im 9 + 0im; 4 + 0im 10 + 0im 18 + 0im]

    replica_codes = [[1, 1, -1], [-1, -1, 1]]
    correlated_signals = map(replica -> Tracking._correlate(downconverted_signals, replica), replica_codes)
    @test @inferred correlated_signals[2][1] == 4 + 0im
end
@testset "Tracking_loop" begin
    test_signal = cis.(2 * π * 50 / 4e6 * (1:4000) + 1 / 3 * π)
    incoming_signals = [1,1,1,1] .* test_signal'
    tracking_loop = Tracking.init_tracking(1/3 * π, 50, 2.0, 1023e3, 1e-3, 4e6, beamform, 18.0, 1.0, 1)
    next_tracking_loop, code_phase, prompts_correlated_signals, prompt_beamformed_signal = tracking_loop(incoming_signals)
    @test @inferred code_phase ≈ 2.0
    @test @inferred prompt_beamformed_signal == beamform(prompts_correlated_signals)[1]

    next_tracking_loop, code_phase, prompts_correlated_signals, prompt_beamformed_signal = next_tracking_loop(incoming_signals)
    @test @inferred code_phase ≈ 2.00096
    @test @inferred prompt_beamformed_signal == beamform(prompts_correlated_signals)[1]
end