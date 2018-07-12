#For now
function beamform(x)
[0.5 0.5 0.5 0.5] * x
end
scale_factor = 1.023e6/1575.43e6

@testset "Correlate" begin
    x = [1im 2im 3im; 4im 5im 6im]
    replica_signal = [1im 2im 3im]
    downconverted_signals = @inferred Tracking.downconvert(x, replica_signal)
    @test downconverted_signals == [1 + 0im 4 + 0im 9 + 0im; 4 + 0im 10 + 0im 18 + 0im]
end
@testset "Downconvert" begin
    x = [1im 2im 3im; 4im 5im 6im]
    replica_signal = [1im 2im 3im]
    downconverted_signals = @inferred Tracking.downconvert(x, replica_signal)
    replica_codes = [1, 1, -1]
    correlated_signals = @inferred Tracking.correlate(downconverted_signals, replica_codes)
    @test correlated_signals[1] == -4 + 0im
end
@testset "Tracking loop" begin
    carrier = cis.(2 * π * 50 / 4e6 * (1:32000) + 1 / 3 * π)
    sampled_code = GNSSSignals.gen_sat_code(1:32000, 1.023e6, 2.0, 4e6, SATELLITE_1_CODE) * 10^(-20/20)
    incoming_signals = [1, 1, 1, 1] .* (carrier[1:4000] .* sampled_code[1:4000])'
    gen_sampled_code, get_code_phase = init_gpsl1_codes()
    track = Tracking.init_tracking(1/3 * π, 1575.43e6, 50, 0.0, 2.0, 1.023e6, 4e6, 18.0, 1.0, 1, gen_sampled_code, get_code_phase)
    next_track, code_phase, prompts_correlated_signals, carrier_freq = track(incoming_signals, beamform)
    @test code_phase ≈ 2.0 atol = 1e-3
    @test carrier_freq ≈ 50.0 atol = 1e-13
    for i = 0:6
        incoming_signals = [1, 1, 1, 1] .* (carrier[(4001 + 4000 * i):(8000 + 4000 * i)] .* sampled_code[(4001 + 4000 * i):(8000 + 4000 * i)])'
        next_track, code_phase, prompts_correlated_signals, carrier_freq = next_track(incoming_signals, beamform)
        @test code_phase ≈ 2.0 atol = 1e-3
        @test carrier_freq ≈ 50.0 atol = 1e-13       
    end

end