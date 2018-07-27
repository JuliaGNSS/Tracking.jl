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
    init_carrier_freq = 50
    carrier = cis.(2 * π * 60 / 4e6 * (1:32000000) + 1 / 3 * π)
    sampled_code = GNSSSignals.gen_sat_code(1:32000000, 1.023e6, 2.1, 4e6, SATELLITE_1_CODE)
    incoming_signals = [1, 1, 1, 1] .* (carrier[1:4000] .* sampled_code[1:4000])'
    gen_sampled_code, get_code_phase = init_gpsl1_codes()
    track = Tracking.init_tracking(1/3 * π, 1575.43e6, init_carrier_freq, 0.0, 2.0, 1.023e6, 4e6, 18.0, 1.0, 1, gen_sampled_code, get_code_phase)
    next_track, results = track(incoming_signals, beamform)
    prompts_correlated_signals = results.prompts_correlated_signals
    println("prompts corr signals", prompts_correlated_signals)
    for i = 0:300
        incoming_signals = [1, 1, 1, 1] .* (carrier[(4001 + 4000 * i):(8000 + 4000 * i)] .* sampled_code[(4001 + 4000 * i):(8000 + 4000 * i)])'
        next_track, results = next_track(incoming_signals, beamform) 
    end
    for i = 301:600
        incoming_signals = [1, 1, 1, 1] .* (carrier[(4001 + 4000 * i):(8000 + 4000 * i)] .* sampled_code[(4001 + 4000 * i):(8000 + 4000 * i)])'
        next_track, results = next_track(incoming_signals, beamform)

        @test results.code_phase ≈ 2.1 atol = 0.02
        @test results.carrier_doppler + init_carrier_freq ≈ 60.0 atol = 0.2  
    end
    println("prompts corr signals after 600", prompts_correlated_signals, " carrier freq", results.carrier_doppler + init_carrier_freq, "code_phase", results.code_phase)
end
