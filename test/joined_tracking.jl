#For now
function beamform(x)
    [0.5 0.5 0.5 0.5] * x
end

@testset "Joined TrackingLoop" begin
    l1_carrier = cis.(2 * π * 51 / 4e6 * (1:32000) + 1 / 3 * π)
    l5_carrier = cis.(2 * π * 101 / 4e7 * (1:320000) + 1 / 4 * π)
    l1_sampled_code = GNSSSignals.gen_sat_code(1:32000, 1.023e6, 2.1, 4e6, SATELLITE_1_CODE)
    l5_sampled_code = GNSSSignals.gen_sat_code(1:320000, 1.023e7, 3.2, 4e7, L5_SAT1_CODE)

    l1_incoming_signals = [1, 1] .* (l1_carrier[1:4000] .* l1_sampled_code[1:4000])'
    l5_incoming_signals = [1, 1] .* (l5_carrier[1:40000] .* l5_sampled_code[1:40000])'

    l1_gen_sampled_code, l1_get_code_phase = init_gpsl1_codes()
    l5_gen_sampled_code, l5_get_code_phase = init_gpsl5_codes()
    l1_system = Tracking.GNSSSystem(1575.42e6, 1023, 4e6, l1_gen_sampled_code, l1_get_code_phase)
    l5_system = Tracking.GNSSSystem(1176.45e6, 10230, 4e7, l1_gen_sampled_code, l1_get_code_phase)
    l1_results = Tracking.TrackingResults(1, 1 / 3 * π, 0.1, 2.0, [1.0 + 3.0 * im] )
    l5_results = Tracking.TrackingResults(1, 1 / 4 * π, 0.2, 3.0, [1.0 + 3.0 * im] )

    joined_tracking = Tracking.init_joined_tracking(l1_system, l5_system, l1_results, l5_results, 50, 100, 18.0, 1.0, 1  )
    joined_results = joined_tracking(l1_incoming_signals, l5_incoming_signals, beamform);

    #= track = Tracking.init_tracking(1/3 * π, 1575.43e6, 50, 0.0, 2.0, 1.023e6, 4e6, 18.0, 1.0, 1, l1_gen_sampled_code, l1_get_code_phase)
    next_track, code_phase, prompts_correlated_signals, carrier_freq = track(incoming_signals, beamform)
    @test code_phase ≈ 2.0 atol = 1e-3
    @test carrier_freq ≈ 50.0 atol = 1e-13
    for i = 0:6
        incoming_signals = [1, 1, 1, 1] .* (carrier[(4001 + 4000 * i):(8000 + 4000 * i)] .* sampled_code[(4001 + 4000 * i):(8000 + 4000 * i)])'
        next_track, code_phase, prompts_correlated_signals, carrier_freq = next_track(incoming_signals, beamform)
        @test code_phase ≈ 2.0 atol = 1e-3
        @test carrier_freq ≈ 50.0 atol = 1e-13       
    end  =#
end