#For now
function beamform(x)
    [0.5 0.5] * x
end

l1_scale_factor = 1.023e6/1575.43e6
l5_scale_factor = 1.023e6/1176.45e6

@testset "joined_tracking_loop" begin

    joined_tracking_loop = init_joined_tracking(
        #     l5_init_carrier_phase, l5_init_carrier_freq, l5_init_code_phase, l5_init_code_freq, l5_f_s, l5_beamform, l5_pll_disc_bandwidth, l5_dll_disc_bandwidth, l5_scale_factor
        1/3 * π, 50, 2.0, 1023e3, 4e6, beamform, 18.0, 1.0, l1_scale_factor,
        2/3 * π, 50, 1.0, 1023e4, 4e7, beamform, 18.0, 1.0, l5_scale_factor,
        1e-3, 1 # Δt, sat_prn
    )

    l1_test_signal = cis.(2 * π * 50 / 4e6 * (1:32000) + 1 / 3 * π)
    l1_samples_code = GNSSSignals.gen_sat_code(1:32000, 1.023e6, 2.0, 4e6, SATELLITE_1_CODE) * 10^(-20/20)
    l1_incoming_signals = [1, 1] .* (l1_test_signal[1:4000] .* l1_samples_code[1:4000])'

    l5_test_signal = cis.(2 * π * 50 / 4e7 * (1:102300) + 2 / 3 * π)
    l5_samples_code = GNSSSignals.gen_sat_code(1:102300, 1.023e7, 1.0, 4e7, L5_SAT1_CODE) * 10^(-20/20)
    l5_incoming_signals = [1, 1] .* (l5_test_signal[1:40000] .* l5_samples_code[1:40000])'

    next_joined_tracking_loop, l1_code_phase, l1_carrier_freq_update, l5_code_phase, l5_carrier_freq_update = joined_tracking_loop(l1_incoming_signals, l5_incoming_signals)
    
    @test l1_code_phase == 2.0
    @test l1_carrier_freq_update ≈ 0.0 atol = 1e-13
    
    @test l5_code_phase ≈ 1.0 atol = 1e-3
    @test l5_carrier_freq_update ≈ 0.0 atol = 1e-13
    

end
