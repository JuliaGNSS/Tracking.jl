#For now
function beamform(x)
    [0.5 0.5 0.5 0.5] * x
end

@testset "Joined TrackingLoop" begin
    doppler_l1 = 5
    doppler_l5 = 5 * 1176.45 / 1575.42
    
    l1_carrier = cis.(2 * π * (50 + doppler_l1) / 4e6 * (1:5000000) + 1 / 3 * π)
    l5_carrier = cis.(2 * π * (100 + doppler_l5) / 4e7 * (1:50000000) + 1 / 4 * π)
    l1_sampled_code = GNSSSignals.gen_sat_code(1:5000000, 1 / 1540 * doppler_l1 + 1.023e6, 2.0, 4e6, SATELLITE_1_CODE)
    l5_sampled_code = GNSSSignals.gen_sat_code(1:50000000, 1 / 1150 * doppler_l5 +  1.023e7, 0.0, 4e7, L5_SAT1_CODE)

    l1_incoming_signals = [1, 1] .* (l1_carrier[1:4000] .* l1_sampled_code[1:4000])'
    l5_incoming_signals = [1, 1] .* (l5_carrier[1:40000] .* l5_sampled_code[1:40000])'

    l1_gen_sampled_code, l1_get_code_phase = GNSSSignals.init_gpsl1_codes()
    l5_gen_sampled_code, l5_get_code_phase = GNSSSignals.init_gpsl5_codes()
    l1_system = Tracking.GNSSSystem(1575.42e6, 1023e3, 4e6, l1_gen_sampled_code, l1_get_code_phase)
    l5_system = Tracking.GNSSSystem(1176.45e6, 1023e4, 4e7, l5_gen_sampled_code, l5_get_code_phase)
    l1_results = Tracking.TrackingResults(0, 1 / 3 * π, 0.0, 2.0, [0.0 + 0.0 * im] )
    l5_results = Tracking.TrackingResults(0, 1 / 4 * π, 0.0, 0.0, [0.0 + 0.0 * im] )

    joined_tracking = Tracking.init_joined_tracking(l1_system, l5_system, l1_results, l5_results, 50, 100, 18.0, 1.0, 1)
    next_joined_tracking, joined_results = joined_tracking(l1_incoming_signals, l5_incoming_signals, beamform)
    @ test joined_results.l1_results.code_doppler == 0
    @ test joined_results.l5_results.code_doppler == 0
    @ test joined_results.l1_results.carrier_doppler == 0
    @ test joined_results.l5_results.carrier_doppler == 0
    @ test joined_results.l1_results.code_phase == l1_get_code_phase(4000, 1023e3, 2.0, 4e6)
    @ test joined_results.l5_results.code_phase == l5_get_code_phase(40000, 1023e4, 0.0, 4e7)
    
    l1_code_phases = zeros(1000)
    l5_code_phases = zeros(1000)
    l1_calculated_code_phases = zeros(1000)
    l5_calculated_code_phases = zeros(1000)
    l1_carrier_dopplers = zeros(1000)
    l5_carrier_dopplers = zeros(1000)
    l1_code_dopplers = zeros(1000)
    l5_code_dopplers = zeros(1000)
    for i = 1:300
        l1_incoming_signals = [1, 1] .* (l1_carrier[(1 + 4000 * i):(4000 * (i+1))] .* l1_sampled_code[(1 + 4000 * i):(4000 * (i+1))])'
        l5_incoming_signals = [1, 1] .* (l5_carrier[(1 + 40000 * i):(40000 * (i+1))] .* l5_sampled_code[(1 + 40000 * i):(40000 * (i+1))])'
        next_joined_tracking, joined_results = next_joined_tracking(l1_incoming_signals, l5_incoming_signals, beamform)
        l1_code_phases[i] = joined_results.l1_results.code_phase
        l5_code_phases[i] = joined_results.l5_results.code_phase
        l1_calculated_code_phases[i] = l1_get_code_phase(4000 * (i + 2), 1 / 1540 * doppler_l1 + 1023e3, 2.0, 4e6)
        l5_calculated_code_phases[i] = l5_get_code_phase(40000 * (i + 2), 1 / 1150 * doppler_l5 + 1023e4, 0.0, 4e7)
        l1_carrier_dopplers[i] = joined_results.l1_results.carrier_doppler
        l5_carrier_dopplers[i] = joined_results.l5_results.carrier_doppler
        l1_code_dopplers[i] = joined_results.l1_results.code_doppler
        l5_code_dopplers[i] = joined_results.l5_results.code_doppler

    end
    for i = 301:1000
        l1_incoming_signals = [1, 1] .* (l1_carrier[(1 + 4000 * i):(4000 * (i+1))] .* l1_sampled_code[(1 + 4000 * i):(4000 * (i+1))])'
        l5_incoming_signals = [1, 1] .* (l5_carrier[(1 + 40000 * i):(40000 * (i+1))] .* l5_sampled_code[(1 + 40000 * i):(40000 * (i+1))])'
        next_joined_tracking, joined_results = next_joined_tracking(l1_incoming_signals, l5_incoming_signals, beamform)
        

        #println("l1 ", joined_results.l1_results.code_phase, " ", l1_get_code_phase(4000 * (i + 2), 1 / 1540 * doppler_l1 + 1023e3, 2.0, 4e6))
        #println("l5 ", joined_results.l5_results.code_phase, " ", l5_get_code_phase(40000 * (i + 2), 1 / 1150 * doppler_l5 + 1023e4, 0.0, 4e7))
        #println("doppler l1 ", joined_results.l1_results.carrier_doppler, " ", doppler_l1)
        #println("doppler l5 ", joined_results.l5_results.carrier_doppler, " ", doppler_l5)
        
        #@test joined_results.l1_results.carrier_doppler ≈ 0 atol = 1e-10
        #@test joined_results.l5_results.carrier_doppler ≈ 0 atol = 1e-10
        #@test results.carrier_doppler + init_carrier_freq ≈ 60.0 atol = 0.2  
        l1_code_phases[i] = joined_results.l1_results.code_phase
        l5_code_phases[i] = joined_results.l5_results.code_phase
        l1_calculated_code_phases[i] = l1_get_code_phase(4000 * (i + 2), 1 / 1540 * doppler_l1 + 1023e3, 2.0, 4e6)
        l5_calculated_code_phases[i] = l5_get_code_phase(40000 * (i + 2), 1 / 1150 * doppler_l5 + 1023e4, 0.0, 4e7)
        l1_carrier_dopplers[i] = joined_results.l1_results.carrier_doppler
        l5_carrier_dopplers[i] = joined_results.l5_results.carrier_doppler
        l1_code_dopplers[i] = joined_results.l1_results.code_doppler
        l5_code_dopplers[i] = joined_results.l5_results.code_doppler

    end
    println(joined_results.l5_results)

    figure("L1 code_phases(b) and calculated (r)")
    plot(l1_code_phases, color = "blue")
    plot(l1_calculated_code_phases, color = "red")
    figure("L5 code_phases")
    plot(mod.(l5_calculated_code_phases, 10230), color = "red")
    plot(mod.(l5_code_phases, 10230), color = "blue")
    figure("L1 carrier_dopplers")
    plot(l1_carrier_dopplers)
    figure("L5 carrier_dopplers")
    plot(l5_carrier_dopplers)
    figure("L1 code dopplers")
    plot(l1_code_dopplers)
    figure("L5 code dopplers")
    plot(l5_code_dopplers)
end