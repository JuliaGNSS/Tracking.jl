@testset "Joined tracking" begin

    function beamform(x)
        [0.5 0.5 0.5 0.5] * x
    end

    l1_center_freq = 1.57542e9
    l5_center_freq = 1.17645e9
    l1_code_freq = 1023e3
    l5_code_freq = 1023e4
    l1_doppler = 10
    l5_doppler = l1_doppler * l5_center_freq / l1_center_freq
    l1_interm_freq = 50
    l5_interm_freq = 100
    l1_carrier_phase = π / 3
    l5_carrier_phase = π / 4
    l1_sample_freq = 4e6
    l5_sample_freq = 4e7
    l1_code_doppler = l1_doppler * l1_code_freq / l1_center_freq
    l5_code_doppler = l5_doppler * l5_code_freq / l5_center_freq
    l1_code_phase = 2.0
    l5_code_phase = 0.0

    run_time = 1
    integration_time = 1e-3
    num_integrations = Int(run_time / integration_time)
    l1_num_samples = Int(run_time * l1_sample_freq)
    l5_num_samples = Int(run_time * l5_sample_freq)
    l1_integration_samples = Int(integration_time * l1_sample_freq)
    l5_integration_samples = Int(integration_time * l5_sample_freq)
    
    l1_carrier = cis.(2 * π * (l1_interm_freq + l1_doppler) / l1_sample_freq * (1:l1_num_samples) + l1_carrier_phase)
    l5_carrier = cis.(2 * π * (l5_interm_freq + l5_doppler) / l5_sample_freq * (1:l5_num_samples) + l5_carrier_phase)
    l1_sampled_code = gen_code(1:l1_num_samples, l1_code_doppler + l1_code_freq, l1_code_phase, l1_sample_freq, SATELLITE_1_CODE)
    l5_sampled_code = gen_code(1:l5_num_samples, l5_code_doppler + l5_code_freq, l5_code_phase, l5_sample_freq, L5_SAT1_CODE)
    l1_signal = l1_carrier .* l1_sampled_code
    l5_signal = l5_carrier .* l5_sampled_code

    l1_results = Tracking.TrackingResults(0, l1_carrier_phase, 0.0, l1_code_phase, [0.0 + 0.0im])
    l5_results = Tracking.TrackingResults(0, l5_carrier_phase, 0.0, l5_code_phase, [0.0 + 0.0im])

    track = Tracking.init_joined_tracking(GPSL1(), GPSL5(), l1_results, l5_results, l1_sample_freq, l5_sample_freq, l1_interm_freq, l5_interm_freq, 10.0, 1.0, 1)
    
    l1_code_phases = zeros(num_integrations)
    l5_code_phases = zeros(num_integrations)
    l1_calculated_code_phases =  mod.((1:num_integrations) * l1_integration_samples * (l1_code_doppler + l1_code_freq) / l1_sample_freq + l1_code_phase, 1023)
    l5_calculated_code_phases = mod.((1:num_integrations) * l5_integration_samples * (l5_code_doppler + l5_code_freq) / l5_sample_freq + l5_code_phase, 10230)
    l1_carrier_dopplers = zeros(num_integrations)
    l5_carrier_dopplers = zeros(num_integrations)
    l1_code_dopplers = zeros(num_integrations)
    l5_code_dopplers = zeros(num_integrations)

    joined_results = nothing
    for i = 1:num_integrations
        l1_current_signal = [1, 1].' .* l1_signal[l1_integration_samples * (i - 1) + 1:l1_integration_samples * i]# .+ complex.(randn(l1_integration_samples, 2), randn(l1_integration_samples, 2)) .* 10^(15 / 20)
        l5_current_signal = [1, 1].' .* l5_signal[l5_integration_samples * (i - 1) + 1:l5_integration_samples * i]# .+ complex.(randn(l5_integration_samples, 2), randn(l5_integration_samples, 2)) .* 10^(15 / 20)
        track, joined_results = track(l1_current_signal, l5_current_signal, beamform)
        l1_code_phases[i] = joined_results.l1_results.code_phase
        l5_code_phases[i] = joined_results.l5_results.code_phase
        l1_carrier_dopplers[i] = joined_results.l1_results.carrier_doppler
        l5_carrier_dopplers[i] = joined_results.l5_results.carrier_doppler
        l1_code_dopplers[i] = joined_results.l1_results.code_doppler
        l5_code_dopplers[i] = joined_results.l5_results.code_doppler
    end

    @test joined_results.l1_results.carrier_doppler ≈ l1_doppler atol = 5e-2
    @test joined_results.l5_results.carrier_doppler ≈ l5_doppler atol = 5e-2
    @test joined_results.l1_results.code_phase ≈ l1_calculated_code_phases[end] atol = 5e-6
    @test joined_results.l5_results.code_phase ≈ l5_calculated_code_phases[end] atol = 5e-5
    @test joined_results.l1_results.code_doppler ≈ l1_code_doppler atol = 5e-5
    @test joined_results.l5_results.code_doppler ≈ l5_code_doppler atol = 5e-4

    #= figure("L1 code_phases(b) and calculated (r)")
    plot(l1_calculated_code_phases, color = "red")
    plot(l1_code_phases, color = "blue")
    figure("L5 code_phases")
    plot(l5_calculated_code_phases, color = "red")
    plot(mod.(l5_code_phases, 10230), color = "blue")
    figure("L1 carrier_dopplers")
    plot(l1_carrier_dopplers)
    figure("L5 carrier_dopplers")
    plot(l5_carrier_dopplers)
    figure("L1 code dopplers")
    plot(l1_code_dopplers)
    figure("L5 code dopplers")
    plot(l5_code_dopplers) =#
end