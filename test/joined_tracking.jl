@testset "Joined tracking" begin

    function beamform(x)
        dot([0.5, 0.5, 0.5, 0.5], x)
    end

    l1_center_freq = 1.57542e9Hz
    l5_center_freq = 1.17645e9Hz
    l1_code_freq = 1023e3Hz
    l5_code_freq = 1023e4Hz
    l1_doppler = 10Hz
    l5_doppler = l1_doppler * l5_center_freq / l1_center_freq
    l1_interm_freq = 50Hz
    l5_interm_freq = 100Hz
    l1_carrier_phase = π / 3
    l5_carrier_phase = π / 4
    l1_sample_freq = 4e6Hz
    l5_sample_freq = 4e7Hz
    l1_code_doppler = l1_doppler * l1_code_freq / l1_center_freq
    l5_code_doppler = l5_doppler * l5_code_freq / l5_center_freq
    l1_code_phase = 2.0
    l5_code_phase = 4.0

    run_time = 1s
    integration_time = 1e-3s
    num_integrations = Int(run_time / integration_time)
    l1_num_samples = convert(Int, run_time * l1_sample_freq)
    l5_num_samples = convert(Int, run_time * l5_sample_freq)
    l1_integration_samples = convert(Int, integration_time * l1_sample_freq)
    l5_integration_samples = convert(Int, integration_time * l5_sample_freq)

    l1_carrier = cis.(2 * π * (l1_interm_freq + l1_doppler) / l1_sample_freq * (1:l1_num_samples) .+ l1_carrier_phase)
    l5_carrier = cis.(2 * π * (l5_interm_freq + l5_doppler) / l5_sample_freq * (1:l5_num_samples) .+ l5_carrier_phase)
    l1_sampled_code = gen_code.(Ref(GPSL1()), 1:l1_num_samples, l1_code_doppler + l1_code_freq, l1_code_phase, l1_sample_freq, 1)
    l5_sampled_code = gen_code.(Ref(GPSL5()), 1:l5_num_samples, l5_code_doppler + l5_code_freq, l5_code_phase, l5_sample_freq, 1)
    l1_signal = l1_carrier .* l1_sampled_code
    l5_signal = l5_carrier .* l5_sampled_code

    l1_inits = Tracking.Initials(0Hz, l1_carrier_phase, 0.0Hz, l1_code_phase)
    l5_inits = Tracking.Initials(0Hz, l5_carrier_phase, 0.0Hz, l5_code_phase)

    track = Tracking.init_tracking([GPSL1(), GPSL5()], [l1_inits, l5_inits], integration_time, [l1_sample_freq, l5_sample_freq], [l1_interm_freq, l5_interm_freq], 10.0Hz, 1.0Hz, 1)

    l1_code_phases = zeros(num_integrations)
    l5_code_phases = zeros(num_integrations)
    # Float64 because of error in Unitful: https://github.com/ajkeller34/Unitful.jl/issues/160
    l1_calculated_code_phases =  mod.((1:num_integrations) * l1_integration_samples * Float64((l1_code_doppler + l1_code_freq) / l1_sample_freq) .+ l1_code_phase, 1023)
    l5_calculated_code_phases = mod.((1:num_integrations) * l5_integration_samples * Float64((l5_code_doppler + l5_code_freq) / l5_sample_freq) .+ l5_code_phase, 10230)
    l1_carrier_dopplers = zeros(num_integrations)
    l5_carrier_dopplers = zeros(num_integrations)
    l1_code_dopplers = zeros(num_integrations)
    l5_code_dopplers = zeros(num_integrations)

    joined_results = nothing
    for i = 1:num_integrations
        l1_current_signal = [1, 1]' .* l1_signal[l1_integration_samples * (i - 1) + 1:l1_integration_samples * i]# .+ complex.(randn(l1_integration_samples, 2), randn(l1_integration_samples, 2)) .* 10^(15 / 20)
        l5_current_signal = [1, 1]' .* l5_signal[l5_integration_samples * (i - 1) + 1:l5_integration_samples * i]# .+ complex.(randn(l5_integration_samples, 2), randn(l5_integration_samples, 2)) .* 10^(15 / 20)
        track, joined_results = track([l1_current_signal, l5_current_signal], beamform)
        l1_code_phases[i] = joined_results[1].code_phase[1]
        l5_code_phases[i] = joined_results[2].code_phase[1]
        l1_carrier_dopplers[i] = joined_results[1].carrier_doppler[1] / 1.0Hz
        l5_carrier_dopplers[i] = joined_results[2].carrier_doppler[1] / 1.0Hz
        l1_code_dopplers[i] = joined_results[1].code_doppler[1] / 1.0Hz
        l5_code_dopplers[i] = joined_results[2].code_doppler[1] / 1.0Hz
    end

    @test joined_results[1].carrier_doppler[1] ≈ l1_doppler atol = 5e-2Hz
    @test joined_results[2].carrier_doppler[1] ≈ l5_doppler atol = 5e-2Hz
    @test joined_results[1].code_phase[1] ≈ l1_calculated_code_phases[end] atol = 5e-6
    @test joined_results[2].code_phase[1] ≈ l5_calculated_code_phases[end] atol = 5e-5
    @test joined_results[1].code_doppler[1] ≈ l1_code_doppler atol = 5e-5Hz
    @test joined_results[2].code_doppler[1] ≈ l5_code_doppler atol = 5e-4Hz

    #=
     figure("L1 code_phases(b) and calculated (r)")
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
    plot(l5_code_dopplers)
    =#
end
