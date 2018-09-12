@testset "Correlate" begin
    signal = [1,1]' .* gen_code(GPSL1(), 1:4000, 1023e3, 10, 4e6, 1)
    replica = gen_code(GPSL1(), 1:4000, 1023e3, 10, 4e6, 1)
    @test @inferred(Tracking.correlate(signal, replica)) == [1.0 1.0]
end

@testset "Downconvert" begin
    signal = cis.(2π * 1e6 / 4e6 * (1:10))
    replica = cis.(2π * 1e6 / 4e6 * (1:10))
    downconverted_signals = @inferred Tracking.downconvert(signal, replica)
    @test downconverted_signals ≈ ones(10)
end

@testset "Track L1" begin

    function beamform(x)
        [0.5 0.5 0.5 0.5] * x
    end

    center_freq = 1.57542e9Hz
    code_freq = 1023e3Hz
    doppler = 10Hz
    interm_freq = 50Hz
    carrier_phase = π / 3
    sample_freq = 4e6Hz
    code_doppler = doppler * code_freq / center_freq
    code_phase = 2.0

    run_time = 500e-3s
    integration_time = 1e-3s
    num_integrations = Int(run_time / integration_time)
    num_samples = Int(run_time * sample_freq)
    integration_samples = Int(integration_time * sample_freq)

    gps_l1 = GPSL1()

    carrier = cis.(2 * π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
    sampled_code = gen_code(gps_l1, 1:num_samples, code_doppler + code_freq, code_phase, sample_freq, 1)
    signal = carrier .* sampled_code

    inits = Initials(0.0Hz, carrier_phase, 0.0Hz, code_phase)
    track = init_tracking(gps_l1, inits, interm_freq, sample_freq, 18.0Hz, 1.0Hz, 1)

    code_dopplers = zeros(num_integrations)
    code_phases = zeros(num_integrations)
    calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 1023)
    carrier_dopplers = zeros(num_integrations)

    results = nothing
    for i = 1:num_integrations
        current_signal = [1, 1, 1, 1]' .* signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,4), randn(integration_samples,4 )) .* 10^(5/20)
        track, results = track(current_signal, beamform)
        code_phases[i] = results.code_phase
        carrier_dopplers[i] = results.carrier_doppler / 1.0Hz
        code_dopplers[i] = results.code_doppler / 1.0Hz
    end

    @test results.carrier_doppler ≈ doppler atol = 5e-2Hz
    @test results.code_phase ≈ calculated_code_phases[end] atol = 5e-6
    @test results.code_doppler ≈ code_doppler atol = 5e-5Hz

    #figure("Tracking code phases")
    #plot(code_phases, color = "blue")
    #plot(calculated_code_phases, color = "red")
    #figure("Tracking carrier_dopplers")
    #plot(carrier_dopplers)
    #figure("Tracking code dopplers")
    #plot(code_dopplers)

end
