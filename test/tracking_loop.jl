@testset "Track L1" begin

    function beamform(x)
        dot([0.5, 0.5, 0.5, 0.5], x)
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
    num_integrations = convert(Int, run_time / integration_time)
    num_samples = convert(Int, run_time * sample_freq)
    integration_samples = convert(Int, integration_time * sample_freq)

    gps_l1 = GPSL1()

    carrier = cis.(2 * π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
    sampled_code = gen_code.(Ref(gps_l1), 1:num_samples, code_doppler + code_freq, code_phase, sample_freq, 1)
    signal = carrier .* sampled_code

    inits = Initials(0.0Hz, carrier_phase, 0.0Hz, code_phase)
    track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)

    code_dopplers = zeros(num_integrations)
    code_phases = zeros(num_integrations)
    calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 1023)
    carrier_dopplers = zeros(num_integrations)

    results = nothing
    for i = 1:num_integrations
        current_signal = [1, 1, 1, 1]' .* signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,4), randn(integration_samples,4 )) .* 10^(5/20)
        track, results = track(current_signal, beamform)
        code_phases[i] = results.code_phase[1]
        carrier_dopplers[i] = results.carrier_doppler[1] / Hz
        code_dopplers[i] = results.code_doppler[1] / Hz
    end

    @test results.carrier_doppler[1] ≈ doppler atol = 5e-2Hz
    @test mod(results.code_phase[1], 1023) ≈ calculated_code_phases[end] atol = 5e-6
    @test results.code_doppler[1] ≈ code_doppler atol = 5e-5Hz

    #=
    figure("Tracking code phases")
    plot(code_phases, color = "blue")
    plot(calculated_code_phases, color = "red")
    figure("Tracking carrier_dopplers")
    plot(carrier_dopplers)
    figure("Tracking code dopplers")
    plot(code_dopplers)
    =#

    # Track all at once
    track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
    track, results = track([1, 1, 1, 1]' .* signal, beamform)

    @test code_phases == results.code_phase
    @test carrier_dopplers == results.carrier_doppler ./ Hz
    @test code_dopplers == results.code_doppler ./ Hz

    track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
    track, results = track([1, 1, 1, 1]' .* signal[1:40000], beamform)

    code_dopplers2 = CircularBuffer{Float64}(num_integrations)
    code_phases2 = CircularBuffer{Float64}(num_integrations)
    carrier_dopplers2 = CircularBuffer{Float64}(num_integrations)
    track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)

    # Track with time smaller than integration_time
    for i = 1:2000:size(signal, 1)
        track, results = track([1, 1, 1, 1]' .* signal[i:i + 2000 - 1], beamform)
        append!(code_dopplers2, results.code_doppler ./ Hz)
        append!(code_phases2, results.code_phase)
        append!(carrier_dopplers2, results.carrier_doppler ./ Hz)
    end

    @test code_phases ≈ code_phases2
    @test carrier_dopplers ≈ carrier_dopplers2
    @test code_dopplers ≈ code_dopplers2
end
