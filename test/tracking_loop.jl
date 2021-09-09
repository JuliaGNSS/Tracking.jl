@testset "Default post correlation Filter" begin
    post_corr_filter = Tracking.get_default_post_corr_filter(EarlyPromptLateCorrelator())
    @test @inferred(post_corr_filter(1.0 + 0.0im)) == 1.0 + 0.0im

    post_corr_filter = Tracking.get_default_post_corr_filter(
        EarlyPromptLateCorrelator(NumAnts(2))
    )
    @test @inferred(post_corr_filter(SVector(1.0 + 0.0im, 2.0 + 0.0im))) == 1.0 + 0.0im
end

@testset "Resize downconverted signal ($type) for multiple ants" for type in (Int16, Int32, Int64, Float32, Float64)
    downconverted_signal_temp = Tracking.DownconvertedSignalCPU(NumAnts(4))
    signal = StructArray(ones(Complex{type}, 2500, 4))
    downconverted_signal = Tracking.resize!(downconverted_signal_temp, 2500, signal)
    if type == Float64
        @test size(downconverted_signal.downconverted_signal_f32) == (0, 4)
        @test size(downconverted_signal.downconverted_signal_f64) == (2500, 4)
    else
        @test size(downconverted_signal.downconverted_signal_f32) == (2500, 4)
        @test size(downconverted_signal.downconverted_signal_f64) == (0, 4)
    end
end

@testset "Resize downconverted signal ($type) for single ant" for type in (Int16, Int32, Int64, Float32, Float64)
    downconverted_signal_temp = Tracking.DownconvertedSignalCPU(NumAnts(1))
    signal = StructArray(ones(Complex{type}, 2500))
    downconverted_signal = Tracking.resize!(downconverted_signal_temp, 2500, signal)
    if type == Float64
        @test size(downconverted_signal.downconverted_signal_f32) == (0,)
        @test size(downconverted_signal.downconverted_signal_f64) == (2500,)
    else
        @test size(downconverted_signal.downconverted_signal_f32) == (2500,)
        @test size(downconverted_signal.downconverted_signal_f64) == (0,)

    end
end

@testset "Integration time" begin
    gpsl1 = GPSL1()
    galileo_e1b = GalileoE1B()
    max_integration_time = 2ms
    bit_found = false
    time = @inferred Tracking.get_integration_time(gpsl1, max_integration_time, bit_found)
    @test time == 1ms

    max_integration_time = 2ms
    bit_found = true
    time = @inferred Tracking.get_integration_time(gpsl1, max_integration_time, bit_found)
    @test time == 2ms

    max_integration_time = 2ms
    bit_found = false
    time = @inferred Tracking.get_integration_time(
        galileo_e1b,
        max_integration_time,
        bit_found
    )
    @test time == 2ms

    max_integration_time = 10ms
    bit_found = false
    time = @inferred Tracking.get_integration_time(
        galileo_e1b,
        max_integration_time,
        bit_found
    )
    @test time == 4ms
end

@testset "Number of chips to integrate" begin
    gpsl1 = GPSL1()
    max_integration_time = 1ms
    code_phase = 200
    bit_found = true
    chips = @inferred Tracking.get_num_chips_to_integrate(
        gpsl1,
        max_integration_time,
        code_phase,
        bit_found
    )
    @test chips == 1023 - 200

    code_phase = 1200
    chips = @inferred Tracking.get_num_chips_to_integrate(
        gpsl1,
        max_integration_time,
        code_phase,
        bit_found
    )
    @test chips == 2046 - 1200
end

@testset "Number of samples to integrate" begin
    gpsl1 = GPSL1()
    max_integration_time = 1ms
    sampling_frequency = 4e6Hz
    current_code_doppler = 0Hz
    current_code_phase = 0
    bit_found = true
    samples = @inferred Tracking.get_num_samples_left_to_integrate(
        gpsl1,
        max_integration_time,
        sampling_frequency,
        current_code_doppler,
        current_code_phase,
        bit_found
    )
    @test samples == 4000
end

@testset "Carrier phase delta" begin
    intermediate_frequency = 100Hz
    carrier_doppler = 50Hz
    frequency = @inferred Tracking.get_current_carrier_frequency(
        intermediate_frequency,
        carrier_doppler
    )
    @test frequency == intermediate_frequency + carrier_doppler
end

@testset "Code phase delta" begin
    gpsl1 = GPSL1()
    code_doppler = 1Hz
    frequency = @inferred Tracking.get_current_code_frequency(gpsl1, code_doppler)
    @test frequency == 1023e3Hz + 1Hz
end

@testset "Doppler aiding" begin
    gpsl1 = GPSL1()
    init_carrier_doppler = 10Hz
    init_code_doppler = 1Hz
    carrier_freq_update = 2Hz
    code_freq_update = -0.5Hz
    velocity_aiding = 3Hz

    carrier_freq, code_freq = @inferred Tracking.aid_dopplers(
        gpsl1,
        init_carrier_doppler,
        init_code_doppler,
        carrier_freq_update,
        code_freq_update,
        velocity_aiding
    )

    @test carrier_freq == 10Hz + 2Hz + 3Hz
    @test code_freq == 1Hz + (2Hz + 3Hz) / 1540 - 0.5Hz
end

@testset "Tracking with signal of type $type" for type in (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + 1023e3Hz
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    state = TrackingState(gpsl1, carrier_doppler - 20Hz, start_code_phase)

    signal_temp = cis.(
            2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase
        ) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn
        )
    scaling = 512
    signal = type <: Integer ? complex.(floor.(type, real.(signal_temp) * scaling), floor.(type, imag.(signal_temp) * scaling)) : Complex{type}.(signal_temp)

    track_result = @inferred track(signal, state, prn, sampling_frequency)

    iterations = 2000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    tracked_prompts = zeros(ComplexF64, iterations)
    for i = 1:iterations
        carrier_phase = mod2pi(2π * carrier_doppler * 4000 * i / sampling_frequency +
            start_carrier_phase + π) - π
        code_phase = mod(
            code_frequency * 4000 * i / sampling_frequency + start_code_phase,
            1023
        )
        signal = cis.(
                2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase
            ) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn
            )
        track_result = @inferred track(
            signal,
            get_state(track_result),
            prn,
            sampling_frequency
        )
        comp_carrier_phase = mod2pi(2π * carrier_doppler * 4000 * (i + 1) /
            sampling_frequency + start_carrier_phase + π) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023
        )
        tracked_code_phases[i] = get_code_phase(track_result)
        tracked_carrier_phases[i] = get_carrier_phase(track_result)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_result)/Hz
        tracked_code_dopplers[i] = get_code_doppler(track_result)/Hz
        tracked_prompts[i] = get_prompt(track_result)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 2e-4
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 5e-2

#    using PyPlot
#    pygui(true)
#    figure("carrier_phases")
#    plot(tracked_carrier_phases)
#    plot(carrier_phases)
#    grid(true)
#    figure("Code Phases")
#    plot(300 * (tracked_code_phases .- code_phases))
#    figure("Carrier Doppler")
#    plot(tracked_carrier_dopplers)
#    figure("Code Doppler")
#    plot(tracked_code_dopplers)
#    figure("Prompt")
#    plot(real.(tracked_prompts))
#    plot(imag.(tracked_prompts))
end

@testset "Track multiple signals with signal of type $type" for type in (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + 1023e3Hz
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    state = TrackingState(gpsl1, carrier_doppler - 20Hz, start_code_phase, num_ants = NumAnts(3))

    signal = cis.(
            2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase
        ) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn
        )
    signal_mat_temp = repeat(signal, outer = (1,3))
    scaling = 512
    signal_mat = type <: Integer ? complex.(floor.(type, real.(signal_mat_temp) * scaling), floor.(type, imag.(signal_mat_temp) * scaling)) : Complex{type}.(signal_mat_temp)

    @test_throws ArgumentError track(signal_mat', state, prn, sampling_frequency)

    track_result = @inferred track(signal_mat, state, prn, sampling_frequency)

    iterations = 2000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    for i = 1:iterations
        carrier_phase = mod2pi(2π * carrier_doppler * 4000 * i / sampling_frequency +
            start_carrier_phase + π) - π
        code_phase = mod(
            code_frequency * 4000 * i / sampling_frequency + start_code_phase,
            1023
        )
        signal = cis.(
                2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase
            ) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn
            )
        signal_mat = repeat(signal, outer = (1,3))
        track_result = @inferred track(
            signal_mat,
            get_state(track_result),
            prn,
            sampling_frequency
        )
        comp_carrier_phase = mod2pi(2π * carrier_doppler * 4000 * (i + 1) /
            sampling_frequency + start_carrier_phase + π) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023
        )
        tracked_code_phases[i] = get_code_phase(track_result)
        tracked_carrier_phases[i] = get_carrier_phase(track_result)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_result)/Hz
        tracked_code_dopplers[i] = get_code_doppler(track_result)/Hz
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 2e-4
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 5e-2

#    using PyPlot
#    pygui(true)
#    figure("carrier_phases")
#    plot(tracked_carrier_phases)
#    plot(carrier_phases)
#    grid(true)
#    figure("Code Phases")
#    plot(300 * (tracked_code_phases .- code_phases))
#    figure("Carrier Doppler")
#    plot(tracked_carrier_dopplers)
#    figure("Code Doppler")
#    plot(tracked_code_dopplers)
end
