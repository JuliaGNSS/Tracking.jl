@testset "Tracking with signal of type $type" for type in
                                                  (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    num_samples = 4000

    track_state = @inferred TrackState(
        gpsl1,
        [SatState(gpsl1, 1, sampling_frequency, start_code_phase, carrier_doppler - 20Hz)];
        num_samples,
        maximum_expected_sampling_frequency = Val(sampling_frequency),
    )

    signal_temp =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, start_code_phase)
    scaling = 512
    signal =
        type <: Integer ?
        complex.(
            round.(type, real.(signal_temp) * scaling),
            round.(type, imag.(signal_temp) * scaling),
        ) : Complex{type}.(signal_temp)

    track_state = @inferred track(signal, track_state, sampling_frequency)

    iterations = 20000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    tracked_prompts = zeros(ComplexF64, iterations)
    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase =
            mod(code_frequency * 4000 * i / sampling_frequency + start_code_phase, 1023)
        signal_temp =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, code_phase)
        signal =
            type <: Integer ?
            complex.(
                round.(type, real.(signal_temp) * scaling),
                round.(type, imag.(signal_temp) * scaling),
            ) : Complex{type}.(signal_temp)
        track_state = @inferred track(signal, track_state, sampling_frequency)
        comp_carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * (i + 1) / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023,
        )
        tracked_code_phases[i] = get_code_phase(track_state)
        tracked_carrier_phases[i] = get_carrier_phase(track_state)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_state) / Hz
        tracked_code_dopplers[i] = get_code_doppler(track_state) / Hz
        @test get_signal_start_sample(track_state) == 4001
        tracked_prompts[i] = get_last_fully_integrated_filtered_prompt(track_state)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 5e-5
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 1e-3

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

@testset "Track multiple systems of type $type" for type in
                                                    (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1()
    galileo_e1b = GalileoE1B()
    carrier_doppler_gps = 200Hz
    carrier_doppler_gal = 1200Hz
    start_code_phase = 100
    code_frequency_gps =
        carrier_doppler_gps * get_code_center_frequency_ratio(gpsl1) +
        get_code_frequency(gpsl1)
    code_frequency_gal =
        carrier_doppler_gal * get_code_center_frequency_ratio(galileo_e1b) +
        get_code_frequency(galileo_e1b)
    sampling_frequency = 15e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    num_samples = 4000

    track_state = @inferred TrackState(
        (
            gps = SystemSatsState(
                gpsl1,
                [
                    SatState(
                        gpsl1,
                        prn,
                        sampling_frequency,
                        start_code_phase,
                        carrier_doppler_gps,
                    ),
                ],
            ),
            gal = SystemSatsState(
                galileo_e1b,
                [
                    SatState(
                        galileo_e1b,
                        prn,
                        sampling_frequency,
                        start_code_phase,
                        carrier_doppler_gal,
                    ),
                ],
            ),
        );
        num_samples,
        maximum_expected_sampling_frequency = Val(sampling_frequency),
    )

    signal_temp =
        cis.(
            2π .* carrier_doppler_gps .* range ./ sampling_frequency .+ start_carrier_phase
        ) .* gen_code(
            4000,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency_gps,
            start_code_phase,
        ) .+
        cis.(
            2π .* carrier_doppler_gal .* range ./ sampling_frequency .+ start_carrier_phase
        ) .* gen_code(
            4000,
            galileo_e1b,
            prn,
            sampling_frequency,
            code_frequency_gal,
            start_code_phase,
        )
    scaling = 5000
    signal =
        type <: Integer ?
        complex.(
            round.(type, real.(signal_temp) * scaling),
            round.(type, imag.(signal_temp) * scaling),
        ) : Complex{type}.(signal_temp)

    track_state = @inferred track(signal, track_state, sampling_frequency)

    iterations = 2000
    for i = 1:iterations
        carrier_phase_gps =
            mod2pi(
                2π * carrier_doppler_gps * 4000 * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase_gps = mod(
            code_frequency_gps * 4000 * i / sampling_frequency + start_code_phase,
            get_code_length(gpsl1),
        )
        carrier_phase_gal =
            mod2pi(
                2π * carrier_doppler_gal * 4000 * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase_gal = mod(
            code_frequency_gal * 4000 * i / sampling_frequency + start_code_phase,
            get_code_length(galileo_e1b),
        )
        signal_temp =
            cis.(
                2π .* carrier_doppler_gps .* range ./ sampling_frequency .+
                carrier_phase_gps
            ) .* gen_code(
                4000,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency_gps,
                code_phase_gps,
            ) .+
            cis.(
                2π .* carrier_doppler_gal .* range ./ sampling_frequency .+
                carrier_phase_gal
            ) .* gen_code(
                4000,
                galileo_e1b,
                prn,
                sampling_frequency,
                code_frequency_gal,
                code_phase_gal,
            )
        signal =
            type <: Integer ?
            complex.(
                round.(type, real.(signal_temp) * scaling),
                round.(type, imag.(signal_temp) * scaling),
            ) : Complex{type}.(signal_temp)
        track_state = @inferred track(signal, track_state, sampling_frequency)
    end
    comp_carrier_phase_gps =
        mod2pi(
            2π * carrier_doppler_gps * 4000 * (iterations + 1) / sampling_frequency +
            start_carrier_phase +
            π,
        ) - π
    comp_code_phase_gps = mod(
        code_frequency_gps * 4000 * (iterations + 1) / sampling_frequency +
        start_code_phase,
        get_code_length(gpsl1),
    )
    comp_carrier_phase_gal =
        mod2pi(
            2π * carrier_doppler_gal * 4000 * (iterations + 1) / sampling_frequency +
            start_carrier_phase +
            π,
        ) - π
    comp_code_phase_gal = mod(
        code_frequency_gal * 4000 * (iterations + 1) / sampling_frequency +
        start_code_phase,
        get_code_length(galileo_e1b),
    )
    @test get_code_phase(track_state, :gps, 1) ≈ comp_code_phase_gps atol = 5e-3
    @test mod(get_carrier_phase(track_state, :gps, 1), π) ≈ mod(comp_carrier_phase_gps, π) atol =
        2e-2
    @test get_code_phase(track_state, :gal, 1) ≈ comp_code_phase_gal atol = 5e-3
    @test mod(get_carrier_phase(track_state, :gal, 1), π) ≈ mod(comp_carrier_phase_gal, π) atol =
        3e-3
end

@testset "Tracking with intermediate frequency of $intermediate_frequency" for intermediate_frequency in
                                                                               (
    0.0Hz,
    -10000.0Hz,
    10000.0Hz,
    -30000.0Hz,
    30000.0Hz,
)
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    num_samples = 4000

    track_state = @inferred TrackState(
        gpsl1,
        [SatState(gpsl1, 1, sampling_frequency, start_code_phase, carrier_doppler - 20Hz)];
        num_samples,
        maximum_expected_sampling_frequency = Val(sampling_frequency),
    )

    signal =
        cis.(
            2π .* (carrier_doppler + intermediate_frequency) .* range ./
            sampling_frequency .+ start_carrier_phase
        ) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn,
        )

    track_state =
        @inferred track(signal, track_state, sampling_frequency; intermediate_frequency)

    iterations = 2000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    tracked_prompts = zeros(ComplexF64, iterations)
    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * (carrier_doppler + intermediate_frequency) * 4000 * i /
                sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase =
            mod(code_frequency * 4000 * i / sampling_frequency + start_code_phase, 1023)
        signal =
            cis.(
                2π .* (carrier_doppler + intermediate_frequency) .* range ./
                sampling_frequency .+ carrier_phase
            ) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn,
            )
        track_state =
            @inferred track(signal, track_state, sampling_frequency; intermediate_frequency)
        comp_carrier_phase =
            mod2pi(
                2π * (carrier_doppler + intermediate_frequency) * 4000 * (i + 1) /
                sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023,
        )
        tracked_code_phases[i] = get_code_phase(track_state)
        tracked_carrier_phases[i] = get_carrier_phase(track_state)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_state) / Hz
        tracked_code_dopplers[i] = get_code_doppler(track_state) / Hz
        tracked_prompts[i] = get_last_fully_integrated_filtered_prompt(track_state)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 5e-5
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 5e-5

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

@testset "Track multiple signals with signal of type $type" for type in (
    Int16,
    Int32,
    Int64,
    Float32,
    Float64,
)
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    num_samples = 4000

    track_state = @inferred TrackState(
        gpsl1,
        [
            SatState(
                gpsl1,
                1,
                sampling_frequency,
                start_code_phase,
                carrier_doppler - 20Hz;
                num_ants = NumAnts(3),
            ),
        ];
        num_samples,
        maximum_expected_sampling_frequency = Val(sampling_frequency),
    )

    @test get_num_ants(track_state) == 3

    signal =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn,
        )
    signal_mat_temp = repeat(signal; outer = (1, 3))
    scaling = 512
    signal_mat =
        type <: Integer ?
        complex.(
            floor.(type, real.(signal_mat_temp) * scaling),
            floor.(type, imag.(signal_mat_temp) * scaling),
        ) : Complex{type}.(signal_mat_temp)

    track_state = @inferred track(signal_mat, track_state, sampling_frequency)

    iterations = 2000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    tracked_prompts = zeros(ComplexF64, iterations)
    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase =
            mod(code_frequency * 4000 * i / sampling_frequency + start_code_phase, 1023)
        signal =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn,
            )
        signal_mat = repeat(signal; outer = (1, 3))
        track_state = @inferred track(signal_mat, track_state, sampling_frequency)
        comp_carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * (i + 1) / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023,
        )
        tracked_code_phases[i] = get_code_phase(track_state)
        tracked_carrier_phases[i] = get_carrier_phase(track_state)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_state) / Hz
        tracked_code_dopplers[i] = get_code_doppler(track_state) / Hz
        tracked_prompts[i] = get_last_fully_integrated_filtered_prompt(track_state)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 5e-5
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 5e-5

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

@testset "Track multiple signals with GPU" begin
    !CUDA.functional() && return
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    num_samples = 4000
    num_ants = NumAnts(3)

    correlator = get_default_correlator(gpsl1, sampling_frequency, num_ants)
    sat_states = [
        SatState(
            gpsl1,
            1,
            sampling_frequency,
            start_code_phase,
            carrier_doppler - 20Hz;
            num_ants,
            correlator,
        ),
    ]

    system_sats_state = SystemSatsState(gpsl1, sat_states)

    track_state = @inferred TrackState(
        system_sats_state;
        num_samples,
        maximum_expected_sampling_frequency = Val(sampling_frequency),
        downconvert_and_correlator = GPUDownconvertAndCorrelator(
            (system_sats_state,),
            num_samples,
        ),
    )

    signal =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn,
        )
    signal_mat = cu(repeat(signal; outer = (1, 3)))

    track_state = @inferred track(signal_mat, track_state, sampling_frequency)

    iterations = 2000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    tracked_prompts = zeros(ComplexF64, iterations)
    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase =
            mod(code_frequency * 4000 * i / sampling_frequency + start_code_phase, 1023)
        signal =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn,
            )
        signal_mat = cu(repeat(signal; outer = (1, 3)))
        track_state = @inferred track(signal_mat, track_state, sampling_frequency)
        comp_carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * (i + 1) / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023,
        )
        tracked_code_phases[i] = get_code_phase(track_state)
        tracked_carrier_phases[i] = get_carrier_phase(track_state)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_state) / Hz
        tracked_code_dopplers[i] = get_code_doppler(track_state) / Hz
        tracked_prompts[i] = get_last_fully_integrated_filtered_prompt(track_state)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 1e-2
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 5e-5

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
