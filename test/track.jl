module TrackTest

using Test: @test, @testset, @inferred
using Random: Random
using Unitful: Hz, kHz, MHz, dBHz
using Statistics: mean
using Dictionaries: dictionary
using GNSSSignals:
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    GPSL5I,
    GalileoE1B,
    gen_code,
    get_code_center_frequency_ratio,
    get_code_frequency,
    get_code

using Tracking:
    TrackedSat,
    TrackState,
    track,
    track!,
    add_satellite!,
    get_code_phase,
    get_carrier_phase,
    get_code_doppler,
    get_carrier_doppler,
    get_prompt,
    estimate_cn0,
    get_signal_start_sample,
    get_last_fully_integrated_filtered_prompt,
    get_last_fully_integrated_correlator,
    get_filtered_prompts,
    get_integrated_samples,
    has_bit_or_secondary_code_been_found,
    get_sat_state,
    get_code_length,
    NumAnts,
    get_num_ants,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL

@testset "Tracking with signal of type $type" for type in
                                                  (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    track_state = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)],
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

@testset "Tracking with large initial Doppler offset" begin
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999 # tracking 1 ms
    start_carrier_phase = π / 2
    iterations = 1000
    tolerance = 1.0  # Hz, tight tolerance since no noise

    function test_convergence(carrier_error, use_assisted)
        doppler_estimator =
            use_assisted ? ConventionalAssistedPLLAndDLL() : ConventionalPLLAndDLL()
        sat_states = [
            TrackedSat(
                gpsl1,
                1,
                start_code_phase,
                carrier_doppler + carrier_error;
                doppler_estimator,
            ),
        ]
        track_state = TrackState(gpsl1, sat_states; doppler_estimator)

        # Initial tracking step
        signal =
            cis.(
                2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase,
            ) .*
            gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, start_code_phase)
        track_state = track(signal, track_state, sampling_frequency)

        # Run tracking iterations
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
                cis.(
                    2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase,
                ) .*
                gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, code_phase)
            track_state = track(signal, track_state, sampling_frequency)
        end

        final_doppler = get_carrier_doppler(track_state) / Hz
        error = abs(final_doppler - carrier_doppler / Hz)
        return error <= tolerance
    end

    # ConventionalPLLAndDLL: converges at 90Hz offset, fails at 100Hz
    @test test_convergence(90Hz, false) == true
    @test test_convergence(100Hz, false) == false

    # ConventionalAssistedPLLAndDLL: converges at 240Hz offset, fails at 250Hz
    @test test_convergence(240Hz, true) == true
    @test test_convergence(250Hz, true) == false
end

@testset "Track multiple systems of type $type" for type in
                                                    (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1CA()
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

    # Pin the loop bandwidths explicitly: this is a convergence test, so it
    # fixes the bandwidth it was calibrated for rather than depending on the
    # per-signal auto defaults (which would give Galileo E1B a tighter 4.5 Hz
    # loop that converges more slowly over this chunk schedule).
    estimator = ConventionalAssistedPLLAndDLL(;
        carrier_loop_filter_bandwidth = 18.0Hz,
        code_loop_filter_bandwidth = 1.0Hz,
    )
    gps_sat = TrackedSat(
        gpsl1,
        prn,
        start_code_phase,
        carrier_doppler_gps;
        doppler_estimator = estimator,
    )
    gal_sat = TrackedSat(
        galileo_e1b,
        prn,
        start_code_phase,
        carrier_doppler_gal;
        doppler_estimator = estimator,
    )
    track_state = @inferred TrackState(
        (gps = dictionary([prn => gps_sat]), gal = dictionary([prn => gal_sat]));
        doppler_estimator = estimator,
    )

    signal_temp =
        cis.(
            2π .* carrier_doppler_gps .* range ./ sampling_frequency .+ start_carrier_phase,
        ) .* gen_code(
            4000,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency_gps,
            start_code_phase,
        ) .+
        cis.(
            2π .* carrier_doppler_gal .* range ./ sampling_frequency .+ start_carrier_phase,
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
                carrier_phase_gps,
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
                carrier_phase_gal,
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
        5e-3
end

@testset "Tracking with intermediate frequency of $intermediate_frequency" for intermediate_frequency in
                                                                               (
    0.0Hz,
    -10000.0Hz,
    10000.0Hz,
    -30000.0Hz,
    30000.0Hz,
)
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    track_state = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)];
    )

    signal =
        cis.(
            2π .* (carrier_doppler + intermediate_frequency) .* range ./
            sampling_frequency .+ start_carrier_phase,
        ) .* get_code.(
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
                sampling_frequency .+ carrier_phase,
            ) .* get_code.(
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
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    track_state = @inferred TrackState(
        gpsl1,
        [
            TrackedSat(
                gpsl1,
                1,
                start_code_phase,
                carrier_doppler - 20Hz;
                num_ants = NumAnts(3),
            ),
        ],
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

@testset "Collect filtered prompts per track call" begin
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    samples_per_code = 4000  # one GPS L1 code period at 4 MHz

    # Generates a signal of `num_codes` full code periods of noise-free GPS L1.
    function make_signal(num_codes, start_carrier_phase, start_code_phase)
        n = samples_per_code * num_codes
        signal = zeros(ComplexF64, n)
        for code_idx = 0:(num_codes-1)
            carrier_phase =
                mod2pi(
                    2π * carrier_doppler * samples_per_code * code_idx /
                    sampling_frequency +
                    start_carrier_phase +
                    π,
                ) - π
            code_phase = mod(
                code_frequency * samples_per_code * code_idx / sampling_frequency +
                start_code_phase,
                1023,
            )
            r = (code_idx*samples_per_code):((code_idx+1)*samples_per_code-1)
            local_range = 0:(samples_per_code-1)
            signal[(code_idx*samples_per_code+1):((code_idx+1)*samples_per_code)] .=
                cis.(
                    2π .* carrier_doppler .* local_range ./ sampling_frequency .+
                    carrier_phase,
                ) .* gen_code(
                    samples_per_code,
                    gpsl1,
                    prn,
                    sampling_frequency,
                    code_frequency,
                    code_phase,
                )
        end
        signal
    end

    track_state =
        TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, carrier_doppler)])

    # First call: 10 code periods -> expect 10 filtered prompts
    signal1 = make_signal(10, π / 2, start_code_phase)
    track_state = track(signal1, track_state, sampling_frequency)
    sat = get_sat_state(track_state, prn)
    prompts1 = get_filtered_prompts(sat)
    @test prompts1 isa Vector{ComplexF64}
    @test length(prompts1) == 10
    # Last entry is exposed by the existing scalar accessor
    @test prompts1[end] == get_last_fully_integrated_filtered_prompt(track_state)
    # Save the buffer reference to verify it is reused, not reallocated.
    buffer_ref = prompts1

    # Second call: 5 code periods -> buffer is reset, so length is 5, not 15.
    signal2 = make_signal(5, π / 2, start_code_phase)
    track_state = track(signal2, track_state, sampling_frequency)
    sat = get_sat_state(track_state, prn)
    prompts2 = get_filtered_prompts(sat)
    @test length(prompts2) == 5
    @test prompts2[end] == get_last_fully_integrated_filtered_prompt(track_state)
    # Same underlying Vector object across calls -> capacity is retained and
    # no reallocation happens on a call that fits in the existing buffer.
    @test prompts2 === buffer_ref

    # Third call with 0 completed integrations: buffer empties.
    signal3 = zeros(ComplexF64, samples_per_code ÷ 4)  # less than one code period
    track_state = track(signal3, track_state, sampling_frequency)
    prompts3 = get_filtered_prompts(get_sat_state(track_state, prn))
    @test isempty(prompts3)
    @test prompts3 === buffer_ref
end

# Integration tests for L1C-D and L1C-P. Their 10-ms primary code period
# (vs 1 ms for L1 C/A) means each `track` call needs enough samples to land
# at least one primary-code boundary; the parameterized loop uses one full
# 10-ms primary period per iteration so each call completes exactly one
# integration. L1C-D has data bits (50 Hz), L1C-P is the pilot — both still
# need a working signal-path with a closed PLL/DLL.
#
# These tests rely on the per-signal default loop bandwidths
# (`default_carrier_loop_filter_bandwidth(::GPSL1C_D)` etc.) which size BL
# at ~0.018/T (Hz). For L1C-D / L1C-P that gives ~1.8 Hz carrier / ~0.1 Hz
# code — well inside the `BL * T < 0.4` stability bound for 10 ms integration.
@testset "Tracking single signal $name with $type samples" for (name, sig_type) in (
        ("GPSL1C_D", GPSL1C_D),
        ("GPSL1C_P", GPSL1C_P),
    ),
    type in (Float32, Float64)

    signal = sig_type()
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(signal) +
        get_code_frequency(signal)
    # L1C-P's TMBOC modulation requires `fs > 2 × code_freq × subcarrier_factor`
    # (~12.28 MHz minimum); use 15 MHz for both signals so the test fits the
    # same loop. L1C-D BPSK has no such floor but happily accepts 15 MHz too.
    sampling_frequency = 15e6Hz
    prn = 1
    primary_period_samples = 150000  # 10 ms at 15 MHz
    range = 0:(primary_period_samples-1)
    start_carrier_phase = π / 2

    track_state = @inferred TrackState(; signal)
    track_state = add_satellite!(
        track_state;
        prn,
        code_phase = start_code_phase,
        carrier_doppler = carrier_doppler - 5Hz,
    )

    function build_signal(carrier_phase, code_phase)
        s =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            gen_code(
                primary_period_samples,
                signal,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase,
            )
        Complex{type}.(s)
    end

    # First iteration to seed the loop filters.
    track!(
        build_signal(start_carrier_phase, start_code_phase),
        track_state,
        sampling_frequency,
    )

    iterations = 100
    primary_code_len = get_code_length(signal)
    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * carrier_doppler * primary_period_samples * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase = mod(
            code_frequency * primary_period_samples * i / sampling_frequency +
            start_code_phase,
            primary_code_len,
        )
        track!(build_signal(carrier_phase, code_phase), track_state, sampling_frequency)
    end
    # By now the PLL/DLL has had 100 × 10 ms = 1 s to converge on a clean
    # signal seeded 5 Hz off — should be well inside 1 Hz.
    final_doppler = get_carrier_doppler(track_state, :default, prn)
    @test abs(final_doppler - carrier_doppler) < 1.0Hz
end

# Regression test for the GPS L1C-D BOC(1,1) DLL early-late spacing.
#
# L1C-D's autocorrelation has nulls at ~±0.29 chip and side-lobes beyond, so
# the C/A-style 0.5-chip early/late taps land on the side-lobes and bias the
# DLL discriminator — the code loop then walks off the main peak as soon as the
# code is perturbed off-center (as it always is after an acquisition handoff).
# `get_default_correlator(::GPSL1C_D)` uses a narrow 0.1-chip spacing that keeps
# the taps on the main peak. This is exactly what the earlier single-signal
# L1C-D test above does NOT catch: it generates the signal at the tracker's own
# code phase (zero offset), where the symmetric autocorrelation gives E == L and
# the bias never shows, and it only asserts carrier-Doppler convergence (the
# bias corrupts the *code* loop, not the carrier).
#
# Here we feed a noisy 45 dB-Hz signal, seed the satellite 0.2 chip off the true
# code phase, and assert the estimated C/N0 stays high. With the narrow default
# spacing the loop pulls in and holds ~45 dB-Hz; with the buggy 0.5-chip spacing
# the code walks off and C/N0 collapses to ~31 dB-Hz. A wider-than-default code
# loop (1 Hz vs ~0.1 Hz) just makes the divergence show within ~1.3 s.
@testset "GPS L1C-D BOC code tracking holds lock under a code offset" begin
    Random.seed!(1234)
    signal = GPSL1C_D()
    prn = 1
    sampling_frequency = 10e6Hz
    period = 100000  # 10 ms at 10 MHz — one L1C-D primary code period
    # Baseband (carrier_doppler = 0), so the code frequency is nominal.
    code_frequency = get_code_frequency(signal)
    primary_len = get_code_length(signal)
    amplitude = 10^(45 / 20)            # 45 dB-Hz with noise_std below
    noise_std = sqrt(10e6)
    range = 0:(period-1)

    estimator = ConventionalAssistedPLLAndDLL(;
        carrier_loop_filter_bandwidth = 1.8Hz,
        code_loop_filter_bandwidth = 1.0Hz,
    )
    track_state = TrackState(; signal, doppler_estimator = estimator)
    # 0.2-chip seed offset — a realistic acquisition handoff error that puts the
    # DLL discriminator onto the BOC autocorrelation slope. Uses the group's
    # default correlator, so reverting the default spacing breaks this test.
    track_state =
        add_satellite!(track_state; prn, code_phase = 0.2, carrier_doppler = 0.0Hz)

    build(code_phase) =
        gen_code(period, signal, prn, sampling_frequency, code_frequency, code_phase) .*
        amplitude .+ randn(ComplexF64, period) .* noise_std

    for i = 0:130
        code_phase = mod(code_frequency * period * i / sampling_frequency, primary_len)
        track!(build(code_phase), track_state, sampling_frequency)
    end

    # Narrow spacing holds ~45 dB-Hz; the 0.5-chip bug collapses to ~31 dB-Hz.
    @test estimate_cn0(track_state, prn) > 40dBHz
end

# Multi-signal integration test for the README's flagship use case: a single
# satellite tracked on the modern GPS L1 signal trio (L1C_P pilot + L1C_D
# data + L1 C/A legacy). Verifies all three signals' correlators run, and
# that the PLL/DLL — driven by `signals[1]` = L1C_P (the pilot) — converges.
@testset "Tracking multi-signal (L1C_P, L1C_D, L1CA) on one sat" begin
    # L1C-P's TMBOC modulation requires fs > ~12.28 MHz; use 15 MHz.
    sampling_frequency = 15e6Hz
    n_samples = 150000  # 10 ms — one L1C-P/D primary period; ten L1CA periods
    prn = 11
    carrier_doppler = 1234.0Hz
    init_offset = 5.0Hz

    # Driver is signals[1] = L1C_P at 10 ms; the per-signal default loop
    # bandwidth picks the right (~1.8 Hz) value automatically.
    track_state =
        TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
    track_state = add_satellite!(
        track_state;
        prn,
        group = :modern_gps,
        code_phase = 0.0,
        carrier_doppler = carrier_doppler - init_offset,
    )

    # Synthesize a signal that contains all three signals on the same carrier
    # (modeling what a real GPS satellite broadcasts on L1).
    range = 0:(n_samples-1)
    function build_signal(code_phase, carrier_phase)
        carrier =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase)
        signal_cp = let s = GPSL1C_P()
            cf =
                carrier_doppler * get_code_center_frequency_ratio(s) + get_code_frequency(s)
            carrier .* gen_code(n_samples, s, prn, sampling_frequency, cf, code_phase)
        end
        signal_cd = let s = GPSL1C_D()
            cf =
                carrier_doppler * get_code_center_frequency_ratio(s) + get_code_frequency(s)
            carrier .* gen_code(n_samples, s, prn, sampling_frequency, cf, code_phase)
        end
        signal_ca = let s = GPSL1CA()
            cf =
                carrier_doppler * get_code_center_frequency_ratio(s) + get_code_frequency(s)
            # L1 C/A code phase wraps every 1023 chips; map the shared phase
            # to its primary period.
            carrier .* gen_code(
                n_samples,
                s,
                prn,
                sampling_frequency,
                cf,
                mod(code_phase, get_code_length(s)),
            )
        end
        signal_cp .+ signal_cd .+ signal_ca
    end

    # Seed + 200 iterations of clean signal at 10 ms per call.
    track!(build_signal(0.0, 0.0), track_state, sampling_frequency)

    iterations = 200
    cp_primary_len = get_code_length(GPSL1C_P())
    for i = 1:iterations
        carrier_phase =
            mod2pi(2π * carrier_doppler * n_samples * i / sampling_frequency + π) - π
        # Drive code_phase off the L1C_P/L1C_D primary length (10230 chips);
        # L1CA mods inside build_signal.
        cf_cp =
            carrier_doppler * get_code_center_frequency_ratio(GPSL1C_P()) +
            get_code_frequency(GPSL1C_P())
        code_phase = mod(cf_cp * n_samples * i / sampling_frequency, cp_primary_len)
        track!(build_signal(code_phase, carrier_phase), track_state, sampling_frequency)
    end

    # All three signals must have completed integrations.
    sat = get_sat_state(track_state, :modern_gps, prn)
    @test length(get_filtered_prompts(sat.signals[1])) > 0  # L1C_P
    @test length(get_filtered_prompts(sat.signals[2])) > 0  # L1C_D
    @test length(get_filtered_prompts(sat.signals[3])) > 0  # L1CA

    # PLL/DLL has had ~2 s to converge on a clean superposition.
    final_doppler = get_carrier_doppler(track_state, :modern_gps, prn)
    @test abs(final_doppler - carrier_doppler) < 5.0Hz
end

# Regression test for issue #117: a satellite must not deadlock when fed
# chunks of exactly one code period.
#
# The trigger is a *continuous* signal sliced into fixed one-primary-code-
# period chunks (the most natural streaming pattern — 1 ms of GPS L5I at
# 20 MHz = 20000 samples = one code period) together with a *negative*
# code Doppler. With negative Doppler a full code block spans slightly
# more than one chunk (`ceil(10230·fs/(f_code+Δ)) = 20001 > 20000`), so a
# block can never complete inside a single chunk; it must carry forward.
#
# Before the fix, once GPS L5I's NH10 secondary code synced, the per-call
# code-phase snap re-anchored `code_phase` to a primary-block boundary on
# *every* call — discarding the partial integration's within-block phase
# and pinning `code_phase` to a value from which no future chunk could
# ever complete the block. The satellite wedged forever: `code_phase`
# frozen, `integrated_samples` growing without bound, the correlator
# output frozen. Crucially this needs a real secondary-code sync, so the
# signal is generated continuously (not re-aligned per call) and fed in
# fixed slices.
@testset "Does not deadlock on one-code-period chunks (issue #117)" begin
    signal = GPSL5I()
    carrier_doppler = -200Hz                 # negative ⇒ negative code Doppler
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(signal) +
        get_code_frequency(signal)
    sampling_frequency = 20e6Hz
    prn = 1
    chunk = 20000                            # 1 ms = one L5I code period at 20 MHz
    num_calls = 60                           # well past the 10-block NH10 sync

    # One continuous signal, sliced into fixed `chunk`-sample pieces.
    total_samples = chunk * (num_calls + 1)
    t = 0:(total_samples-1)
    long_signal =
        cis.(2π .* carrier_doppler .* t ./ sampling_frequency) .*
        gen_code(total_samples, signal, prn, sampling_frequency, code_frequency, 0.0)

    track_state = TrackState(; signal)
    track_state = add_satellite!(
        track_state;
        prn,
        code_phase = 0.0,
        carrier_doppler = carrier_doppler - 20Hz,
    )

    track!((@view long_signal[1:chunk]), track_state, sampling_frequency)
    prompts = ComplexF64[]
    for i = 1:num_calls
        chunk_signal = @view long_signal[(i*chunk+1):((i+1)*chunk)]
        track!(chunk_signal, track_state, sampling_frequency)
        push!(prompts, get_prompt(get_last_fully_integrated_correlator(track_state, prn)))
    end

    sat = get_sat_state(track_state, prn)

    # The secondary code must actually have synced — that's the regime
    # the bug lived in.
    @test has_bit_or_secondary_code_been_found(sat)

    # Deadlock signature: `integrated_samples` grows without bound because
    # no integration ever completes. After the fix it completes every
    # block, so it stays within one chunk's worth of samples.
    @test get_integrated_samples(track_state, prn) <= 2 * chunk

    # Deadlock signature: the correlator output freezes. After the fix the
    # prompt keeps updating, so the last several prompts are not all equal.
    @test length(unique(prompts[(end-9):end])) > 1

    # And it actually tracks: the loops stay locked on the true Doppler
    # rather than diverging or freezing.
    @test abs(get_carrier_doppler(track_state, prn) - carrier_doppler) < 5.0Hz
end

end
