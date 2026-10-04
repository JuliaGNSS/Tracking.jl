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
    get_noise_density,
    get_weights,
    AbstractPostCorrFilter,
    DefaultPostCorrFilter,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL
import Tracking
using StaticArrays: SVector
using Unitful: ustrip, uconvert

# Normalise a code replica to unit peak amplitude, so Galileo E1B's integer CBOC
# amplitudes (±25/±13) don't swamp ±1 BPSK codes when summed into one sample buffer.
_unit_code(c) = (m = maximum(abs, c); m == 0 ? float(c) : c ./ m)

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

    # Pull-in ranges of the default 18 Hz loop (90 / 240 Hz before #244).
    # ConventionalPLLAndDLL: converges at 40Hz offset, fails at 60Hz. In between
    # the outcome depends on where cycle slips leave the loop.
    @test test_convergence(40Hz, false) == true
    @test test_convergence(60Hz, false) == false

    # ConventionalAssistedPLLAndDLL: converges at 170Hz offset, fails at 180Hz,
    # matching a classic per-code-period FLL update.
    @test test_convergence(170Hz, true) == true
    @test test_convergence(180Hz, true) == false
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
        ) .* _unit_code(
            gen_code(
                4000,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency_gps,
                start_code_phase,
            ),
        ) .+
        cis.(
            2π .* carrier_doppler_gal .* range ./ sampling_frequency .+ start_carrier_phase,
        ) .* _unit_code(
            gen_code(
                4000,
                galileo_e1b,
                prn,
                sampling_frequency,
                code_frequency_gal,
                start_code_phase,
            ),
        )
    # Scale to fill (half) the integer type's range from the (unit-amplitude) signal.
    scaling =
        type <: Integer ? fld(typemax(type), 2 * ceil(Int, maximum(abs, signal_temp))) : 1
    signal =
        type <: Integer ?
        complex.(
            round.(type, real.(signal_temp) * scaling),
            round.(type, imag.(signal_temp) * scaling),
        ) : Complex{type}.(signal_temp)

    track_state = @inferred track(signal, track_state, sampling_frequency)

    # 1.33 s: the 18 Hz carrier loop settles the 90° start phase in about 1.2 s.
    iterations = 5000
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
            ) .* _unit_code(
                gen_code(
                    4000,
                    gpsl1,
                    prn,
                    sampling_frequency,
                    code_frequency_gps,
                    code_phase_gps,
                ),
            ) .+
            cis.(
                2π .* carrier_doppler_gal .* range ./ sampling_frequency .+
                carrier_phase_gal,
            ) .* _unit_code(
                gen_code(
                    4000,
                    galileo_e1b,
                    prn,
                    sampling_frequency,
                    code_frequency_gal,
                    code_phase_gal,
                ),
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
    # Galileo E1B: GNSSSignals' Int8 CBOC approximation leaves a deterministic
    # ~0.006 rad residual for every sample type, hence 1e-2 (GNSSSignals #90).
    @test mod(get_carrier_phase(track_state, :gal, 1), π) ≈ mod(comp_carrier_phase_gal, π) atol =
        1e-2
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

    # 3.5 s: the 18 Hz carrier loop settles the 20 Hz start offset in about 3.2 s.
    iterations = 3500
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

    # 3.5 s: the 18 Hz carrier loop settles the 20 Hz start offset in about 3.2 s.
    iterations = 3500
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
end

# Multi-antenna C/N₀ with independent noise at a different level per antenna, so
# reading the wrong antenna's floor is observable (unlike the identical columns above).
@testset "Multi-antenna C/N₀ under the default filter matches its own column" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    prn = 1
    samples_per_code = 5000
    num_ants = 3
    # Distinct per-antenna noise levels: antenna m's floor is `scaleₘ²·σ²/f_s`.
    scales = (1.0, 1.5, 2.0)

    function make_array_signal(rng, code_idx)
        code = get_code.(
            gpsl1,
            get_code_frequency(gpsl1) .* (0:(samples_per_code-1)) ./ sampling_frequency .+ code_idx * get_code_length(gpsl1),
            prn,
        )
        clean = 30.0 .* ComplexF64.(code)
        sigma = sqrt(ustrip(Hz, sampling_frequency))
        reduce(
            hcat,
            (
                clean .+ (scale * sigma) .* randn(rng, ComplexF64, samples_per_code) for
                scale in scales
            ),
        )
    end

    # One array run, and one single-antenna control fed the *same* last column.
    array_state = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, prn, 0.0, 0.0Hz; num_ants = NumAnts(num_ants))],
    )
    single_state = @inferred TrackState(gpsl1, [TrackedSat(gpsl1, prn, 0.0, 0.0Hz)])
    rng_array = Random.MersenneTwister(2718)
    rng_single = Random.MersenneTwister(2718)
    for i = 0:199
        array_state = @inferred track(
            make_array_signal(rng_array, i),
            array_state,
            sampling_frequency,
        )
        single_state = @inferred track(
            view(make_array_signal(rng_single, i), :, num_ants),
            single_state,
            sampling_frequency,
        )
    end

    # The array's window measures a covariance, one entry per antenna pair.
    R = get_noise_density(array_state.noise_estimators.GPSL1CA)
    @test size(R) == (num_ants, num_ants)
    # `DefaultPostCorrFilter` selects the last antenna, so the covariance reduced
    # through its weights is that antenna's floor: identical C/N₀ to the control.
    @test estimate_cn0(array_state, 1) === estimate_cn0(single_state, 1)
    @test isfinite(ustrip(uconvert(dBHz, estimate_cn0(array_state, 1))))
end

@testset "A beamformer's C/N₀ follows its own weights" begin
    # An equal-weight combiner over independent noise sees `‖w‖²` = 1/3 of the floor
    # while the signal adds coherently: C/N₀ ~10log10(3) dB above one antenna.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    prn = 1
    samples_per_code = 5000
    num_ants = 3

    struct MeanBeamformer <: AbstractPostCorrFilter end
    Tracking.update(f::MeanBeamformer, prompt) = f
    Tracking.get_weights(::MeanBeamformer, ::NumAnts{1}) = 1.0 + 0.0im
    Tracking.get_weights(::MeanBeamformer, ::NumAnts{M}) where {M} =
        SVector{M,ComplexF64}(ntuple(_ -> 1 / M + 0.0im, M))

    function make_array_signal(rng)
        code = get_code.(
            gpsl1,
            get_code_frequency(gpsl1) .* (0:(samples_per_code-1)) ./ sampling_frequency,
            prn,
        )
        clean = 30.0 .* ComplexF64.(code)
        sigma = sqrt(ustrip(Hz, sampling_frequency))
        reduce(
            hcat,
            (clean .+ sigma .* randn(rng, ComplexF64, samples_per_code) for _ = 1:num_ants),
        )
    end

    function run(filter)
        state = TrackState(
            gpsl1,
            [
                TrackedSat(
                    gpsl1,
                    prn,
                    0.0,
                    0.0Hz;
                    num_ants = NumAnts(num_ants),
                    post_corr_filter = filter,
                ),
            ],
        )
        rng = Random.MersenneTwister(31415)
        for _ = 1:200
            state = track(make_array_signal(rng), state, sampling_frequency)
        end
        state
    end

    beamformed = ustrip(uconvert(dBHz, estimate_cn0(run(MeanBeamformer()), 1)))
    one_antenna = ustrip(uconvert(dBHz, estimate_cn0(run(DefaultPostCorrFilter()), 1)))
    @test isfinite(beamformed)
    @test beamformed - one_antenna ≈ 10 * log10(num_ants) atol = 1.0
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
    # Same Vector object across calls -> no reallocation.
    @test prompts2 === buffer_ref

    # Third call with 0 completed integrations: buffer empties.
    signal3 = zeros(ComplexF64, samples_per_code ÷ 4)  # less than one code period
    track_state = track(signal3, track_state, sampling_frequency)
    prompts3 = get_filtered_prompts(get_sat_state(track_state, prn))
    @test isempty(prompts3)
    @test prompts3 === buffer_ref
end

# L1C-D (data) and L1C-P (pilot) closed-loop tracking. Each call feeds one full 10 ms
# primary period, so it completes exactly one integration. Uses the default loop
# bandwidths (carrier capped to 9 Hz at 10 ms, code 1 Hz).
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
    # L1C-P's TMBOC needs fs > ~12.28 MHz; 15 MHz for both signals.
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

    iterations = 200
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
    # 2 s to converge from 5 Hz off to within 1 Hz.
    final_doppler = get_carrier_doppler(track_state, :default, prn)
    @test abs(final_doppler - carrier_doppler) < 1.0Hz
end

# Regression test for GPS L1C-D BOC(1,1) code tracking: with C/A-style early/late taps
# the DLL is biased by the BOC side lobes and walks off the peak once the code starts
# off-center (see `get_default_correlator` in gps/l1c_d.jl). The test above misses
# this, as it starts at zero code offset and checks only the carrier.
#
# Noisy 45 dB-Hz signal, seeded 0.2 chip off; C/N0 must stay high. The 4 Hz code loop
# makes a divergence show within ~1.3 s.
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
        code_loop_filter_bandwidth = 4.0Hz,
    )
    track_state = TrackState(; signal, doppler_estimator = estimator)
    # Realistic acquisition handoff error; uses the default correlator under test.
    track_state =
        add_satellite!(track_state; prn, code_phase = 0.2, carrier_doppler = 0.0Hz)

    build(code_phase) =
        gen_code(period, signal, prn, sampling_frequency, code_frequency, code_phase) .*
        amplitude .+ randn(ComplexF64, period) .* noise_std

    for i = 0:130
        code_phase = mod(code_frequency * period * i / sampling_frequency, primary_len)
        track!(build(code_phase), track_state, sampling_frequency)
    end

    # Holds ~45 dB-Hz when locked; a walked-off code loop collapses to ~31 dB-Hz.
    @test estimate_cn0(track_state, prn) > 40dBHz
end

# The README's use case: one satellite on L1C_P + L1C_D + L1 C/A. All three
# correlators must run and the loops, driven by `signals[1]` = L1C_P, must converge.
@testset "Tracking multi-signal (L1C_P, L1C_D, L1CA) on one sat" begin
    # L1C-P's TMBOC modulation requires fs > ~12.28 MHz; use 15 MHz.
    sampling_frequency = 15e6Hz
    n_samples = 150000  # 10 ms — one L1C-P/D primary period; ten L1CA periods
    prn = 11
    carrier_doppler = 1234.0Hz
    init_offset = 5.0Hz

    # Driver is signals[1] = L1C_P at 10 ms, with the default loop bandwidths.
    track_state =
        TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
    track_state = add_satellite!(
        track_state;
        prn,
        group = :modern_gps,
        code_phase = 0.0,
        carrier_doppler = carrier_doppler - init_offset,
    )

    # All three signals on the same carrier.
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
            # Map the shared code phase into L1 C/A's 1023-chip period.
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

# Regression test for issue #117: no deadlock on chunks of exactly one code period.
# With negative code Doppler a code block spans slightly more than one chunk (20001 >
# 20000 samples), so every block must carry across calls; after NH10 secondary-code
# sync this once wedged the satellite. Needs a real sync, so one continuous signal is
# sliced into fixed chunks.
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

    # The regime the bug lived in.
    @test has_bit_or_secondary_code_been_found(sat)

    # Deadlock signatures: `integrated_samples` grows without bound ...
    @test get_integrated_samples(track_state, prn) <= 2 * chunk

    # ... and the prompt freezes.
    @test length(unique(prompts[(end-9):end])) > 1

    # And the loops stay locked.
    @test abs(get_carrier_doppler(track_state, prn) - carrier_doppler) < 5.0Hz
end

end
