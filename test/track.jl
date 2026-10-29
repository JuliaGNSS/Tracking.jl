module TrackTest

using Test: @test, @testset, @inferred
using Random: Random
using Unitful: Hz, dBHz
using Dictionaries: dictionary
using GNSSSignals:
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    GPSL5I,
    GalileoE1B,
    gen_code,
    get_chip_rate_to_center_frequency_ratio,
    get_chip_rate,
    get_code,
    get_code_length,
    get_band,
    get_band_id

using Tracking:
    TrackedSat,
    TrackState,
    track,
    track!,
    add_satellite!,
    get_code_phase,
    get_carrier_phase,
    get_carrier_phase_correction,
    get_code_doppler,
    get_carrier_doppler,
    get_resume_index,
    get_last_fully_integrated_filtered_prompt,
    get_last_fully_integrated_correlator,
    get_filtered_prompts,
    get_integrated_samples,
    get_sat_state,
    get_doppler_estimator_state,
    reset_loop_filters!
import TrackingLoops
using TrackingLoops:
    get_prompt,
    estimate_cn0,
    has_bit_or_secondary_code_been_found,
    NumAnts,
    get_num_ants,
    get_noise_density,
    AbstractPostCorrFilter,
    DefaultPostCorrFilter,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    NCOReferencedPLLAndDLL,
    SatNCOReferencedPLLAndDLL,
    VectorPLLAndDLL,
    SatVectorPLLAndDLL,
    LoopRecord
import Tracking
using StaticArrays: SVector
using Unitful: ustrip, uconvert

# Normalise a code replica to unit peak amplitude. The CBOC code (Galileo E1B) has
# integer sub-carrier amplitudes (±25/±13) against ±1 for BPSK, and would swamp
# GPS in a shared fixed-point buffer.
_unit_code(c) = (m = maximum(abs, c); m == 0 ? float(c) : c ./ m)

@testset "Tracking with samples of type $type" for type in
                                                   (Int16, Int32, Int64, Float32, Float64)
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency =
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1) +
        get_chip_rate(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    track_state = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)],
    )

    samples_temp =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, start_code_phase)
    scaling = 512
    samples =
        type <: Integer ?
        complex.(
            round.(type, real.(samples_temp) * scaling),
            round.(type, imag.(samples_temp) * scaling),
        ) : Complex{type}.(samples_temp)

    track_state = @inferred track(samples, track_state, sampling_frequency)

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
        samples_temp =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, code_phase)
        samples =
            type <: Integer ?
            complex.(
                round.(type, real.(samples_temp) * scaling),
                round.(type, imag.(samples_temp) * scaling),
            ) : Complex{type}.(samples_temp)
        track_state = @inferred track(samples, track_state, sampling_frequency)
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
        @test get_resume_index(track_state) == 4001
        tracked_prompts[i] = get_last_fully_integrated_filtered_prompt(track_state)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase / 2π
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 5e-5
    @test tracked_carrier_phases[end] + 0.5 ≈ carrier_phases[end] atol = 1e-3 / 2π
end

# GPS L1 C/A PRN 1 at a 1 kHz Doppler in noise, cut into 1 ms chunks at 4 MHz.
function _noisy_gpsl1_chunks(num_chunks)
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    num_samples = 4000 * num_chunks
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    t = (0:(num_samples-1)) ./ ustrip(Hz, fs)
    code = gen_code(num_samples, gpsl1, 1, fs, get_chip_rate(gpsl1) + code_doppler, 100.0)
    rng = Random.Xoshiro(1)
    samples =
        ComplexF32.(code .* cis.(2π .* ustrip(Hz, carrier_doppler) .* t .+ 0.3)) .+
        0.5f0 .* randn(rng, ComplexF32, num_samples)
    gpsl1, fs, carrier_doppler, [samples[((k-1)*4000+1):(k*4000)] for k = 1:num_chunks]
end

@testset "Through track, the NCO-referenced and the vector loops are the conventional assisted loop" begin
    # `track` hands over the chunk's own replica word and no landing sample, so
    # the delay-aware estimator must equal the conventional one at the same
    # bandwidths (it pulls in at 18 Hz by default), also after a reset. Outside
    # a vector loop the vector estimator is the scalar loop it wraps, bit for bit.
    gpsl1, fs, carrier_doppler, chunks = _noisy_gpsl1_chunks(200)
    function run(estimator)
        track_state = TrackState(; signal = gpsl1, doppler_estimator = estimator)
        track_state = add_satellite!(
            track_state;
            prn = 1,
            code_phase = 100.2,
            carrier_doppler = 1005.0Hz,
        )
        history = Tuple{typeof(1.0Hz),Float64}[]
        for chunk in chunks
            track_state = track(chunk, track_state, fs)
            sat_state = get_sat_state(track_state, 1)
            push!(history, (get_carrier_doppler(sat_state), get_code_phase(sat_state)))
        end
        track_state = track(first(chunks), reset_loop_filters!(track_state), fs)
        sat_state = get_sat_state(track_state, 1)
        push!(history, (get_carrier_doppler(sat_state), get_code_phase(sat_state)))
        # A scalar loop's replica is its carrier estimate; so is a satellite's
        # out of the vector loop.
        @test get_carrier_phase_correction(track_state, 1, fs) == 0.0
        @test get_carrier_phase_correction(track_state, :default, 1, fs) == 0.0
        history, track_state
    end
    conventional, _ = run(ConventionalAssistedPLLAndDLL())
    narrow_conventional, _ =
        run(ConventionalAssistedPLLAndDLL(; wide_carrier_loop_filter_bandwidth = 18.0Hz))
    referenced, track_state = run(NCOReferencedPLLAndDLL())
    @test referenced == narrow_conventional
    @test get_doppler_estimator_state(get_sat_state(track_state, 1)) isa
          SatNCOReferencedPLLAndDLL
    @test abs(conventional[200][1] - carrier_doppler) < 5Hz
    for (estimator, scalar) in (
        (VectorPLLAndDLL(gpsl1), conventional),
        (VectorPLLAndDLL(gpsl1; inner = NCOReferencedPLLAndDLL()), narrow_conventional),
    )
        vector, track_state = run(estimator)
        @test vector == scalar
        @test get_doppler_estimator_state(get_sat_state(track_state, 1)) isa
              SatVectorPLLAndDLL
    end
end

# An estimator implementing only TrackingLoops' interface, counting the records
# it hands to the loop it wraps.
struct CountingEstimator{E<:TrackingLoops.AbstractDopplerEstimator} <:
       TrackingLoops.AbstractDopplerEstimator
    inner::E
end
struct SatCountingEstimator{S}
    inner::S
    num_records::Int
end
TrackingLoops.init_estimator_state(
    estimator::CountingEstimator,
    driver_signal,
    carrier_doppler,
    code_doppler,
) = SatCountingEstimator(
    TrackingLoops.init_estimator_state(
        estimator.inner,
        driver_signal,
        carrier_doppler,
        code_doppler,
    ),
    0,
)
function TrackingLoops.step_loop(
    estimator::CountingEstimator,
    state::SatCountingEstimator,
    record,
    words,
    landing_sample::Int64,
)
    inner, carrier_doppler, code_doppler =
        TrackingLoops.step_loop(estimator.inner, state.inner, record, words, landing_sample)
    SatCountingEstimator(inner, state.num_records + 1), carrier_doppler, code_doppler
end

@testset "reset_loop_filters! keeps the stage the carrier loop has reached" begin
    # A reset is a handoff, not a re-acquisition. A new satellite starts
    # FLL-assisted.
    gpsl1, fs, carrier_doppler, chunks = _noisy_gpsl1_chunks(800)
    track_state = TrackState(; signal = gpsl1)
    track_state =
        add_satellite!(track_state; prn = 1, code_phase = 100.2, carrier_doppler = 1005.0Hz)
    stage(prn) = TrackingLoops.carrier_loop_stage(
        get_doppler_estimator_state(get_sat_state(track_state, prn)),
    )
    @test stage(1) == TrackingLoops.FLL_ASSISTED_PLL
    for chunk in chunks
        track!(chunk, track_state, fs)
    end
    @test stage(1) == TrackingLoops.NARROW_PLL
    reset_loop_filters!(track_state)
    @test stage(1) == TrackingLoops.NARROW_PLL
    track_state = add_satellite!(track_state; prn = 2, carrier_doppler = 0.0Hz)
    @test stage(2) == TrackingLoops.FLL_ASSISTED_PLL
end

@testset "Through track, an estimator implementing step_loop runs without adaptation" begin
    gpsl1, fs, carrier_doppler, chunks = _noisy_gpsl1_chunks(100)
    function run(estimator)
        track_state = TrackState(; signal = gpsl1, doppler_estimator = estimator)
        track_state = add_satellite!(
            track_state;
            prn = 1,
            code_phase = 100.2,
            carrier_doppler = 1005.0Hz,
        )
        for chunk in chunks
            track!(chunk, track_state, fs)
        end
        # No `reset_estimator_state` method: the reset re-initializes.
        reset_loop_filters!(track_state)
        get_sat_state(track_state, 1)
    end
    reference = run(ConventionalAssistedPLLAndDLL())
    counted = run(CountingEstimator(ConventionalAssistedPLLAndDLL()))
    @test get_carrier_doppler(counted) == get_carrier_doppler(reference)
    @test get_code_phase(counted) == get_code_phase(reference)
    @test get_doppler_estimator_state(counted).num_records == 0
end

@testset "A scalar loop has no navigation solution to report" begin
    track_state = TrackState(; signal = GPSL1CA())
    @test TrackingLoops.navigation_solution(track_state) === nothing
    @test TrackingLoops.navigation_status(track_state) === nothing
    @test TrackingLoops.navigation_cycle(track_state) === nothing
    @test TrackingLoops.navigation_epoch(track_state) === nothing
    @test TrackingLoops.satellite_report(track_state, GPSL1CA(), 1) === nothing
    vector =
        TrackState(; signal = GPSL1CA(), doppler_estimator = VectorPLLAndDLL(GPSL1CA()))
    @test TrackingLoops.navigation_cycle(vector) == 0
end

# An estimator that keeps every record it is handed, around the loop it wraps.
struct RecordingEstimator{E<:TrackingLoops.AbstractDopplerEstimator} <:
       TrackingLoops.AbstractDopplerEstimator
    inner::E
    records::Vector{Tuple{Int,Int,Int,Float64}}
end
TrackingLoops.init_estimator_state(
    estimator::RecordingEstimator,
    driver_signal,
    carrier_doppler,
    code_doppler,
) = TrackingLoops.init_estimator_state(
    estimator.inner,
    driver_signal,
    carrier_doppler,
    code_doppler,
)
function TrackingLoops.step_loop(
    estimator::RecordingEstimator,
    state,
    record::LoopRecord,
    words,
    landing_sample::Int64,
)
    push!(
        estimator.records,
        (record.prn, record.sample_index, record.fold_end, record.code_phase),
    )
    TrackingLoops.step_loop(estimator.inner, state, record, words, landing_sample)
end

@testset "Through track, the records name their satellite and share one sample grid" begin
    # What the vector loop needs: the satellite, the replica's code phase at the
    # record's end, and a sample index that runs on across `track` calls.
    gpsl1, fs, carrier_doppler, chunks = _noisy_gpsl1_chunks(20)
    estimator =
        RecordingEstimator(ConventionalAssistedPLLAndDLL(), Tuple{Int,Int,Int,Float64}[])
    track_state = TrackState(; signal = gpsl1, doppler_estimator = estimator)
    track_state =
        add_satellite!(track_state; prn = 7, code_phase = 100.2, carrier_doppler = 1005.0Hz)
    for chunk in chunks
        track_state = track(chunk, track_state, fs)
    end
    records = estimator.records
    @test length(records) >= 18
    @test all(r -> r[1] == 7, records)
    sample_indices = [r[2] for r in records]
    @test issorted(sample_indices) && allunique(sample_indices)
    @test last(sample_indices) > 4000 * 18
    @test all(r -> r[3] >= r[2], records)
    # One code period per record: each ends within a sample of a block boundary.
    code_rate = ustrip(Hz, get_chip_rate(gpsl1)) / ustrip(Hz, fs)
    @test all(r -> abs(rem(r[4], 1023, RoundNearest)) <= code_rate, records)
    @test track_state.sample_offsets[get_band_id(get_band(gpsl1))] == 4000 * 20
    # `track` leaves its input's sample grid where it was.
    before = copy(track_state.sample_offsets)
    track(first(chunks), track_state, fs)
    @test track_state.sample_offsets == before
    track!(first(chunks), track_state, fs)
    @test track_state.sample_offsets[get_band_id(get_band(gpsl1))] == 4000 * 21
    # `track` shares the estimator with its input; a `deepcopy` branches it.
    num_records = length(estimator.records)
    branch = track(first(chunks), deepcopy(track_state), fs)
    @test length(estimator.records) == num_records
    @test length(branch.doppler_estimator.records) > num_records
    track(first(chunks), track_state, fs)
    @test length(estimator.records) > num_records
end

@testset "Tracking with large initial Doppler offset" begin
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1) +
        get_chip_rate(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999 # tracking 1 ms
    start_carrier_phase = π / 2
    # 3 s: these test the pull-in range, not its speed. The default 18 Hz loop
    # takes up to about 1.1 s with the FLL assist and 1.9 s without (65 Hz).
    iterations = 3000
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

        samples =
            cis.(
                2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase,
            ) .*
            gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, start_code_phase)
        track_state = track(samples, track_state, sampling_frequency)

        for i = 1:iterations
            carrier_phase =
                mod2pi(
                    2π * carrier_doppler * 4000 * i / sampling_frequency +
                    start_carrier_phase +
                    π,
                ) - π
            code_phase =
                mod(code_frequency * 4000 * i / sampling_frequency + start_code_phase, 1023)
            samples =
                cis.(
                    2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase,
                ) .*
                gen_code(4000, gpsl1, prn, sampling_frequency, code_frequency, code_phase)
            track_state = track(samples, track_state, sampling_frequency)
        end

        final_doppler = get_carrier_doppler(track_state) / Hz
        error = abs(final_doppler - carrier_doppler / Hz)
        return error <= tolerance
    end

    # Without an FLL the loop pulls in by cycle slipping at the wide 50 Hz: it
    # converged up to 68 Hz and not from 94 Hz on; in between it depends on
    # where the slips leave it.
    @test test_convergence(65Hz, false) == true
    @test test_convergence(100Hz, false) == false

    # The FLL's unambiguous range, ±1/(4T) = 250 Hz at 1 ms, independent of the
    # bandwidth: past it the loop locks onto the alias 500 Hz away.
    @test test_convergence(245Hz, true) == true
    @test test_convergence(255Hz, true) == false
end

@testset "Track multiple constellations with samples of type $type" for type in (
    Int16,
    Int32,
    Int64,
    Float32,
    Float64,
)
    gpsl1 = GPSL1CA()
    galileo_e1b = GalileoE1B()
    carrier_doppler_gps = 200Hz
    carrier_doppler_gal = 1200Hz
    start_code_phase = 100
    code_frequency_gps =
        carrier_doppler_gps * get_chip_rate_to_center_frequency_ratio(gpsl1) +
        get_chip_rate(gpsl1)
    code_frequency_gal =
        carrier_doppler_gal * get_chip_rate_to_center_frequency_ratio(galileo_e1b) +
        get_chip_rate(galileo_e1b)
    sampling_frequency = 15e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    # Pinned to the bandwidths this convergence test was calibrated for.
    estimator = ConventionalAssistedPLLAndDLL(;
        wide_carrier_loop_filter_bandwidth = 18.0Hz,
        narrow_carrier_loop_filter_bandwidth = 18.0Hz,
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

    samples_temp =
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
    # Fill half the integer type's range.
    scaling =
        type <: Integer ? fld(typemax(type), 2 * ceil(Int, maximum(abs, samples_temp))) : 1
    samples =
        type <: Integer ?
        complex.(
            round.(type, real.(samples_temp) * scaling),
            round.(type, imag.(samples_temp) * scaling),
        ) : Complex{type}.(samples_temp)

    track_state = @inferred track(samples, track_state, sampling_frequency)

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
        samples_temp =
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
        samples =
            type <: Integer ?
            complex.(
                round.(type, real.(samples_temp) * scaling),
                round.(type, imag.(samples_temp) * scaling),
            ) : Complex{type}.(samples_temp)
        track_state = @inferred track(samples, track_state, sampling_frequency)
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
    @test mod(get_carrier_phase(track_state, :gps, 1), 0.5) ≈
          mod(comp_carrier_phase_gps / 2π, 0.5) atol = 2e-2 / 2π
    @test get_code_phase(track_state, :gal, 1) ≈ comp_code_phase_gal atol = 5e-3
    # The Int8 CBOC code approximates the BOC(1,1)/BOC(6,1) mix and leaves a
    # deterministic ~0.006 rad residual, the same for every sample type
    # (GNSSSignals #90).
    @test mod(get_carrier_phase(track_state, :gal, 1), 0.5) ≈
          mod(comp_carrier_phase_gal / 2π, 0.5) atol = 1e-2 / 2π
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
    code_frequency = carrier_doppler / 1540 + get_chip_rate(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    track_state = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)];
    )

    samples =
        cis.(
            2π .* (carrier_doppler + intermediate_frequency) .* range ./
            sampling_frequency .+ start_carrier_phase,
        ) .* get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn,
        )

    track_state =
        @inferred track(samples, track_state, sampling_frequency; intermediate_frequency)

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
        samples =
            cis.(
                2π .* (carrier_doppler + intermediate_frequency) .* range ./
                sampling_frequency .+ carrier_phase,
            ) .* get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn,
            )
        track_state = @inferred track(
            samples,
            track_state,
            sampling_frequency;
            intermediate_frequency,
        )
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
        carrier_phases[i] = comp_carrier_phase / 2π
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 5e-5
    @test tracked_carrier_phases[end] + 0.5 ≈ carrier_phases[end] atol = 5e-5 / 2π
end

@testset "Track multiple signals with samples of type $type" for type in (
    Int16,
    Int32,
    Int64,
    Float32,
    Float64,
)
    gpsl1 = GPSL1CA()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_chip_rate(gpsl1)
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

    samples =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn,
        )
    signal_mat_temp = repeat(samples; outer = (1, 3))
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
        samples =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn,
            )
        signal_mat = repeat(samples; outer = (1, 3))
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
        carrier_phases[i] = comp_carrier_phase / 2π
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 5e-5
    @test tracked_carrier_phases[end] + 0.5 ≈ carrier_phases[end] atol = 5e-5 / 2π
end

# The testset above has identical, noise-free columns. Here each antenna has
# independent noise at its own level, so reading the wrong floor shows.
@testset "Multi-antenna C/N₀ under the default filter matches its own column" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    prn = 1
    samples_per_code = 5000
    num_ants = 3
    # Antenna m's floor is `scaleₘ²·σ²/f_s`.
    scales = (1.0, 1.5, 2.0)

    function make_array_signal(rng, code_idx)
        code = get_code.(
            gpsl1,
            get_chip_rate(gpsl1) .* (0:(samples_per_code-1)) ./ sampling_frequency .+
            code_idx * get_code_length(gpsl1),
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

    # A single-antenna control is fed the array's last column.
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

    R = get_noise_density(array_state.noise_estimators.GPSL1CA)
    @test size(R) == (num_ants, num_ants)
    # `DefaultPostCorrFilter` picks the last antenna, so the covariance reduces
    # to its floor and the C/N₀ is identical, not merely close.
    @test estimate_cn0(array_state, 1) === estimate_cn0(single_state, 1)
    @test isfinite(ustrip(uconvert(dBHz, estimate_cn0(array_state, 1))))
end

@testset "A beamformer's C/N₀ follows its own weights" begin
    # An equal-weight combiner over independent noise gets `‖w‖²` = 1/3 of the
    # floor while the signal adds coherently: ~10log10(3) dB of gain.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    prn = 1
    samples_per_code = 5000
    num_ants = 3

    struct MeanBeamformer <: AbstractPostCorrFilter end
    TrackingLoops.update(f::MeanBeamformer, prompt) = f
    TrackingLoops.get_weights(::MeanBeamformer, ::NumAnts{1}) = 1.0 + 0.0im
    TrackingLoops.get_weights(::MeanBeamformer, ::NumAnts{M}) where {M} =
        SVector{M,ComplexF64}(ntuple(_ -> 1 / M + 0.0im, M))

    function make_array_signal(rng)
        code = get_code.(
            gpsl1,
            get_chip_rate(gpsl1) .* (0:(samples_per_code-1)) ./ sampling_frequency,
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
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1) +
        get_chip_rate(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    samples_per_code = 4000  # one GPS L1 code period at 4 MHz

    function make_signal(num_codes, start_carrier_phase, start_code_phase)
        n = samples_per_code * num_codes
        samples = zeros(ComplexF64, n)
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
            samples[(code_idx*samples_per_code+1):((code_idx+1)*samples_per_code)] .=
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
        samples
    end

    track_state =
        TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, carrier_doppler)])

    signal1 = make_signal(10, π / 2, start_code_phase)
    track_state = track(signal1, track_state, sampling_frequency)
    sat = get_sat_state(track_state, prn)
    prompts1 = get_filtered_prompts(sat)
    @test prompts1 isa Vector{ComplexF64}
    @test length(prompts1) == 10
    @test prompts1[end] == get_last_fully_integrated_filtered_prompt(track_state)
    buffer_ref = prompts1

    # The buffer is reset per call (5, not 15) and reused, not reallocated.
    signal2 = make_signal(5, π / 2, start_code_phase)
    track_state = track(signal2, track_state, sampling_frequency)
    sat = get_sat_state(track_state, prn)
    prompts2 = get_filtered_prompts(sat)
    @test length(prompts2) == 5
    @test prompts2[end] == get_last_fully_integrated_filtered_prompt(track_state)
    @test prompts2 === buffer_ref

    # No completed integration: the buffer empties.
    signal3 = zeros(ComplexF64, samples_per_code ÷ 4)
    track_state = track(signal3, track_state, sampling_frequency)
    prompts3 = get_filtered_prompts(get_sat_state(track_state, prn))
    @test isempty(prompts3)
    @test prompts3 === buffer_ref
end

# L1C-D and L1C-P, one 10 ms primary period per call, at the default bandwidths
# (18 Hz carrier capped to 9 Hz at 10 ms, 1 Hz code).
@testset "Tracking single signal $name with $type samples" for (name, sig_type) in (
        ("GPSL1C_D", GPSL1C_D),
        ("GPSL1C_P", GPSL1C_P),
    ),
    type in (Float32, Float64)

    signal = sig_type()
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(signal) +
        get_chip_rate(signal)
    # L1C-P's TMBOC needs fs > ~12.28 MHz.
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
    # 1 s from a 5 Hz offset, noise-free.
    final_doppler = get_carrier_doppler(track_state, :default, prn)
    @test abs(final_doppler - carrier_doppler) < 1.0Hz
end

# GPS L1C-D's BOC(1,1) autocorrelation has nulls at ~±0.29 chip, so 0.5-chip
# early/late taps land on the side lobes and the code loop walks off the peak
# once the code starts off-center. The test above starts on the code phase,
# where E == L, so it cannot catch this. Seeded 0.2 chip off at 45 dB-Hz, the
# default correlator holds ~45 dB-Hz and 0.5-chip taps collapse to ~31 dB-Hz;
# the 4 Hz code loop makes that show within ~1.3 s.
@testset "GPS L1C-D BOC code tracking holds lock under a code offset" begin
    Random.seed!(1234)
    signal = GPSL1C_D()
    prn = 1
    sampling_frequency = 10e6Hz
    period = 100000  # 10 ms at 10 MHz — one L1C-D primary code period
    # Baseband (carrier_doppler = 0), so the code frequency is nominal.
    code_frequency = get_chip_rate(signal)
    primary_len = get_code_length(signal)
    amplitude = 10^(45 / 20)            # 45 dB-Hz with noise_std below
    noise_std = sqrt(10e6)
    range = 0:(period-1)

    estimator = ConventionalAssistedPLLAndDLL(;
        wide_carrier_loop_filter_bandwidth = 1.8Hz,
        narrow_carrier_loop_filter_bandwidth = 1.8Hz,
        code_loop_filter_bandwidth = 4.0Hz,
    )
    track_state = TrackState(; signal, doppler_estimator = estimator)
    # A realistic handoff error, on the group's default correlator.
    track_state =
        add_satellite!(track_state; prn, code_phase = 0.2, carrier_doppler = 0.0Hz)

    build(code_phase) =
        gen_code(period, signal, prn, sampling_frequency, code_frequency, code_phase) .*
        amplitude .+ randn(ComplexF64, period) .* noise_std

    for i = 0:130
        code_phase = mod(code_frequency * period * i / sampling_frequency, primary_len)
        track!(build(code_phase), track_state, sampling_frequency)
    end

    @test estimate_cn0(track_state, prn) > 40dBHz
end

# The README's use case: one satellite on L1C_P, L1C_D and L1 C/A, with the loop
# driven by `signals[1]`, the L1C_P pilot.
@testset "Tracking multi-signal (L1C_P, L1C_D, L1CA) on one sat" begin
    # L1C-P's TMBOC needs fs > ~12.28 MHz.
    sampling_frequency = 15e6Hz
    n_samples = 150000  # 10 ms — one L1C-P/D primary period; ten L1CA periods
    prn = 11
    carrier_doppler = 1234.0Hz
    init_offset = 5.0Hz

    track_state =
        TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
    track_state = add_satellite!(
        track_state;
        prn,
        group_key = :modern_gps,
        code_phase = 0.0,
        carrier_doppler = carrier_doppler - init_offset,
    )

    # All three signals on one carrier.
    range = 0:(n_samples-1)
    function build_signal(code_phase, carrier_phase)
        carrier =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase)
        signal_cp = let s = GPSL1C_P()
            cf =
                carrier_doppler * get_chip_rate_to_center_frequency_ratio(s) +
                get_chip_rate(s)
            carrier .* gen_code(n_samples, s, prn, sampling_frequency, cf, code_phase)
        end
        signal_cd = let s = GPSL1C_D()
            cf =
                carrier_doppler * get_chip_rate_to_center_frequency_ratio(s) +
                get_chip_rate(s)
            carrier .* gen_code(n_samples, s, prn, sampling_frequency, cf, code_phase)
        end
        signal_ca = let s = GPSL1CA()
            cf =
                carrier_doppler * get_chip_rate_to_center_frequency_ratio(s) +
                get_chip_rate(s)
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

    track!(build_signal(0.0, 0.0), track_state, sampling_frequency)

    iterations = 200
    cp_primary_len = get_code_length(GPSL1C_P())
    for i = 1:iterations
        carrier_phase =
            mod2pi(2π * carrier_doppler * n_samples * i / sampling_frequency + π) - π
        cf_cp =
            carrier_doppler * get_chip_rate_to_center_frequency_ratio(GPSL1C_P()) +
            get_chip_rate(GPSL1C_P())
        code_phase = mod(cf_cp * n_samples * i / sampling_frequency, cp_primary_len)
        track!(build_signal(code_phase, carrier_phase), track_state, sampling_frequency)
    end

    sat = get_sat_state(track_state, :modern_gps, prn)
    @test length(get_filtered_prompts(sat.signals[1])) > 0  # L1C_P
    @test length(get_filtered_prompts(sat.signals[2])) > 0  # L1C_D
    @test length(get_filtered_prompts(sat.signals[3])) > 0  # L1CA

    final_doppler = get_carrier_doppler(track_state, :modern_gps, prn)
    @test abs(final_doppler - carrier_doppler) < 5.0Hz
end

# A continuous signal in one-code-period chunks with a negative code Doppler: a
# block spans 20001 > 20000 samples, so it always carries into the next chunk.
# A code-phase snap on every call after the NH10 sync would discard that partial
# block and wedge the satellite (#117).
@testset "Does not deadlock on one-code-period chunks (issue #117)" begin
    signal = GPSL5I()
    carrier_doppler = -200Hz                 # negative ⇒ negative code Doppler
    code_frequency =
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(signal) +
        get_chip_rate(signal)
    sampling_frequency = 20e6Hz
    prn = 1
    chunk = 20000                            # 1 ms = one L5I code period at 20 MHz
    num_calls = 60                           # well past the 10-block NH10 sync

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

    @test has_bit_or_secondary_code_been_found(sat)
    # A wedged satellite grows `integrated_samples` without bound and freezes
    # its prompt.
    @test get_integrated_samples(track_state, prn) <= 2 * chunk
    @test length(unique(prompts[(end-9):end])) > 1
    @test abs(get_carrier_doppler(track_state, prn) - carrier_doppler) < 5.0Hz
end

end
