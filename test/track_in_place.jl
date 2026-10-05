module TrackInPlaceTest

using Test: @test, @testset
using Random: MersenneTwister
using Unitful: Hz
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    GalileoE1B,
    gen_code,
    get_code_center_frequency_ratio,
    get_code_frequency

using Tracking:
    TrackedSat,
    TrackState,
    BandMeasurement,
    track,
    track!,
    reset_start_sample_and_bit_buffer!,
    downconvert_and_correlate!,
    estimate_dopplers_and_filter_prompt!,
    get_code_phase,
    get_carrier_phase,
    get_code_doppler,
    get_carrier_doppler,
    get_signal_start_sample,
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_sat_state,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    NumAnts,
    ConventionalAssistedPLLAndDLL,
    NoiseRefCN0Estimator,
    NWPRCN0Estimator,
    append_noise_observation!,
    noise_observation_from_samples,
    CorrelatorOutput,
    append_correlator_output!,
    get_default_correlator,
    update_accumulator

# Build a simple 4 ms GPS-L1 PRN-1 signal with known carrier doppler & code phase.
function make_signal(sampling_frequency)
    gpsl1 = GPSL1CA()
    carrier_doppler = 200.0Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    range = 0:3999
    start_carrier_phase = π / 2
    signal_template =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        gen_code(4000, gpsl1, 1, sampling_frequency, code_frequency, start_code_phase)
    ComplexF32.(signal_template), gpsl1, carrier_doppler, start_code_phase
end

@testset "track! produces same result as track for $DC" for DC in (
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
)
    sampling_frequency = 4e6Hz
    signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)

    sat_immutable = TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)
    sat_mutable = TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)

    ts_immutable = TrackState(gpsl1, [sat_immutable])
    ts_mutable = TrackState(gpsl1, [sat_mutable])

    dc = DC()

    ts_immutable =
        track(signal, ts_immutable, sampling_frequency; downconvert_and_correlator = dc)
    track!(signal, ts_mutable, sampling_frequency; downconvert_and_correlator = dc)

    @test get_code_phase(ts_immutable) == get_code_phase(ts_mutable)
    @test get_carrier_phase(ts_immutable) == get_carrier_phase(ts_mutable)
    @test get_code_doppler(ts_immutable) == get_code_doppler(ts_mutable)
    @test get_carrier_doppler(ts_immutable) == get_carrier_doppler(ts_mutable)
    @test get_last_fully_integrated_filtered_prompt(ts_immutable) ==
          get_last_fully_integrated_filtered_prompt(ts_mutable)
    @test get_filtered_prompts(get_sat_state(ts_immutable)) ==
          get_filtered_prompts(get_sat_state(ts_mutable))
end

@testset "track! returns same TrackState identity" begin
    sampling_frequency = 4e6Hz
    signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)
    track_state =
        TrackState(gpsl1, [TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)])
    returned = track!(signal, track_state, sampling_frequency)
    @test returned === track_state
end

# Per-stage allocation measurement. Measured inside typed helper functions, since
# `@allocated` at module scope picks up boxing from untyped globals. Warmup calls page in
# Bumper.jl's slab buffer and settle `filtered_prompts`' capacity.

function measure_reset!(track_state)
    for _ = 1:8
        reset_start_sample_and_bit_buffer!(track_state)
    end
    @allocated reset_start_sample_and_bit_buffer!(track_state)
end

function measure_dc!(dc, signal, track_state, sampling_frequency)
    measurements = (L1 = BandMeasurement(signal, sampling_frequency),)
    for _ = 1:8
        downconvert_and_correlate!(dc, measurements, track_state)
    end
    @allocated downconvert_and_correlate!(dc, measurements, track_state)
end

function measure_est!(track_state, sampling_frequency)
    # Estimator only reads `sampling_frequency` off the measurement.
    measurements = (L1 = BandMeasurement(ComplexF64[], sampling_frequency),)
    for _ = 1:8
        estimate_dopplers_and_filter_prompt!(track_state, measurements)
    end
    @allocated estimate_dopplers_and_filter_prompt!(track_state, measurements)
end

@testset "track! per-stage is allocation-free in steady state ($DC)" for DC in (
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
)
    sampling_frequency = 4e6Hz
    signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)

    # Pinned to an estimator that reads no noise density, to measure the correlate and
    # estimate stages alone; the noise-referenced default is tested below (Julia >= 1.11).
    track_state = TrackState(
        gpsl1,
        [
            TrackedSat(
                gpsl1,
                1,
                start_code_phase,
                carrier_doppler - 20Hz;
                cn0_estimator = NWPRCN0Estimator(gpsl1),
            ),
        ],
    )
    dc = DC()

    # Warm up: compile, reach steady-state bit buffer and `filtered_prompts` capacity.
    for _ = 1:8
        track!(signal, track_state, sampling_frequency; downconvert_and_correlator = dc)
    end

    @test measure_reset!(track_state) == 0
    # The threaded backend keeps a small Polyester `@batch` launch residual; the loose
    # cap still catches per-sat allocations.
    @test measure_est!(track_state, sampling_frequency) == 0
    if DC === CPUDownconvertAndCorrelator
        @test measure_dc!(dc, signal, track_state, sampling_frequency) == 0
    else
        @test measure_dc!(dc, signal, track_state, sampling_frequency) <= 1024
    end
end

# Same measurement with a per-signal noise reference: `_update_signal_noise!` (correlate
# step) and the per-signal density tuple (estimate step) must not allocate either.
# Julia >= 1.11 only, as in `test/noise_estimators/correlator.jl` (1.10 leaves a small
# per-sub-integration allocation in the despread).
@static if VERSION >= v"1.11"
    @testset "a noise-referenced signal adds no allocation ($DC)" for DC in (
        CPUDownconvertAndCorrelator,
        CPUThreadedDownconvertAndCorrelator,
    )
        sampling_frequency = 4e6Hz
        signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)
        track_state = TrackState(
            gpsl1,
            [
                TrackedSat(
                    gpsl1,
                    1,
                    start_code_phase,
                    carrier_doppler - 20Hz;
                    cn0_estimator = NoiseRefCN0Estimator(),
                ),
            ],
        )
        dc = DC()
        @test keys(track_state.noise_estimators) == (:GPSL1CA,)
        observation = noise_observation_from_samples(4000.0, 4000, sampling_frequency)
        for _ = 1:8
            append_noise_observation!(track_state, observation, :GPSL1CA)
            track!(signal, track_state, sampling_frequency; downconvert_and_correlator = dc)
        end

        @test measure_est!(track_state, sampling_frequency) == 0
        if DC === CPUDownconvertAndCorrelator
            @test measure_dc!(dc, signal, track_state, sampling_frequency) == 0
        else
            @test measure_dc!(dc, signal, track_state, sampling_frequency) <= 1024
        end
    end
end

# The noise reference must add nothing to the threaded `@batch` launch residual, which
# the loose cap above cannot see (Polyester copies by-value captures into its argument
# tuple, so the descriptors are reached through one box; see `_park_noise_items!`).
# Compared against an identical state without one rather than a byte count, which is
# Polyester's to change. Two satellites so the loop has more than one item.
@static if VERSION >= v"1.11"
    @testset "the noise reference adds nothing to the threaded launch residual" begin
        sampling_frequency = 4e6Hz
        signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)
        make_state(make_cn0_estimator) = TrackState(
            gpsl1,
            [
                TrackedSat(
                    gpsl1,
                    prn,
                    start_code_phase,
                    carrier_doppler - 20Hz;
                    cn0_estimator = make_cn0_estimator(),
                ) for prn = 1:2
            ],
        )
        plain_state = make_state(() -> NWPRCN0Estimator(gpsl1))
        noise_state = make_state(NoiseRefCN0Estimator)
        @test keys(plain_state.noise_estimators) == ()
        @test keys(noise_state.noise_estimators) == (:GPSL1CA,)

        plain_dc = CPUThreadedDownconvertAndCorrelator()
        noise_dc = CPUThreadedDownconvertAndCorrelator()
        for _ = 1:8
            track!(
                signal,
                plain_state,
                sampling_frequency;
                downconvert_and_correlator = plain_dc,
            )
            track!(
                signal,
                noise_state,
                sampling_frequency;
                downconvert_and_correlator = noise_dc,
            )
        end

        @test measure_dc!(noise_dc, signal, noise_state, sampling_frequency) ==
              measure_dc!(plain_dc, signal, plain_state, sampling_frequency)

        # The box is provisioned once and must survive the `TrackState` copies made by
        # `track` and `track!`.
        box = noise_state.noise_descriptor[]
        @test box isa Base.RefValue
        track!(
            signal,
            noise_state,
            sampling_frequency;
            downconvert_and_correlator = noise_dc,
        )
        @test noise_state.noise_descriptor[] === box
        @test track(
            signal,
            noise_state,
            sampling_frequency;
            downconvert_and_correlator = noise_dc,
        ).noise_descriptor === noise_state.noise_descriptor

        # No density-reading estimator ⇒ no box.
        @test isnothing(plain_state.noise_descriptor[])
    end
end

# Pre-sync allocation guard (issue #198): before bit sync, `track!` searches for the bit
# edge on every code block, so a leak there scales with signal length. The per-stage
# test above misses it, since its real signal syncs during warmup. Noise never syncs, so
# the search runs on every block; a warm call must allocate nothing at any length.
#
# `_detect_bit_edge_cfar` scores only once `num_blocks >= 2 * blocks_per_bit` (40 for
# L1 C/A), so 2 and 20 blocks cover the search alone and 200 / 900 also cover the CFAR
# statistic and Student-t threshold. The measured (second) call sweeps blocks n+1 … 2n,
# i.e. `dof ≈ n/20 … n/10`; 200 and 900 hit the `dof` ranges where an allocating
# `beta_inc_inv` once slipped past this guard.
#
# Julia >= 1.11 only: 1.10 leaves a per-block allocation in the pre-sync soft-bit path.
if VERSION >= v"1.11"
    @testset "track! is allocation-free during acquisition (pre-sync bit search)" begin
        sampling_frequency = 5e6Hz
        samples_per_block = 5000               # 1 ms GPS L1CA code period @ 5 MHz
        # Fresh pre-sync tracker, warmed once; `nblocks` sets the chunk length.
        function measure_track_alloc(nblocks)
            gpsl1 = GPSL1CA()
            track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 10.5, 1000.0Hz)])
            dc = CPUDownconvertAndCorrelator()
            signal = rand(MersenneTwister(nblocks), ComplexF32, samples_per_block * nblocks)
            call() = track!(
                signal,
                track_state,
                sampling_frequency;
                downconvert_and_correlator = dc,
            )
            call()                             # warmup: seat buffers, compile
            @allocated call()
        end
        @test measure_track_alloc(2) == 0
        @test measure_track_alloc(20) == 0
        # Past `2 * blocks_per_bit`: the CFAR statistic is scored (see above).
        @test measure_track_alloc(200) == 0
        @test measure_track_alloc(900) == 0
    end

    # Same guard for the secondary-code detector (GPS L5I, NH10), whose own accumulator
    # update `_update_secondary_accumulators!` is scored past `2 × 10` blocks; noise
    # never locks, and 200 / 900 blocks hit the same `dof` ranges as above.
    @testset "track! is allocation-free during acquisition (pre-sync secondary-code search)" begin
        sampling_frequency = 25e6Hz
        samples_per_block = 25000              # 1 ms GPS L5I code period @ 25 MHz
        function measure_track_alloc_l5i(nblocks)
            gpsl5i = GPSL5I()
            track_state = TrackState(gpsl5i, [TrackedSat(gpsl5i, 1, 10.5, 1000.0Hz)])
            dc = CPUDownconvertAndCorrelator()
            signal = rand(MersenneTwister(nblocks), ComplexF32, samples_per_block * nblocks)
            call() = track!(
                signal,
                track_state,
                sampling_frequency;
                downconvert_and_correlator = dc,
            )
            call()                             # warmup: seat buffers, compile
            @allocated call()
        end
        @test measure_track_alloc_l5i(200) == 0
        @test measure_track_alloc_l5i(900) == 0
    end
end

# Folding a record normalizes its correlator through a closure. Not specialized on
# it, `apply` dispatched dynamically and allocated per record for a
# five-tap (VEML) correlator, as Galileo E1's default is.
@testset "folding a record allocates nothing, whatever the correlator ($name)" for (
    name,
    signal,
) in (
    ("EarlyPromptLate", GPSL1CA()),
    ("VeryEarlyPromptLate", GalileoE1B()),
)
    sampling_frequency = 16.368e6Hz
    n = 65472
    track_state = TrackState(
        signal,
        [TrackedSat(signal, 11, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator(signal))],
    )
    correlator = get_default_correlator(signal)
    record = CorrelatorOutput(
        update_accumulator(
            correlator,
            map(_ -> (1.0 + 0.1im) * n, correlator.accumulators),
        ),
        n,
        n,
    )
    measurements = (L1 = BandMeasurement(ComplexF64[], sampling_frequency),)
    allocated = map(1:10) do _
        reset_start_sample_and_bit_buffer!(track_state)
        append_correlator_output!(track_state, record, 11)
        @allocated estimate_dopplers_and_filter_prompt!(track_state, measurements)
    end
    @test last(allocated) == 0
end

end
