module TrackInPlaceTest

using Test: @test, @testset
using Random: MersenneTwister
using Unitful: Hz
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    GalileoE1C,
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
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_sat_state,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator
using TrackingLoops:
    NoiseRefCN0Estimator,
    NWPRCN0Estimator,
    append_noise_observation!,
    noise_observation_from_samples

# 1 ms of GPS L1 C/A PRN 1 with a known carrier Doppler and code phase.
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

# Per-stage allocation measurement, inside functions: `@allocated` at module
# scope counts the boxing of untyped globals. The warmup calls page in Bumper's
# slab buffer and settle the capacity of `filtered_prompts`.

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
    # Only `sampling_frequency` is read.
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

    # An estimator that reads no noise density: the default's per-signal despread
    # has its own test below, gated to Julia >= 1.11.
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

    # Warm up: compile, and settle the bit buffer and `filtered_prompts`.
    for _ = 1:8
        track!(signal, track_state, sampling_frequency; downconvert_and_correlator = dc)
    end

    @test measure_reset!(track_state) == 0
    # The threaded backend keeps a small residual (`@batch`'s
    # ManualMemory.Reference and Bumper checkpoints inside Polyester's closure),
    # capped loosely enough to still catch per-satellite allocations.
    @test measure_est!(track_state, sampling_frequency) == 0
    if DC === CPUDownconvertAndCorrelator
        @test measure_dc!(dc, signal, track_state, sampling_frequency) == 0
    else
        @test measure_dc!(dc, signal, track_state, sampling_frequency) <= 1024
    end
end

# The same with a per-signal noise reference: the walk over `noise_estimators`
# in the correlate step and the per-signal densities in the estimate step. On
# Julia 1.10 the despread leaves a few hundred bytes per sub-integration that
# 1.11 elides.
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

# What the noise reference costs the threaded launch, below the loose cap above.
# Polyester copies what a `@batch` region references into its argument tuple by
# value, so the noise descriptor is reached through one box (a pointer). Asserted
# as equality against a state that reads no density, since the byte count is
# Polyester's to change; two satellites, so the loop has a residual at all.
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

        # The box is provisioned once and survives the copies `track` and
        # `track!` make.
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

        # A state that reads no density never provisions one.
        @test isnothing(plain_state.noise_descriptor[])
    end
end

# Pre-sync allocation guard (#198). The per-stage test above measures its
# estimate step with no records to fold, so never the pre-sync bit-edge search
# (`TrackingLoops._buffer_find_bit`). Noise never locks, so the search runs on
# every block, and a warm call must allocate nothing at any length: a per-block
# leak grows with the chunk.
#
# `_detect_bit_edge_cfar` scores nothing below `2 * blocks_per_bit` (40 blocks),
# so 2 and 20 blocks cannot see a leak in the Student-t threshold. Its
# incomplete-beta inverse once allocated ~600 B per block, at `dof` around 11 and
# from ~42 up. The measured second call sweeps `dof ≈ n/20 … n/10`, so 200 blocks
# covers `dof ≈ 9…19` and 900 blocks `dof ≈ 44…89`.
#
# Gated to Julia ≥ 1.11: 1.10 leaves ~770 B per block in the pre-sync soft-bit
# path.
if VERSION >= v"1.11"
    @testset "track! is allocation-free during acquisition (pre-sync bit search)" begin
        sampling_frequency = 5e6Hz
        samples_per_block = 5000               # 1 ms GPS L1CA code period @ 5 MHz
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
            call()                             # warmup
            @allocated call()
        end
        @test measure_track_alloc(2) == 0
        @test measure_track_alloc(20) == 0
        # Past `2 * blocks_per_bit`: the threshold is scored.
        @test measure_track_alloc(200) == 0
        @test measure_track_alloc(900) == 0
    end

    # The same for the secondary-code detector (GPS L5I, NH10), which shares the
    # threshold but has its own per-rotation accumulators. Both lengths score
    # every block and sweep the same `dof` stretches as above.
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
            call()                             # warmup
            @allocated call()
        end
        @test measure_track_alloc_l5i(200) == 0
        @test measure_track_alloc_l5i(900) == 0
    end

    # And for a pilot whose secondary code starts with a negative chip (Galileo
    # E1C, CS25): its pre-sync prompts are multiplied by that chip before the
    # fold. Both lengths are past `2 × 25` blocks.
    @testset "track! is allocation-free during acquisition (pre-sync pilot)" begin
        sampling_frequency = 16.368e6Hz
        samples_per_block = 65472              # 4 ms Galileo E1C code period
        function measure_track_alloc_e1c(nblocks)
            e1c = GalileoE1C()
            track_state = TrackState(e1c, [TrackedSat(e1c, 1, 10.5, 1000.0Hz)])
            dc = CPUDownconvertAndCorrelator()
            signal = rand(MersenneTwister(nblocks), ComplexF32, samples_per_block * nblocks)
            call() = track!(
                signal,
                track_state,
                sampling_frequency;
                downconvert_and_correlator = dc,
            )
            call()                             # warmup
            @allocated call()
        end
        @test measure_track_alloc_e1c(60) == 0
        @test measure_track_alloc_e1c(120) == 0
    end
end

end
