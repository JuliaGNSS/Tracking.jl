module TrackInPlaceTest

using Test: @test, @testset
using Random: MersenneTwister
using Unitful: Hz
using GNSSSignals: GPSL1CA, gen_code, get_code_center_frequency_ratio, get_code_frequency

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
    ConventionalAssistedPLLAndDLL

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

# Per-stage allocation measurement. `@allocated` in module scope will
# pick up boxing overhead from non-typed local lookups, so the work is
# done inside typed helper functions — that matches what `track!` looks
# like when called from real user code (a function with concrete
# argument types) and what BenchmarkTools' `@benchmark $signal $ts ...`
# measures.
#
# Each helper does a few warmup calls first so Bumper.jl's slab buffer
# is paged in and `push!` to `filtered_prompts` has its capacity settled.

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
    # Estimator only reads `sampling_frequency` off the measurement;
    # samples are unused.
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

    track_state =
        TrackState(gpsl1, [TrackedSat(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)])
    dc = DC()

    # Run the full track! several times so all stages compile, the bit
    # buffer is in its steady-state shape, and `filtered_prompts`'
    # capacity is settled (the first call grows it via `push!`).
    for _ = 1:8
        track!(signal, track_state, sampling_frequency; downconvert_and_correlator = dc)
    end

    @test measure_reset!(track_state) == 0
    # Single-threaded backend: every per-stage call is now genuinely
    # allocation-free. The threaded backend still pays a small residual
    # for `@batch`'s ManualMemory.Reference + Bumper SlabCheckpoints
    # that don't elide through Polyester's task closure. Cap loosely so
    # a regression to genuine per-sat allocations would still fire.
    @test measure_est!(track_state, sampling_frequency) == 0
    if DC === CPUDownconvertAndCorrelator
        @test measure_dc!(dc, signal, track_state, sampling_frequency) == 0
    else
        @test measure_dc!(dc, signal, track_state, sampling_frequency) <= 1024
    end
end

# Acquisition (pre-sync) allocation guard. Before bit sync, `track!` runs the
# per-code-block bit-edge search (`_buffer_find_bit`) once per code block. A
# regression there — e.g. a closure-capture `Core.Box` + boxed `SyncResult`
# (~80 B/block) — makes `track!` allocate in proportion to the signal length
# instead of staying allocation-free after the first (buffer-seating) call.
# This is exactly the "allocates per completed integration, scales with signal
# length" symptom reported in #198: measuring at two chunk lengths and asserting
# both are 0 pins the allocation flat across completion counts (a per-block leak
# would grow 10× from 2 to 20 blocks). The box shows up at any runtime thread
# count on the single-threaded backend, so this single-threaded `== 0` catches
# it on CI without needing a multi-threaded run. The
# per-stage test above does NOT catch it: it tracks a real signal that reaches
# bit sync within the warmup, so by the time it measures it exercises only the
# post-sync path, where `_buffer_find_bit` is no longer called.
#
# Feed noise instead — its prompts never form a consistent energy peak, so the
# CFAR bit-edge detector never locks and the pre-sync search runs on every code
# block. After one warmup call (which seats the per-satellite buffers), a warm
# call must allocate nothing, at any signal length. Measuring both a short and a
# long chunk makes a per-block leak fail loudly: with the box it allocated ~80 B
# per block (160 B at 2 blocks, 1600 B at 20), so the 20-block assertion below
# would have caught it. Single-threaded backend so the assertion is a clean
# `== 0` (the threaded backend keeps a small Polyester residual — see above).
#
# Gated to Julia ≥ 1.11. On 1.10 the compiler leaves a per-block allocation in
# the pre-sync soft-bit path (~770 B/block — measured 2144 B at 2 blocks vs
# 15968 B at 20; unrelated to the box, which was ~80 B/block) that 1.11+ elides,
# so `== 0` only holds on 1.11+. The box regression this guards against still
# shows up on 1.11+ (and in the benchmark suite's memory table, which runs on
# the release Julia), so the guard is not lost.
if VERSION >= v"1.11"
    @testset "track! is allocation-free during acquisition (pre-sync bit search)" begin
        sampling_frequency = 5e6Hz
        samples_per_block = 5000               # 1 ms GPS L1CA code period @ 5 MHz
        # Build a fresh pre-sync tracker, warm it once, then measure one more call.
        # `nblocks` sets the chunk length (and thus the pre-sync loop-iteration
        # count); the RNG is seeded per length so the signal is deterministic.
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
    end
end

end
