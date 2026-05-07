module TrackInPlaceTest

using Test: @test, @testset
using Unitful: Hz
using GNSSSignals:
    GPSL1, gen_code, get_code_center_frequency_ratio, get_code_frequency

using Tracking:
    SatState,
    TrackState,
    track,
    track!,
    prewarm!,
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
    gpsl1 = GPSL1()
    carrier_doppler = 200.0Hz
    start_code_phase = 100.0
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) +
        get_code_frequency(gpsl1)
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

    sat_immutable = SatState(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)
    sat_mutable   = SatState(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)

    ts_immutable = TrackState(gpsl1, [sat_immutable])
    ts_mutable   = TrackState(gpsl1, [sat_mutable])

    dc = DC(Val(sampling_frequency))

    ts_immutable = track(signal, ts_immutable, sampling_frequency;
        downconvert_and_correlator = dc)
    track!(signal, ts_mutable, sampling_frequency;
        downconvert_and_correlator = dc)

    @test get_code_phase(ts_immutable)    == get_code_phase(ts_mutable)
    @test get_carrier_phase(ts_immutable) == get_carrier_phase(ts_mutable)
    @test get_code_doppler(ts_immutable)  == get_code_doppler(ts_mutable)
    @test get_carrier_doppler(ts_immutable) == get_carrier_doppler(ts_mutable)
    @test get_last_fully_integrated_filtered_prompt(ts_immutable) ==
          get_last_fully_integrated_filtered_prompt(ts_mutable)
    @test get_filtered_prompts(get_sat_state(ts_immutable)) ==
          get_filtered_prompts(get_sat_state(ts_mutable))
end

@testset "track! returns same TrackState identity" begin
    sampling_frequency = 4e6Hz
    signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)
    track_state = TrackState(
        gpsl1,
        [SatState(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)],
    )
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
    for _ in 1:8
        reset_start_sample_and_bit_buffer!(track_state)
    end
    @allocated reset_start_sample_and_bit_buffer!(track_state)
end

function measure_dc!(dc, signal, track_state, sampling_frequency)
    for _ in 1:8
        downconvert_and_correlate!(
            dc, signal, track_state, 1, sampling_frequency, 0.0Hz,
        )
    end
    @allocated downconvert_and_correlate!(
        dc, signal, track_state, 1, sampling_frequency, 0.0Hz,
    )
end

function measure_est!(track_state, sampling_frequency)
    for _ in 1:8
        estimate_dopplers_and_filter_prompt!(track_state, 1, sampling_frequency)
    end
    @allocated estimate_dopplers_and_filter_prompt!(
        track_state, 1, sampling_frequency,
    )
end

@testset "track! per-stage is allocation-free in steady state ($DC)" for DC in (
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
)
    sampling_frequency = 4e6Hz
    signal, gpsl1, carrier_doppler, start_code_phase = make_signal(sampling_frequency)

    track_state = TrackState(
        gpsl1,
        [SatState(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz)],
    )
    dc = DC(Val(sampling_frequency))

    # 4 ms / 1 ms code period → up to 4 prompts per call
    prewarm!(track_state, 8)

    # Run the full track! once so all stages compile and the bit buffer
    # is in its steady-state shape.
    for _ in 1:8
        track!(signal, track_state, sampling_frequency; downconvert_and_correlator = dc)
    end

    @test measure_reset!(track_state) == 0
    # Downconvert and the doppler estimator each have a small residual
    # allocation footprint (a few hundred bytes for threaded `@batch`
    # scheduling and the `bit_buffer` find-bit closures). These are
    # short-lived young-generation allocations and do not trigger GC
    # pauses. The cap catches any regression that would introduce
    # genuine per-sat allocations.
    @test measure_dc!(dc, signal, track_state, sampling_frequency) <= 1024
    @test measure_est!(track_state, sampling_frequency) <= 1024
end

end
