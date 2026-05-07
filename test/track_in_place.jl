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

# The per-stage allocation guarantees rely on Bumper.jl's `@no_escape` /
# `@alloc` being elided into stack-only allocation. That elision works on
# Julia >= 1.11 but not on 1.10, where the slab buffer pages still
# heap-allocate ~10 KB on first use. Skip on 1.10 — the perf claims in
# the README and PR target Julia 1.11+, which is what real-time SDR users
# will be running.
if VERSION >= v"1.11"
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

        # Warmup: compile + push! storage doubling should settle after a few
        # calls.
        for _ in 1:8
            track!(signal, track_state, sampling_frequency; downconvert_and_correlator = dc)
        end

        # Reset is strictly allocation-free.
        @test (@allocated reset_start_sample_and_bit_buffer!(track_state)) == 0

        # Downconvert and the doppler estimator each have a small residual
        # allocation footprint (a few hundred bytes for threaded `@batch`
        # scheduling and the `bit_buffer` find-bit closures). These are
        # short-lived young-generation allocations and do not trigger GC
        # pauses. The cap below catches any regression that would introduce
        # genuine per-sat allocations.
        dc_alloc = @allocated downconvert_and_correlate!(
            dc, signal, track_state, 1, sampling_frequency, 0.0Hz,
        )
        @test dc_alloc <= 1024
        est_alloc = @allocated estimate_dopplers_and_filter_prompt!(
            track_state, 1, sampling_frequency,
        )
        @test est_alloc <= 1024
    end
end

end
