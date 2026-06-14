module MultiSignalTest

# Smoke test for the multi-signal code path. Builds a TrackedSat with a
# 2-tuple of TrackedSignals (both GPSL1CA, identical config) and runs one
# track() call on a synthetic L1 C/A signal. Verifies:
#
#  1. The tuple-walking code path in `_update_tracked_sat_correlator` and
#     `_update_tracked_sat_doppler` handles N>1 signals without error.
#  2. Both signals' correlators accumulate.
#  3. The PLL/DLL only fires for `signals[1]` (the estimator-driver signal) —
#     `signals[2]`'s `last_fully_integrated_correlator` updates but the
#     sat-level Doppler reflects only `signals[1]`'s discriminator.
#
# This exercises the mechanics. Real multi-signal scenarios (one sat
# carrying L1 C/A + L1C-D + L1C-P with different primary periods) come
# online with the full user-facing capabilities API in a later step.

using Test: @test, @testset, @test_throws
using Unitful: Hz, dBHz
using Dictionaries: dictionary
using GNSSSignals: GPSL1CA, gen_code, get_code_frequency, get_code_center_frequency_ratio
using Tracking:
    TrackedSat,
    TrackedSignal,
    TrackState,
    add_satellite!,
    track,
    estimate_cn0,
    get_filtered_prompts,
    get_sat_state,
    get_signal,
    get_correlator,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_filtered_prompt,
    get_post_corr_filter,
    get_cn0_estimator,
    get_bit_buffer,
    get_bits,
    get_num_bits,
    get_integrated_samples,
    get_carrier_doppler,
    get_code_doppler,
    has_bit_or_secondary_code_been_found,
    ConventionalAssistedPLLAndDLL,
    EarlyPromptLateCorrelator,
    NumAnts,
    DefaultPostCorrFilter,
    MomentsCN0Estimator,
    BitBuffer,
    init_estimator_state
using GNSSSignals: GPSL1C_P, GPSL1C_D, GalileoE1B

@testset "Multi-signal track on a single sat" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)
    start_code_phase = 0.0
    prn = 1
    num_samples = 5000  # 1 ms at 5 MHz

    # Build a 2-signal sat: both signals are GPSL1CA. The library doesn't
    # yet expose a public multi-signal constructor (that's step 5), so we
    # assemble the TrackedSat manually.
    estimator = ConventionalAssistedPLLAndDLL()
    signal_a = TrackedSignal(
        gpsl1;
        num_ants = NumAnts(1),
        correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
        post_corr_filter = DefaultPostCorrFilter(),
    )
    signal_b = TrackedSignal(
        gpsl1;
        num_ants = NumAnts(1),
        correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
        post_corr_filter = DefaultPostCorrFilter(),
    )
    bare_sat = TrackedSat(
        prn,
        float(start_code_phase),
        float(code_doppler),
        0.0,
        float(carrier_doppler),
        1,
        (signal_a, signal_b),
        nothing,
    )
    doppler_state = init_estimator_state(estimator, bare_sat)
    sat = TrackedSat(
        bare_sat.prn, bare_sat.code_phase, bare_sat.code_doppler,
        bare_sat.carrier_phase, bare_sat.carrier_doppler,
        bare_sat.signal_start_sample, bare_sat.signals, doppler_state,
    )

    track_state = TrackState(gpsl1, sat; doppler_estimator = estimator)

    # Synthetic L1 C/A signal — same as track tests use.
    range_ = 0:(num_samples - 1)
    signal_buf =
        cis.(2π .* carrier_doppler .* range_ ./ sampling_frequency) .*
        gen_code(num_samples, gpsl1, prn, sampling_frequency, code_frequency, start_code_phase)

    new_track_state = track(signal_buf, track_state, sampling_frequency)

    # Both signals saw a completed integration this call (1 ms = 1 code
    # period of L1 C/A).
    new_sat = get_sat_state(new_track_state, prn)
    @test length(get_filtered_prompts(new_sat.signals[1])) == 1
    @test length(get_filtered_prompts(new_sat.signals[2])) == 1

    # Both signals' last_fully_integrated_correlator received non-zero
    # values (they actually correlated against the input).
    corr_1 = get_last_fully_integrated_correlator(new_sat.signals[1])
    corr_2 = get_last_fully_integrated_correlator(new_sat.signals[2])
    @test corr_1.accumulators != zero(corr_1.accumulators)
    @test corr_2.accumulators != zero(corr_2.accumulators)

    # Identical signals + identical correlator config => identical results.
    @test corr_1.accumulators ≈ corr_2.accumulators
end

@testset "Multi-signal track over a multi-code-period chunk" begin
    # Regression test for the tile-share tuple kernel's code-replica
    # indexing. `gen_code_replica!` writes its first sample at the
    # absolute buffer index `start_sample`; the in-register single-signal
    # kernel reads the replica at that same absolute offset, but the
    # multi-signal `downconvert_and_correlate_fused_tuple!` kernel used to
    # read it from index 1. The two agree only for the *first* integration
    # of a chunk (`start_sample == 1`) — every later integration within
    # the same chunk correlated against the wrong code, decorrelating the
    # loop. The single-code-period chunk used by the smoke test above
    # never exercises `start_sample > 1`, so it missed this entirely.
    #
    # Here we feed a chunk spanning ten L1 C/A code periods (so nine
    # integrations run at `start_sample > 1`) and require the multi-signal
    # estimator-driver signal to match a single-signal reference bit for
    # bit.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)
    start_code_phase = 0.0
    prn = 1
    num_samples = 50_000  # 10 ms at 5 MHz => 10 L1 C/A code periods

    range_ = 0:(num_samples - 1)
    signal_buf =
        cis.(2π .* carrier_doppler .* range_ ./ sampling_frequency) .*
        gen_code(num_samples, gpsl1, prn, sampling_frequency, code_frequency, start_code_phase)

    # Single-signal reference.
    ref_state = TrackState(; signal = gpsl1)
    ref_state = add_satellite!(ref_state; prn, code_phase = start_code_phase, carrier_doppler)
    ref_state = track(signal_buf, ref_state, sampling_frequency)
    ref_sat = get_sat_state(ref_state, prn)

    # Two-signal sat, both GPSL1CA with identical config.
    multi_state = TrackState(; signals = (default = (gpsl1, gpsl1),))
    multi_state = add_satellite!(multi_state; prn, code_phase = start_code_phase, carrier_doppler)
    multi_state = track(signal_buf, multi_state, sampling_frequency)
    multi_sat = get_sat_state(multi_state, prn)

    # Ten completed integrations for the driver signal.
    @test length(get_filtered_prompts(multi_sat.signals[1])) == 10

    # The estimator-driver signal must match the single-signal reference:
    # same correlator, same shared carrier/code Doppler, same prompts.
    ref_corr = get_last_fully_integrated_correlator(ref_sat).accumulators
    multi_corr = get_last_fully_integrated_correlator(multi_sat.signals[1]).accumulators
    @test multi_corr ≈ ref_corr rtol = 1e-4
    @test get_carrier_doppler(multi_sat) ≈ get_carrier_doppler(ref_sat) rtol = 1e-6
    @test get_code_doppler(multi_sat) ≈ get_code_doppler(ref_sat) rtol = 1e-6
    @test get_filtered_prompts(multi_sat.signals[1]) ≈ get_filtered_prompts(ref_sat) rtol = 1e-4

    # Both signals in the multi sat are identical => identical correlators.
    @test get_last_fully_integrated_correlator(multi_sat.signals[2]).accumulators ≈ multi_corr
end

@testset "Per-signal accessors on multi-signal TrackedSat / TrackState" begin
    # Build a TrackState with one 3-signal sat (GPSL1C_P, GPSL1C_D, GPSL1CA).
    # The signal types are distinct, so the type-based selector form has no
    # collisions.
    track_state = TrackState(;
        signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),),
    )
    track_state = add_satellite!(track_state;
        prn = 11, group = :modern_gps,
        code_phase = 0.0, carrier_doppler = 1234.0Hz,
    )
    sat = get_sat_state(track_state, :modern_gps, 11)

    @testset "TrackedSat: integer index" begin
        @test get_signal(sat, 1) isa GPSL1C_P
        @test get_signal(sat, 2) isa GPSL1C_D
        @test get_signal(sat, 3) isa GPSL1CA
        # Per-signal getters: the value happens to match the no-selector
        # form for signal[1] (the `only`-equivalent slot when there's one
        # signal), but here we exercise all three indices to prove the
        # dispatch works on multi-signal sats.
        @test get_num_bits(sat, 1) == 0
        @test get_num_bits(sat, 2) == 0
        @test get_num_bits(sat, 3) == 0
        @test get_integrated_samples(sat, 1) == 0
        @test has_bit_or_secondary_code_been_found(sat, 1) == false
        @test get_cn0_estimator(sat, 1) isa MomentsCN0Estimator
        @test get_bit_buffer(sat, 1) isa BitBuffer
        @test get_post_corr_filter(sat, 1) isa DefaultPostCorrFilter
    end

    @testset "TrackedSat: signal type" begin
        @test get_signal(sat, GPSL1C_P) isa GPSL1C_P
        @test get_signal(sat, GPSL1C_D) isa GPSL1C_D
        @test get_signal(sat, GPSL1CA)  isa GPSL1CA
        @test get_num_bits(sat, GPSL1CA) == 0
        @test get_cn0_estimator(sat, GPSL1C_P) isa MomentsCN0Estimator
    end

    @testset "TrackedSat: error on missing signal type" begin
        @test_throws ArgumentError get_signal(sat, GalileoE1B)
        @test_throws ArgumentError get_num_bits(sat, GalileoE1B)
    end

    @testset "TrackedSat: error on duplicate signal type" begin
        # Build a sat that tracks GPSL1CA twice (different correlator slots).
        # The type-based selector should refuse to disambiguate.
        gpsl1 = GPSL1CA()
        s_a = TrackedSignal(gpsl1; num_ants = NumAnts(1),
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            post_corr_filter = DefaultPostCorrFilter())
        s_b = TrackedSignal(gpsl1; num_ants = NumAnts(1),
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            post_corr_filter = DefaultPostCorrFilter())
        bare = TrackedSat(1, 0.0, 0.0Hz, 0.0, 0.0Hz, 1, (s_a, s_b), nothing)
        est = ConventionalAssistedPLLAndDLL()
        de  = init_estimator_state(est, bare)
        dup_sat = TrackedSat(bare.prn, bare.code_phase, bare.code_doppler,
            bare.carrier_phase, bare.carrier_doppler, bare.signal_start_sample,
            bare.signals, de)
        # Integer index still works.
        @test get_signal(dup_sat, 1) isa GPSL1CA
        @test get_signal(dup_sat, 2) isa GPSL1CA
        # Type selector errors — pointing the user at integer indexing.
        @test_throws ArgumentError get_signal(dup_sat, GPSL1CA)
        @test_throws ArgumentError get_num_bits(dup_sat, GPSL1CA)
    end

    @testset "TrackState forwarding: (group, prn, sig)" begin
        @test get_num_bits(track_state, :modern_gps, 11, 1) == 0
        @test get_num_bits(track_state, :modern_gps, 11, GPSL1CA) == 0
        @test get_cn0_estimator(track_state, :modern_gps, 11, 2) isa MomentsCN0Estimator
        @test get_cn0_estimator(track_state, :modern_gps, 11, GPSL1C_D) isa MomentsCN0Estimator
        @test get_bit_buffer(track_state, :modern_gps, 11, 3) isa BitBuffer
        @test get_bits(track_state, :modern_gps, 11, GPSL1CA) == 0
        @test estimate_cn0(track_state, :modern_gps, 11, 1) == 0.0dBHz
        @test estimate_cn0(track_state, :modern_gps, 11, GPSL1CA) == 0.0dBHz
    end

    @testset "estimate_cn0(::TrackedSignal) direct dispatch" begin
        # The per-signal entry point that didn't exist before.
        sig = sat.signals[1]
        @test estimate_cn0(sig) == 0.0dBHz
    end
end

end
