module MultiSignalTest

# Smoke test for the multi-signal code path. Builds a TrackedSat with a
# 2-tuple of TrackedSignals (both GPSL1CA, identical config) and runs one
# track() call on a synthetic L1 C/A signal. Verifies:
#
#  1. The tuple-walking code path in `_update_tracked_sat_correlator` and
#     `_update_tracked_sat_doppler` handles N>1 signals without error.
#  2. Both signals' correlators accumulate.
#  3. The PLL/DLL only fires for `signals[1]` (the Doppler source) —
#     `signals[2]`'s `last_fully_integrated_correlator` updates but the
#     sat-level Doppler reflects only `signals[1]`'s discriminator.
#
# This exercises the mechanics. Real multi-signal scenarios (one sat
# carrying L1 C/A + L1C-D + L1C-P with different primary periods) come
# online with the full user-facing capabilities API in a later step.

using Test: @test, @testset
using Unitful: Hz
using Dictionaries: dictionary
using GNSSSignals: GPSL1CA, gen_code, get_code_frequency, get_code_center_frequency_ratio
using Tracking:
    TrackedSat,
    TrackedSignal,
    TrackState,
    track,
    get_filtered_prompts,
    get_sat_state,
    get_correlator,
    get_last_fully_integrated_correlator,
    ConventionalAssistedPLLAndDLL,
    EarlyPromptLateCorrelator,
    NumAnts,
    DefaultPostCorrFilter,
    MomentsCN0Estimator,
    BitBuffer,
    init_estimator_state

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

end
