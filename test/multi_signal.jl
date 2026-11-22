module MultiSignalTest

# The multi-signal code path, on a satellite tracking GPS L1 C/A and Galileo E1B:
# both signals correlate, and only `signals[1]` drives the loop.

using Test: @test, @testset, @test_throws
using GNSSSignals: TMBOC
using TrackingLoopFilters: ThirdOrderAssistedBilinearLF
using Unitful: Hz, dBHz
using GNSSSignals: GPSL1CA, gen_code, get_chip_rate, get_chip_rate_to_center_frequency_ratio
using Tracking:
    TrackedSat,
    TrackedSignal,
    TrackState,
    add_satellite!,
    track,
    get_filtered_prompts,
    get_sat_state,
    get_signal,
    get_correlator_bank,
    get_last_fully_integrated_correlator,
    get_post_corr_filter,
    get_cn0_estimator,
    get_bit_buffer,
    get_num_bits,
    get_integrated_samples,
    get_carrier_doppler,
    get_code_doppler
using TrackingLoops:
    estimate_cn0,
    get_soft_bits,
    is_bit_or_secondary_code_synced,
    ScalarTrackingLoop,
    EarlyPromptLateCorrelatorBank,
    NumAnts,
    DefaultPostCorrFilter,
    NoiseReferencedCN0Estimator,
    BitSynchronizer
using GNSSSignals: GPSL1C_P, GPSL1C_D, GalileoE1B

@testset "Multi-signal track on a single sat" begin
    gpsl1 = GPSL1CA()
    galileo_e1b = GalileoE1B()
    sampling_frequency = 5e6Hz
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)
    start_code_phase = 0.0
    prn = 1
    num_samples = 20_000  # 4 ms at 5 MHz: four L1 C/A code periods, one E1B

    estimator = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    signal_a = TrackedSignal(
        gpsl1;
        num_ants = NumAnts(1),
        correlator_bank = EarlyPromptLateCorrelatorBank(; num_ants = NumAnts(1)),
        post_corr_filter = DefaultPostCorrFilter(),
    )
    signal_b = TrackedSignal(
        galileo_e1b;
        num_ants = NumAnts(1),
        correlator_bank = EarlyPromptLateCorrelatorBank(; num_ants = NumAnts(1)),
        post_corr_filter = DefaultPostCorrFilter(),
    )
    sat = TrackedSat(
        (signal_a, signal_b),
        prn,
        start_code_phase,
        carrier_doppler;
        doppler_estimator = estimator,
    )

    track_state = TrackState(gpsl1, sat; doppler_estimator = estimator)

    range_ = 0:(num_samples-1)
    signal_buf =
        cis.(2π .* carrier_doppler .* range_ ./ sampling_frequency) .* (
            gen_code(
                num_samples,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                start_code_phase,
            ) .+ gen_code(
                num_samples,
                galileo_e1b,
                prn,
                sampling_frequency,
                code_frequency,
                start_code_phase,
            )
        )

    new_track_state = track(signal_buf, track_state, sampling_frequency)

    new_sat = get_sat_state(new_track_state, prn)
    @test length(get_filtered_prompts(new_sat.signals[1])) == 4
    @test length(get_filtered_prompts(new_sat.signals[2])) == 1

    corr_1 = get_last_fully_integrated_correlator(new_sat.signals[1])
    corr_2 = get_last_fully_integrated_correlator(new_sat.signals[2])
    @test corr_1.correlator_outputs != zero(corr_1.correlator_outputs)
    @test corr_2.correlator_outputs != zero(corr_2.correlator_outputs)

    # Each signal's last integration correlates aligned with its own code.
    @test maximum(abs, corr_1.correlator_outputs) > num_samples / 8
    @test maximum(abs, corr_2.correlator_outputs) > num_samples / 2

    # A per-signal accessor needs a selector on a multi-signal sat.
    @test_throws "pass a signal selector" get_correlator_bank(new_sat)
    @test_throws "pass a signal selector" get_cn0_estimator(new_sat)
    @test_throws "pass a signal selector" get_correlator_bank(new_track_state, prn)
    @test get_correlator_bank(new_sat, 1) isa EarlyPromptLateCorrelatorBank
    @test get_correlator_bank(new_sat, 2) isa EarlyPromptLateCorrelatorBank
end

@testset "Multi-signal track over a multi-code-period chunk" begin
    # `gen_code_replica!` writes from `first(sample_range)`, and the multi-signal
    # `downconvert_and_correlate_fused_tuple!` must read from there too, not
    # from 1. Ten code periods put nine integrations past index 1; the
    # driver must match a single-signal reference bit for bit.
    gpsl1 = GPSL1CA()
    galileo_e1b = GalileoE1B()
    sampling_frequency = 5e6Hz
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)
    start_code_phase = 0.0
    prn = 1
    num_samples = 50_000  # 10 ms at 5 MHz => 10 L1 C/A code periods

    range_ = 0:(num_samples-1)
    signal_buf =
        cis.(2π .* carrier_doppler .* range_ ./ sampling_frequency) .* (
            gen_code(
                num_samples,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                start_code_phase,
            ) .+ gen_code(
                num_samples,
                galileo_e1b,
                prn,
                sampling_frequency,
                code_frequency,
                start_code_phase,
            )
        )

    ref_state = TrackState(; signal = gpsl1)
    ref_state =
        add_satellite!(ref_state; prn, code_phase = start_code_phase, carrier_doppler)
    ref_state = track(signal_buf, ref_state, sampling_frequency)
    ref_sat = get_sat_state(ref_state, prn)

    multi_state = TrackState(; signals = (default = (gpsl1, galileo_e1b),))
    multi_state =
        add_satellite!(multi_state; prn, code_phase = start_code_phase, carrier_doppler)
    multi_state = track(signal_buf, multi_state, sampling_frequency)
    multi_sat = get_sat_state(multi_state, prn)

    @test length(get_filtered_prompts(multi_sat.signals[1])) == 10

    ref_corr = get_last_fully_integrated_correlator(ref_sat).correlator_outputs
    multi_corr =
        get_last_fully_integrated_correlator(multi_sat.signals[1]).correlator_outputs
    @test multi_corr ≈ ref_corr rtol = 1e-4
    @test get_carrier_doppler(multi_sat) ≈ get_carrier_doppler(ref_sat) rtol = 1e-6
    @test get_code_doppler(multi_sat) ≈ get_code_doppler(ref_sat) rtol = 1e-6
    @test get_filtered_prompts(multi_sat.signals[1]) ≈ get_filtered_prompts(ref_sat) rtol =
        1e-4

    # The passenger reads the same samples: its 4 ms integrations correlate aligned.
    @test maximum(
        abs,
        get_last_fully_integrated_correlator(multi_sat.signals[2]).correlator_outputs,
    ) > num_samples / 5
end

@testset "Per-signal accessors on multi-signal TrackedSat / TrackState" begin
    # Distinct signal types, so a type selector is unambiguous.
    track_state =
        TrackState(; signals = (modern_gps = (GPSL1C_P(TMBOC), GPSL1C_D(), GPSL1CA()),))
    track_state = add_satellite!(
        track_state;
        prn = 11,
        group_key = :modern_gps,
        code_phase = 0.0,
        carrier_doppler = 1234.0Hz,
    )
    sat = get_sat_state(track_state, :modern_gps, 11)

    @testset "TrackedSat: integer index" begin
        @test get_signal(sat, 1) isa GPSL1C_P
        @test get_signal(sat, 2) isa GPSL1C_D
        @test get_signal(sat, 3) isa GPSL1CA
        @test get_num_bits(sat, 1) == 0
        @test get_num_bits(sat, 2) == 0
        @test get_num_bits(sat, 3) == 0
        @test get_integrated_samples(sat, 1) == 0
        @test is_bit_or_secondary_code_synced(sat, 1) == false
        @test get_cn0_estimator(sat, 1) isa NoiseReferencedCN0Estimator
        @test get_bit_buffer(sat, 1) isa BitSynchronizer
        @test get_post_corr_filter(sat, 1) isa DefaultPostCorrFilter
    end

    @testset "TrackedSat: signal type" begin
        @test get_signal(sat, GPSL1C_P) isa GPSL1C_P
        @test get_signal(sat, GPSL1C_D) isa GPSL1C_D
        @test get_signal(sat, GPSL1CA) isa GPSL1CA
        @test get_num_bits(sat, GPSL1CA) == 0
        @test get_cn0_estimator(sat, GPSL1C_P) isa NoiseReferencedCN0Estimator
    end

    @testset "TrackedSat: error on missing signal type" begin
        @test_throws ArgumentError get_signal(sat, GalileoE1B)
        @test_throws ArgumentError get_num_bits(sat, GalileoE1B)
    end

    @testset "TrackedSat: a signal given twice is rejected" begin
        # GPSL1CA twice: TrackingLoops could not tell the driver from the other.
        gpsl1 = GPSL1CA()
        s_a = TrackedSignal(
            gpsl1;
            num_ants = NumAnts(1),
            correlator_bank = EarlyPromptLateCorrelatorBank(; num_ants = NumAnts(1)),
            post_corr_filter = DefaultPostCorrFilter(),
        )
        s_b = TrackedSignal(
            gpsl1;
            num_ants = NumAnts(1),
            correlator_bank = EarlyPromptLateCorrelatorBank(; num_ants = NumAnts(1)),
            post_corr_filter = DefaultPostCorrFilter(),
        )
        @test_throws "once per satellite" TrackedSat(
            (s_a, s_b),
            1,
            0.0,
            0.0Hz;
            doppler_estimator = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF),
        )
    end

    @testset "TrackState forwarding: (group_key, prn, sig)" begin
        @test get_num_bits(track_state, :modern_gps, 11, 1) == 0
        @test get_num_bits(track_state, :modern_gps, 11, GPSL1CA) == 0
        @test get_cn0_estimator(track_state, :modern_gps, 11, 2) isa
              NoiseReferencedCN0Estimator
        @test get_cn0_estimator(track_state, :modern_gps, 11, GPSL1C_D) isa
              NoiseReferencedCN0Estimator
        @test get_bit_buffer(track_state, :modern_gps, 11, 3) isa BitSynchronizer
        @test isempty(get_soft_bits(track_state, :modern_gps, 11, GPSL1CA))
        @test estimate_cn0(track_state, :modern_gps, 11, 1) == -Inf * dBHz
        @test estimate_cn0(track_state, :modern_gps, 11, GPSL1CA) == -Inf * dBHz
    end

    @testset "estimate_cn0(::TrackedSignal) direct dispatch" begin
        sig = sat.signals[1]
        @test estimate_cn0(sig) == -Inf * dBHz
    end
end

end
