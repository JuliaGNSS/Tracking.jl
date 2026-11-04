module MultiSignalTest

# The multi-signal code path, on a satellite tracking GPS L1 C/A twice: both
# signals correlate, and only `signals[1]` drives the loop.

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
    get_beamformed_prompts,
    get_tracked_sat,
    get_signal,
    get_correlator_bank,
    get_last_integration_correlator_bank,
    get_beamformer,
    get_cn0_estimator,
    get_bit_synchronizer,
    get_num_bits,
    get_num_integrated_samples,
    get_carrier_doppler,
    get_code_doppler
using TrackingLoops:
    estimate_cn0,
    get_soft_bits,
    is_bit_or_secondary_code_synced,
    ScalarTrackingLoop,
    EarlyPromptLateCorrelatorBank,
    NumAntennas,
    DefaultBeamformer,
    NoiseReferencedCN0Estimator,
    BitSynchronizer
using GNSSSignals: GPSL1C_P, GPSL1C_D, GalileoE1B

@testset "Multi-signal track on a single sat" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)
    start_code_phase = 0.0
    prn = 1
    num_samples = 5000  # 1 ms at 5 MHz

    tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    signal_a = TrackedSignal(
        gpsl1;
        num_antennas = NumAntennas(1),
        correlator_bank = EarlyPromptLateCorrelatorBank(; num_antennas = NumAntennas(1)),
        beamformer = DefaultBeamformer(),
    )
    signal_b = TrackedSignal(
        gpsl1;
        num_antennas = NumAntennas(1),
        correlator_bank = EarlyPromptLateCorrelatorBank(; num_antennas = NumAntennas(1)),
        beamformer = DefaultBeamformer(),
    )
    sat = TrackedSat(
        (signal_a, signal_b),
        prn,
        start_code_phase,
        carrier_doppler;
        tracking_loop,
    )

    track_state = TrackState(gpsl1, sat; tracking_loop)

    range_ = 0:(num_samples-1)
    signal_buf =
        cis.(2π .* carrier_doppler .* range_ ./ sampling_frequency) .* gen_code(
            num_samples,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            start_code_phase,
        )

    new_track_state = track(signal_buf, track_state, sampling_frequency)

    new_sat = get_tracked_sat(new_track_state, prn)
    @test length(get_beamformed_prompts(new_sat.signals[1])) == 1
    @test length(get_beamformed_prompts(new_sat.signals[2])) == 1

    corr_1 = get_last_integration_correlator_bank(new_sat.signals[1])
    corr_2 = get_last_integration_correlator_bank(new_sat.signals[2])
    @test corr_1.correlator_outputs != zero(corr_1.correlator_outputs)
    @test corr_2.correlator_outputs != zero(corr_2.correlator_outputs)

    @test corr_1.correlator_outputs ≈ corr_2.correlator_outputs

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
    sampling_frequency = 5e6Hz
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)
    start_code_phase = 0.0
    prn = 1
    num_samples = 50_000  # 10 ms at 5 MHz => 10 L1 C/A code periods

    range_ = 0:(num_samples-1)
    signal_buf =
        cis.(2π .* carrier_doppler .* range_ ./ sampling_frequency) .* gen_code(
            num_samples,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            start_code_phase,
        )

    ref_state = TrackState(; signal = gpsl1)
    ref_state =
        add_satellite!(ref_state; prn, code_phase = start_code_phase, carrier_doppler)
    ref_state = track(signal_buf, ref_state, sampling_frequency)
    ref_sat = get_tracked_sat(ref_state, prn)

    multi_state = TrackState(; signals = (default = (gpsl1, gpsl1),))
    multi_state =
        add_satellite!(multi_state; prn, code_phase = start_code_phase, carrier_doppler)
    multi_state = track(signal_buf, multi_state, sampling_frequency)
    multi_sat = get_tracked_sat(multi_state, prn)

    @test length(get_beamformed_prompts(multi_sat.signals[1])) == 10

    ref_corr = get_last_integration_correlator_bank(ref_sat).correlator_outputs
    multi_corr =
        get_last_integration_correlator_bank(multi_sat.signals[1]).correlator_outputs
    @test multi_corr ≈ ref_corr rtol = 1e-4
    @test get_carrier_doppler(multi_sat) ≈ get_carrier_doppler(ref_sat) rtol = 1e-6
    @test get_code_doppler(multi_sat) ≈ get_code_doppler(ref_sat) rtol = 1e-6
    @test get_beamformed_prompts(multi_sat.signals[1]) ≈ get_beamformed_prompts(ref_sat) rtol =
        1e-4

    @test get_last_integration_correlator_bank(multi_sat.signals[2]).correlator_outputs ≈
          multi_corr
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
    sat = get_tracked_sat(track_state, :modern_gps, 11)

    @testset "TrackedSat: integer index" begin
        @test get_signal(sat, 1) isa GPSL1C_P
        @test get_signal(sat, 2) isa GPSL1C_D
        @test get_signal(sat, 3) isa GPSL1CA
        @test get_num_bits(sat, 1) == 0
        @test get_num_bits(sat, 2) == 0
        @test get_num_bits(sat, 3) == 0
        @test get_num_integrated_samples(sat, 1) == 0
        @test is_bit_or_secondary_code_synced(sat, 1) == false
        @test get_cn0_estimator(sat, 1) isa NoiseReferencedCN0Estimator
        @test get_bit_synchronizer(sat, 1) isa BitSynchronizer
        @test get_beamformer(sat, 1) isa DefaultBeamformer
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

    @testset "TrackedSat: error on duplicate signal type" begin
        # GPSL1CA twice: a type selector is ambiguous.
        gpsl1 = GPSL1CA()
        s_a = TrackedSignal(
            gpsl1;
            num_antennas = NumAntennas(1),
            correlator_bank = EarlyPromptLateCorrelatorBank(;
                num_antennas = NumAntennas(1),
            ),
            beamformer = DefaultBeamformer(),
        )
        s_b = TrackedSignal(
            gpsl1;
            num_antennas = NumAntennas(1),
            correlator_bank = EarlyPromptLateCorrelatorBank(;
                num_antennas = NumAntennas(1),
            ),
            beamformer = DefaultBeamformer(),
        )
        dup_sat = TrackedSat(
            (s_a, s_b),
            1,
            0.0,
            0.0Hz;
            tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF),
        )
        # Integer index still works.
        @test get_signal(dup_sat, 1) isa GPSL1CA
        @test get_signal(dup_sat, 2) isa GPSL1CA
        @test_throws ArgumentError get_signal(dup_sat, GPSL1CA)
        @test_throws ArgumentError get_num_bits(dup_sat, GPSL1CA)
    end

    @testset "TrackState forwarding: (group_key, prn, sig)" begin
        @test get_num_bits(track_state, :modern_gps, 11, 1) == 0
        @test get_num_bits(track_state, :modern_gps, 11, GPSL1CA) == 0
        @test get_cn0_estimator(track_state, :modern_gps, 11, 2) isa
              NoiseReferencedCN0Estimator
        @test get_cn0_estimator(track_state, :modern_gps, 11, GPSL1C_D) isa
              NoiseReferencedCN0Estimator
        @test get_bit_synchronizer(track_state, :modern_gps, 11, 3) isa BitSynchronizer
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
