@testset "Doppler aiding" begin
    gpsl1 = GPSL1()
    init_carrier_doppler = 10Hz
    init_code_doppler = 1Hz
    carrier_freq_update = 2Hz
    code_freq_update = -0.5Hz

    carrier_freq, code_freq = @inferred Tracking.aid_dopplers(
        gpsl1,
        init_carrier_doppler,
        init_code_doppler,
        carrier_freq_update,
        code_freq_update,
    )

    @test carrier_freq == 10Hz + 2Hz
    @test code_freq == 1Hz + 2Hz / 1540 - 0.5Hz
end

@testset "Conventional PLL and DLL" begin
    gpsl1 = GPSL1()
    sat_state = SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz)

    pll_and_dll = @inferred ConventionalPLLAndDLL(500.0Hz, 100.0Hz)

    @test pll_and_dll.init_carrier_doppler == 500.0Hz
    @test pll_and_dll.init_code_doppler == 100.0Hz
    @test pll_and_dll.carrier_loop_filter == ThirdOrderBilinearLF()
    @test pll_and_dll.code_loop_filter == SecondOrderBilinearLF()
    @test pll_and_dll.carrier_loop_filter_bandwidth == 18.0Hz
    @test pll_and_dll.code_loop_filter_bandwidth == 1.0Hz
    @test pll_and_dll.post_corr_filter == DefaultPostCorrFilter()

    sampling_frequency = 5e6Hz

    gpsl1 = GPSL1()

    code_doppler = 1.0Hz
    carrier_doppler = 100.0Hz
    correlator = EarlyPromptLateCorrelator(complex.([1000, 2000, 1000], 0.0), -1:1, 3, 2, 1)
    num_samples = 5000
    prn = 1

    doppler_estimator = ConventionalPLLAndDLL(carrier_doppler, code_doppler)

    sat_state = SatState(
        prn,
        0.5,
        code_doppler,
        0.01,
        carrier_doppler,
        150,
        correlator,
        EarlyPromptLateCorrelator(complex.(zeros(3), zeros(3)), -1:1, 3, 2, 1),
        complex(0.0, 0.0),
        14,
        Tracking.SecondaryCodeOrBitDetector(),
        Tracking.MomentsCN0Estimator(20),
        Tracking.BitBuffer(),
        doppler_estimator,
        nothing,
        NoSatPostProcess(),
    )

    track_state = TrackState(gpsl1, sat_state; num_samples)

    system_sats_sample_params =
        Tracking.init_sample_params(track_state.multiple_system_sats_state, 1)
    next_system_sats_sample_params = Tracking.calc_sample_params(
        track_state.multiple_system_sats_state,
        system_sats_sample_params,
        num_samples,
        sampling_frequency,
        1,
    )

    new_track_state = Tracking.estimate_dopplers_and_filter_prompt(
        track_state,
        next_system_sats_sample_params,
        sampling_frequency,
        0.75ms,
    )

    @test get_carrier_doppler(get_sat_state(new_track_state, prn)) == 100.0Hz
    @test get_code_doppler(get_sat_state(new_track_state, prn)) == 1.0Hz
    @test get_last_fully_integrated_filtered_prompt(get_sat_state(new_track_state, prn)) ==
          0.0

    sat_state = SatState(
        prn,
        0.5,
        code_doppler,
        0.01,
        carrier_doppler,
        4500,
        correlator,
        EarlyPromptLateCorrelator(complex.(zeros(3), zeros(3)), -1:1, 3, 2, 1),
        complex(0.0, 0.0),
        14,
        Tracking.SecondaryCodeOrBitDetector(),
        Tracking.MomentsCN0Estimator(20),
        Tracking.BitBuffer(),
        doppler_estimator,
        nothing,
        NoSatPostProcess(),
    )

    track_state = TrackState(gpsl1, sat_state; num_samples)

    new_track_state = Tracking.estimate_dopplers_and_filter_prompt(
        track_state,
        next_system_sats_sample_params,
        sampling_frequency,
        0.75ms,
    )

    @test get_carrier_doppler(get_sat_state(new_track_state, prn)) == 100.0Hz
    @test get_code_doppler(get_sat_state(new_track_state, prn)) == 1.0Hz
    @test get_last_fully_integrated_filtered_prompt(get_sat_state(new_track_state, prn)) ==
          2000 / 4500
end
