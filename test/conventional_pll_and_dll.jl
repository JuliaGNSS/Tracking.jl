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
    pll_and_dll = @inferred ConventionalPLLAndDLL(500.0Hz, 100.0Hz)

    @test pll_and_dll.init_carrier_doppler == 500.0Hz
    @test pll_and_dll.init_code_doppler == 100.0Hz
    @test pll_and_dll.carrier_loop_filter == ThirdOrderBilinearLF()
    @test pll_and_dll.code_loop_filter == SecondOrderBilinearLF()
    @test pll_and_dll.carrier_loop_filter_bandwidth == 18.0Hz
    @test pll_and_dll.code_loop_filter_bandwidth == 1.0Hz

    sampling_frequency = 5e6Hz

    gpsl1 = GPSL1()

    carrier_doppler = 100.0Hz
    correlator = EarlyPromptLateCorrelator(complex.([1000, 2000, 1000], 0.0), -1:1, 3, 2, 1)
    prn = 1

    doppler_estimator = ConventionalPLLAndDLL(carrier_doppler, 10Hz)
    code_phase = 0.5
    preferred_num_code_blocks_to_integrate = 1

    sat_state = SatState(
        gpsl1,
        prn,
        sampling_frequency,
        code_phase,
        carrier_doppler;
        doppler_estimator,
    )

    # Number of samples too small to generate a new estimate for phases and dopplers
    num_samples = 2000
    track_state = TrackState(gpsl1, sat_state; num_samples)

    new_track_state = Tracking.estimate_dopplers_and_filter_prompt(
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    # Since number of samples is too small the state doesn't change
    @test get_carrier_doppler(get_sat_state(new_track_state, prn)) == carrier_doppler
    @test get_code_doppler(get_sat_state(new_track_state, prn)) ==
          get_code_center_frequency_ratio(gpsl1) * carrier_doppler
    @test get_last_fully_integrated_filtered_prompt(get_sat_state(new_track_state, prn)) ==
          0.0

    # This time it is large enough to produce new dopplers and phases
    num_samples = 5000
    track_state = TrackState(gpsl1, sat_state; num_samples)

    new_track_state = Tracking.estimate_dopplers_and_filter_prompt(
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    @test get_carrier_doppler(get_sat_state(new_track_state, prn)) == 100.0Hz
    @test get_code_doppler(get_sat_state(new_track_state, prn)) ==
          get_code_center_frequency_ratio(gpsl1) * carrier_doppler
    @test get_last_fully_integrated_filtered_prompt(get_sat_state(new_track_state, prn)) ==
          0.0
end
