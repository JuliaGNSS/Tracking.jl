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

    pll_and_dll = @inferred ConventionalPLLAndDLL(
        sat_state,
    )

    @test pll_and_dll.init_carrier_doppler == 500.0Hz
    @test pll_and_dll.init_code_doppler == 500.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test pll_and_dll.carrier_loop_filter == ThirdOrderBilinearLF()
    @test pll_and_dll.code_loop_filter == SecondOrderBilinearLF()
    @test pll_and_dll.carrier_loop_filter_bandwidth == 18.0Hz
    @test pll_and_dll.code_loop_filter_bandwidth == 1.0Hz
    @test pll_and_dll.post_corr_filter == DefaultPostCorrFilter()

    sampling_frequency = 5e6Hz

    gpsl1 = GPSL1()
    system_sats_state = SystemSatsState(
        gpsl1, 
        [
            SatState(
                1,
                0.5,
                1.0Hz,
                0.01,
                100.0Hz,
                150,
                EarlyPromptLateCorrelator(complex.([1000, 2000, 1000], 0.0), -1:1, 3, 2, 1),
                EarlyPromptLateCorrelator(complex.(zeros(3), zeros(3)), -1:1, 3, 2, 1),
                complex(0.0, 0.0),
                14,
                Tracking.SecondaryCodeOrBitDetector(),
                Tracking.PromptsBuffer(20),
                Tracking.BitBuffer()
            ),
        ]
    )

    doppler_estimators = ConventionalPLLsAndDLLs(map(sat_states -> map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states.states), (system_sats_state,)))

    system_sats_sample_params = Tracking.init_sample_params(
        (system_sats_state,),
        1,
    )
    next_system_sats_sample_params = Tracking.calc_sample_params(
        (system_sats_state,),
        system_sats_sample_params,
        5000,
        sampling_frequency,
        1,
    )

    next_doppler_estimators, dopplers_and_filtered_prompts = Tracking.estimate_dopplers_and_filter_prompt(
        doppler_estimators,
        (system_sats_state,),
        next_system_sats_sample_params,
        sampling_frequency
    )

    @test dopplers_and_filtered_prompts[1][1].carrier_doppler == 100.0Hz
    @test dopplers_and_filtered_prompts[1][1].code_doppler == 1.0Hz
    @test dopplers_and_filtered_prompts[1][1].filtered_prompt == 2000 / 150
end