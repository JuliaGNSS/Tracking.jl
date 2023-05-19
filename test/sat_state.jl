@testset "Satellite state" begin

    gpsl1 = GPSL1()
    sat_state = @inferred SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz)
    @test get_prn(sat_state) == 1
    @test get_code_phase(sat_state) == 10.0
    @test get_code_doppler(sat_state) == 500.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 500.0Hz
    @test get_integrated_samples(sat_state) == 0.0
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_filtered_prompt(sat_state) == complex(0.0, 0.0)
    @test get_sample_of_last_fully_integrated_correlator(sat_state) == -1
    @test Tracking.found(get_secondary_code_or_bit_detector(sat_state)) == false
    @test length(get_prompts_buffer(sat_state)) == 0
    @test length(get_bit_buffer(sat_state)) == 0

    acq = AcquisitionResults(gpsl1, 5, 5e6Hz, 100.0Hz, 524.6, 45.0, 1.0, randn(100,100), -500:100.0:500)
    sat_state = @inferred SatState(acq)
    @test get_prn(sat_state) == 5
    @test get_code_phase(sat_state) == 524.6
    @test get_code_doppler(sat_state) == 100.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 100.0Hz
    @test get_integrated_samples(sat_state) == 0.0
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_filtered_prompt(sat_state) == complex(0.0, 0.0)
    @test get_sample_of_last_fully_integrated_correlator(sat_state) == -1
    @test Tracking.found(get_secondary_code_or_bit_detector(sat_state)) == false
    @test length(get_prompts_buffer(sat_state)) == 0
    @test length(get_bit_buffer(sat_state)) == 0

    sampling_frequency = 5e6Hz

    system_sats_state = @inferred SystemSatsState(
        gpsl1,
        [
            SatState(
                1,
                10.5,
                0.0Hz,
                0.49,
                0.0Hz,
                150,
                EarlyPromptLateCorrelator(complex.(ones(3), ones(3)), -1:1, 3, 2, 1),
                EarlyPromptLateCorrelator(complex.(zeros(3), zeros(3)), -1:1, 3, 2, 1),
                complex(0.0, 0.0),
                14,
                Tracking.SecondaryCodeOrBitDetector(),
                Tracking.PromptsBuffer(20),
                Tracking.BitBuffer()
            ),
            SatState(
                2,
                11.5,
                1.0Hz,
                0.8,
                10.0Hz,
                100,
                EarlyPromptLateCorrelator(complex.(ones(3) * 2, ones(3) * 3), -1:1, 3, 2, 1),
                EarlyPromptLateCorrelator(complex.(zeros(3), zeros(3)), -1:1, 3, 2, 1),
                complex(0.0, 0.0),
                15,
                Tracking.SecondaryCodeOrBitDetector(),
                Tracking.PromptsBuffer(20),
                Tracking.BitBuffer()
            ),
            SatState(
                3,
                12.5,
                2.0Hz,
                -0.5,
                21.0Hz,
                10,
                EarlyPromptLateCorrelator(complex.(ones(3) * 3, ones(3) * 4), -1:1, 3, 2, 1),
                EarlyPromptLateCorrelator(complex.(ones(3), zeros(3)), -1:1, 3, 2, 1),
                complex(0.0, 0.0),
                -1,
                Tracking.SecondaryCodeOrBitDetector(),
                Tracking.PromptsBuffer(20),
                Tracking.BitBuffer()
            )
        ]
    )

    dopplers_and_filtered_prompts = [
        Tracking.DopplersAndFilteredPrompt(200.0Hz, 2.0Hz, complex(1.0, 2.0)),
        Tracking.DopplersAndFilteredPrompt(150.0Hz, 1.5Hz, complex(2.0, 5.0)),
        Tracking.DopplersAndFilteredPrompt(0.0Hz, 0.0Hz, complex(0.0, 0.0)),
    ]

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

    next_system_sats_states = @inferred Tracking.update(
        (system_sats_state,),
        (dopplers_and_filtered_prompts,),
        next_system_sats_sample_params,
    )

    state1 = next_system_sats_states[1].states[1]
    @test get_carrier_doppler(state1) == 200.0Hz
    @test get_code_doppler(state1) == 2.0Hz
    @test get_last_fully_integrated_filtered_prompt(state1) == complex(1.0, 2.0)
    @test get_integrated_samples(state1) == 0
    @test get_correlator(state1).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(state1).accumulators == complex.(ones(3), ones(3))
    @test get_sample_of_last_fully_integrated_correlator(state1) == 4949

    state2 = next_system_sats_states[1].states[2]
    @test get_carrier_doppler(state2) == 150.0Hz
    @test get_code_doppler(state2) == 1.5Hz
    @test get_last_fully_integrated_filtered_prompt(state2) == complex(2.0, 5.0)
    @test get_integrated_samples(state2) == 0
    @test get_correlator(state2).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(state2).accumulators == complex.(ones(3) * 2, ones(3) * 3)
    @test get_sample_of_last_fully_integrated_correlator(state2) == 4944

    state3 = next_system_sats_states[1].states[3]
    @test get_carrier_doppler(state3) == 21.0Hz
    @test get_code_doppler(state3) == 2.0Hz
    @test get_last_fully_integrated_filtered_prompt(state3) == complex(0.0, 0.0)
    @test get_integrated_samples(state3) == 10
    @test get_correlator(state3).accumulators == complex.(ones(3) * 3, ones(3) * 4)
    @test get_last_fully_integrated_correlator(state3).accumulators == complex.(ones(3), zeros(3))
    @test get_sample_of_last_fully_integrated_correlator(state3) == -1

    correlators = [
        EarlyPromptLateCorrelator(complex.(ones(3) * 3, ones(3) * 4), -1:1, 3, 2, 1)
        EarlyPromptLateCorrelator(complex.(ones(3) * 2, ones(3) * 3), -1:1, 3, 2, 1)
        EarlyPromptLateCorrelator(complex.(ones(3) * 1, ones(3) * 2), -1:1, 3, 2, 1)
    ]

    next_system_sats_states = @inferred Tracking.update(
        (system_sats_state,),
        sampling_frequency,
        0.0Hz,
        (correlators,),
        next_system_sats_sample_params,
        5000
    )

    state1 = next_system_sats_states[1].states[1]
    state2 = next_system_sats_states[1].states[2]
    state3 = next_system_sats_states[1].states[3]
    @test get_code_phase(state1) == mod(10.5 + next_system_sats_sample_params[1][1].signal_samples_to_integrate * (get_code_frequency(gpsl1) + 0.0Hz) / sampling_frequency, 1023)
    @test get_code_phase(state2) == mod(11.5 + next_system_sats_sample_params[1][2].signal_samples_to_integrate * (get_code_frequency(gpsl1) + 1.0Hz) / sampling_frequency, 1023)
    @test get_code_phase(state3) == mod(12.5 + next_system_sats_sample_params[1][3].signal_samples_to_integrate * (get_code_frequency(gpsl1) + 2.0Hz) / sampling_frequency, 1023)

    @test get_carrier_phase(state1) ≈ mod2pi(2π * (0.49 + next_system_sats_sample_params[1][1].signal_samples_to_integrate * 0.0Hz / sampling_frequency) + π) - π
    @test get_carrier_phase(state2) ≈ mod2pi(2π * (0.8 + next_system_sats_sample_params[1][2].signal_samples_to_integrate * 10.0Hz / sampling_frequency) + π) - π
    @test get_carrier_phase(state3) ≈ mod2pi(2π * (-0.5 + next_system_sats_sample_params[1][3].signal_samples_to_integrate * 21.0Hz / sampling_frequency) + π) - π

    @test get_correlator(state1) == correlators[1]
    @test get_correlator(state2) == correlators[2]
    @test get_correlator(state3) == correlators[3]

    @test get_sample_of_last_fully_integrated_correlator(state1) == 14 - 5000
    @test get_sample_of_last_fully_integrated_correlator(state2) == 15 - 5000
    @test get_sample_of_last_fully_integrated_correlator(state3) == -1 - 5000

    next_next_system_sats_sample_params = @inferred Tracking.calc_sample_params(
        (system_sats_state,),
        next_system_sats_sample_params,
        5000,
        sampling_frequency,
        1,
    )

    next_system_sats_states = @inferred Tracking.update(
        (system_sats_state,),
        sampling_frequency,
        0.0Hz,
        (correlators,),
        next_next_system_sats_sample_params,
        5000
    )

    state1 = next_system_sats_states[1].states[1]
    state2 = next_system_sats_states[1].states[2]
    state3 = next_system_sats_states[1].states[3]
    @test get_sample_of_last_fully_integrated_correlator(state1) == 14
    @test get_sample_of_last_fully_integrated_correlator(state2) == 15
    @test get_sample_of_last_fully_integrated_correlator(state3) == -1
end