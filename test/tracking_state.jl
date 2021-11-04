@testset "Tracking state" begin

    carrier_doppler = 100Hz
    code_phase = 100
    gpsl1 = GPSL1()
    state = TrackingState(1, gpsl1, carrier_doppler, code_phase)

    @test @inferred(Tracking.get_prn(state)) == 1
    @test @inferred(Tracking.get_code_phase(state)) == 100
    @test @inferred(Tracking.get_carrier_phase(state)) == 0.0
    @test @inferred(Tracking.get_init_code_doppler(state)) == 0Hz
    @test @inferred(Tracking.get_init_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_code_doppler(state)) == 0Hz
    @test @inferred(Tracking.get_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_sc_bit_detector(state)) == SecondaryCodeOrBitDetector()
    @test @inferred(Tracking.get_carrier_doppler_estimator_state(state)) == CostasLoopState(ThirdOrderBilinearLF())
    @test @inferred(Tracking.get_code_doppler_estimator_state(state)) == EarlyPromptLateLoopState(SecondOrderBilinearLF())
    @test @inferred(Tracking.get_prompt_accumulator(state)) == 0.0
    @test @inferred(Tracking.get_integrated_samples(state)) == 0

    state = TrackingState(
        1,
        gpsl1,
        carrier_doppler,
        code_phase;
        code_doppler = 1.0Hz,
        carrier_phase = 0.0,
        carrier_doppler_estimator_state = CostasLoopState(ThirdOrderBilinearLF()),
        code_doppler_estimator_state = EarlyPromptLateLoopState(SecondOrderBilinearLF()),
        sc_bit_detector = SecondaryCodeOrBitDetector(),
        correlator = Tracking.get_default_correlator(gpsl1, NumAnts(1)),
        integrated_samples = 0,
        prompt_accumulator = zero(ComplexF64)
    )

    @test @inferred(Tracking.get_prn(state)) == 1
    @test @inferred(Tracking.get_code_phase(state)) == 100
    @test @inferred(Tracking.get_carrier_phase(state)) == 0.0
    @test @inferred(Tracking.get_init_code_doppler(state)) == 1.0Hz
    @test @inferred(Tracking.get_init_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_code_doppler(state)) == 1.0Hz
    @test @inferred(Tracking.get_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_sc_bit_detector(state)) == SecondaryCodeOrBitDetector()
    @test @inferred(Tracking.get_carrier_doppler_estimator_state(state)) == CostasLoopState(ThirdOrderBilinearLF())
    @test @inferred(Tracking.get_code_doppler_estimator_state(state)) == EarlyPromptLateLoopState(SecondOrderBilinearLF())
    @test @inferred(Tracking.get_prompt_accumulator(state)) == 0.0
    @test @inferred(Tracking.get_integrated_samples(state)) == 0

    state = TrackingState(1, gpsl1, carrier_doppler, code_phase, num_ants = NumAnts(2))
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(2))

end
