@testset "Tracking state" begin

    carrier_doppler = 100Hz
    code_phase = 100
    gpsl1 = GPSL1()
    state = TrackingState(1, gpsl1, carrier_doppler, code_phase)

    @test @inferred(Tracking.get_prn(state)) == 1
    @test @inferred(Tracking.get_code_phase(state)) == 100
    @test @inferred(Tracking.get_carrier_phase(state)) == 0.0
    @test @inferred(Tracking.get_init_code_doppler(state)) == 100Hz / 1540
    @test @inferred(Tracking.get_init_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_code_doppler(state)) == 100Hz / 1540
    @test @inferred(Tracking.get_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_sc_bit_detector(state)) == SecondaryCodeOrBitDetector()
    @test @inferred(Tracking.get_carrier_loop_filter(state)) == ThirdOrderBilinearLF()
    @test @inferred(Tracking.get_code_loop_filter(state)) == SecondOrderBilinearLF()
    @test @inferred(Tracking.get_prompt_accumulator(state)) == 0.0
    @test @inferred(Tracking.get_integrated_samples(state)) == 0
    @test @inferred(get_post_corr_filter(state)) == Tracking.DefaultPostCorrFilter()

    state = TrackingState(
        1,
        gpsl1,
        carrier_doppler,
        code_phase;
        code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1),
        carrier_phase = 0.0,
        carrier_loop_filter = ThirdOrderBilinearLF(),
        code_loop_filter = SecondOrderBilinearLF(),
        sc_bit_detector = SecondaryCodeOrBitDetector(),
        correlator = Tracking.get_default_correlator(gpsl1, NumAnts(1)),
        integrated_samples = 0,
        prompt_accumulator = zero(ComplexF64)
    )

    @test @inferred(Tracking.get_prn(state)) == 1
    @test @inferred(Tracking.get_code_phase(state)) == 100
    @test @inferred(Tracking.get_carrier_phase(state)) == 0.0
    @test @inferred(Tracking.get_init_code_doppler(state)) == 100Hz / 1540
    @test @inferred(Tracking.get_init_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_code_doppler(state)) == 100Hz / 1540
    @test @inferred(Tracking.get_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_sc_bit_detector(state)) == SecondaryCodeOrBitDetector()
    @test @inferred(Tracking.get_carrier_loop_filter(state)) == ThirdOrderBilinearLF()
    @test @inferred(Tracking.get_code_loop_filter(state)) == SecondOrderBilinearLF()
    @test @inferred(Tracking.get_prompt_accumulator(state)) == 0.0
    @test @inferred(Tracking.get_integrated_samples(state)) == 0

    state = TrackingState(1, gpsl1, carrier_doppler, code_phase, num_ants = NumAnts(2))
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(2))

    acq = AcquisitionResults(GPSL1(), 5, 5e6Hz, 100.0Hz, 524.6, 45.0, randn(100,100), -500:100.0:500)
    state = TrackingState(acq)
    @test @inferred(Tracking.get_prn(state)) == 5
    @test @inferred(Tracking.get_code_phase(state)) == 524.6
    @test @inferred(Tracking.get_carrier_phase(state)) == 0.0
    @test @inferred(Tracking.get_init_code_doppler(state)) == 100Hz / 1540
    @test @inferred(Tracking.get_init_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_code_doppler(state)) == 100Hz / 1540
    @test @inferred(Tracking.get_carrier_doppler(state)) == 100Hz
    @test @inferred(Tracking.get_correlator(state)) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_sc_bit_detector(state)) == SecondaryCodeOrBitDetector()
    @test @inferred(Tracking.get_carrier_loop_filter(state)) == ThirdOrderBilinearLF()
    @test @inferred(Tracking.get_code_loop_filter(state)) == SecondOrderBilinearLF()
    @test @inferred(Tracking.get_prompt_accumulator(state)) == 0.0
    @test @inferred(Tracking.get_integrated_samples(state)) == 0
end
