@testset "Galileo E1B" begin
    galileo_e1b = GalileoE1B()
    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf, 29)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf, 6)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xc, 40)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf, 8)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf0, 8)) == true

    sampling_frequency = 5e6Hz

    @test @inferred(
        Tracking.get_default_correlator(galileo_e1b, sampling_frequency, NumAnts(1))
    ) == VeryEarlyPromptLateCorrelator(
        galileo_e1b,
        sampling_frequency;
        num_ants = NumAnts(1),
    )
    @test @inferred(
        Tracking.get_default_correlator(galileo_e1b, sampling_frequency, NumAnts(3))
    ) == VeryEarlyPromptLateCorrelator(
        galileo_e1b,
        sampling_frequency;
        num_ants = NumAnts(3),
    )
end
