@testset "CUDA: Galileo E1B" begin
    galileo_e1b = GalileoE1B(use_gpu = Val(true))
    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf, 29)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf, 6)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xc, 40)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf, 8)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(galileo_e1b, 0xf0, 8)) == true

    @test @inferred(Tracking.get_default_correlator(galileo_e1b, NumAnts(1))) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(galileo_e1b, NumAnts(3))) == EarlyPromptLateCorrelator(NumAnts(3))
end
