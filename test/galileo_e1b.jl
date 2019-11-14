@testset "Galileo E1B" begin
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GalileoE1B, 0xf, 29)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GalileoE1B, 0xf, 6)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GalileoE1B, 0xc, 40)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GalileoE1B, 0xf, 8)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GalileoE1B, 0xf0, 8)) == true

    @test @inferred(Tracking.get_default_correlator(GalileoE1B, NumAnts(1))) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(GalileoE1B, NumAnts(3))) == EarlyPromptLateCorrelator(NumAnts(3))
end
