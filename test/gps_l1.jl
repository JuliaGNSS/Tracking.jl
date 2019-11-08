@testset "GPS L1" begin
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 0xfffff00000, 50)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 0xfffff, 10)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 0xfffff, 40)) == true

    @test @inferred(Tracking.get_default_correlator(GPSL1, NumAnts(1))) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(GPSL1, NumAnts(3))) == EarlyPromptLateCorrelator(NumAnts(3))
end
