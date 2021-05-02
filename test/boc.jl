@testset "BOCcos" begin
    system = BOCcos(GPSL1(), 1, 1)
    @test @inferred(Tracking.is_upcoming_integration_new_bit(system, 0xfffff00000, 50)) == true
    @test @inferred(Tracking.is_upcoming_integration_new_bit(system, 0xfffff, 10)) == false
    @test @inferred(Tracking.is_upcoming_integration_new_bit(system, 0xfffff, 40)) == true

    @test @inferred(Tracking.get_default_correlator(system, NumAnts(1))) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(system, NumAnts(3))) == EarlyPromptLateCorrelator(NumAnts(3))

    system = BOCcos(GPSL5(), 1, 1)
    @test @inferred(Tracking.is_upcoming_integration_new_bit(system, 0x35, 50)) == true
    @test @inferred(Tracking.is_upcoming_integration_new_bit(system, 0x35, 5)) == false
    @test @inferred(Tracking.is_upcoming_integration_new_bit(system, 0x3ca, 10)) == true # 0x3ca == 1111001010

    @test @inferred(Tracking.get_default_correlator(system, NumAnts(1))) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(system, NumAnts(3))) == EarlyPromptLateCorrelator(NumAnts(3))

end
