@testset "GPS L5" begin
    gpsl5 = GPSL5()
    @test @inferred(Tracking.is_upcoming_integration_new_bit(gpsl5, 0x35, 50)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(gpsl5, 0x35, 5)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(gpsl5, 0x3ca, 10)) == true # 0x3ca == 1111001010

    @test @inferred(Tracking.get_default_correlator(gpsl5, NumAnts(1))) == EarlyPromptLateCorrelator(NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(gpsl5, NumAnts(3))) == EarlyPromptLateCorrelator(NumAnts(3))
end
