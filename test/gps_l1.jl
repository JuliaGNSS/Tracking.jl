@testset "GPS L1" begin
    gpsl1 = GPSL1()
    @test @inferred(Tracking.is_upcoming_integration_new_bit(gpsl1, 0xfffff00000, 50)) == true

    @test @inferred(Tracking.is_upcoming_integration_new_bit(gpsl1, 0xfffff, 10)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(gpsl1, 0xfffff, 40)) == true

    sampling_frequency = 5e6Hz

    @test @inferred(Tracking.get_default_correlator(gpsl1, sampling_frequency, NumAnts(1))) ==
        EarlyPromptLateCorrelator(gpsl1, sampling_frequency, num_ants = NumAnts(1))
    @test @inferred(Tracking.get_default_correlator(gpsl1, sampling_frequency, NumAnts(3))) ==
        EarlyPromptLateCorrelator(gpsl1, sampling_frequency, num_ants = NumAnts(3))
end
