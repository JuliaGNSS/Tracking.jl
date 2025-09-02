module GPSL1Test

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1" begin
    gpsl1 = GPSL1()
    @test @inferred(is_upcoming_integration_new_bit(gpsl1, 0xfffff00000, 50)) == true

    @test @inferred(is_upcoming_integration_new_bit(gpsl1, 0xfffff, 10)) == false

    @test @inferred(is_upcoming_integration_new_bit(gpsl1, 0xfffff, 40)) == true

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(gpsl1, sampling_frequency, NumAnts(1))) ==
          EarlyPromptLateCorrelator(gpsl1, sampling_frequency; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1, sampling_frequency, NumAnts(3))) ==
          EarlyPromptLateCorrelator(gpsl1, sampling_frequency; num_ants = NumAnts(3))
end

end
