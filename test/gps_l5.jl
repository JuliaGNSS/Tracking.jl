module GPSL5Test

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL5
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L5" begin
    gpsl5 = GPSL5()
    @test @inferred(is_upcoming_integration_new_bit(gpsl5, 0x35, 50)) == true

    @test @inferred(is_upcoming_integration_new_bit(gpsl5, 0x35, 5)) == false

    @test @inferred(is_upcoming_integration_new_bit(gpsl5, 0x3ca, 10)) == true # 0x3ca == 1111001010

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(gpsl5, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl5, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))
end

end
