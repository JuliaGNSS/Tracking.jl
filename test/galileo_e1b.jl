module GalileoE1BTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE1B
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    VeryEarlyPromptLateCorrelator,
    NumAnts

@testset "Galileo E1B" begin
    galileo_e1b = GalileoE1B()
    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, 0xf, 29)) == true

    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, 0xf, 6)) == false

    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, 0xc, 40)) == false

    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, 0xf, 8)) == true

    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, 0xf0, 8)) == true

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))
end

end
