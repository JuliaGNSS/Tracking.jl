module GalileoE1BTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE1B
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

@testset "Galileo E1B" begin
    galileo_e1b = GalileoE1B()
    prn = 1
    # 4-on-4-off bit-edge template — matched at positive polarity.
    res = @inferred(is_upcoming_integration_new_bit(galileo_e1b, prn, UInt8(0x0f), 29))
    @test res.found == true
    @test res.polarity == +1

    # Buffer not yet at 8 blocks.
    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, prn, UInt8(0x0f), 6)).found == false

    # Pattern is `0xc = 1100` — only 2 same/2 same, not 4+4.
    @test @inferred(is_upcoming_integration_new_bit(galileo_e1b, prn, UInt8(0xc), 40)).found == false

    res = @inferred(is_upcoming_integration_new_bit(galileo_e1b, prn, UInt8(0x0f), 8))
    @test res.found == true
    @test res.polarity == +1

    res = @inferred(is_upcoming_integration_new_bit(galileo_e1b, prn, UInt8(0xf0), 8))
    @test res.found == true
    @test res.polarity == -1

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 4 ms primary period (4092 chips at 1.023 Mcps) → BL·T ≈ 0.018 gives
    # 4.5 Hz carrier / 0.25 Hz code.
    @test @inferred(default_carrier_loop_filter_bandwidth(galileo_e1b)) ≈ 4.5Hz
    @test @inferred(default_code_loop_filter_bandwidth(galileo_e1b)) ≈ 0.25Hz

    # 8-block sync window fits in a UInt8.
    @test @inferred(get_code_block_buffer_type(galileo_e1b)) === UInt8
end

end
