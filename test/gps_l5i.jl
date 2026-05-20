module GPSL5ITest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL5I
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L5" begin
    gpsl5 = GPSL5I()
    @test @inferred(is_upcoming_integration_new_bit(gpsl5, 0x35, 50)) == true

    @test @inferred(is_upcoming_integration_new_bit(gpsl5, 0x35, 5)) == false

    @test @inferred(is_upcoming_integration_new_bit(gpsl5, 0x3ca, 10)) == true # 0x3ca == 1111001010

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(gpsl5, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl5, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # L5I has a 1 ms primary code period (10230 chips at 10.23 Mcps), same
    # as L1 C/A — so the per-signal default BL falls out the same: 18 Hz /
    # 1 Hz. Longer coherent integration is unlocked by secondary-code sync
    # at runtime, not changed by this default.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl5)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl5)) ≈ 1.0Hz

    # 20-block sync window (2 × NH10) fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(gpsl5)) === UInt32
end

end
