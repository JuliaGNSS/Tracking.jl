module GPSL1CPTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1C_P
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1C-P" begin
    gpsl1c_p = GPSL1C_P()

    # L1C-P is the pilot: no navigation data. Step 4 of the sync-detection
    # redesign brings the 1800-chip overlay search backed by
    # `BitIntegers.UInt1800`; until then the detector always reports
    # `found = false` so the inner loop stays at 10 ms primary-code
    # boundaries.
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, UInt64(0x0), 0)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, UInt64(0x1), 1)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, UInt64(0xffffffff), 32)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, typemax(UInt64), 128)).found == false

    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 10 ms primary period at BL·T ≈ 0.018 → 1.8 Hz carrier / 0.1 Hz code.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1c_p)) ≈ 1.8Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_p)) ≈ 0.1Hz

    # Placeholder until step 4 of the sync-detection redesign brings
    # `UInt1800` via BitIntegers. Step 4 will replace this assertion with
    # `=== BitIntegers.UInt1800`.
    @test @inferred(get_code_block_buffer_type(gpsl1c_p)) === UInt64
end

end
