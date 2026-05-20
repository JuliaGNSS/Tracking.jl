module GPSL1CDTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1C_D
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1C-D" begin
    gpsl1c_d = GPSL1C_D()

    # L1C-D broadcasts one channel symbol per primary code period (100 sps,
    # 10 ms primary period). Step 5 of the sync-detection redesign will
    # collapse `is_upcoming_integration_new_bit` to return
    # `SyncResult(true, 0, +1)` unconditionally; for now it stays
    # `found = false` and the inner loop runs at 10 ms boundaries.
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_d, UInt8(0x0), 0)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_d, UInt8(0x1), 1)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_d, UInt8(0xff), 32)).found == false

    @test @inferred(get_default_correlator(gpsl1c_d, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1c_d, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 10 ms primary period at BL·T ≈ 0.018 → 1.8 Hz carrier / 0.1 Hz code.
    # 10× tighter than the L1 C/A default; required for stable 10 ms tracking.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1c_d)) ≈ 1.8Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_d)) ≈ 0.1Hz

    # 1 symbol = 1 primary period; sync buffer is dead state, but a
    # concrete type is still required.
    @test @inferred(get_code_block_buffer_type(gpsl1c_d)) === UInt8
end

end
