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

    # L1C-D broadcasts one channel symbol per primary code period
    # (100 sps, 10 ms primary period) — no sub-symbol boundary to find,
    # so the detector reports `found = true` from the start. Polarity
    # ambiguity is resolved downstream by GNSSDecoder.jl via the CNAV-2
    # preamble.
    prn = 1
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred is_upcoming_integration_new_bit(gpsl1c_d, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

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
