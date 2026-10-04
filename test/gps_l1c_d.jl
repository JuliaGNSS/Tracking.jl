module GPSL1CDTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms
using GNSSSignals: GPSL1C_D
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    effective_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1C-D" begin
    gpsl1c_d = GPSL1C_D()

    # One symbol per 10 ms primary period: `found = true` from the start (see
    # src/gps/l1c_d.jl).
    prn = 1
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl1c_d, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    # BOC(1,1) → VeryEarlyPromptLate (VEML; see src/gps/l1c_d.jl).
    @test @inferred(get_default_correlator(gpsl1c_d, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1c_d, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Flat defaults; at 10 ms the carrier is capped to 9 Hz, the DLL not yet.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1c_d)) ≈ 18.0Hz
    @test effective_carrier_loop_filter_bandwidth(18.0Hz, 10ms) ≈ 9.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_d)) ≈ 1.0Hz

    # Sync buffer is dead state but needs a concrete type.
    @test @inferred(get_code_block_buffer_type(gpsl1c_d)) === UInt8
end

end
