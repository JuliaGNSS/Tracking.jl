module GPSL1CDTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1C_D
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1C-D" begin
    gpsl1c_d = GPSL1C_D()

    # One symbol per 10 ms block: synced from the start (polarity is left to
    # the decoder).
    prn = 1
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl1c_d, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    # BOC(1,1): VEML against side-peak locks.
    @test @inferred(get_default_correlator(gpsl1c_d, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1c_d, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 18 Hz, capped to 9 Hz at filter time for 10 ms; 1 Hz unclamped.
    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(gpsl1c_d)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_d)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(gpsl1c_d)) === UInt8
end

end
