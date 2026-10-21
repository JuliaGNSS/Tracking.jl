module BeiDouB2bTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: BeiDouB2bI, get_band_id, get_secondary_code_length
import Tracking
import TrackingLoops
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "BeiDou B2b-I" begin
    b2b = BeiDouB2bI()
    prn = 6   # the ICD defines ranging codes for PRN 6-58 only

    # One symbol per 1 ms block and no overlay: immediate.
    @test get_secondary_code_length(b2b) == 1
    for num_blocks in (1, 2, 33)
        res = @inferred detect_bit_or_secondary_code_sync(b2b, prn, UInt8(0x0), num_blocks)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    # 1207.14 MHz, named `E5b` for both constellations; measurements are keyed
    # by it.
    @test get_band_id(b2b) === :E5b

    @test @inferred(get_default_correlator(b2b, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b2b, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(b2b)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b2b)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(b2b)) === UInt8

    @test TrackingLoops.uses_soft_secondary_code_detection(b2b) == false
    @test TrackingLoops.uses_soft_bit_edge_detection(b2b) == false
end

end
