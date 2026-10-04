module BeiDouB2bTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: BeiDouB2bI, get_band_id, get_secondary_code_length
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "BeiDou B2b-I" begin
    b2b = BeiDouB2bI()
    prn = 6   # the ICD defines ranging codes for PRN 6-58 only

    # One symbol per 1 ms primary period, no secondary code: immediate sync.
    @test get_secondary_code_length(b2b) == 1
    for num_blocks in (1, 2, 33)
        res = @inferred detect_bit_or_secondary_code_sync(b2b, prn, UInt8(0x0), num_blocks)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    # Shares the E5b carrier with Galileo; the band id keys per-band antenna
    # counts and measurements.
    @test get_band_id(b2b) === :E5b

    # BPSK(10) → EarlyPromptLate.
    @test @inferred(get_default_correlator(b2b, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b2b, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Flat defaults, uncapped at 1 ms.
    @test @inferred(default_carrier_loop_filter_bandwidth(b2b)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b2b)) ≈ 1.0Hz

    # Sync buffer is dead state.
    @test @inferred(get_code_block_buffer_type(b2b)) === UInt8

    @test Tracking.uses_soft_secondary_code_detection(b2b) == false
    @test Tracking.uses_soft_bit_edge_detection(b2b) == false
end

end
