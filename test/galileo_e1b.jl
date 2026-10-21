module GalileoE1BTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE1B, GalileoE1B_BOC11
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

# `GalileoE1B` (CBOC) and `GalileoE1B_BOC11` differ only in modulation.
@testset "Galileo E1B ($(nameof(typeof(galileo_e1b))))" for galileo_e1b in (
    GalileoE1B(),
    GalileoE1B_BOC11(),
)

    # One symbol per 4 ms block (OS SIS ICD Tables 11 & 15): synced from the
    # start (polarity is left to the decoder).
    prn = 1
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred detect_bit_or_secondary_code_sync(galileo_e1b, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 18 Hz is inside the 0.09 / T cap at 4 ms; 1 Hz unclamped.
    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(galileo_e1b)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(galileo_e1b)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(galileo_e1b)) === UInt8
end

end
