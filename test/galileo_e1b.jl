module GalileoE1BTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms
using GNSSSignals: GalileoE1B, GalileoE1B_BOC11
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    effective_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

# CBOC and its BOC(1,1) approximation differ only in modulation, so all
# tracking-side traits must match.
@testset "Galileo E1B ($(nameof(typeof(galileo_e1b))))" for galileo_e1b in (
    GalileoE1B(),
    GalileoE1B_BOC11(),
)

    # One symbol per 4 ms primary period: `found = true` from the start (see
    # src/galileo/e1b.jl).
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

    # Flat defaults; a 4 ms integration is below the carrier cap.
    @test @inferred(default_carrier_loop_filter_bandwidth(galileo_e1b)) ≈ 18.0Hz
    @test effective_carrier_loop_filter_bandwidth(18.0Hz, 4ms) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(galileo_e1b)) ≈ 1.0Hz

    # Sync buffer is dead state but needs a concrete type.
    @test @inferred(get_code_block_buffer_type(galileo_e1b)) === UInt8
end

end
