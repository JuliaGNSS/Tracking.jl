module GalileoE1BTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE1B, GalileoE1B_BOC11
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

# Both `GalileoE1B` (full CBOC) and `GalileoE1B_BOC11` (BOC(1,1)
# approximation) share primary code, code length, code rate, data rate,
# and band — only the modulation differs. All tracking-side traits should
# therefore behave identically; the testset is parameterised over both.
@testset "Galileo E1B ($(nameof(typeof(galileo_e1b))))" for galileo_e1b in
                                                            (GalileoE1B(), GalileoE1B_BOC11())

    # E1B broadcasts one I/NAV channel symbol per primary code period
    # (250 sym/s, 4 ms primary period; Galileo OS SIS ICD Tables 11 & 15)
    # — no sub-symbol boundary to find, so the detector reports
    # `found = true` from the start. Polarity ambiguity is resolved
    # downstream by GNSSDecoder.jl via the I/NAV preamble.
    prn = 1
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred is_upcoming_integration_new_bit(galileo_e1b, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(galileo_e1b, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 4 ms primary period (4092 chips at 1.023 Mcps) → BL·T ≈ 0.018 gives
    # 4.5 Hz carrier / 0.25 Hz code.
    @test @inferred(default_carrier_loop_filter_bandwidth(galileo_e1b)) ≈ 4.5Hz
    @test @inferred(default_code_loop_filter_bandwidth(galileo_e1b)) ≈ 0.25Hz

    # 1 symbol = 1 primary period; sync buffer is dead state, but a
    # concrete type is still required.
    @test @inferred(get_code_block_buffer_type(galileo_e1b)) === UInt8
end

end
