module GPSL2CTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms, s
using GNSSSignals: GPSL2CM, GPSL2CL, get_band, get_band_id
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    effective_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L2CM" begin
    gpsl2cm = GPSL2CM()
    prn = 1

    # One symbol per 20 ms block: synced from the start (polarity is left to
    # the decoder).
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl2cm, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    @test @inferred(get_default_correlator(gpsl2cm, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # Flat 18 Hz / 1 Hz, capped at filter time to 4.5 Hz / 0.9 Hz for 20 ms.
    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(gpsl2cm)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl2cm)) ≈ 1.0Hz
    @test @inferred(
        effective_code_loop_filter_bandwidth(
            default_code_loop_filter_bandwidth(gpsl2cm),
            20ms,
        )
    ) ≈ 0.9Hz

    @test @inferred(get_code_block_buffer_type(gpsl2cm)) === UInt8

    @test get_band_id(get_band(gpsl2cm)) == :L2
    @test get_band_id(get_band(GPSL2CL())) == :L2
end

@testset "GPS L2CL" begin
    gpsl2cl = GPSL2CL()
    prn = 1

    # A dataless 1.5 s code with no overlay: nothing to sync, one block per
    # integration.
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 10), (UInt8(0xff), 1000))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl2cl, prn, bits, n)
        @test res.found == false
    end

    @test @inferred(get_default_correlator(gpsl2cl, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # Capped at filter time to 0.06 Hz / 0.012 Hz for 1.5 s.
    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(gpsl2cl)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl2cl)) ≈ 1.0Hz
    @test @inferred(
        effective_code_loop_filter_bandwidth(
            default_code_loop_filter_bandwidth(gpsl2cl),
            1.5s,
        )
    ) ≈ 0.012Hz

    @test @inferred(get_code_block_buffer_type(gpsl2cl)) === UInt8
end

end
