module GPSL2CTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms, s
using GNSSSignals: GPSL2CM, GPSL2CL, get_band, get_band_id
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    effective_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    effective_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L2CM" begin
    gpsl2cm = GPSL2CM()
    prn = 1

    # L2CM broadcasts one CNAV symbol per 20 ms L2CM code period (50 sps) —
    # one block per symbol, no sub-symbol boundary to find, so the detector
    # reports `found = true` from the start. Polarity ambiguity is resolved
    # downstream by GNSSDecoder.jl.
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 32))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl2cm, prn, bits, n)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    # L2CM is BPSK — EarlyPromptLate default.
    @test @inferred(get_default_correlator(gpsl2cm, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # Flat defaults, capped at the 20 ms primary period.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl2cm)) ≈ 18.0Hz
    @test @inferred(
        effective_carrier_loop_filter_bandwidth(
            default_carrier_loop_filter_bandwidth(gpsl2cm),
            20ms,
        )
    ) ≈ 4.5Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl2cm)) ≈ 1.0Hz
    @test @inferred(
        effective_code_loop_filter_bandwidth(
            default_code_loop_filter_bandwidth(gpsl2cm),
            20ms,
        )
    ) ≈ 0.9Hz

    # 1 symbol = 1 primary period; sync buffer is dead state, but a concrete
    # type is still required.
    @test @inferred(get_code_block_buffer_type(gpsl2cm)) === UInt8

    # L2C introduces the L2 band; get_band_id maps it into the multi-band
    # `track` measurement keys.
    @test get_band_id(get_band(gpsl2cm)) == :L2
    @test get_band_id(get_band(GPSL2CL())) == :L2
end

@testset "GPS L2CL" begin
    gpsl2cl = GPSL2CL()
    prn = 1

    # L2CL is a dataless pilot with no secondary/overlay code and a single
    # 767250-chip code (1.5 s period). There is no bit and no secondary code
    # to lock, so the detector never reports `found` — the tracker keeps one
    # code block per integration and simply tracks.
    for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 10), (UInt8(0xff), 1000))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl2cl, prn, bits, n)
        @test res.found == false
    end

    @test @inferred(get_default_correlator(gpsl2cl, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # Flat defaults, capped at the 1.5 s primary period.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl2cl)) ≈ 18.0Hz
    @test @inferred(
        effective_carrier_loop_filter_bandwidth(
            default_carrier_loop_filter_bandwidth(gpsl2cl),
            1.5s,
        )
    ) ≈ 0.06Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl2cl)) ≈ 1.0Hz
    @test @inferred(
        effective_code_loop_filter_bandwidth(
            default_code_loop_filter_bandwidth(gpsl2cl),
            1.5s,
        )
    ) ≈ 0.012Hz

    # No sync feature; the search buffer is dead state but a concrete type is
    # still required.
    @test @inferred(get_code_block_buffer_type(gpsl2cl)) === UInt8
end

end
