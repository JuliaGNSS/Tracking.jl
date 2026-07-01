module GalileoE5aTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE5aI, GalileoE5aQ, get_secondary_code_length
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

rotl(x::T, r, N) where {T} =
    r == 0 ? x : ((x << r) | (x >> (N - r))) & ((one(T) << N) - one(T))

@testset "Galileo E5a-I" begin
    e5a_i = GalileoE5aI()
    prn = 1
    N = get_secondary_code_length(e5a_i)  # 20 (CS20)
    @test N == 20

    # E5a-I carries F/NAV data at 50 sps; one CS20 period (20 primary blocks)
    # is one symbol, and the CS20 secondary code is the sync feature. Below
    # one full period the detector returns `found = false`.
    @test @inferred(detect_bit_or_secondary_code_sync(e5a_i, prn, UInt32(0x0), N - 1)).found ==
          false

    @testset "CS20 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt32, e5a_i, prn)
        for r in (0, 9, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(e5a_i, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(e5a_i, prn, negated, N)
        @test res.found == true
        @test res.polarity == -1
    end

    # BPSK on L5 → EarlyPromptLate default.
    @test @inferred(get_default_correlator(e5a_i, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # 1 ms primary period (10230 chips at 10.23 Mcps) → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(e5a_i)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e5a_i)) ≈ 1.0Hz

    # 20-block CS20 window fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(e5a_i)) === UInt32
end

@testset "Galileo E5a-Q" begin
    e5a_q = GalileoE5aQ()
    prn = 1
    N = get_secondary_code_length(e5a_q)  # 100 (per-PRN CS100)
    @test N == 100

    @test @inferred(detect_bit_or_secondary_code_sync(e5a_q, prn, UInt128(0x0), N - 1)).found ==
          false

    @testset "CS100 per-PRN search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt128, e5a_q, prn)
        for r in (0, 37, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(e5a_q, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt128) << N) - one(UInt128))
        res = @inferred detect_bit_or_secondary_code_sync(e5a_q, prn, negated, N)
        @test res.found == true
        @test res.polarity == -1
    end

    # A different PRN uses a different CS100 column, so the reference differs.
    @test Tracking._packed_secondary_code(UInt128, e5a_q, 1) !=
          Tracking._packed_secondary_code(UInt128, e5a_q, 2)

    @test @inferred(get_default_correlator(e5a_q, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    @test @inferred(default_carrier_loop_filter_bandwidth(e5a_q)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e5a_q)) ≈ 1.0Hz

    # 100-block CS100 window needs UInt128.
    @test @inferred(get_code_block_buffer_type(e5a_q)) === UInt128
end

end
