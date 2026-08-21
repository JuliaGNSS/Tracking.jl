module BeiDouB2aTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: BeiDouB2aI, BeiDouB2aQ, get_band_id, get_secondary_code_length
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

@testset "BeiDou B2a data" begin
    b2a_i = BeiDouB2aI()
    prn = 1
    N = get_secondary_code_length(b2a_i)  # 5
    @test N == 5

    # B2a carries B-CNAV2 at 200 sym/s; one 5-chip secondary period is one
    # symbol. Below one full period the detector returns `found = false`.
    @test @inferred(detect_bit_or_secondary_code_sync(b2a_i, prn, UInt32(0x0), N - 1)).found ==
          false

    @testset "Secondary-code search — clean lock at known phase / polarity" begin
        # The 5-chip code is `00010` (BDS-SIS-ICD-B2a-1.0 §5.2.1), shared
        # across PRNs, packed newest-first.
        reference = Tracking._packed_secondary_code(UInt32, b2a_i, prn)
        @test reference == UInt32(0b00010)
        for r = 0:(N-1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(b2a_i, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(b2a_i, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    # Shares the 1176.45 MHz carrier with GPS L5 / Galileo E5a → band `L5`.
    @test get_band_id(b2a_i) === :L5

    @test @inferred(get_default_correlator(b2a_i, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b2a_i, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 1 ms primary period (10230 chips at 10.23 Mcps) → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(b2a_i)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b2a_i)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(b2a_i)) === UInt32

    # 5 chips is inside the soft detector's 100-chip cap.
    @test Tracking.uses_soft_secondary_code_detection(b2a_i) == true
end

@testset "BeiDou B2a pilot" begin
    b2a_q = BeiDouB2aQ()
    prn = 1
    N = get_secondary_code_length(b2a_q)  # 100
    @test N == 100

    @test @inferred(detect_bit_or_secondary_code_sync(b2a_q, prn, UInt128(0x0), N - 1)).found ==
          false

    @testset "Secondary-code search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt128, b2a_q, prn)
        for r in (0, 44, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(b2a_q, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt128) << N) - one(UInt128))
        res = @inferred detect_bit_or_secondary_code_sync(b2a_q, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @test get_band_id(b2a_q) === :L5

    @test @inferred(get_default_correlator(b2a_q, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    @test @inferred(default_carrier_loop_filter_bandwidth(b2a_q)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b2a_q)) ≈ 1.0Hz

    # 100-block overlay window needs UInt128.
    @test @inferred(get_code_block_buffer_type(b2a_q)) === UInt128

    @test Tracking.uses_soft_secondary_code_detection(b2a_q) == true
end

end
