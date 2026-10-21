module GalileoE5bTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE5aQ, GalileoE5bI, GalileoE5bQ, get_secondary_code_length
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

# Rotate the low `N` bits of `x` left by `r`: a buffer whose upcoming
# integration is `r` secondary chips into the period.
rotl(x::T, r, N) where {T} =
    r == 0 ? x : ((x << r) | (x >> (N - r))) & ((one(T) << N) - one(T))

@testset "Galileo E5b-I" begin
    e5b_i = GalileoE5bI()
    prn = 1
    N = get_secondary_code_length(e5b_i)  # 4 (CS4)
    @test N == 4

    # I/NAV at 250 sym/s: one CS4 period is one symbol.
    @test @inferred(detect_bit_or_secondary_code_sync(e5b_i, prn, UInt32(0x0), N - 1)).found ==
          false

    @testset "CS4 search — clean lock at known phase / polarity" begin
        # CS4 is `1110` (ICD §3.5.1) for every SVID.
        reference = TrackingLoops._packed_secondary_code(UInt32, e5b_i, prn)
        @test reference == UInt32(0b1110)
        for r = 0:(N-1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(e5b_i, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(e5b_i, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Hamming tolerance" begin
        # floor(0.025 × 4) = 0 errors allowed, and no one-flip neighbour of
        # `1110` is a rotation of it or its negation. Bits above `N` are masked.
        reference = TrackingLoops._packed_secondary_code(UInt32, e5b_i, prn)
        @test detect_bit_or_secondary_code_sync(e5b_i, prn, reference, N).found == true
        for bit = 0:(N-1)
            @test detect_bit_or_secondary_code_sync(
                e5b_i,
                prn,
                reference ⊻ (UInt32(1) << bit),
                N,
            ).found == false
        end
    end

    # A standalone BPSK(10) sideband.
    @test @inferred(get_default_correlator(e5b_i, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(e5b_i, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(e5b_i)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e5b_i)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(e5b_i)) === UInt32

    @test TrackingLoops.uses_soft_secondary_code_detection(e5b_i) == true
end

@testset "Galileo E5b-Q" begin
    e5b_q = GalileoE5bQ()
    prn = 1
    N = get_secondary_code_length(e5b_q)  # 100 (CS100_{n+50})
    @test N == 100

    @test @inferred(detect_bit_or_secondary_code_sync(e5b_q, prn, UInt128(0x0), N - 1)).found ==
          false

    @testset "CS100 search — clean lock at known phase / polarity" begin
        reference = TrackingLoops._packed_secondary_code(UInt128, e5b_q, prn)
        for r in (0, 37, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(e5b_q, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt128) << N) - one(UInt128))
        res = @inferred detect_bit_or_secondary_code_sync(e5b_q, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    # E5b-Q uses CS100_51..100, E5a-Q the lower half: different overlays.
    @test TrackingLoops._packed_secondary_code(UInt128, e5b_q, prn) !=
          TrackingLoops._packed_secondary_code(UInt128, GalileoE5aQ(), prn)

    @test @inferred(get_default_correlator(e5b_q, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(e5b_q)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e5b_q)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(e5b_q)) === UInt128

    # At the soft detector's 100-chip cap.
    @test TrackingLoops.uses_soft_secondary_code_detection(e5b_q) == true
end

end
