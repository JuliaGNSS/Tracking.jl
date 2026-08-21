module GalileoE5bTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE5aQ, GalileoE5bI, GalileoE5bQ, get_secondary_code_length
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

# Rotate the low `N` bits of `x` left by `r` (emulates a prompt buffer whose
# upcoming integration sits `r` secondary chips into the period).
rotl(x::T, r, N) where {T} =
    r == 0 ? x : ((x << r) | (x >> (N - r))) & ((one(T) << N) - one(T))

@testset "Galileo E5b-I" begin
    e5b_i = GalileoE5bI()
    prn = 1
    N = get_secondary_code_length(e5b_i)  # 4 (CS4)
    @test N == 4

    # E5b-I carries I/NAV at 250 sym/s; one CS4 period (4 primary blocks) is
    # one symbol, and CS4 is the sync feature. Below one full period the
    # detector returns `found = false`.
    @test @inferred(detect_bit_or_secondary_code_sync(e5b_i, prn, UInt32(0x0), N - 1)).found ==
          false

    @testset "CS4 search — clean lock at known phase / polarity" begin
        # CS4 is `1110` (hex E, ICD §3.5.1), shared across SVIDs — so the
        # newest-first packed reference is exactly 0b1110.
        reference = Tracking._packed_secondary_code(UInt32, e5b_i, prn)
        @test reference == UInt32(0b1110)
        # All four rotations of `1110` are distinct (and distinct from their
        # negations), so the sweep recovers a unique phase.
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
        # 2.5 % over a 4-block window discretizes to exact match
        # (floor(0.025 × 4) = 0) — any single bit-flip *inside* the 4-bit window
        # rejects (none of the four one-flip neighbours of `1110` is a rotation
        # of it or of its negation). Bits above the window are masked off by the
        # search, so the flip has to land in the low `N`.
        reference = Tracking._packed_secondary_code(UInt32, e5b_i, prn)
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

    # Modelled as a standalone BPSK(10) sideband → EarlyPromptLate default.
    @test @inferred(get_default_correlator(e5b_i, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(e5b_i, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 1 ms primary period (10230 chips at 10.23 Mcps) → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(e5b_i)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e5b_i)) ≈ 1.0Hz

    # 4-block CS4 window; UInt32 matches the other short-secondary signals.
    @test @inferred(get_code_block_buffer_type(e5b_i)) === UInt32

    # CS4 (4 chips) is short enough for the soft, CFAR secondary-code detector.
    @test Tracking.uses_soft_secondary_code_detection(e5b_i) == true
end

@testset "Galileo E5b-Q" begin
    e5b_q = GalileoE5bQ()
    prn = 1
    N = get_secondary_code_length(e5b_q)  # 100 (CS100_{n+50})
    @test N == 100

    @test @inferred(detect_bit_or_secondary_code_sync(e5b_q, prn, UInt128(0x0), N - 1)).found ==
          false

    @testset "CS100 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt128, e5b_q, prn)
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

    # E5b-Q draws the upper half (CS100_51..100) of the same CS100 table
    # Galileo E5a-Q draws the lower half of, so the same SVID must get a
    # different overlay on the two components.
    @test Tracking._packed_secondary_code(UInt128, e5b_q, prn) !=
          Tracking._packed_secondary_code(UInt128, GalileoE5aQ(), prn)

    @test @inferred(get_default_correlator(e5b_q, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # 1 ms primary period, same as E5b-I → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(e5b_q)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e5b_q)) ≈ 1.0Hz

    # 100-block CS100 window needs UInt128.
    @test @inferred(get_code_block_buffer_type(e5b_q)) === UInt128

    # CS100 (100 chips) is at the soft detector's length cap, so still soft.
    @test Tracking.uses_soft_secondary_code_detection(e5b_q) == true
end

end
