module GPSL5QTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL5Q, get_secondary_code_length
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

@testset "GPS L5Q" begin
    gpsl5q = GPSL5Q()
    prn = 1
    N = get_secondary_code_length(gpsl5q)  # 20 (NH20)
    @test N == 20

    # Below the NH20 horizon the detector returns `found = false` without
    # running the sweep; from one full period on it locks.
    @test @inferred(detect_bit_or_secondary_code_sync(gpsl5q, prn, UInt32(0x0), N - 1)).found ==
          false

    @testset "NH20 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt32, gpsl5q, prn)
        for r in (0, 7, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(gpsl5q, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        # Negated polarity = complement within the N-bit window.
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(gpsl5q, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Hamming tolerance" begin
        # 2.5 % over a 20-block window discretizes to exact match
        # (floor(0.025 × 20) = 0) — any single bit-flip rejects.
        reference = Tracking._packed_secondary_code(UInt32, gpsl5q, prn)
        @test detect_bit_or_secondary_code_sync(gpsl5q, prn, reference, N).found == true
        @test detect_bit_or_secondary_code_sync(
            gpsl5q,
            prn,
            reference ⊻ UInt32(0x1),
            N,
        ).found == false
    end

    # L5Q is BPSK on L5 — EarlyPromptLate default, same as L5I.
    @test @inferred(get_default_correlator(gpsl5q, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl5q, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 1 ms primary period (10230 chips at 10.23 Mcps), same as L5I → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl5q)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl5q)) ≈ 1.0Hz

    # 20-block NH20 window fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(gpsl5q)) === UInt32
end

end
