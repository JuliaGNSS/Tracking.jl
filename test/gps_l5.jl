module GPSL5Test

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL5I, GPSL5Q, get_secondary_code_length
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

@testset "GPS L5I" begin
    gpsl5 = GPSL5I()
    prn = 1
    # Packed NH10 is `0x3ca`, the complement of the ICD's `0000110101`, since
    # `_packed_secondary_code` sets a bit per `+1` chip (see its docstring).
    @test Tracking._packed_secondary_code(UInt32, gpsl5, prn) == UInt32(0x3ca)
    res = @inferred(detect_bit_or_secondary_code_sync(gpsl5, prn, UInt32(0x3ca), 50))
    @test res.found == true
    @test res.polarity == +1

    # Buffer not yet at 10 blocks.
    @test @inferred(detect_bit_or_secondary_code_sync(gpsl5, prn, UInt32(0x3ca), 5)).found ==
          false

    # 0x035 is the negated NH10 (matches at negative polarity).
    res = @inferred(detect_bit_or_secondary_code_sync(gpsl5, prn, UInt32(0x035), 10))
    @test res.found == true
    @test res.polarity == -1

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(gpsl5, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl5, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Flat defaults, uncapped at 1 ms. A longer integration after the
    # secondary-code sync is capped at runtime, not by these defaults.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl5)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl5)) ≈ 1.0Hz

    # 20-block sync window (2 × NH10) fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(gpsl5)) === UInt32

    @test Tracking.uses_soft_secondary_code_detection(gpsl5) == true

    @testset "Hamming tolerance" begin
        # floor(0.025 × 10) = 0: exact match only.
        template = Tracking._packed_secondary_code(UInt32, gpsl5, prn)
        @test detect_bit_or_secondary_code_sync(gpsl5, prn, template, 10).found == true
        @test detect_bit_or_secondary_code_sync(gpsl5, prn, template ⊻ UInt32(0x1), 10).found ==
              false
    end
end

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
        # floor(0.025 × 20) = 0: exact match only.
        reference = Tracking._packed_secondary_code(UInt32, gpsl5q, prn)
        @test detect_bit_or_secondary_code_sync(gpsl5q, prn, reference, N).found == true
        @test detect_bit_or_secondary_code_sync(gpsl5q, prn, reference ⊻ UInt32(0x1), N).found ==
              false
    end

    # BPSK → EarlyPromptLate.
    @test @inferred(get_default_correlator(gpsl5q, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl5q, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Flat defaults, uncapped at 1 ms.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl5q)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl5q)) ≈ 1.0Hz

    # 20-block NH20 window fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(gpsl5q)) === UInt32

    @test Tracking.uses_soft_secondary_code_detection(gpsl5q) == true
end

end
