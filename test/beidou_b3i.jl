module BeiDouB3ITest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: BeiDouB1I, BeiDouB3I, get_secondary_code_length
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

const MEO_PRN = 6   # MEO/IGSO (D1) — carries NH20
const GEO_PRN = 1   # GEO (D2) — carries no overlay

@testset "BeiDou B3I" begin
    b3i = BeiDouB3I()
    N = get_secondary_code_length(b3i)  # 20 (NH20)
    @test N == 20

    @test @inferred(detect_bit_or_secondary_code_sync(b3i, MEO_PRN, UInt32(0x0), N - 1)).found ==
          false

    @testset "NH20 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt32, b3i, MEO_PRN)
        # B3I uses the same NH20 overlay as B1I (BDS-SIS-ICD-B3I-1.0 §5.2.1).
        @test reference == Tracking._packed_secondary_code(UInt32, BeiDouB1I(), MEO_PRN)
        for r in (0, 13, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(b3i, MEO_PRN, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(b3i, MEO_PRN, negated, N)
        @test res.found == true
        @test res.polarity == -1
    end

    # As on B1I, the GEO satellites carry no overlay — modelled as an all-ones
    # column, which the soft detector handles as a data-bit-edge search.
    @test Tracking._packed_secondary_code(UInt32, b3i, GEO_PRN) ==
          (one(UInt32) << N) - one(UInt32)
    @test Tracking.uses_soft_secondary_code_detection(b3i) == true
    @test Tracking.uses_soft_bit_edge_detection(b3i) == false

    # Plain BPSK (`LOC`) → EarlyPromptLate default.
    @test @inferred(get_default_correlator(b3i, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b3i, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 1 ms primary period (10230 chips at 10.23 Mcps) → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(b3i)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b3i)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(b3i)) === UInt32
end

end
