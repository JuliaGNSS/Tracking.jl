module GalileoE1CTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE1C, GalileoE1C_BOC11, get_secondary_code_length
import Tracking
import TrackingLoops
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

rotl(x::T, r, N) where {T} =
    r == 0 ? x : ((x << r) | (x >> (N - r))) & ((one(T) << N) - one(T))

# `GalileoE1C` (CBOC) and `GalileoE1C_BOC11` differ only in modulation.
@testset "Galileo E1C ($(nameof(typeof(galileo_e1c))))" for galileo_e1c in (
    GalileoE1C(),
    GalileoE1C_BOC11(),
)
    prn = 1
    N = get_secondary_code_length(galileo_e1c)  # 25 (CS25)
    @test N == 25

    # The E1 pilot syncs on its CS25 secondary code.
    @test @inferred(
        detect_bit_or_secondary_code_sync(galileo_e1c, prn, UInt32(0x0), N - 1)
    ).found == false

    @testset "CS25 search — clean lock at known phase / polarity" begin
        reference = TrackingLoops._packed_secondary_code(UInt32, galileo_e1c, prn)
        for r in (0, 11, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(galileo_e1c, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(galileo_e1c, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @test @inferred(get_default_correlator(galileo_e1c, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(galileo_e1c, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(galileo_e1c)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(galileo_e1c)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(galileo_e1c)) === UInt32

    @test TrackingLoops.uses_soft_secondary_code_detection(galileo_e1c) == true
end

end
