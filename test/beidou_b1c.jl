module BeiDouB1CTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using Random: MersenneTwister, randperm
using GNSSSignals: BeiDouB1C_D, BeiDouB1C_P, get_band_id, get_secondary_code_length
import Tracking
import TrackingLoops
using TrackingLoops:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    get_bit_edge_or_secondary_code_tolerance,
    VeryEarlyPromptLateCorrelator,
    NumAnts

const B1C_P_MAX_ERRORS =
    floor(Int, get_bit_edge_or_secondary_code_tolerance(BeiDouB1C_P()) * 1800)

@testset "BeiDou B1C data" begin
    b1c_d = BeiDouB1C_D()
    prn = 1

    # One symbol per 10 ms block and no overlay: immediate, as GPS L1C-D.
    @test get_secondary_code_length(b1c_d) == 1
    for num_blocks in (1, 2, 7)
        res =
            @inferred detect_bit_or_secondary_code_sync(b1c_d, prn, UInt8(0x0), num_blocks)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    @test get_band_id(b1c_d) === :L1

    # BOC(1,1): VEML against side-peak locks.
    @test @inferred(get_default_correlator(b1c_d, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b1c_d, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 18 Hz, capped to 9 Hz at filter time for 10 ms; 1 Hz unclamped.
    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(b1c_d)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b1c_d)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(b1c_d)) === UInt8

    @test TrackingLoops.uses_soft_secondary_code_detection(b1c_d) == false
    @test TrackingLoops.uses_soft_bit_edge_detection(b1c_d) == false
end

@testset "BeiDou B1C pilot" begin
    b1c_p = BeiDouB1C_P()
    prn = 1
    N = get_secondary_code_length(b1c_p)  # 1800
    @test N == 1800

    # An 18 s overlay, as GPS L1C-P: the hard rotation sweep.
    @test TrackingLoops.uses_soft_secondary_code_detection(b1c_p) == false

    for num_blocks in (0, 1, 1799)
        @test @inferred(
            detect_bit_or_secondary_code_sync(
                b1c_p,
                prn,
                TrackingLoops.UInt1800(0x1),
                num_blocks,
            )
        ).found == false
    end

    @test get_band_id(b1c_p) === :L1

    # QMBOC's principal axis is BOC(1,1): VEML.
    @test @inferred(get_default_correlator(b1c_p, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b1c_p, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(b1c_p)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b1c_p)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(b1c_p)) === TrackingLoops.UInt1800

    @testset "Overlay search — clean lock at known phase / polarity" begin
        # A `PerPRNSecondaryCode`, so the generic packer applies.
        reference = TrackingLoops._packed_secondary_code(TrackingLoops.UInt1800, b1c_p, prn)
        rotl(x, r) = r == 0 ? x : ((x << r) | (x >> (1800 - r)))
        for r in (0, 137, 1799)
            received = rotl(reference, r)
            res = @inferred detect_bit_or_secondary_code_sync(b1c_p, prn, received, 1800)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end

        all_ones =
            (TrackingLoops.UInt1800(1) << 1799) |
            ((TrackingLoops.UInt1800(1) << 1799) - one(TrackingLoops.UInt1800))
        negated = reference ⊻ all_ones
        res = @inferred detect_bit_or_secondary_code_sync(b1c_p, prn, negated, 1800)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Overlay search — tolerance" begin
        # 2.5 % of 1800 discretizes to 45 errors, as for GPS L1C-P.
        @test B1C_P_MAX_ERRORS == 45
        overlay = TrackingLoops._packed_secondary_code(TrackingLoops.UInt1800, b1c_p, prn)
        rng = MersenneTwister(42)

        flip(x, n) = begin
            corrupted = x
            for idx in randperm(rng, 1800)[1:n]
                corrupted ⊻= TrackingLoops.UInt1800(1) << (idx - 1)
            end
            corrupted
        end

        for n_errors in (0, 1, 10, B1C_P_MAX_ERRORS)
            res =
                detect_bit_or_secondary_code_sync(b1c_p, prn, flip(overlay, n_errors), 1800)
            @test res.found == true
            @test res.phase == 0
        end

        # One error more rejects (random flips could match another rotation;
        # this seed does not).
        res = detect_bit_or_secondary_code_sync(
            b1c_p,
            prn,
            flip(overlay, B1C_P_MAX_ERRORS + 1),
            1800,
        )
        @test res.found == false
    end
end

end
