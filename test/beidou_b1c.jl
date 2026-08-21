module BeiDouB1CTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using Random: MersenneTwister, randperm
using GNSSSignals: BeiDouB1C_D, BeiDouB1C_P, get_band_id, get_secondary_code_length
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    get_bit_edge_or_secondary_code_tolerance,
    VeryEarlyPromptLateCorrelator,
    NumAnts

const B1C_P_MAX_ERRORS =
    floor(Int, get_bit_edge_or_secondary_code_tolerance(BeiDouB1C_P()) * 1800)

@testset "BeiDou B1C data" begin
    b1c_d = BeiDouB1C_D()
    prn = 1

    # One B-CNAV1 symbol per 10 ms primary code period (100 sym/s) and no
    # secondary code, so the detector fires immediately — the GPS L1C-D case.
    @test get_secondary_code_length(b1c_d) == 1
    for num_blocks in (1, 2, 7)
        res =
            @inferred detect_bit_or_secondary_code_sync(b1c_d, prn, UInt8(0x0), num_blocks)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    @test get_band_id(b1c_d) === :L1

    # BOC(1,1): the VeryEarlyPromptLate default, whose very-early/very-late
    # taps feed the VEML discriminator that mitigates the BOC side-peak false
    # locks — same as GPS L1C and Galileo E1.
    @test @inferred(get_default_correlator(b1c_d, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b1c_d, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 10 ms primary period (10230 chips at 1.023 Mcps) → 1.8 Hz carrier, and
    # short enough to leave the DLL at its unclamped 1 Hz reference.
    @test @inferred(default_carrier_loop_filter_bandwidth(b1c_d)) ≈ 1.8Hz
    @test @inferred(default_code_loop_filter_bandwidth(b1c_d)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(b1c_d)) === UInt8

    @test Tracking.uses_soft_secondary_code_detection(b1c_d) == false
    @test Tracking.uses_soft_bit_edge_detection(b1c_d) == false
end

@testset "BeiDou B1C pilot" begin
    b1c_p = BeiDouB1C_P()
    prn = 1
    N = get_secondary_code_length(b1c_p)  # 1800
    @test N == 1800

    # Like GPS L1C-P's, an 1800-chip / 18 s overlay is far too long to
    # integrate coherently per bin, so B1C pilot stays on the hard-decision
    # rotation sweep rather than the soft CFAR secondary-code detector.
    @test Tracking.uses_soft_secondary_code_detection(b1c_p) == false

    # Below the 1800-block horizon the detector returns `found = false`
    # without running the sweep.
    for num_blocks in (0, 1, 1799)
        @test @inferred(
            detect_bit_or_secondary_code_sync(
                b1c_p,
                prn,
                Tracking.UInt1800(0x1),
                num_blocks,
            )
        ).found == false
    end

    @test get_band_id(b1c_p) === :L1

    # BOC(1,1)-class (the QMBOC pilot's principal axis *is* BOC(1,1)) →
    # VeryEarlyPromptLate default, same as the data component.
    @test @inferred(get_default_correlator(b1c_p, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b1c_p, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 10 ms primary period → 1.8 Hz / 1 Hz, same as the data component.
    @test @inferred(default_carrier_loop_filter_bandwidth(b1c_p)) ≈ 1.8Hz
    @test @inferred(default_code_loop_filter_bandwidth(b1c_p)) ≈ 1.0Hz

    # 1800-chip per-PRN overlay → the exact-width UInt1800 that GPS L1C-P's
    # identically sized overlay already introduced.
    @test @inferred(get_code_block_buffer_type(b1c_p)) === Tracking.UInt1800

    @testset "Overlay search — clean lock at known phase / polarity" begin
        # Unlike GPS L1C-P, B1C pilot needs no bespoke packer: GNSSSignals
        # exposes the overlay as a `PerPRNSecondaryCode`, which the generic
        # `_packed_secondary_code` reads.
        reference = Tracking._packed_secondary_code(Tracking.UInt1800, b1c_p, prn)
        rotl(x, r) = r == 0 ? x : ((x << r) | (x >> (1800 - r)))
        for r in (0, 137, 1799)
            received = rotl(reference, r)
            res = @inferred detect_bit_or_secondary_code_sync(b1c_p, prn, received, 1800)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end

        all_ones =
            (Tracking.UInt1800(1) << 1799) |
            ((Tracking.UInt1800(1) << 1799) - one(Tracking.UInt1800))
        negated = reference ⊻ all_ones
        res = @inferred detect_bit_or_secondary_code_sync(b1c_p, prn, negated, 1800)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Overlay search — tolerance" begin
        # 2.5 % of 1800 discretizes to 45 errors, as for GPS L1C-P.
        @test B1C_P_MAX_ERRORS == 45
        overlay = Tracking._packed_secondary_code(Tracking.UInt1800, b1c_p, prn)
        rng = MersenneTwister(42)

        flip(x, n) = begin
            corrupted = x
            for idx in randperm(rng, 1800)[1:n]
                corrupted ⊻= Tracking.UInt1800(1) << (idx - 1)
            end
            corrupted
        end

        # Up to the error budget: still locks, at the unrotated phase.
        for n_errors in (0, 1, 10, B1C_P_MAX_ERRORS)
            res =
                detect_bit_or_secondary_code_sync(b1c_p, prn, flip(overlay, n_errors), 1800)
            @test res.found == true
            @test res.phase == 0
        end

        # One error past the budget must reject. With random flips it is in
        # principle possible (vanishingly unlikely over an 1800-chip code) for
        # the corrupted buffer to land within tolerance of some *other*
        # rotation; the fixed seed pins that it does not here.
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
