module GPSL1CPTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1C_P
using Random: MersenneTwister, randperm
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
    NumAnts,
    SyncResult

const L1C_P_MAX_ERRORS =
    floor(Int, get_bit_edge_or_secondary_code_tolerance(GPSL1C_P()) * 1800)

@testset "GPS L1C-P" begin
    gpsl1c_p = GPSL1C_P()

    # An 18 s overlay is too long for coherent bins: the hard rotation sweep.
    @test TrackingLoops.uses_soft_secondary_code_detection(gpsl1c_p) == false

    prn = 1
    @test @inferred(
        detect_bit_or_secondary_code_sync(gpsl1c_p, prn, TrackingLoops.UInt1800(0x0), 0)
    ).found == false
    @test @inferred(
        detect_bit_or_secondary_code_sync(gpsl1c_p, prn, TrackingLoops.UInt1800(0x1), 1)
    ).found == false
    @test @inferred(
        detect_bit_or_secondary_code_sync(
            gpsl1c_p,
            prn,
            TrackingLoops.UInt1800(0xffffffff),
            1799,
        )
    ).found == false

    # TMBOC: VEML against BOC side-peak locks.
    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 10 ms leaves the 1 Hz DLL unclamped.
    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(gpsl1c_p)) == 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_p)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(gpsl1c_p)) === TrackingLoops.UInt1800

    @testset "Overlay search — clean lock at known phase / polarity" begin
        # A buffer whose upcoming integration is overlay chip `r`.
        reference =
            TrackingLoops._packed_secondary_code(TrackingLoops.UInt1800, gpsl1c_p, prn)
        rotl(x, r) = r == 0 ? x : ((x << r) | (x >> (1800 - r)))
        for r in (0, 137, 1799)
            received = rotl(reference, r)
            res = @inferred detect_bit_or_secondary_code_sync(gpsl1c_p, prn, received, 1800)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end

        # Negative polarity: the complement within the 1800-bit window.
        all_ones =
            (TrackingLoops.UInt1800(1) << 1799) |
            ((TrackingLoops.UInt1800(1) << 1799) - one(TrackingLoops.UInt1800))
        negated = reference ⊻ all_ones
        res = @inferred detect_bit_or_secondary_code_sync(gpsl1c_p, prn, negated, 1800)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Overlay search — tolerance" begin
        overlay =
            TrackingLoops._packed_secondary_code(TrackingLoops.UInt1800, gpsl1c_p, prn)
        rng = MersenneTwister(42)

        for n_errors in (0, 1, 10, L1C_P_MAX_ERRORS)
            corrupted = overlay
            indices = randperm(rng, 1800)[1:n_errors]
            for idx in indices
                corrupted ⊻= TrackingLoops.UInt1800(1) << (idx - 1)
            end
            res = detect_bit_or_secondary_code_sync(gpsl1c_p, prn, corrupted, 1800)
            @test res.found == true
            @test res.phase == 0
        end

        # One error more rejects.
        n_errors = L1C_P_MAX_ERRORS + 1
        corrupted = overlay
        indices = randperm(rng, 1800)[1:n_errors]
        for idx in indices
            corrupted ⊻= TrackingLoops.UInt1800(1) << (idx - 1)
        end
        # Random flips could match another rotation; this seed does not.
        res = detect_bit_or_secondary_code_sync(
            gpsl1c_p,
            prn,
            corrupted,
            n_errors == L1C_P_MAX_ERRORS + 1 ? 1800 : 1800,
        )
        @test res.found == false
    end
end

end
