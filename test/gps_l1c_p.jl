module GPSL1CPTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1C_P
using Random: MersenneTwister, randperm
import Tracking
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts,
    SyncResult,
    L1C_P_MAX_ERRORS

@testset "GPS L1C-P" begin
    gpsl1c_p = GPSL1C_P()

    # Below the 1800-block horizon, the L1C-P detector returns `found =
    # false` without running the sweep. Above it, the sweep runs against
    # the per-PRN overlay.
    prn = 1
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, prn, Tracking.UInt1800(0x0), 0)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, prn, Tracking.UInt1800(0x1), 1)).found == false
    @test @inferred(is_upcoming_integration_new_bit(gpsl1c_p, prn, Tracking.UInt1800(0xffffffff), 1799)).found == false

    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 10 ms primary period at BL·T ≈ 0.018 → 1.8 Hz carrier / 0.1 Hz code.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1c_p)) ≈ 1.8Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_p)) ≈ 0.1Hz

    # 1800-chip per-PRN overlay → exact-width UInt1800.
    @test @inferred(get_code_block_buffer_type(gpsl1c_p)) === Tracking.UInt1800

    @testset "Overlay search — clean lock at known phase / polarity" begin
        # Build PRN 1's overlay-packed UInt1800 by reusing the internal
        # helper, then rotate-left by a known offset and feed it back.
        # The detector should recover the offset and lock at positive
        # polarity (no bit-flips, distance 0 < L1C_P_MAX_ERRORS).
        overlay = Tracking._pack_overlay(gpsl1c_p, prn)
        for k in (0, 137, 1799)
            rotated = k == 0 ? overlay :
                ((overlay >> k) | (overlay << (1800 - k)))
            res = @inferred is_upcoming_integration_new_bit(gpsl1c_p, prn, rotated, 1800)
            @test res.found == true
            @test res.phase == k
            @test res.polarity == +1
        end

        # Negative polarity = bitwise NOT of the overlay within the
        # 1800-bit window. The exact-width UInt1800 makes `~` equivalent
        # to XOR with all-ones; build that explicitly.
        all_ones = (Tracking.UInt1800(1) << 1799) | ((Tracking.UInt1800(1) << 1799) - one(Tracking.UInt1800))
        negated = overlay ⊻ all_ones
        res = @inferred is_upcoming_integration_new_bit(gpsl1c_p, prn, negated, 1800)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Overlay search — tolerance" begin
        overlay = Tracking._pack_overlay(gpsl1c_p, prn)
        rng = MersenneTwister(42)

        # Up to `L1C_P_MAX_ERRORS` bit-flips: still locks.
        for n_errors in (0, 1, 10, L1C_P_MAX_ERRORS)
            corrupted = overlay
            indices = randperm(rng, 1800)[1:n_errors]
            for idx in indices
                corrupted ⊻= Tracking.UInt1800(1) << (idx - 1)
            end
            res = is_upcoming_integration_new_bit(gpsl1c_p, prn, corrupted, 1800)
            @test res.found == true
            @test res.phase == 0
        end

        # One extra error — must reject.
        n_errors = L1C_P_MAX_ERRORS + 1
        corrupted = overlay
        indices = randperm(rng, 1800)[1:n_errors]
        for idx in indices
            corrupted ⊻= Tracking.UInt1800(1) << (idx - 1)
        end
        # Note: with random flips it's *possible* (very small probability)
        # for the corrupted buffer to coincide with the overlay rotated
        # by some other phase within tolerance. With this fixed seed we
        # verified that doesn't happen for PRN 1.
        res = is_upcoming_integration_new_bit(gpsl1c_p, prn, corrupted, n_errors == L1C_P_MAX_ERRORS + 1 ? 1800 : 1800)
        @test res.found == false
    end
end

end
