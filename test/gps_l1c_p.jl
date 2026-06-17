module GPSL1CPTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1C_P
using Random: MersenneTwister, randperm
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    get_bit_edge_or_secondary_code_tolerance,
    EarlyPromptLateCorrelator,
    NumAnts,
    SyncResult

const L1C_P_MAX_ERRORS =
    floor(Int, get_bit_edge_or_secondary_code_tolerance(GPSL1C_P()) * 1800)

@testset "GPS L1C-P" begin
    gpsl1c_p = GPSL1C_P()

    # Below the 1800-block horizon, the L1C-P detector returns `found =
    # false` without running the sweep. Above it, the sweep runs against
    # the per-PRN overlay.
    prn = 1
    @test @inferred(
        detect_bit_or_secondary_code_sync(gpsl1c_p, prn, Tracking.UInt1800(0x0), 0)
    ).found == false
    @test @inferred(
        detect_bit_or_secondary_code_sync(gpsl1c_p, prn, Tracking.UInt1800(0x1), 1)
    ).found == false
    @test @inferred(
        detect_bit_or_secondary_code_sync(
            gpsl1c_p,
            prn,
            Tracking.UInt1800(0xffffffff),
            1799,
        )
    ).found == false

    # TMBOC(6,1,4/33): narrow 0.1-chip early-late spacing keeps the taps on
    # the BOC main peak rather than the side-lobes (see get_default_correlator).
    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(1))) ==
          EarlyPromptLateCorrelator(;
        num_ants = NumAnts(1),
        preferred_early_late_to_prompt_code_shift = 0.1,
    )
    @test @inferred(get_default_correlator(gpsl1c_p, NumAnts(3))) ==
          EarlyPromptLateCorrelator(;
        num_ants = NumAnts(3),
        preferred_early_late_to_prompt_code_shift = 0.1,
    )

    # 10 ms primary period at BL·T ≈ 0.018 → 1.8 Hz carrier / 0.1 Hz code.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1c_p)) ≈ 1.8Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1c_p)) ≈ 0.1Hz

    # 1800-chip per-PRN overlay → exact-width UInt1800.
    @test @inferred(get_code_block_buffer_type(gpsl1c_p)) === Tracking.UInt1800

    @testset "Overlay search — clean lock at known phase / polarity" begin
        # Build PRN 1's newest-first overlay reference, then rotate it *left*
        # by `r` to emulate a prompt buffer whose upcoming integration is
        # overlay chip `r`. The rotation search recovers `phase == r` (the
        # upcoming chip) at positive polarity (distance 0 ≤ max_errors).
        reference = Tracking._packed_secondary_code(Tracking.UInt1800, gpsl1c_p, prn)
        rotl(x, r) = r == 0 ? x : ((x << r) | (x >> (1800 - r)))
        for r in (0, 137, 1799)
            received = rotl(reference, r)
            res = @inferred detect_bit_or_secondary_code_sync(gpsl1c_p, prn, received, 1800)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end

        # Negative polarity = bitwise NOT of the reference within the
        # 1800-bit window. The exact-width UInt1800 makes `~` equivalent
        # to XOR with all-ones; build that explicitly. No rotation, so the
        # recovered upcoming chip is 0.
        all_ones =
            (Tracking.UInt1800(1) << 1799) |
            ((Tracking.UInt1800(1) << 1799) - one(Tracking.UInt1800))
        negated = reference ⊻ all_ones
        res = @inferred detect_bit_or_secondary_code_sync(gpsl1c_p, prn, negated, 1800)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "Overlay search — tolerance" begin
        overlay = Tracking._packed_secondary_code(Tracking.UInt1800, gpsl1c_p, prn)
        rng = MersenneTwister(42)

        # Up to `L1C_P_MAX_ERRORS` bit-flips: still locks.
        for n_errors in (0, 1, 10, L1C_P_MAX_ERRORS)
            corrupted = overlay
            indices = randperm(rng, 1800)[1:n_errors]
            for idx in indices
                corrupted ⊻= Tracking.UInt1800(1) << (idx - 1)
            end
            res = detect_bit_or_secondary_code_sync(gpsl1c_p, prn, corrupted, 1800)
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
