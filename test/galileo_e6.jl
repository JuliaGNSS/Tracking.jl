module GalileoE6Test

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals:
    GalileoE5aQ, GalileoE6B, GalileoE6C, get_carrier_phase_offset, get_secondary_code_length
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

@testset "Galileo E6-B" begin
    e6b = GalileoE6B()
    prn = 1

    # One C/NAV symbol per 1 ms primary code period (1000 sym/s) and no
    # secondary code, so the detector fires immediately and unconditionally —
    # same shape as Galileo E1B / GPS L1C-D.
    @test get_secondary_code_length(e6b) == 1
    for num_blocks in (1, 2, 17)
        res = @inferred detect_bit_or_secondary_code_sync(e6b, prn, UInt8(0x0), num_blocks)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1
    end

    # Plain BPSK(5) (`LOC`) — no subcarrier, so the C/A-style EarlyPromptLate
    # default applies, not the BOC VeryEarlyPromptLate one E1B takes.
    @test @inferred(get_default_correlator(e6b, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(e6b, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 1 ms primary period (5115 chips at 5.115 Mcps) → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(e6b)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e6b)) ≈ 1.0Hz

    # No sub-symbol boundary to search — the buffer is dead state, UInt8.
    @test @inferred(get_code_block_buffer_type(e6b)) === UInt8

    # No secondary code and one block per symbol, so neither soft detector
    # applies: the trivial hard detector above is the whole story.
    @test Tracking.uses_soft_secondary_code_detection(e6b) == false
    @test Tracking.uses_soft_bit_edge_detection(e6b) == false
end

@testset "Galileo E6-C" begin
    e6c = GalileoE6C()
    prn = 1
    N = get_secondary_code_length(e6c)  # 100 (CS100_n)
    @test N == 100

    @test @inferred(detect_bit_or_secondary_code_sync(e6c, prn, UInt128(0x0), N - 1)).found ==
          false

    @testset "CS100 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt128, e6c, prn)
        for r in (0, 61, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(e6c, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt128) << N) - one(UInt128))
        res = @inferred detect_bit_or_secondary_code_sync(e6c, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    # E6-C draws the *same* CS100_1..50 half of the table as Galileo E5a-Q, so
    # for a given SVID the two overlays coincide — pinned here because the
    # detector's per-PRN reference is derived generically from the signal.
    @test Tracking._packed_secondary_code(UInt128, e6c, prn) ==
          Tracking._packed_secondary_code(UInt128, GalileoE5aQ(), prn)

    @test @inferred(get_default_correlator(e6c, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))

    # 1 ms primary period, same as E6-B → 18 Hz / 1 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(e6c)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(e6c)) ≈ 1.0Hz

    # 100-block CS100 window needs UInt128.
    @test @inferred(get_code_block_buffer_type(e6c)) === UInt128

    @test Tracking.uses_soft_secondary_code_detection(e6c) == true

    # E6-C is the anti-phase arm of the E6 composite ((e_B − e_C)/√2, OS SIS
    # ICD Eq. 10). Tracking reads that offset generically when both components
    # share one driver loop, so pin the value the loop relies on.
    @test get_carrier_phase_offset(e6c) ≈ π
    @test get_carrier_phase_offset(GalileoE6B()) == 0.0
end

end
