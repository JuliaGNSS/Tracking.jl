module GPSL1CATest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    get_bit_edge_or_secondary_code_tolerance,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1" begin
    gpsl1 = GPSL1CA()
    # PRN argument is ignored by the L1 C/A detector but required by the
    # uniform `detect_bit_or_secondary_code_sync(signal, prn, ...)` API.
    prn = 1

    # Matched at negative polarity (20 ones followed by 20 zeros).
    res = @inferred(detect_bit_or_secondary_code_sync(gpsl1, prn, UInt64(0xfffff00000), 50))
    @test res.found == true
    @test res.polarity == -1

    # Not enough integrations yet — buffer hasn't filled to 40 blocks.
    @test @inferred(detect_bit_or_secondary_code_sync(gpsl1, prn, UInt64(0xfffff), 10)).found == false

    # Matched at positive polarity (20 zeros followed by 20 ones).
    res = @inferred(detect_bit_or_secondary_code_sync(gpsl1, prn, UInt64(0xfffff), 40))
    @test res.found == true
    @test res.polarity == +1

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(gpsl1, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Per-signal default loop bandwidths: sized for the 1 ms primary code
    # period at BL·T ≈ 0.018. Historically the package shipped 18 Hz / 1 Hz
    # as the universal default — that value falls out of the formula here.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1)) ≈ 1.0Hz

    # 40-block sync window fits in a UInt64.
    @test @inferred(get_code_block_buffer_type(gpsl1)) === UInt64

    @testset "Hamming tolerance" begin
        # Default tolerance is 2.5 % → floor(0.025 × 40) = 1 error allowed.
        @test @inferred(get_bit_edge_or_secondary_code_tolerance(gpsl1)) ≈ 0.025
        template = UInt64(0xfffff)
        @test detect_bit_or_secondary_code_sync(gpsl1, prn, template ⊻ UInt64(0x1), 40).found == true
        # 2 errors → reject (above the 2.5 % ceiling).
        @test detect_bit_or_secondary_code_sync(gpsl1, prn, template ⊻ UInt64(0x3), 40).found == false
    end
end

# User override of the tolerance via dispatch on
# `get_bit_edge_or_secondary_code_tolerance`. The detector picks up the
# override immediately — no TrackState rebuild needed.
#
# `Core.eval` is used so the override and its rollback execute at test
# time rather than at module-parse time (literal method-definition
# expressions get hoisted to module scope and the last one would win
# unconditionally).
@testset "GPS L1 — tolerance override" begin
    gpsl1 = GPSL1CA()
    template = UInt64(0xfffff)
    Core.eval(Tracking, :(get_bit_edge_or_secondary_code_tolerance(::$GPSL1CA) = 0.10))
    try
        # 10 % over a 40-block window = 4 errors allowed.
        @test detect_bit_or_secondary_code_sync(gpsl1, 1, template ⊻ UInt64(0xf), 40).found == true
        # 5 errors → still reject — override raised the ceiling, not removed it.
        @test detect_bit_or_secondary_code_sync(gpsl1, 1, template ⊻ UInt64(0x1f), 40).found == false
    finally
        # Restore the package-wide default for any tests that run after this one.
        Core.eval(Tracking, :(get_bit_edge_or_secondary_code_tolerance(::$GPSL1CA) = 0.025))
    end
end

end
