module GPSL1CATest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA
import Tracking
using Tracking:
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    uses_soft_bit_edge_detection,
    get_bit_edge_detection_confidence,
    _detect_bit_edge_cfar,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1" begin
    gpsl1 = GPSL1CA()

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

    @testset "Soft bit-edge detection traits" begin
        # GPS L1 C/A uses the soft, maximum-energy CFAR bit-edge detector.
        @test @inferred(uses_soft_bit_edge_detection(gpsl1)) == true
        # Default confidence target.
        @test @inferred(get_bit_edge_detection_confidence(gpsl1)) ≈ 0.999
    end
end

# User override of the bit-edge detection confidence via dispatch on
# `get_bit_edge_detection_confidence`. The detector picks up the override
# immediately — no TrackState rebuild needed.
#
# `Core.eval` is used so the override and its rollback execute at test time
# rather than at module-parse time (literal method-definition expressions get
# hoisted to module scope and the last one would win unconditionally).
@testset "GPS L1 — confidence override" begin
    gpsl1 = GPSL1CA()
    @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.999
    Core.eval(Tracking, :(get_bit_edge_detection_confidence(::$GPSL1CA) = 0.95))
    try
        @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.95
    finally
        # Restore the package-wide default for any tests that run after this one.
        Core.eval(Tracking, :(get_bit_edge_detection_confidence(::$GPSL1CA) = 0.999))
    end
    @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.999
end

end
