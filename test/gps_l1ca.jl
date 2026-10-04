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

    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1)) ≈ 1.0Hz

    # 40-block sync window fits in a UInt64.
    @test @inferred(get_code_block_buffer_type(gpsl1)) === UInt64

    @testset "Soft bit-edge detection traits" begin
        @test @inferred(uses_soft_bit_edge_detection(gpsl1)) == true
        @test @inferred(get_bit_edge_detection_confidence(gpsl1)) ≈ 0.999
    end
end

# Users override the confidence by dispatch; the detector picks it up without a
# TrackState rebuild. `Core.eval` runs the override and its rollback at test time
# (literal method definitions are hoisted, so the last one would always win).
@testset "GPS L1 — confidence override" begin
    gpsl1 = GPSL1CA()
    @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.999
    Core.eval(Tracking, :(get_bit_edge_detection_confidence(::$GPSL1CA) = 0.95))
    try
        @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.95
    finally
        # Restore the default for later tests.
        Core.eval(Tracking, :(get_bit_edge_detection_confidence(::$GPSL1CA) = 0.999))
    end
    @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.999
end

end
