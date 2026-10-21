module GPSL1CATest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA
import Tracking
import TrackingLoops
using TrackingLoops:
    get_default_correlator,
    get_code_block_buffer_type,
    default_narrow_carrier_loop_filter_bandwidth,
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

    @test @inferred(default_narrow_carrier_loop_filter_bandwidth(gpsl1)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1)) ≈ 1.0Hz

    @test @inferred(get_code_block_buffer_type(gpsl1)) === UInt64

    @testset "Soft bit-edge detection traits" begin
        @test @inferred(uses_soft_bit_edge_detection(gpsl1)) == true
        @test @inferred(get_bit_edge_detection_confidence(gpsl1)) ≈ 0.999
    end
end

# A user override of `get_bit_edge_detection_confidence` takes effect without a
# TrackState rebuild. `Core.eval` defines it at test time, not for the whole
# file, and in this module, where a user's would live: in `Tracking` it would be
# type piracy, which Aqua refuses.
@testset "GPS L1 — confidence override" begin
    gpsl1 = GPSL1CA()
    @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.999
    Core.eval(
        @__MODULE__,
        :(TrackingLoops.get_bit_edge_detection_confidence(::$GPSL1CA) = 0.95),
    )
    try
        @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.95
    finally
        # Remove the override for the tests after this one.
        Base.delete_method(which(get_bit_edge_detection_confidence, (GPSL1CA,)))
    end
    @test get_bit_edge_detection_confidence(gpsl1) ≈ 0.999
end

end
