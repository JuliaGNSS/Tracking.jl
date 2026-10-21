module GPSL1CATest

# GPS L1 C/A's sync traits are tested in TrackingLoops (`test/signals.jl` there).

using Test: @test, @testset
using GNSSSignals: GPSL1CA
import Tracking
import TrackingLoops
using TrackingLoops: get_bit_edge_detection_confidence

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
