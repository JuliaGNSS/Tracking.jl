module NoCN0EstimatorTest

# `NoCN0Estimator` itself is tested in TrackingLoops.

using Test: @test, @testset
using Unitful: Hz, dBHz
using GNSSSignals: GPSL1CA, GPSL1C_P
import Tracking
using Tracking: TrackedSat, get_cn0_estimator
using TrackingLoops: NoCN0Estimator, NWPRCN0Estimator, estimate_cn0

@testset "NoCN0Estimator plugs into one signal of a multi-signal satellite" begin
    sat = TrackedSat(
        (GPSL1C_P(), GPSL1CA()),
        1,
        0.0,
        0.0Hz;
        cn0_estimator = (NWPRCN0Estimator(), NoCN0Estimator()),
    )
    @test get_cn0_estimator(sat, GPSL1CA) isa NoCN0Estimator
    @test estimate_cn0(sat, GPSL1CA) == -Inf * dBHz
end

end
