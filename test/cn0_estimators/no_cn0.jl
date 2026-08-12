module NoCN0EstimatorTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms, dBHz
using GNSSSignals: GPSL1CA, GPSL1C_P
import Tracking
using Tracking:
    NoCN0Estimator,
    NWPRCN0Estimator,
    CN0UpdateContext,
    BitBuffer,
    update,
    estimate_cn0,
    TrackedSat,
    get_cn0_estimator

@testset "NoCN0Estimator measures nothing and says so" begin
    estimator = NoCN0Estimator()
    # Stateless: the same instance comes back, whichever `update` form is used.
    @test @inferred(update(estimator, 1.0 + 0.0im)) === estimator
    context = CN0UpdateContext(GPSL1CA(), BitBuffer{UInt64}(), 1)
    @test @inferred(update(estimator, 1.0 + 0.0im, context)) === estimator
    fold(e, n) = (for _ = 1:n
        e = update(e, 1.0 + 0.0im)
    end; e)
    fold(estimator, 10)
    @test @allocated(fold(estimator, 1000)) == 0

    # `-Inf dB-Hz`, not `NaN dB-Hz`: `NaN dB-Hz >= threshold` is `true` for every
    # threshold with Unitful's `Level` comparison, so a NaN would clear every lock
    # detector it met. `-Inf` is the safe answer to "is this signal locked?".
    @test @inferred(estimate_cn0(estimator, 1ms)) == -Inf * dBHz
    @test !(estimate_cn0(estimator, 1ms) >= 20dBHz)
    @test NaN * dBHz >= 20dBHz          # ... which is why NaN is not used

    # It is a legal `fallback`, which is the point: it replaces the moment ratio's
    # noise floor with "no estimate" for the phases that admit no coherent window.
    honest = NWPRCN0Estimator(; num_narrowband_code_blocks = 20, fallback = estimator)
    @test @inferred(estimate_cn0(honest, 1ms)) == -Inf * dBHz
    for _ = 1:19
        honest = update(honest, 1.0 + 0.0im)
    end
    @test estimate_cn0(honest, 1ms) == -Inf * dBHz     # window still open
    honest = update(honest, 1.0 + 0.0im)
    @test estimate_cn0(honest, 1ms) == Inf * dBHz      # first window closed

    # And it plugs into a satellite like any other estimator, including for a
    # single signal of a multi-signal sat.
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
