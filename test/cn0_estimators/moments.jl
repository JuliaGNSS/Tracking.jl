module MomentsCN0EstimatorTest

using Test: @test, @testset, @inferred
using Random: Random
using Unitful: kHz, MHz, Hz, ms, dBHz
using GNSSSignals: GPSL1CA, get_code
import Tracking
using TrackingLoops:
    MomentsCN0Estimator,
    update,
    EarlyPromptLateCorrelator,
    get_prompt,
    get_correlator_sample_shifts,
    estimate_cn0

@testset "CN0 estimation" begin
    Random.seed!(1234)
    carrier_doppler = 0Hz
    start_code_phase = 0
    code_frequency = 1023kHz
    sampling_frequency = 4MHz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    cn0_estimator = MomentsCN0Estimator(100)
    num_samples = 4000
    gpsl1 = GPSL1CA()

    for i = 1:100
        samples =
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ start_code_phase,
                prn,
            ) .* 10^(45 / 20) .+ randn(ComplexF64, length(range)) .* sqrt(4e6)
        code = get_code.(
            gpsl1,
            code_frequency .* (-2:4001) ./ sampling_frequency .+ start_code_phase,
            prn,
        )
        correlator = EarlyPromptLateCorrelator()
        sample_shifts =
            get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
        correlator = Tracking.downconvert_and_correlate_fused!(
            correlator,
            samples,
            code,
            sample_shifts,
            0.0Hz,
            sampling_frequency,
            0.0,
            1:num_samples,
        )
        cn0_estimator = update(cn0_estimator, get_prompt(correlator))
    end
    @test @inferred(Base.length(cn0_estimator)) == 100
    cn0_estimate = @inferred estimate_cn0(cn0_estimator, 1ms)

    @test cn0_estimate ≈ 45dBHz atol = 1.0dBHz
end

end
