@testset "Update CN0 Estimator" begin

    cn0_estimator = MomentsCN0Estimator(20)
    @test @inferred(Tracking.get_prompt_buffer(cn0_estimator)) == zero(SVector{20, ComplexF64})
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 0
    @test @inferred(Tracking.length(cn0_estimator)) == 0

    next_cn0_estimator = @inferred Tracking.update(cn0_estimator, 1 + 2im)
    @test @inferred(Tracking.get_prompt_buffer(next_cn0_estimator))[1] == 1 + 2im
    @test @inferred(Tracking.get_current_index(next_cn0_estimator)) == 1
    @test @inferred(Tracking.length(next_cn0_estimator)) == 1

    cn0_estimator = MomentsCN0Estimator(ones(SVector{20, ComplexF64}), 20, 20)
    next_cn0_estimator = @inferred Tracking.update(cn0_estimator, 1 + 2im)
    @test @inferred(Tracking.get_prompt_buffer(next_cn0_estimator))[1] == 1 + 2im
    @test @inferred(Tracking.get_current_index(next_cn0_estimator)) == 1
    @test @inferred(Tracking.length(next_cn0_estimator)) == 20

    cn0_estimator = MomentsCN0Estimator(ones(SVector{20, ComplexF64}), 19, 20)
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 19
    @test @inferred(Tracking.length(cn0_estimator)) == 20

    @test @allocated(Tracking.update(MomentsCN0Estimator(20), 1 + 2im)) == 0
end

@testset "CN0 estimation" begin

    Random.seed!(1234)
    carrier_doppler = 0Hz
    start_code_phase = 0
    code_frequency = 1023kHz
    sampling_frequency = 4MHz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    cn0_estimator = MomentsCN0Estimator(20)
    correlator_sample_shifts = SVector(-2, 0, 2)
    start_sample = 1
    num_samples = 4000
    gpsl1 = GPSL1()

    for i = 1:20
        signal = get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ start_code_phase,
                prn
            ) .* 10^(45 / 20) .+
            randn(ComplexF64, length(range)) .* sqrt(4e6)
        code = get_code.(
            gpsl1,
            code_frequency .* (-2:4001) ./ sampling_frequency .+ start_code_phase,
            prn
        )
        correlator = EarlyPromptLateCorrelator()
        correlator = Tracking.correlate(
            correlator,
            signal,
            code,
            correlator_sample_shifts,
            start_sample,
            num_samples,
        )
        cn0_estimator = Tracking.update(
            cn0_estimator,
            get_prompt(correlator, correlator_sample_shifts)
        )
    end
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 20
    @test @inferred(Tracking.length(cn0_estimator)) == 20
    cn0_estimate = @inferred Tracking.estimate_cn0(cn0_estimator, 1ms)

    @test cn0_estimate ≈ 45dBHz atol = 1.05dBHz

end
