@testset "CUDA: CN0 estimation" begin
    CUDA.allowscalar(true)
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
    gpsl1 = GPSL1(use_gpu = Val(true))

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
        signal_struct = StructArray(signal)
        correlator = Tracking.correlate(
            correlator,
            signal_struct,
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
