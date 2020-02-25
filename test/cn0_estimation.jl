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
    sample_frequency = 4MHz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    cn0_estimator = MomentsCN0Estimator(20)
    early_late_sample_shift = 2
    start_sample = 1
    num_samples = 4000
    agc_bits = 5

    for i = 1:20
        signal = get_code.(
                GPSL1,
                code_frequency .* range ./ sample_frequency .+ start_code_phase,
                prn
            ) .* 10^(45 / 20) .+
            randn(ComplexF64, length(range)) .* sqrt(4e6)
        attenuation = Tracking.find_max(signal)
        agc_signal = StructArray{Complex{Int16}}((
            floor.(Int16, real.(signal * 1 << agc_bits / attenuation)),
            floor.(Int16, imag.(signal * 1 << agc_bits / attenuation))
        ))
        code = get_code.(
            GPSL1,
            code_frequency .* (-2:4001) ./ sample_frequency .+ start_code_phase,
            prn
        )
        correlator = EarlyPromptLateCorrelator()
        correlator = Tracking.correlate(
            correlator,
            agc_signal,
            code,
            early_late_sample_shift,
            start_sample,
            num_samples,
            1.0,
            agc_bits,
            Val(1)
        )
        cn0_estimator = Tracking.update(cn0_estimator, get_prompt(correlator))
    end
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 20
    @test @inferred(Tracking.length(cn0_estimator)) == 20
    cn0_estimate = @inferred Tracking.estimate_cn0(cn0_estimator, 1ms)

    @test cn0_estimate ≈ 45dBHz atol = 0.30dBHz

end
