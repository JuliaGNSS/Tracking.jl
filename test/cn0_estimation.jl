@testset "Prompts buffer" begin

    prompts_buffer = Tracking.PromptsBuffer(20)
    @test @inferred(Tracking.get_prompt_buffer(prompts_buffer)) == zero(SVector{20, ComplexF64})
    @test @inferred(Tracking.get_current_index(prompts_buffer)) == 0
    @test @inferred(Tracking.length(prompts_buffer)) == 0

    next_prompts_buffer = @inferred Tracking.update(prompts_buffer, 1 + 2im)
    @test @inferred(Tracking.get_prompt_buffer(next_prompts_buffer))[1] == 1 + 2im
    @test @inferred(Tracking.get_current_index(next_prompts_buffer)) == 1
    @test @inferred(Tracking.length(next_prompts_buffer)) == 1

    prompts_buffer = Tracking.PromptsBuffer(ones(SVector{20, ComplexF64}), 20, 20)
    next_prompts_buffer = @inferred Tracking.update(prompts_buffer, 1 + 2im)
    @test @inferred(Tracking.get_prompt_buffer(next_prompts_buffer))[1] == 1 + 2im
    @test @inferred(Tracking.get_current_index(next_prompts_buffer)) == 1
    @test @inferred(Tracking.length(next_prompts_buffer)) == 20

    prompts_buffer = Tracking.PromptsBuffer(ones(SVector{20, ComplexF64}), 19, 20)
    @test @inferred(Tracking.get_current_index(prompts_buffer)) == 19
    @test @inferred(Tracking.length(prompts_buffer)) == 20
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
    cn0_estimator = MomentsCN0Estimator()
    prompts_buffer = Tracking.PromptsBuffer(20)
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
        correlator = EarlyPromptLateCorrelator(gpsl1, sampling_frequency)
        signal_struct = StructArray(signal)
        correlator = Tracking.correlate(
            correlator,
            signal_struct,
            code,
            start_sample,
            num_samples,
        )
        prompts_buffer = Tracking.update(
            prompts_buffer,
            get_prompt(correlator)
        )
    end
    @test @inferred(Tracking.get_current_index(prompts_buffer)) == 20
    @test @inferred(Tracking.length(prompts_buffer)) == 20
    cn0_estimate = @inferred Tracking.estimate_cn0(prompts_buffer, cn0_estimator, 1ms)

    @test cn0_estimate ≈ 45dBHz atol = 2.05dBHz

end
