@testset "Moments CN0 estimator" begin
    cn0_estimator = Tracking.MomentsCN0Estimator(20)
    @test @inferred(Tracking.get_prompt_buffer(cn0_estimator)) ==
          zero(SVector{20,ComplexF64})
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 0
    @test @inferred(Tracking.length(cn0_estimator)) == 0

    next_cn0_estimator = @inferred Tracking.update(cn0_estimator, 1 + 2im)
    @test @inferred(Tracking.get_prompt_buffer(next_cn0_estimator))[1] == 1 + 2im
    @test @inferred(Tracking.get_current_index(next_cn0_estimator)) == 1
    @test @inferred(Tracking.length(next_cn0_estimator)) == 1

    cn0_estimator = Tracking.MomentsCN0Estimator(ones(SVector{20,ComplexF64}), 20, 20)
    next_cn0_estimator = @inferred Tracking.update(cn0_estimator, 1 + 2im)
    @test @inferred(Tracking.get_prompt_buffer(next_cn0_estimator))[1] == 1 + 2im
    @test @inferred(Tracking.get_current_index(next_cn0_estimator)) == 1
    @test @inferred(Tracking.length(next_cn0_estimator)) == 20

    cn0_estimator = Tracking.MomentsCN0Estimator(ones(SVector{20,ComplexF64}), 19, 20)
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 19
    @test @inferred(Tracking.length(cn0_estimator)) == 20
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
    cn0_estimator = Tracking.MomentsCN0Estimator(100)
    start_sample = 1
    num_samples = 4000
    gpsl1 = GPSL1()

    for i = 1:100
        signal =
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ start_code_phase,
                prn,
            ) .* 10^(45 / 20) .+ randn(ComplexF64, length(range)) .* sqrt(4e6)
        code =
            get_code.(
                gpsl1,
                code_frequency .* (-2:4001) ./ sampling_frequency .+ start_code_phase,
                prn,
            )
        correlator = EarlyPromptLateCorrelator(gpsl1, sampling_frequency)
        signal_struct = StructArray(signal)
        correlator =
            Tracking.correlate(correlator, signal_struct, code, start_sample, num_samples)
        cn0_estimator = Tracking.update(cn0_estimator, get_prompt(correlator))
    end
    @test @inferred(Tracking.get_current_index(cn0_estimator)) == 100
    @test @inferred(Tracking.length(cn0_estimator)) == 100
    cn0_estimate = @inferred estimate_cn0(cn0_estimator, 1ms)

    @test cn0_estimate ≈ 45dBHz atol = 1.0dBHz
end

@testset "CN0 estimation integration test" begin
    Random.seed!(1234)
    carrier_doppler = 0Hz
    start_code_phase = 0
    code_frequency = 1023kHz
    sampling_frequency = 4MHz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    cn0_estimator = Tracking.MomentsCN0Estimator(100)
    start_sample = 1
    num_samples = 4000
    gpsl1 = GPSL1()

    track_state = @inferred TrackState(
        gpsl1,
        [SatState(gpsl1, prn, sampling_frequency, start_code_phase, carrier_doppler)];
        num_samples,
    )

    for i = 1:100
        signal =
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ start_code_phase,
                prn,
            ) .* 10^(45 / 20) .+ randn(ComplexF64, length(range)) .* sqrt(4e6)
        start_code_phase =
            code_frequency * length(range) ./ sampling_frequency + start_code_phase
        track_state = @inferred track(signal, track_state, sampling_frequency)
    end
    cn0_estimate = @inferred estimate_cn0(track_state, 1)
    @test cn0_estimate ≈ 45dBHz atol = 1.0dBHz
end