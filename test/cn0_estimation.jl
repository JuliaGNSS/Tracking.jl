@testset "CN0 estimation" begin

    prompt_correlator_output = complex(2,2)
    cn0_state = Tracking.CN0State(0.0, 0.0, NaN, 20ms, 10ms)
    new_cn0_state = Tracking.estimate_CN0(cn0_state, 1ms, prompt_correlator_output)
    @test new_cn0_state.summed_total_power == 8.0
    @test new_cn0_state.summed_abs_inphase_ampl == 2.0
    @test isnan(new_cn0_state.last_valid_cn0)
    @test new_cn0_state.update_time == 20ms
    @test new_cn0_state.current_updated_time == 11ms

    signal = complex.(ones(19), ones(19))
    summed_total_power = sum(abs2.(signal))
    summed_abs_inphase_ampl = sum(abs.(real.(signal)))
    prompt_correlator_output = complex(1,1)
    cn0_state = Tracking.CN0State(summed_total_power, summed_abs_inphase_ampl, NaN, 20ms, 19ms)
    new_cn0_state = Tracking.estimate_CN0(cn0_state, 1ms, prompt_correlator_output)
    @test new_cn0_state.summed_total_power == 0.0
    @test new_cn0_state.summed_abs_inphase_ampl == 0.0
    @test new_cn0_state.last_valid_cn0 == 1000.0
    @test new_cn0_state.update_time == 20ms
    @test new_cn0_state.current_updated_time == 0ms
end

@testset "Tracking: CN0 estimation" begin
    Random.seed!(1234)
    num_samples = 55000
    sample_freq = 2.5e6Hz
    carrier = cis.(2π * (1:num_samples) .* 50Hz ./ sample_freq .+ 1.2)
    code = get_code.(GPSL1, (1:num_samples) .* 1023e3Hz ./ sample_freq .+ 2.0, 1)
    signal = carrier .* code .+ randn(ComplexF64, num_samples) .* 10 .^ (20 ./ 20)
    correlator_outputs = zeros(SVector{3,ComplexF64})
    code_shift = Tracking.CodeShift(GPSL1, sample_freq, 0.5)
    inits = TrackingInitials(20Hz, 1.2, 0.0Hz, 2.0)
    dopplers = Tracking.Dopplers(inits)
    phases = Tracking.Phases(inits)
    carrier_loop = Tracking.init_3rd_order_bilinear_loop_filter(18Hz)
    code_loop = Tracking.init_2nd_order_bilinear_loop_filter(1Hz)
    last_valid_correlator_outputs = zeros(typeof(correlator_outputs))
    last_valid_filtered_correlator_outputs = zeros(typeof(correlator_outputs))
    data_bits = Tracking.DataBits(GPSL1)
    cn0_state = Tracking.CN0State(20ms)
    results = @inferred Tracking._tracking(GPSL1, correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, signal, sample_freq, 30Hz, inits, dopplers, phases, code_shift, carrier_loop, code_loop, 1, x -> x, 0.5ms, 1ms, 1, 0, 0, data_bits, cn0_state, 0.0Hz)
    @test 10 * log10(results[2].cn0) ≈ 45 atol = 2.6 #1?
end
