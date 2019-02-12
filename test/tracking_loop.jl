@testset "Code shift" begin
    gpsl1 = GPSL1()
    code_shift = Tracking.CodeShift{3}(gpsl1, 4e6Hz, 0.5)
    @test code_shift.samples == 2
    @test code_shift.actual_shift == 0.5115

    @test Tracking.init_shifts(code_shift) == [-2,0,2]
end

@testset "Downconvert and correlate" begin

    gpsl1 = GPSL1()

    @test @inferred(Tracking.wrap_code_idx(gpsl1, 0, 0)) == (0, 0)
    @test @inferred(Tracking.wrap_code_idx(gpsl1, 1023, 0)) == (0, -1023)
    @test @inferred(Tracking.wrap_code_idx(gpsl1, -1, 0)) == (1022, 0)
    @test @inferred(Tracking.wrap_code_idx(gpsl1, -4, 0)) == (1019, 0)
    @test @inferred(Tracking.wrap_code_idx(gpsl1, 1026, -1023)) == (3, -1023 * 2)

    carrier = gen_carrier.(1:4000, 50Hz, 1.2, 4e6Hz)
    code = gen_code.(Ref(gpsl1), 1:4000, 1023e3Hz, 2.0, 4e6Hz, 1)
    signal = carrier .* code
    gen_carrier_replica(x) = gen_carrier(x, -50Hz, -1.2, 4e6Hz)
    calc_code_replica_phase_unsafe(x) = calc_code_phase_unsafe(x, 1023e3Hz, 2.0, 4e6Hz)
    prev_correlator_outputs = zeros(SVector{3,ComplexF64})
    code_shift = Tracking.CodeShift{3}(gpsl1, 4e6Hz, 0.5)
    outputs = @inferred Tracking.downconvert_and_correlate(prev_correlator_outputs, signal, gpsl1, 1, 4000, gen_carrier_replica, calc_code_replica_phase_unsafe, code_shift, 1)
    @test outputs ≈ [1952, 4000, 1952]

    dopplers = Tracking.Dopplers(10Hz, 0Hz)
    phases = Tracking.Phases(1.2, 2.0)
    outputs = @inferred Tracking.correlate_and_dump(prev_correlator_outputs, signal, gpsl1, 4e6Hz, 40Hz, dopplers, phases, code_shift, 1, 4000, 1)
    @test outputs ≈ [1952, 4000, 1952]
end

@testset "Aiding" begin
    gpsl1 = GPSL1()
    inits = Initials(20.0Hz, 0.0, 1.0Hz, 0)

    dopplers = @inferred Tracking.aid_dopplers(gpsl1, inits, 2.5Hz, 0.1Hz, 0.0Hz)
    @test dopplers.carrier == 22.5Hz
    @test dopplers.code == 1.1Hz + 2.5Hz / 1540
end

@testset "Correlator output" begin
    gpsl1 = GPSL1()
    code_shift = Tracking.CodeShift{3}(gpsl1, 4e6Hz, 0.5)
    output = @inferred Tracking.init_correlator_outputs(code_shift)
    @test output === zeros(SVector{3, ComplexF64})
end

@testset "Integration time" begin
    @test @inferred(Tracking.calc_actual_integration_time(4000, 4e6)) == 1e-3

    gpsl1 = GPSL1()
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, Tracking.Phases(0.0, 100.0), 1ms)
    @test num_samples == ceil(Int, (1023 - 100) * 4e6Hz / 1023e3Hz)
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, Tracking.Phases(0.0, 0.0), 1ms)
    @test num_samples == 4000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, Tracking.Phases(0.0, 1023), 1ms)
    @test num_samples == 4000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, Tracking.Phases(0.0, 1022.8), 1ms)
    @test num_samples == 1
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, Tracking.Phases(0.0, 1023.1), 1ms)
    @test num_samples == 4000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, Tracking.Phases(0.0, 1024.0), 1ms)
    @test num_samples == ceil(Int, (1023 - 1) * 4e6Hz / 1023e3Hz)

    signal = zeros(1000)
    @test @inferred(Tracking.calc_num_samples_signal_bound(signal, 100)) == 901
    @test @inferred(Tracking.calc_num_samples_signal_bound(signal, 1)) == 1000
    @test @inferred(Tracking.calc_num_samples_signal_bound(signal, 1000)) == 1

    data_bits = Tracking.DataBits(gpsl1)
    @test @inferred(Tracking.calc_integration_time(data_bits, 2ms)) == 1ms
    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 0, 0)
    @test @inferred(Tracking.calc_integration_time(data_bits, 2ms)) == 2ms
end

@testset "Phases" begin
    gpsl1 = GPSL1()
    dopplers = Tracking.Dopplers(20Hz, 1Hz)
    phases = Tracking.Phases(0.0, 100.0)

    data_bits = Tracking.DataBits(gpsl1)
    adjusted_code_phase = @inferred Tracking.adjust_code_phase(gpsl1, data_bits, 100)
    @test adjusted_code_phase == 100

    next_phases = @inferred Tracking.calc_next_phases(gpsl1, 50Hz, 4e6Hz, dopplers, phases, 200, data_bits)
    @test next_phases.carrier ≈ 2π * 200 * (50Hz + 20Hz) / 4e6Hz + 0.0
    @test next_phases.code ≈ 200 * (1023e3Hz + 1Hz) / 4e6Hz + 100.0

    gpsl5 = GPSL5()
    data_bits = Tracking.DataBits(gpsl5)
    adjusted_code_phase = @inferred Tracking.adjust_code_phase(gpsl5, data_bits, 10430)
    @test adjusted_code_phase == 200

    next_phases = @inferred Tracking.calc_next_phases(gpsl5, 50Hz, 4e7Hz, dopplers, phases, 50000, data_bits)
    @test next_phases.carrier ≈ 2π * 50000 * (50Hz + 20Hz) / 4e7Hz + 0.0
    @test next_phases.code ≈ 50000 * (1023e4Hz + 1Hz) / 4e7Hz + 100.0 - 10230

    data_bits = Tracking.DataBits{GPSL5}(0, 0, 10, 0.0, 0, 0)
    adjusted_code_phase = @inferred Tracking.adjust_code_phase(gpsl5, data_bits, 10430)
    @test adjusted_code_phase == 10430
    next_phases = @inferred Tracking.calc_next_phases(gpsl5, 50Hz, 4e7Hz, dopplers, phases, 50000, data_bits)
    @test next_phases.code ≈ 50000 * (1023e4Hz + 1Hz) / 4e7Hz + 100.0
end


@testset "Tracking" begin
    gpsl1 = GPSL1()
    carrier = gen_carrier.(1:24000, 50Hz, 1.2, 4e6Hz)
    code = gen_code.(Ref(gpsl1), 1:24000, 1023e3Hz, 2.0, 4e6Hz, 1)
    signal = carrier .* code
    correlator_outputs = zeros(SVector{3,ComplexF64})
    code_shift = Tracking.CodeShift{3}(gpsl1, 4e6Hz, 0.5)
    inits = Initials(20Hz, 1.2, 0.0Hz, 2.0)
    dopplers = Tracking.Dopplers(inits)
    phases = Tracking.Phases(inits)
    carrier_loop = Tracking.init_3rd_order_bilinear_loop_filter(18Hz)
    code_loop = Tracking.init_2nd_order_bilinear_loop_filter(1Hz)
    last_valid_correlator_outputs = zeros(typeof(correlator_outputs))
    data_bits = Tracking.DataBits(gpsl1)
    cn0_buffer = Tracking.CN0Buffer(20ms)
    results = @inferred Tracking._tracking(correlator_outputs, last_valid_correlator_outputs, signal, gpsl1, 4e6Hz, 30Hz, inits, dopplers, phases, code_shift, carrier_loop, code_loop, 1, x -> x, 0.5ms, 1ms, 1, 0, 0, data_bits, cn0_buffer, 0.0Hz)
    @test results[2].carrier_doppler ≈ 20Hz
    @test results[2].code_doppler ≈ 0Hz atol = 3e-3Hz #??
    @test results[2].code_phase ≈ 2 atol = 2e-5
end


@testset "Track L1" begin

     center_freq = 1.57542e9Hz
     code_freq = 1023e3Hz
     doppler = 10Hz
     interm_freq = 50Hz
     carrier_phase = π / 3
     sample_freq = 4e6Hz
     code_doppler = doppler * code_freq / center_freq
     data_doppler = doppler * 50Hz / center_freq
     code_phase = 2.0
     min_integration_time = 0.5ms

     run_time = 1500e-3s
     integration_time = 1e-3s
     num_integrations = convert(Int, run_time / integration_time)
     num_samples = convert(Int, run_time * sample_freq)
     integration_samples = convert(Int, integration_time * sample_freq)

     gps_l1 = GPSL1()

     carrier = cis.(2π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
     sampled_code = gen_code.(Ref(gps_l1), 1:num_samples, code_doppler + code_freq, code_phase, sample_freq, 1)
     signal = carrier .* sampled_code

     inits = Initials(0.0Hz, carrier_phase, 0.0Hz, code_phase)
     track = @inferred init_tracking(gps_l1, inits, sample_freq, interm_freq, 1, min_integration_time = min_integration_time, max_integration_time = integration_time)

     code_dopplers = zeros(num_integrations)
     code_phases = zeros(num_integrations)
     calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 1023)
     carrier_dopplers = zeros(num_integrations)

     results = nothing
     for i = 1:num_integrations
         current_signal = signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,2), randn(integration_samples,2)) .* 10^(5/20)
         track, results = @inferred track(current_signal)
         @test results.num_processed_prns == 1
         @test results.num_bits == 0
         @test results.data_bits == 0
         code_phases[i] = results.code_phase
         carrier_dopplers[i] = results.carrier_doppler / Hz
         code_dopplers[i] = results.code_doppler / Hz
     end

     @test results.carrier_doppler ≈ doppler atol = 4e-4Hz
     @test mod(results.code_phase, 1023) ≈ calculated_code_phases[end] atol = 3e-6
     @test results.code_doppler ≈ code_doppler atol = 4e-6Hz

#=
     using PyPlot
     pygui(true)
     figure("Tracking code phase error")
     plot(code_phases - calculated_code_phases)
     figure("Tracking carrier_dopplers")
     plot(carrier_dopplers)
     figure("Tracking code dopplers")
     plot(code_dopplers)
=#

     one_prn_code = gen_code.(Ref(gps_l1), 0:1022, 1Hz, 0.0, 1Hz, 1)
     switching_prn_code = vcat(repeat(-one_prn_code, 20), repeat(one_prn_code, 20))
     sampled_switching_code = switching_prn_code[mod.(floor.(Int, (1:num_samples) .* (code_doppler + code_freq) ./ sample_freq .+ code_phase), 40920) .+ 1]
     signal = carrier .* sampled_switching_code
     track = @inferred init_tracking(gps_l1, inits, sample_freq, interm_freq, 1, pll_bandwidth = 18.0Hz, dll_bandwidth = 1.0Hz, min_integration_time = min_integration_time, max_integration_time = integration_time)
     for i = 1:num_integrations
         current_signal = signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,2), randn(integration_samples,2)) .* 10^(5/20)
         track, results = @inferred track(current_signal)
         @test results.num_processed_prns == 1
         if i-1 % 20 == 0 && i > 1
             @test results.num_bits == 1
             @test results.data_bits == (i % 40 == 0)
         else
             @test results.num_bits == 0
             @test results.data_bits == 0
         end
         code_phases[i] = results.code_phase
         carrier_dopplers[i] = results.carrier_doppler / Hz
         code_dopplers[i] = results.code_doppler / Hz
     end
     @test results.carrier_doppler ≈ doppler atol = 4e-4Hz
     @test mod(results.code_phase, 1023) ≈ calculated_code_phases[end] atol = 5e-5
     @test results.code_doppler ≈ code_doppler atol = 2e-3Hz
end

@testset "Track L5" begin

     center_freq = 1.17645e9Hz
     code_freq = 10230e3Hz
     doppler = 10Hz
     interm_freq = 50Hz
     carrier_phase = π / 3
     sample_freq = 40e6Hz
     code_doppler = doppler * code_freq / center_freq
     code_phase = 65950.0
     min_integration_time = 0.5ms

     run_time = 500e-3s
     integration_time = 1e-3s
     num_integrations = convert(Int, run_time / integration_time)
     num_samples = convert(Int, run_time * sample_freq)
     integration_samples = convert(Int, integration_time * sample_freq)

     gps_l5 = GPSL5()

     carrier = cis.(2 * π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
     sampled_code = gen_code.(Ref(gps_l5), 1:num_samples, code_doppler + code_freq, code_phase, sample_freq, 1)
     signal = carrier .* sampled_code

     inits = Initials(0.0Hz, carrier_phase, 0.0Hz, mod(code_phase, 10230))
     track = @inferred init_tracking(gps_l5, inits, sample_freq, interm_freq, 1, pll_bandwidth = 18.0Hz, dll_bandwidth = 1.0Hz, min_integration_time = min_integration_time, max_integration_time = integration_time)

     code_dopplers = zeros(num_integrations)
     code_phases = zeros(num_integrations)
     calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 102300)
     carrier_dopplers = zeros(num_integrations)
     real_prompts = zeros(num_integrations)

     results = nothing
     for i = 1:num_integrations
         current_signal = signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,2), randn(integration_samples,2)) .* 10^(5/20)
         track, results = @inferred track(current_signal)
         @test results.num_processed_prns == 1
         if i-1 % 10 == 0 && i > 1
             @test results.num_bits == 1
             @test results.data_bits == 1
         else
             @test results.num_bits == 0
             @test results.data_bits == 0
         end
         code_phases[i] = results.code_phase
         carrier_dopplers[i] = results.carrier_doppler / Hz
         code_dopplers[i] = results.code_doppler / Hz
         real_prompts[i] = real(prompt(results.prompt))
     end

     @test results.carrier_doppler[1] ≈ doppler atol = 5e-2Hz
     @test mod(results.code_phase[1], 102300) ≈ calculated_code_phases[end] atol = 5e-5
     @test results.code_doppler[1] ≈ code_doppler atol = 5e-4Hz

     #=
     using PyPlot
     pygui(true)
     figure("Tracking code phases error")
     plot(code_phases - calculated_code_phases)
     figure("Tracking carrier_dopplers")
     plot(carrier_dopplers)
     figure("Tracking code dopplers")
     plot(code_dopplers)
     figure("Prompts")
     plot(real_prompts)
     =#

end
