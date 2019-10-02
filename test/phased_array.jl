@testset "Phased array tracking" begin
    num_samples = 24000
    sample_freq = 4e6Hz
    carrier = cis.(2π * (1:num_samples) .* 50Hz ./ sample_freq .+ 1.2)
    code = get_code.(GPSL1, (1:num_samples) .* 1023e3Hz ./ sample_freq .+ 2.0, 1)
    signal = carrier .* code * [1 1 1 1]
    correlator_outputs = zeros(SMatrix{3,4,ComplexF64})
    code_shift = Tracking.CodeShift(GPSL1, 4e6Hz, 0.5)
    inits = TrackingInitials(20Hz, 1.2, 0.0Hz, 2.0)
    dopplers = Tracking.Dopplers(inits)
    phases = Tracking.Phases(inits)
    carrier_loop = init_3rd_order_bilinear_loop_filter(18Hz)
    code_loop = init_2nd_order_bilinear_loop_filter(1Hz)
    last_valid_correlator_outputs = copy(correlator_outputs)
    last_valid_filtered_correlator_outputs = zeros(SVector{3,ComplexF64})
    data_bits = Tracking.DataBits(GPSL1)
    cn0_state = Tracking.CN0State(20ms)
    results = @inferred Tracking._tracking(GPSL1, correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, signal, 4e6Hz, 30Hz, inits, dopplers, phases, code_shift, carrier_loop, code_loop, 1, x -> x[:,1], 0.5ms, 1ms, 1, 0, 0, data_bits, cn0_state, 0.0Hz)
    @test results[2].carrier_doppler ≈ 20Hz
    @test results[2].code_doppler ≈ 0Hz atol = 3e-3Hz #??
    @test results[2].code_phase ≈ 2 atol = 2e-5
end

@testset "Phased array track L1" begin

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

     carrier = cis.(2π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
     sampled_code = get_code.(GPSL1, (1:num_samples) .* (code_doppler + code_freq) ./ sample_freq .+ code_phase, 1)
     signal = carrier .* sampled_code * [1 1 1 1]

     inits = TrackingInitials(0.0Hz, carrier_phase, 0.0Hz, code_phase)
     track = @inferred Tracking.init_tracking(GPSL1, inits, sample_freq, interm_freq, 1, num_ants = NumAnts(4), min_integration_time = min_integration_time, max_integration_time = integration_time)

     code_dopplers = zeros(num_integrations)
     code_phases = zeros(num_integrations)
     calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 1023)
     carrier_dopplers = zeros(num_integrations)

     results = nothing
     for i = 1:num_integrations
         current_signal = signal[integration_samples * (i - 1) + 1:integration_samples * i, :]# .+ complex.(randn(integration_samples,2), randn(integration_samples,2)) .* 10^(5/20)
         track, results = @inferred track(current_signal, x -> x[:,1])
         @test results.num_processed_prns == 1
         @test results.num_bits == 0
         @test results.data_bits == 0
         @test prompt(results.correlator_outputs)[1] ≈ prompt(results.correlator_outputs)[2]
         @test prompt(results.correlator_outputs)[2] ≈ prompt(results.correlator_outputs)[3]
         @test prompt(results.correlator_outputs)[3] ≈ prompt(results.correlator_outputs)[4]
         code_phases[i] = results.code_phase
         carrier_dopplers[i] = results.carrier_doppler / Hz
         code_dopplers[i] = results.code_doppler / Hz
     end

     @test results.carrier_doppler ≈ doppler atol = 4e-4Hz
     @test mod(results.code_phase, 1023) ≈ calculated_code_phases[end] atol = 3e-6
     @test results.code_doppler ≈ code_doppler atol = 4e-6Hz
end

@testset "Phased array prompt early late" begin
    A = @SMatrix [i*j for i = 1:3, j = 1:4]
    @test Tracking.early(A) == A[3,:]
    @test Tracking.prompt(A) == A[2,:]
    @test Tracking.late(A) == A[1,:]

    A = @SMatrix [i*j for i = 1:5, j = 1:4]
    @test Tracking.veryearly(A) == A[5,:]
    @test Tracking.early(A) == A[4,:]
    @test Tracking.prompt(A) == A[3,:]
    @test Tracking.late(A) == A[2,:]
    @test Tracking.verylate(A) == A[1,:]
end
