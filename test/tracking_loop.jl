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
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 1000, Tracking.Phases(0.0, 100.0), 1ms)
    @test num_samples == ceil(Int, (1023 - 100) * 4e6Hz / 1023e3Hz) - 1000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 0, Tracking.Phases(0.0, 0.0), 1ms)
    @test num_samples == 4000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 0, Tracking.Phases(0.0, 1023), 1ms)
    @test num_samples == 4000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 0, Tracking.Phases(0.0, 1022.8), 1ms)
    @test num_samples == 1
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 0, Tracking.Phases(0.0, 1023.1), 1ms)
    @test num_samples == 4000
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 0, Tracking.Phases(0.0, 1024.0), 1ms)
    @test num_samples == ceil(Int, (1023 - 1) * 4e6Hz / 1023e3Hz)
    num_samples = @inferred Tracking.calc_num_samples_left_to_integrate(gpsl1, 4e6Hz, 1, Tracking.Phases(0.0, 1022.8), 1ms)
    @test num_samples == 0

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
struct DataBits{T<:AbstractGNSSSystem}
    synchronisation_buffer::UInt
    num_bits_in_synchronisation_buffer::UInt
    first_found_after_num_prns::Int
    prompt_accumulator::Float64
    buffer::UInt
    num_bits_in_buffer::UInt
end
@testset "Buffer data bits" begin
    gpsl1 = GPSL1()
    gpsl5 = GPSL5()

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 0, 0)) == false
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 12 << 41 + 2^20-1, 40)) == true
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 12 << 41 + (2^20-1) << 20, 40)) == true
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL1, 12 << 41 + 2^20-1, 20)) == false

    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL5, 12 << 11 + 0x35, 10)) == true # 0x35 == 0000110101
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL5, 12 << 11 + 0x35, 7)) == false
    @test @inferred(Tracking.is_upcoming_integration_new_bit(GPSL5, 12 << 11 + 0x3ca, 10)) == true # 0x3ca == 1111001010

    data_bits = Tracking.DataBits(gpsl1)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl1, 5.0, 1)
    @test next_data_bits.synchronisation_buffer == 1
    @test next_data_bits.num_bits_in_synchronisation_buffer == 1
    @test next_data_bits.first_found_after_num_prns == -1
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL1}(12 << 40 + 2^19-1, 39, -1, 0.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl1, 5.0, 3)
    @test next_data_bits.synchronisation_buffer == 12 << 41 + 2^20-1
    @test next_data_bits.num_bits_in_synchronisation_buffer == 40
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL5}(12 << 10 + 0x1a, 9, -1, 0.0, 0, 0) # 0x1a == 000011010
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 3)
    @test next_data_bits.synchronisation_buffer == 12 << 11 + 0x35 # 0x35 == 0000110101
    @test next_data_bits.num_bits_in_synchronisation_buffer == 10
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 0.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 3)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 5.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 0

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 4)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 1
    @test next_data_bits.num_bits_in_buffer == 1

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 0, 0)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, -10.0, 4)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 0
    @test next_data_bits.num_bits_in_buffer == 1

    data_bits = Tracking.DataBits{GPSL1}(0, 20, 4, 8.0, 1, 1)
    next_data_bits = @inferred Tracking.buffer(data_bits, gpsl5, 5.0, 4)
    @test next_data_bits.synchronisation_buffer == 0
    @test next_data_bits.num_bits_in_synchronisation_buffer == 20
    @test next_data_bits.first_found_after_num_prns == 4
    @test next_data_bits.prompt_accumulator == 0.0
    @test next_data_bits.buffer == 3
    @test next_data_bits.num_bits_in_buffer == 2
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
    filtered_prompt_correlator_buffer = Tracking.init_prompt_correlator_buffer(gpsl1)
    last_valid_correlator_outputs = zeros(typeof(correlator_outputs))
    data_bits = Tracking.DataBits(gpsl1)
    results = @inferred Tracking._tracking(correlator_outputs, last_valid_correlator_outputs, signal, gpsl1, 4e6Hz, 30Hz, inits, dopplers, phases, code_shift, carrier_loop, code_loop, 1, x -> x, 0.5ms, 1ms, 1, 0, UInt(0), data_bits, 0.0Hz)
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
     code_phase = 2.0
     min_integration_time = 0.5ms

     run_time = 2500e-3s
     integration_time = 1e-3s
     num_integrations = convert(Int, run_time / integration_time)
     num_samples = convert(Int, run_time * sample_freq)
     integration_samples = convert(Int, integration_time * sample_freq)

     gps_l1 = GPSL1()

     carrier = cis.(2π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
     sampled_code = gen_code.(Ref(gps_l1), 1:num_samples, code_doppler + code_freq, code_phase, sample_freq, 1)
     signal = carrier .* sampled_code

     inits = Initials(0.0Hz, carrier_phase, 0.0Hz, code_phase)
     track = init_tracking(gps_l1, inits, sample_freq, interm_freq, 18.0Hz, 1.0Hz, min_integration_time, integration_time, 1)

     code_dopplers = zeros(num_integrations)
     code_phases = zeros(num_integrations)
     calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 1023)
     carrier_dopplers = zeros(num_integrations)

     results = nothing
     for i = 1:num_integrations
         current_signal = signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,2), randn(integration_samples,2)) .* 10^(5/20)
         track, results = track(current_signal)
         code_phases[i] = results.code_phase
         carrier_dopplers[i] = results.carrier_doppler / Hz
         code_dopplers[i] = results.code_doppler / Hz
     end

     @test results.carrier_doppler ≈ doppler atol = 5e-2Hz
     @test mod(results.code_phase, 1023) ≈ calculated_code_phases[end] atol = 3e-4
     @test results.code_doppler ≈ code_doppler atol = 2e-3Hz

#=
     figure("Tracking code phases")
     plot(code_phases, color = "blue")
     plot(calculated_code_phases, color = "red")
     figure("Tracking carrier_dopplers")
     plot(carrier_dopplers)
     figure("Tracking code dopplers")
     plot(code_dopplers)
=#

#     # Track all at once
#     track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#     track, results = track([1, 1]' .* signal, beamform)
#
#     @test code_phases == results.code_phase
#     @test carrier_dopplers == results.carrier_doppler ./ Hz
#     @test code_dopplers == results.code_doppler ./ Hz
#
#     track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#     track, results = track([1, 1]' .* signal[1:40000], beamform)
#
#     code_dopplers2 = CircularBuffer{Float64}(num_integrations)
#     code_phases2 = CircularBuffer{Float64}(num_integrations)
#     carrier_dopplers2 = CircularBuffer{Float64}(num_integrations)
#     track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#
#     # Track with time smaller than integration_time
#     for i = 1:2000:size(signal, 1)
#         track, results = track([1, 1]' .* signal[i:i + 2000 - 1], beamform)
#         append!(code_dopplers2, results.code_doppler ./ Hz)
#         append!(code_phases2, results.code_phase)
#         append!(carrier_dopplers2, results.carrier_doppler ./ Hz)
#     end
#
#     @test code_phases ≈ code_phases2
#     @test carrier_dopplers ≈ carrier_dopplers2
#     @test code_dopplers ≈ code_dopplers2
#
#     inits = Initials(0.0Hz, carrier_phase + 1.0, 0.0Hz, code_phase + 0.1)
#     track = init_tracking(gps_l1, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#     track, results = track([1, 1]' .* signal, beamform)
#
#     @test results.code_phase[end] ≈ calculated_code_phases[end] atol = 4e-3
#
#     #=
#     figure("Tracking code phases")
#     plot(results.code_phase, color = "blue")
#     plot(calculated_code_phases, color = "red")
#     figure("Tracking carrier_dopplers")
#     plot(results.carrier_doppler ./ Hz)
#     figure("Tracking code dopplers")
#     plot(results.code_doppler ./ Hz)
#     =#
end
#
# @testset "Track L5" begin
#
#     function beamform(x)
#         dot([0.5, 0.5], x)
#     end
#
#     center_freq = 1.17645e9Hz
#     code_freq = 10230e3Hz
#     doppler = 10Hz
#     interm_freq = 50Hz
#     carrier_phase = π / 3
#     sample_freq = 40e6Hz
#     code_doppler = doppler * code_freq / center_freq
#     code_phase = 65950.0
#
#     run_time = 500e-3s
#     integration_time = 1e-3s
#     num_integrations = convert(Int, run_time / integration_time)
#     num_samples = convert(Int, run_time * sample_freq)
#     integration_samples = convert(Int, integration_time * sample_freq)
#
#     gps_l5 = GPSL5()
#
#     carrier = cis.(2 * π * (interm_freq + doppler) / sample_freq * (1:num_samples) .+ carrier_phase)
#     sampled_code = gen_code.(Ref(gps_l5), 1:num_samples, code_doppler + code_freq, code_phase, sample_freq, 1)
#     signal = carrier .* sampled_code
#
#     inits = Initials(0.0Hz, carrier_phase, 0.0Hz, mod(code_phase, 10230))
#     track = init_tracking(gps_l5, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#
#     code_dopplers = zeros(num_integrations)
#     code_phases = zeros(num_integrations)
#     calculated_code_phases  = mod.((1:num_integrations) * integration_samples * (code_doppler + code_freq) / sample_freq .+ code_phase, 102300)
#     carrier_dopplers = zeros(num_integrations)
#     real_prompts = zeros(num_integrations)
#
#     results = nothing
#     for i = 1:num_integrations
#         current_signal = [1, 1]' .* signal[integration_samples * (i - 1) + 1:integration_samples * i]# .+ complex.(randn(integration_samples,2), randn(integration_samples,2)) .* 10^(5/20)
#         track, results = track(current_signal, beamform)
#         code_phases[i] = results.code_phase[1]
#         carrier_dopplers[i] = results.carrier_doppler[1] / Hz
#         code_dopplers[i] = results.code_doppler[1] / Hz
#         real_prompts[i] = [1, 1]' * real.(results.prompt[1])
#     end
#
#     @test results.carrier_doppler[1] ≈ doppler atol = 5e-2Hz
#     @test mod(results.code_phase[1], 102300) ≈ calculated_code_phases[end] atol = 5e-5
#     @test results.code_doppler[1] ≈ code_doppler atol = 5e-4Hz
#
#     #=
#     figure("Tracking code phases error")
#     plot(code_phases - calculated_code_phases)
#     figure("Tracking carrier_dopplers")
#     plot(carrier_dopplers)
#     figure("Tracking code dopplers")
#     plot(code_dopplers)
#     figure("Prompts")
#     plot(real_prompts)
#     =#
#
#     # Track all at once
#     track = init_tracking(gps_l5, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#     track, results = track([1, 1]' .* signal, beamform)
#
#     @test code_phases == results.code_phase
#     @test carrier_dopplers == results.carrier_doppler ./ Hz
#     @test code_dopplers == results.code_doppler ./ Hz
#
#     track = init_tracking(gps_l5, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#     track, results = track([1, 1]' .* signal[1:40000], beamform)
#
#     code_dopplers2 = CircularBuffer{Float64}(num_integrations)
#     code_phases2 = CircularBuffer{Float64}(num_integrations)
#     carrier_dopplers2 = CircularBuffer{Float64}(num_integrations)
#     track = init_tracking(gps_l5, inits, integration_time, sample_freq, interm_freq, 18.0Hz, 1.0Hz, 1)
#
#     # Track with time smaller than integration_time
#     for i = 1:20000:size(signal, 1)
#         track, results = track([1, 1]' .* signal[i:i + 20000 - 1], beamform)
#         append!(code_dopplers2, results.code_doppler ./ Hz)
#         append!(code_phases2, results.code_phase)
#         append!(carrier_dopplers2, results.carrier_doppler ./ Hz)
#     end
#
#     @test code_phases ≈ code_phases2
#     @test carrier_dopplers ≈ carrier_dopplers2 atol = 5e-3
#     @test code_dopplers ≈ code_dopplers2 atol = 5e-4
#     @test carrier_dopplers[end] ≈ carrier_dopplers2[end]
#     @test code_dopplers[end] ≈ code_dopplers2[end]
# end
