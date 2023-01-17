@testset "Downconvert and Correlator" begin
    gpsl1 = GPSL1()

    system_sats_states = (SystemSatsState(
        gpsl1,
        [
            SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz)
        ]
    ),)

    downconvert_and_correlator = @inferred CPUDownconvertAndCorrelator(system_sats_states, 5000)
    @test size(downconvert_and_correlator.buffers[1][1].code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator.buffers[1][1].code_replica_buffer) == Int16
    @test size(downconvert_and_correlator.buffers[1][1].carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator.buffers[1][1].carrier_replica_buffer) == ComplexF32
    @test size(downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer) == (5000,)
    @test eltype(downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer) == ComplexF32

    downconvert_and_correlator_f64 = @inferred CPUDownconvertAndCorrelator(Float64, system_sats_states, 5000)
    @test size(downconvert_and_correlator_f64.buffers[1][1].code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator_f64.buffers[1][1].code_replica_buffer) == Int16
    @test size(downconvert_and_correlator_f64.buffers[1][1].carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator_f64.buffers[1][1].carrier_replica_buffer) == ComplexF64
    @test size(downconvert_and_correlator_f64.buffers[1][1].downconvert_signal_buffer) == (5000,)
    @test eltype(downconvert_and_correlator_f64.buffers[1][1].downconvert_signal_buffer) == ComplexF64

    system_sats_states_numants2 = (SystemSatsState(
        gpsl1,
        [
            SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz, num_ants = NumAnts(2))
        ]
    ),)

    downconvert_and_correlator_numants2 = @inferred CPUDownconvertAndCorrelator(system_sats_states_numants2, 5000)
    @test size(downconvert_and_correlator_numants2.buffers[1][1].code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator_numants2.buffers[1][1].code_replica_buffer) == Int16
    @test size(downconvert_and_correlator_numants2.buffers[1][1].carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator_numants2.buffers[1][1].carrier_replica_buffer) == ComplexF32
    @test size(downconvert_and_correlator_numants2.buffers[1][1].downconvert_signal_buffer) == (5000, 2)
    @test eltype(downconvert_and_correlator_numants2.buffers[1][1].downconvert_signal_buffer) == ComplexF32

    galileo_e1b = GalileoE1B()
    system_sats_states = (
        SystemSatsState(
            gpsl1,
            [
                SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz)
            ]
        ),
        SystemSatsState(
            galileo_e1b,
            [
                SatState(galileo_e1b, 1, 5e6Hz, 10.0, 500.0Hz)
            ]
        )
    )

    downconvert_and_correlator2 = @inferred CPUDownconvertAndCorrelator(system_sats_states, 5000)
    @test size(downconvert_and_correlator2.buffers[1][1].code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator2.buffers[1][1].code_replica_buffer) == Int16
    @test size(downconvert_and_correlator2.buffers[1][1].carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator2.buffers[1][1].carrier_replica_buffer) == ComplexF32
    @test size(downconvert_and_correlator2.buffers[1][1].downconvert_signal_buffer) == (5000,)
    @test eltype(downconvert_and_correlator2.buffers[1][1].downconvert_signal_buffer) == ComplexF32

    @test size(downconvert_and_correlator2.buffers[2][1].code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator2.buffers[2][1].code_replica_buffer) == Float32
    @test size(downconvert_and_correlator2.buffers[2][1].carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator2.buffers[2][1].carrier_replica_buffer) == ComplexF32
    @test size(downconvert_and_correlator2.buffers[2][1].downconvert_signal_buffer) == (5000,)
    @test eltype(downconvert_and_correlator2.buffers[2][1].downconvert_signal_buffer) == ComplexF32
end

@testset "Downconvert and correlate" begin

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    intermediate_frequency = 0.0Hz

    system_sats_state = SystemSatsState(
        gpsl1,
        [
            SatState(gpsl1, 1, sampling_frequency, code_phase, 1000.0Hz),
            SatState(gpsl1, 2, sampling_frequency, 11.0, 500.0Hz)
        ]
    )
    system_sats_states = (system_sats_state,)

    downconvert_and_correlator = @inferred CPUDownconvertAndCorrelator(system_sats_states, num_samples_signal)

    system_sats_sample_params = Tracking.init_sample_params(
        system_sats_states,
        1,
    )
    next_system_sats_sample_params = Tracking.calc_sample_params(
        system_sats_states,
        system_sats_sample_params,
        num_samples_signal,
        sampling_frequency,
        1,
    )

    signal = gen_code(
            num_samples_signal,
            gpsl1,
            1,
            sampling_frequency,
            get_code_frequency(gpsl1) + 1000Hz * get_code_center_frequency_ratio(gpsl1),
            code_phase
        ) .*
        cis.(2π * (0:num_samples_signal-1) * 1000.0Hz / sampling_frequency)

    correlators = @inferred Tracking.downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        sampling_frequency,
        intermediate_frequency,
        system_sats_states,
        next_system_sats_sample_params,
    )

    @test real.(correlators[1][1].accumulators) ≈ [2921, 4949, 2917]

    signal = gen_code(
            num_samples_signal,
            gpsl1,
            2,
            sampling_frequency,
            get_code_frequency(gpsl1) + 500Hz * get_code_center_frequency_ratio(gpsl1),
            11.0
        ) .*
        cis.(2π * (0:num_samples_signal-1) * 500.0Hz / sampling_frequency)

    correlators = @inferred Tracking.downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        sampling_frequency,
        intermediate_frequency,
        system_sats_states,
        next_system_sats_sample_params,
    )

    @test real.(correlators[1][2].accumulators) ≈ [2919, 4947, 2915]
end