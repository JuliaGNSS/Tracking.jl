@testset "Downconvert and Correlator" begin
    gpsl1 = GPSL1()

    correlator = get_default_correlator(gpsl1, 5e6Hz)

    downconvert_and_correlator =
        @inferred CPUSatDownconvertAndCorrelator(gpsl1, correlator, 5000)

    @test size(downconvert_and_correlator.code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator.code_replica_buffer) == Int16
    @test size(downconvert_and_correlator.carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator.carrier_replica_buffer) == ComplexF32
    @test size(downconvert_and_correlator.downconvert_signal_buffer) == (5000,)
    @test eltype(downconvert_and_correlator.downconvert_signal_buffer) == ComplexF32

    downconvert_and_correlator_f64 =
        @inferred CPUSatDownconvertAndCorrelator(Float64, gpsl1, correlator, 5000)
    @test size(downconvert_and_correlator_f64.code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator_f64.code_replica_buffer) == Int16
    @test size(downconvert_and_correlator_f64.carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator_f64.carrier_replica_buffer) == ComplexF64
    @test size(downconvert_and_correlator_f64.downconvert_signal_buffer) == (5000,)
    @test eltype(downconvert_and_correlator_f64.downconvert_signal_buffer) == ComplexF64

    system_sats_states_numants2 = (
        SystemSatsState(
            gpsl1,
            [SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz; num_ants = NumAnts(2))],
        ),
    )

    correlator_numants2 = get_default_correlator(gpsl1, 5e6Hz, NumAnts(2))

    downconvert_and_correlator_numants2 =
        @inferred CPUSatDownconvertAndCorrelator(gpsl1, correlator_numants2, 5000)

    @test size(downconvert_and_correlator_numants2.code_replica_buffer) == (5004,)
    @test eltype(downconvert_and_correlator_numants2.code_replica_buffer) == Int16
    @test size(downconvert_and_correlator_numants2.carrier_replica_buffer) == (5000,)
    @test eltype(downconvert_and_correlator_numants2.carrier_replica_buffer) == ComplexF32
    @test size(downconvert_and_correlator_numants2.downconvert_signal_buffer) == (5000, 2)
    @test eltype(downconvert_and_correlator_numants2.downconvert_signal_buffer) ==
          ComplexF32

    num_samples = 5000
    galileo_e1b = GalileoE1B()
    multiple_system_sats_state = (
        SystemSatsState(gpsl1, [SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz)]),
        SystemSatsState(galileo_e1b, [SatState(galileo_e1b, 1, 5e6Hz, 10.0, 500.0Hz)]),
    )

    track_state = TrackState(
        multiple_system_sats_state;
        num_samples,
        maximum_expected_sampling_frequency = Val(5e6Hz),
    )

    @test size(track_state.downconvert_and_correlator.buffers[1][1].code_replica_buffer) ==
          (5004,)
    @test eltype(
        track_state.downconvert_and_correlator.buffers[1][1].code_replica_buffer,
    ) == Int16
    @test size(
        track_state.downconvert_and_correlator.buffers[1][1].carrier_replica_buffer,
    ) == (5000,)
    @test eltype(
        track_state.downconvert_and_correlator.buffers[1][1].carrier_replica_buffer,
    ) == ComplexF32
    @test size(
        track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer,
    ) == (5000,)
    @test eltype(
        track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer,
    ) == ComplexF32

    @test size(track_state.downconvert_and_correlator.buffers[2][1].code_replica_buffer) ==
          (5006,)
    @test eltype(
        track_state.downconvert_and_correlator.buffers[2][1].code_replica_buffer,
    ) == Float32
    @test size(
        track_state.downconvert_and_correlator.buffers[2][1].carrier_replica_buffer,
    ) == (5000,)
    @test eltype(
        track_state.downconvert_and_correlator.buffers[2][1].carrier_replica_buffer,
    ) == ComplexF32
    @test size(
        track_state.downconvert_and_correlator.buffers[2][1].downconvert_signal_buffer,
    ) == (5000,)
    @test eltype(
        track_state.downconvert_and_correlator.buffers[2][1].downconvert_signal_buffer,
    ) == ComplexF32
end

@testset "Downconvert and correlate with $type" for type in (:CPU, :GPU)
    type == :GPU && !CUDA.functional() && return

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    intermediate_frequency = 0.0Hz

    system_sats_state = SystemSatsState(
        gpsl1,
        [
            SatState(gpsl1, 1, sampling_frequency, code_phase, 1000.0Hz),
            SatState(gpsl1, 2, sampling_frequency, 11.0, 500.0Hz),
        ];
    )
    multiple_system_sats_state = (system_sats_state,)

    maximum_expected_sampling_frequency = Val(sampling_frequency)

    downconvert_and_correlator =
        type == :CPU ?
        CPUDownconvertAndCorrelator(
            maximum_expected_sampling_frequency,
            multiple_system_sats_state,
            num_samples_signal,
        ) : GPUDownconvertAndCorrelator(multiple_system_sats_state, num_samples_signal)

    array_transform = type == :CPU ? Array : cu

    track_state = TrackState(
        multiple_system_sats_state;
        num_samples = num_samples_signal,
        downconvert_and_correlator,
        maximum_expected_sampling_frequency,
    )

    preferred_num_code_blocks_to_integrate = 1

    signal = array_transform(
        gen_code(
            num_samples_signal,
            gpsl1,
            1,
            sampling_frequency,
            get_code_frequency(gpsl1) + 1000Hz * get_code_center_frequency_ratio(gpsl1),
            code_phase,
        ) .* cis.(2π * (0:num_samples_signal-1) * 1000.0Hz / sampling_frequency),
    )

    next_track_state = @inferred Tracking.downconvert_and_correlate(
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
        num_samples_signal,
    )

    # GPU uses floating point arithmetic and might differ a little with the fixed point arithmetic
    @test real.(
        next_track_state.multiple_system_sats_state[1].states[1].last_fully_integrated_correlator.accumulators
    ) ≈ [2921, 4949, 2917] ||
          real.(
        next_track_state.multiple_system_sats_state[1].states[1].last_fully_integrated_correlator.accumulators
    ) ≈ [2921, 4949, 2921]

    signal = array_transform(
        gen_code(
            num_samples_signal,
            gpsl1,
            2,
            sampling_frequency,
            get_code_frequency(gpsl1) + 500Hz * get_code_center_frequency_ratio(gpsl1),
            11.0,
        ) .* cis.(2π * (0:num_samples_signal-1) * 500.0Hz / sampling_frequency),
    )

    next_track_state = @inferred Tracking.downconvert_and_correlate(
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
        num_samples_signal,
    )

    # GPU uses floating point arithmetic and might differ a little with the fixed point arithmetic
    @test real.( #CPU
        next_track_state.multiple_system_sats_state[1].states[2].last_fully_integrated_correlator.accumulators
    ) ≈ [2919, 4947, 2915] ||
          real.( #GPU
        next_track_state.multiple_system_sats_state[1].states[2].last_fully_integrated_correlator.accumulators
    ) ≈ [2919, 4947, 2919]
end
