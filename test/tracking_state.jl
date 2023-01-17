@testset "Tracking state" begin
    sampling_frequency = 5e6Hz
    num_samples = 5000
    gpsl1 = GPSL1()
    sat_states = [
        SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
        SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz)
    ]

    track_state = @inferred TrackState(
        gpsl1,
        sat_states;
        num_samples,
    )

    @test track_state.system_sats_states[1].system isa GPSL1
    @test track_state.system_sats_states[1].states[1].prn == 1
    @test track_state.system_sats_states[1].states[2].prn == 2

    @test track_state.doppler_estimator.plls_and_dlls[1][1].init_carrier_doppler == 10.0Hz
    @test track_state.doppler_estimator.plls_and_dlls[1][2].init_carrier_doppler == 20.0Hz

    @test track_state.downconvert_and_correlator isa CPUDownconvertAndCorrelator

    system_sats_states = (SystemSatsState(
        gpsl1, 
        [
            SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
            SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz)
        ]
    ),)

    track_state = @inferred TrackState(
        system_sats_states;
        num_samples,
    )

    @test track_state.system_sats_states[1].system isa GPSL1
    @test track_state.system_sats_states[1].states[1].prn == 1
    @test track_state.system_sats_states[1].states[2].prn == 2

    @test track_state.doppler_estimator.plls_and_dlls[1][1].init_carrier_doppler == 10.0Hz
    @test track_state.doppler_estimator.plls_and_dlls[1][2].init_carrier_doppler == 20.0Hz

    @test track_state.downconvert_and_correlator isa CPUDownconvertAndCorrelator

    sat_states_num_ants2 = [
        SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz, num_ants = NumAnts(2)),
        SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz, num_ants = NumAnts(2))
    ]

    track_state = @inferred TrackState(
        gpsl1,
        sat_states_num_ants2;
        num_samples
    )

    @test track_state.downconvert_and_correlator isa CPUDownconvertAndCorrelator
    @test track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer isa AbstractMatrix
    @test size(track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer, 2) == 2
    @test size(track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer, 1) == 5000
    @test track_state.downconvert_and_correlator.buffers[1][2].downconvert_signal_buffer isa AbstractMatrix
    @test size(track_state.downconvert_and_correlator.buffers[1][2].downconvert_signal_buffer, 2) == 2
    @test size(track_state.downconvert_and_correlator.buffers[1][2].downconvert_signal_buffer, 1) == 5000
end

@testset "Add and remove satellite state to track state" begin
    sampling_frequency = 5e6Hz
    num_samples = 5000
    gpsl1 = GPSL1()
    sat_states = [
        SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
        SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz)
    ]

    track_state = @inferred TrackState(
        gpsl1,
        sat_states;
        num_samples,
    )

    @inferred add_sats!(
        track_state,
        gpsl1,
        SatState(gpsl1, 3, sampling_frequency, 5.5, 80.0Hz),
    )

    @test length(track_state.system_sats_states[1].states) == 3
    @test track_state.system_sats_states[1].states[3].prn == 3

    @inferred remove_sats!(
        track_state,
        gpsl1,
        2,
    )
    @test length(track_state.system_sats_states[1].states) == 2
    @test track_state.system_sats_states[1].states[1].prn == 1
    @test track_state.system_sats_states[1].states[2].prn == 3

    @inferred add_sats!(
        track_state,
        gpsl1,
        [
            SatState(gpsl1, 6, sampling_frequency, 12.5, 60.0Hz),
            SatState(gpsl1, 8, sampling_frequency, 67.5, 120.0Hz),
        ]
    )

    @test length(track_state.system_sats_states[1].states) == 4
    @test track_state.system_sats_states[1].states[1].prn == 1
    @test track_state.system_sats_states[1].states[2].prn == 3
    @test track_state.system_sats_states[1].states[3].prn == 6
    @test track_state.system_sats_states[1].states[4].prn == 8

    @inferred remove_sats!(
        track_state,
        gpsl1,
        [1, 6],
    )

    @test length(track_state.system_sats_states[1].states) == 2
    @test track_state.system_sats_states[1].states[1].prn == 3
    @test track_state.system_sats_states[1].states[2].prn == 8
    
end