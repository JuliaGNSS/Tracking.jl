@testset "Tracking state" begin
    sampling_frequency = 5e6Hz
    num_samples = 5000
    gpsl1 = GPSL1()
    sat_states = [
        SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
        SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz),
    ]

    track_state = @inferred TrackState(gpsl1, sat_states; num_samples)

    @test get_system(track_state) isa GPSL1
    @test get_sat_state(track_state, 1).prn == 1
    @test get_sat_state(track_state, 2).prn == 2

    @test track_state.doppler_estimator.plls_and_dlls[1][1].init_carrier_doppler == 10.0Hz
    @test track_state.doppler_estimator.plls_and_dlls[1][2].init_carrier_doppler == 20.0Hz

    @test track_state.downconvert_and_correlator isa CPUDownconvertAndCorrelator

    system_sats_states = (
        SystemSatsState(
            gpsl1,
            [
                SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
                SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz),
            ],
        ),
    )

    track_state = @inferred TrackState(system_sats_states; num_samples)

    @test @inferred(get_system(track_state)) isa GPSL1
    @test @inferred(get_sat_state(track_state, 1)).prn == 1
    @test @inferred(get_sat_state(track_state, 2)).prn == 2

    @test track_state.doppler_estimator.plls_and_dlls[1][1].init_carrier_doppler == 10.0Hz
    @test track_state.doppler_estimator.plls_and_dlls[1][2].init_carrier_doppler == 20.0Hz

    @test track_state.downconvert_and_correlator isa CPUDownconvertAndCorrelator

    sat_states_num_ants2 = [
        SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz; num_ants = NumAnts(2)),
        SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz; num_ants = NumAnts(2)),
    ]

    track_state = @inferred TrackState(gpsl1, sat_states_num_ants2; num_samples)

    @test track_state.downconvert_and_correlator isa CPUDownconvertAndCorrelator
    @test track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer isa
          AbstractMatrix
    @test size(
        track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer,
        2,
    ) == 2
    @test size(
        track_state.downconvert_and_correlator.buffers[1][1].downconvert_signal_buffer,
        1,
    ) == 5000
    @test track_state.downconvert_and_correlator.buffers[1][2].downconvert_signal_buffer isa
          AbstractMatrix
    @test size(
        track_state.downconvert_and_correlator.buffers[1][2].downconvert_signal_buffer,
        2,
    ) == 2
    @test size(
        track_state.downconvert_and_correlator.buffers[1][2].downconvert_signal_buffer,
        1,
    ) == 5000
end

@testset "Add and remove satellite state to and from track state" begin
    sampling_frequency = 5e6Hz
    num_samples = 5000
    gpsl1 = GPSL1()
    sat_states = [
        SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
        SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz),
    ]

    track_state = @inferred TrackState(gpsl1, sat_states; num_samples)

    @inferred add_sats!(
        track_state,
        1,
        gpsl1,
        SatState(gpsl1, 3, sampling_frequency, 5.5, 80.0Hz),
    )

    @test length(@inferred(get_sat_states(track_state))) == 3
    @test @inferred(get_sat_state(track_state, 3)).prn == 3
    @test @inferred(get_sat_state(track_state, 1, 3)).prn == 3

    @inferred remove_sats!(track_state, 1, 2)
    @test length(@inferred(get_sat_states(track_state))) == 2
    @test @inferred(get_sat_state(track_state, 1)).prn == 1
    @test @inferred(get_sat_state(track_state, 2)).prn == 3

    @inferred add_sats!(
        track_state,
        gpsl1,
        [
            SatState(gpsl1, 6, sampling_frequency, 12.5, 60.0Hz),
            SatState(gpsl1, 8, sampling_frequency, 67.5, 120.0Hz),
        ],
    )

    @test length(@inferred(get_sat_states(track_state))) == 4
    @test @inferred(get_sat_state(track_state, 1)).prn == 1
    @test @inferred(get_sat_state(track_state, 2)).prn == 3
    @test @inferred(get_sat_state(track_state, 3)).prn == 6
    @test @inferred(get_sat_state(track_state, 4)).prn == 8

    @inferred remove_sats!(track_state, [1, 6])

    @test length(@inferred(get_sat_states(track_state))) == 2
    @test @inferred(get_sat_state(track_state, 1)).prn == 3
    @test @inferred(get_sat_state(track_state, 2)).prn == 8
end

@testset "Add and remove satellite state to track state with multiple systems" begin
    sampling_frequency = 5e6Hz
    num_samples = 5000
    gpsl1 = GPSL1()
    galileo_e1b = GalileoE1B()

    system_sats_states = (
        gps = SystemSatsState(
            gpsl1,
            [
                SatState(gpsl1, 1, sampling_frequency, 10.5, 10.0Hz),
                SatState(gpsl1, 2, sampling_frequency, 11.5, 20.0Hz),
            ],
        ),
        gal = SystemSatsState(
            galileo_e1b,
            [
                SatState(galileo_e1b, 2, sampling_frequency, 10.5, 10.0Hz),
                SatState(galileo_e1b, 3, sampling_frequency, 11.5, 20.0Hz),
            ],
        ),
    )

    track_state = @inferred TrackState(system_sats_states; num_samples)

    @inferred add_sats!(
        track_state,
        1,
        gpsl1,
        SatState(gpsl1, 3, sampling_frequency, 5.5, 80.0Hz),
    )

    @test length(@inferred(get_sat_states(track_state, Val(:gps)))) == 3
    @test length(@inferred(get_sat_states(track_state, Val(:gal)))) == 2
    @test length(@inferred(get_sat_states(track_state, Val(1)))) == 3
    @test length(@inferred(get_sat_states(track_state, Val(2)))) == 2
    @test @inferred(get_sat_state(track_state, Val(:gps), 3)).prn == 3
    @test @inferred(get_sat_state(track_state, Val(1), 3)).prn == 3

    @inferred remove_sats!(track_state, 1, 2)
    @test length(get_sat_states(track_state, :gps)) == 2
    @test length(get_sat_states(track_state, :gal)) == 2
    @test length(get_sat_states(track_state, 1)) == 2
    @test length(get_sat_states(track_state, 2)) == 2
    @test get_sat_state(track_state, :gps, 1).prn == 1
    @test get_sat_state(track_state, :gps, 2).prn == 3

    @inferred add_sats!(
        track_state,
        :gps,
        gpsl1,
        SatState(gpsl1, 2, sampling_frequency, 5.5, 80.0Hz),
    )
    @test length(get_sat_states(track_state, :gps)) == 3
    @test length(get_sat_states(track_state, :gal)) == 2
    @test get_sat_state(track_state, :gps, 1).prn == 1
    @test get_sat_state(track_state, :gps, 2).prn == 3
    @test get_sat_state(track_state, :gps, 3).prn == 2

    @inferred remove_sats!(track_state, :gps, 2)
    @test length(get_sat_states(track_state, :gps)) == 2
    @test length(get_sat_states(track_state, :gal)) == 2
    @test get_sat_state(track_state, :gps, 1).prn == 1
    @test get_sat_state(track_state, :gps, 2).prn == 3

    @inferred add_sats!(
        track_state,
        :gps,
        gpsl1,
        [
            SatState(gpsl1, 6, sampling_frequency, 12.5, 60.0Hz),
            SatState(gpsl1, 8, sampling_frequency, 67.5, 120.0Hz),
        ],
    )

    @test length(get_sat_states(track_state, :gps)) == 4
    @test length(get_sat_states(track_state, :gal)) == 2
    @test get_sat_state(track_state, :gps, 1).prn == 1
    @test get_sat_state(track_state, :gps, 2).prn == 3
    @test get_sat_state(track_state, :gps, 3).prn == 6
    @test get_sat_state(track_state, :gps, 4).prn == 8

    @inferred remove_sats!(track_state, :gps, [1, 6])

    @test length(get_sat_states(track_state, :gps)) == 2
    @test length(get_sat_states(track_state, :gal)) == 2
    @test get_sat_state(track_state, :gps, 1).prn == 3
    @test get_sat_state(track_state, :gps, 2).prn == 8
end
