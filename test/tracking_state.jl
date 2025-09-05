module TrackingStateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, GalileoE1B
using Tracking:
    SatState,
    TrackState,
    SystemSatsState,
    get_system,
    get_prn,
    get_num_ants,
    get_integrated_samples,
    get_signal_start_sample,
    get_correlator,
    get_last_fully_integrated_correlator,
    get_post_corr_filter,
    get_cn0_estimator,
    get_bit_buffer,
    get_sat_state,
    get_sat_states,
    merge_sats,
    filter_out_sats,
    DefaultPostCorrFilter,
    MomentsCN0Estimator,
    BitBuffer,
    NumAnts

@testset "Tracking state" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1()
    sat_states = [SatState(gpsl1, 1, 10.5, 10.0Hz), SatState(gpsl1, 2, 11.5, 20.0Hz)]

    track_state = @inferred TrackState(gpsl1, sat_states)

    @test get_system(track_state) isa GPSL1
    @test get_prn(track_state, 1) == 1
    @test get_prn(track_state, 2) == 2
    @test get_num_ants(track_state, 1) == 1
    @test get_integrated_samples(track_state, 1) == 0
    @test get_signal_start_sample(track_state, 1) == 1
    @test get_correlator(track_state, 1).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(track_state, 1).accumulators == zeros(3)
    @test get_post_corr_filter(track_state, 1) isa DefaultPostCorrFilter
    @test get_cn0_estimator(track_state, 1) isa MomentsCN0Estimator
    @test get_bit_buffer(track_state, 1) isa BitBuffer

    @test track_state.doppler_estimator.states[1][1].init_carrier_doppler == 10.0Hz
    @test track_state.doppler_estimator.states[1][2].init_carrier_doppler == 20.0Hz

    system_sats_states = (
        SystemSatsState(
            gpsl1,
            [SatState(gpsl1, 1, 10.5, 10.0Hz), SatState(gpsl1, 2, 11.5, 20.0Hz)],
        ),
    )

    track_state = @inferred TrackState(system_sats_states)

    @test @inferred(get_system(track_state)) isa GPSL1
    @test @inferred(get_sat_state(track_state, 1)).prn == 1
    @test @inferred(get_sat_state(track_state, 2)).prn == 2

    @test track_state.doppler_estimator.states[1][1].init_carrier_doppler == 10.0Hz
    @test track_state.doppler_estimator.states[1][2].init_carrier_doppler == 20.0Hz

    sat_states_num_ants2 = [
        SatState(gpsl1, 1, 10.5, 10.0Hz; num_ants = NumAnts(2)),
        SatState(gpsl1, 2, 11.5, 20.0Hz; num_ants = NumAnts(2)),
    ]
end

@testset "Add and remove satellite state to and from track state" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1()
    sat_states = [SatState(gpsl1, 1, 10.5, 10.0Hz), SatState(gpsl1, 2, 11.5, 20.0Hz)]

    track_state = @inferred TrackState(gpsl1, sat_states)

    new_track_state = @inferred merge_sats(track_state, 1, SatState(gpsl1, 3, 5.5, 80.0Hz))

    @test length(@inferred(get_sat_states(new_track_state))) == 3
    @test @inferred(get_prn(new_track_state, 3)) == 3
    @test @inferred(get_prn(new_track_state, 1, 3)) == 3

    filtered_new_track_state = @inferred filter_out_sats(new_track_state, 1, 2)
    @test length(@inferred(get_sat_states(filtered_new_track_state))) == 2
    @test @inferred(get_prn(filtered_new_track_state, 1)) == 1
    @test @inferred(get_prn(filtered_new_track_state, 3)) == 3

    new_new_track_state = @inferred merge_sats(
        filtered_new_track_state,
        [SatState(gpsl1, 6, 12.5, 60.0Hz), SatState(gpsl1, 8, 67.5, 120.0Hz)],
    )

    @test length(@inferred(get_sat_states(new_new_track_state))) == 4
    @test @inferred(get_prn(new_new_track_state, 1)) == 1
    @test @inferred(get_prn(new_new_track_state, 3)) == 3
    @test @inferred(get_prn(new_new_track_state, 6)) == 6
    @test @inferred(get_prn(new_new_track_state, 8)) == 8

    filtered_new_new_track_state = @inferred filter_out_sats(new_new_track_state, [1, 6])

    @test length(@inferred(get_sat_states(filtered_new_new_track_state))) == 2
    @test @inferred(get_prn(filtered_new_new_track_state, 3)) == 3
    @test @inferred(get_prn(filtered_new_new_track_state, 8)) == 8
end

@testset "Add and remove satellite state to track state with multiple systems" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1()
    galileo_e1b = GalileoE1B()

    system_sats_states = (
        gps = SystemSatsState(
            gpsl1,
            [SatState(gpsl1, 1, 10.5, 10.0Hz), SatState(gpsl1, 2, 11.5, 20.0Hz)],
        ),
        gal = SystemSatsState(
            galileo_e1b,
            [
                SatState(galileo_e1b, 2, 10.5, 10.0Hz),
                SatState(galileo_e1b, 3, 11.5, 20.0Hz),
            ],
        ),
    )

    track_state = @inferred TrackState(system_sats_states)

    new_track_state = @inferred merge_sats(track_state, 1, SatState(gpsl1, 3, 5.5, 80.0Hz))

    @test length(@inferred(get_sat_states(new_track_state, Val(:gps)))) == 3
    @test length(@inferred(get_sat_states(new_track_state, Val(:gal)))) == 2
    @test length(@inferred(get_sat_states(new_track_state, Val(1)))) == 3
    @test length(@inferred(get_sat_states(new_track_state, Val(2)))) == 2
    @test @inferred(get_prn(new_track_state, Val(:gps), 3)) == 3
    @test @inferred(get_prn(new_track_state, Val(1), 3)) == 3

    filtered_new_track_state = @inferred filter_out_sats(new_track_state, 1, 2)
    @test length(get_sat_states(filtered_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_track_state, :gal)) == 2
    @test length(get_sat_states(filtered_new_track_state, 1)) == 2
    @test length(get_sat_states(filtered_new_track_state, 2)) == 2
    @test get_prn(filtered_new_track_state, :gps, 1) == 1
    @test get_prn(filtered_new_track_state, :gps, 3) == 3

    new_new_track_state = @inferred merge_sats(
        filtered_new_track_state,
        :gps,
        SatState(gpsl1, 2, 5.5, 80.0Hz),
    )
    @test length(get_sat_states(new_new_track_state, :gps)) == 3
    @test length(get_sat_states(new_new_track_state, :gal)) == 2
    @test get_prn(new_new_track_state, :gps, 1) == 1
    @test get_prn(new_new_track_state, :gps, 3) == 3
    @test get_prn(new_new_track_state, :gps, 2) == 2

    filtered_new_new_track_state = @inferred filter_out_sats(new_new_track_state, :gps, 2)
    @test length(get_sat_states(filtered_new_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_new_track_state, :gal)) == 2
    @test get_prn(filtered_new_new_track_state, :gps, 1) == 1
    @test get_prn(filtered_new_new_track_state, :gps, 3) == 3

    new_new_new_track_state = @inferred merge_sats(
        filtered_new_new_track_state,
        :gps,
        [SatState(gpsl1, 6, 12.5, 60.0Hz), SatState(gpsl1, 8, 67.5, 120.0Hz)],
    )

    @test length(get_sat_states(new_new_new_track_state, :gps)) == 4
    @test length(get_sat_states(new_new_new_track_state, :gal)) == 2
    @test get_prn(new_new_new_track_state, :gps, 1) == 1
    @test get_prn(new_new_new_track_state, :gps, 3) == 3
    @test get_prn(new_new_new_track_state, :gps, 6) == 6
    @test get_prn(new_new_new_track_state, :gps, 8) == 8

    filtered_new_new_new_track_state =
        @inferred filter_out_sats(new_new_new_track_state, :gps, [1, 6])

    @test length(get_sat_states(filtered_new_new_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_new_new_track_state, :gal)) == 2
    @test get_prn(filtered_new_new_new_track_state, :gps, 3) == 3
    @test get_prn(filtered_new_new_new_track_state, :gps, 8) == 8
end

end
