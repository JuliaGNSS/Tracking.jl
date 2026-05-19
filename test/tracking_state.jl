module TrackingStateTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz
using GNSSSignals: GPSL1CA, GalileoE1B
using Dictionaries: Dictionary, dictionary
import Tracking
using Tracking:
    TrackedSat,
    TrackState,
    add_satellite!,
    get_signal,
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
    remove_satellite!,
    remove_satellite,
    DefaultPostCorrFilter,
    MomentsCN0Estimator,
    BitBuffer,
    NumAnts,
    ConventionalAssistedPLLAndDLL,
    ConventionalPLLAndDLL

@testset "Tracking state" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    sat_states = [TrackedSat(gpsl1, 1, 10.5, 10.0Hz), TrackedSat(gpsl1, 2, 11.5, 20.0Hz)]

    track_state = @inferred TrackState(gpsl1, sat_states)

    @test get_signal(track_state) isa GPSL1CA
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

    @test get_sat_states(track_state)[1].doppler_estimator_state.init_carrier_doppler == 10.0Hz
    @test get_sat_states(track_state)[2].doppler_estimator_state.init_carrier_doppler == 20.0Hz

    track_state2 = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz), TrackedSat(gpsl1, 2, 11.5, 20.0Hz)],
    )

    @test @inferred(get_signal(track_state2)) isa GPSL1CA
    @test @inferred(get_sat_state(track_state2, 1)).prn == 1
    @test @inferred(get_sat_state(track_state2, 2)).prn == 2

    @test get_sat_states(track_state2)[1].doppler_estimator_state.init_carrier_doppler == 10.0Hz
    @test get_sat_states(track_state2)[2].doppler_estimator_state.init_carrier_doppler == 20.0Hz

    sat_states_num_ants2 = [
        TrackedSat(gpsl1, 1, 10.5, 10.0Hz; num_ants = NumAnts(2)),
        TrackedSat(gpsl1, 2, 11.5, 20.0Hz; num_ants = NumAnts(2)),
    ]
end

@testset "Legacy TrackState(satellites::SatelliteDicts) infers default estimator" begin
    # `_default_estimator_for_satellite_dicts` picks the most-constraining
    # default across all per-system driver signals when no estimator kwarg
    # is supplied.
    gpsl1 = GPSL1CA()
    galileo = GalileoE1B()
    estimator = ConventionalAssistedPLLAndDLL()
    sats = (
        gps = dictionary([
            1 => TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator),
        ]),
        gal = dictionary([
            2 => TrackedSat(galileo, 2, 11.5, 20.0Hz; doppler_estimator = estimator),
        ]),
    )
    ts = TrackState(sats)
    @test length(get_sat_states(ts, :gps)) == 1
    @test length(get_sat_states(ts, :gal)) == 1
end

@testset "Legacy TrackState(dict) default estimator inference" begin
    # `_default_estimator_for_sats_dict` (non-empty path) sizes a default
    # estimator from the dict's first sat's driver signal.
    gpsl1 = GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL()
    sat = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator)
    ts = TrackState(dictionary([1 => sat]))
    @test length(get_sat_states(ts)) == 1
end

@testset "Legacy TrackState constructor rejects sats with a different estimator" begin
    # `_assert_doppler_estimator_types_match` errors when sat-state and
    # the configured estimator would produce different concrete types.
    gpsl1 = GPSL1CA()
    sat_default = TrackedSat(gpsl1, 1, 10.5, 10.0Hz)  # default estimator
    different = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test_throws ArgumentError TrackState(gpsl1, [sat_default];
                                          doppler_estimator = different)
end

@testset "Internal merge_sats(::SatelliteDicts, group_idx, ::Dictionary)" begin
    # The sat-state-level helper that the legacy TrackState `merge_sats`
    # delegates to. Exercise it directly so the recursion is recorded.
    gpsl1 = GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL()
    a = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator)
    b = TrackedSat(gpsl1, 2, 11.5, 20.0Hz; doppler_estimator = estimator)
    sats = (gps = dictionary([1 => a]),)
    new_sats = Tracking.merge_sats(sats, :gps, dictionary([2 => b]))
    @test length(new_sats.gps) == 2
    @test new_sats.gps[2].prn == 2

    # `_copy_slot_vectors` walks a `SatelliteDicts` and copies each dict's
    # values vector — the result should share keys but detach the slot
    # vector.
    copies = Tracking._copy_slot_vectors(sats)
    @test length(copies.gps) == 1
    @test copies.gps.values !== sats.gps.values
end

@testset "remove_satellite! kwarg form mutates in place" begin
    ts = TrackState(; signal = GPSL1CA())
    add_satellite!(ts; prn = 5, carrier_doppler = 100.0Hz)
    add_satellite!(ts; prn = 6, carrier_doppler = 200.0Hz)
    @test length(get_sat_states(ts, :default)) == 2

    ret = remove_satellite!(ts; prn = 5)
    @test ret === ts
    @test length(get_sat_states(ts, :default)) == 1
    @test get_prn(ts, :default, 6) == 6
end

@testset "Add and remove satellite state to and from track state" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    sat_states = [TrackedSat(gpsl1, 1, 10.5, 10.0Hz), TrackedSat(gpsl1, 2, 11.5, 20.0Hz)]

    track_state = @inferred TrackState(gpsl1, sat_states)

    new_track_state = @inferred merge_sats(track_state, 1, TrackedSat(gpsl1, 3, 5.5, 80.0Hz))

    @test length(@inferred(get_sat_states(new_track_state))) == 3
    @test @inferred(get_prn(new_track_state, 3)) == 3
    @test @inferred(get_prn(new_track_state, 1, 3)) == 3

    filtered_new_track_state = @inferred remove_satellite(new_track_state; prn = 2)
    @test length(@inferred(get_sat_states(filtered_new_track_state))) == 2
    @test @inferred(get_prn(filtered_new_track_state, 1)) == 1
    @test @inferred(get_prn(filtered_new_track_state, 3)) == 3

    new_new_track_state = @inferred merge_sats(
        filtered_new_track_state,
        [TrackedSat(gpsl1, 6, 12.5, 60.0Hz), TrackedSat(gpsl1, 8, 67.5, 120.0Hz)],
    )

    @test length(@inferred(get_sat_states(new_new_track_state))) == 4
    @test @inferred(get_prn(new_new_track_state, 1)) == 1
    @test @inferred(get_prn(new_new_track_state, 3)) == 3
    @test @inferred(get_prn(new_new_track_state, 6)) == 6
    @test @inferred(get_prn(new_new_track_state, 8)) == 8

    filtered_new_new_track_state =
        @inferred remove_satellite(
            @inferred(remove_satellite(new_new_track_state; prn = 1));
            prn = 6,
        )

    @test length(@inferred(get_sat_states(filtered_new_new_track_state))) == 2
    @test @inferred(get_prn(filtered_new_new_track_state, 3)) == 3
    @test @inferred(get_prn(filtered_new_new_track_state, 8)) == 8
end

@testset "Add and remove satellite state to track state with multiple systems" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    galileo_e1b = GalileoE1B()

    estimator = ConventionalAssistedPLLAndDLL()
    gps_sats = dictionary([
        1 => TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator),
        2 => TrackedSat(gpsl1, 2, 11.5, 20.0Hz; doppler_estimator = estimator),
    ])
    gal_sats = dictionary([
        2 => TrackedSat(galileo_e1b, 2, 10.5, 10.0Hz; doppler_estimator = estimator),
        3 => TrackedSat(galileo_e1b, 3, 11.5, 20.0Hz; doppler_estimator = estimator),
    ])
    satellites = (gps = gps_sats, gal = gal_sats)

    track_state = @inferred TrackState(satellites; doppler_estimator = estimator)

    new_track_state = @inferred merge_sats(track_state, 1, TrackedSat(gpsl1, 3, 5.5, 80.0Hz))

    @test length(@inferred(get_sat_states(new_track_state, Val(:gps)))) == 3
    @test length(@inferred(get_sat_states(new_track_state, Val(:gal)))) == 2
    @test length(@inferred(get_sat_states(new_track_state, Val(1)))) == 3
    @test length(@inferred(get_sat_states(new_track_state, Val(2)))) == 2
    @test @inferred(get_prn(new_track_state, Val(:gps), 3)) == 3
    @test @inferred(get_prn(new_track_state, Val(1), 3)) == 3

    filtered_new_track_state =
        @inferred remove_satellite(new_track_state; prn = 2, group = :gps)
    @test length(get_sat_states(filtered_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_track_state, :gal)) == 2
    @test length(get_sat_states(filtered_new_track_state, 1)) == 2
    @test length(get_sat_states(filtered_new_track_state, 2)) == 2
    @test get_prn(filtered_new_track_state, :gps, 1) == 1
    @test get_prn(filtered_new_track_state, :gps, 3) == 3

    new_new_track_state = @inferred merge_sats(
        filtered_new_track_state,
        :gps,
        TrackedSat(gpsl1, 2, 5.5, 80.0Hz),
    )
    @test length(get_sat_states(new_new_track_state, :gps)) == 3
    @test length(get_sat_states(new_new_track_state, :gal)) == 2
    @test get_prn(new_new_track_state, :gps, 1) == 1
    @test get_prn(new_new_track_state, :gps, 3) == 3
    @test get_prn(new_new_track_state, :gps, 2) == 2

    filtered_new_new_track_state =
        @inferred remove_satellite(new_new_track_state; prn = 2, group = :gps)
    @test length(get_sat_states(filtered_new_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_new_track_state, :gal)) == 2
    @test get_prn(filtered_new_new_track_state, :gps, 1) == 1
    @test get_prn(filtered_new_new_track_state, :gps, 3) == 3

    new_new_new_track_state = @inferred merge_sats(
        filtered_new_new_track_state,
        :gps,
        [TrackedSat(gpsl1, 6, 12.5, 60.0Hz), TrackedSat(gpsl1, 8, 67.5, 120.0Hz)],
    )

    @test length(get_sat_states(new_new_new_track_state, :gps)) == 4
    @test length(get_sat_states(new_new_new_track_state, :gal)) == 2
    @test get_prn(new_new_new_track_state, :gps, 1) == 1
    @test get_prn(new_new_new_track_state, :gps, 3) == 3
    @test get_prn(new_new_new_track_state, :gps, 6) == 6
    @test get_prn(new_new_new_track_state, :gps, 8) == 8

    filtered_new_new_new_track_state =
        @inferred remove_satellite(
            @inferred(remove_satellite(new_new_new_track_state; prn = 1, group = :gps));
            prn = 6, group = :gps,
        )

    @test length(get_sat_states(filtered_new_new_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_new_new_track_state, :gal)) == 2
    @test get_prn(filtered_new_new_new_track_state, :gps, 3) == 3
    @test get_prn(filtered_new_new_new_track_state, :gps, 8) == 8
end

@testset "TrackedSat helpers" begin
    using Tracking:
        TrackedSat,
        get_doppler_estimator_state,
        reset_start_sample_and_bit_buffer!,
        SatConventionalPLLAndDLL

    gpsl1 = GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL()

    # Build TrackedSats directly with the doppler estimator attached.
    sat = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator)
    sat2 = TrackedSat(gpsl1, 2, 11.5, 20.0Hz; doppler_estimator = estimator)

    # TrackState accepts vectors, dictionaries, or single sats.
    ts_from_vec = TrackState(gpsl1, [sat, sat2]; doppler_estimator = estimator)
    @test length(get_sat_states(ts_from_vec)) == 2
    ts_from_one = TrackState(gpsl1, sat; doppler_estimator = estimator)
    @test length(get_sat_states(ts_from_one)) == 1
    ts_from_dict = TrackState(
        dictionary([1 => sat, 2 => sat2]);
        doppler_estimator = estimator,
    )
    @test length(get_sat_states(ts_from_dict)) == 2

    # Direct accessors on TrackedSat
    @test sat.prn == 1
    @test get_doppler_estimator_state(sat) isa SatConventionalPLLAndDLL

    # In-place reset on a single TrackedSat (the `!` overload).
    track_state = TrackState(gpsl1, [sat, sat2]; doppler_estimator = estimator)
    track_state2 = reset_start_sample_and_bit_buffer!(track_state)
    @test track_state2 === track_state
    # Each TrackedSat in the dict has signal_start_sample reset to 1.
    for tsat in get_sat_states(track_state).values
        @test tsat.signal_start_sample == 1
    end
end

@testset "merge_sats invokes update_estimator_on_handoff" begin
    import Tracking
    using Tracking: AbstractDopplerEstimator

    # Immutable estimator with growing shared state held in a resizable
    # Vector. update_estimator_on_handoff mutates the vectors in place
    # and returns the same `est` object — concrete type is preserved, no
    # new heap-allocated wrapper per call.
    struct CountingEstimator <: AbstractDopplerEstimator
        sats_added::Vector{Int}
        last_batch_prns::Vector{Int}
    end
    CountingEstimator() = CountingEstimator(Int[0], Int[])

    struct SatCountingEstimator end

    Tracking.init_estimator_state(::CountingEstimator, ::TrackedSat) = SatCountingEstimator()

    function Tracking.update_estimator_on_handoff(est::CountingEstimator, new_sats)
        est.sats_added[1] += length(new_sats)
        empty!(est.last_batch_prns)
        for s in new_sats
            push!(est.last_batch_prns, s.prn)
        end
        return est
    end

    gpsl1 = GPSL1CA()
    estimator = CountingEstimator()
    track_state = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator)];
        doppler_estimator = estimator,
    )
    @test track_state.doppler_estimator.sats_added[1] == 0

    track_state2 = merge_sats(
        track_state,
        1,
        TrackedSat(gpsl1, 2, 11.5, 20.0Hz; doppler_estimator = estimator),
    )
    @test track_state2.doppler_estimator === estimator  # same object, mutated in place
    @test track_state2.doppler_estimator.sats_added[1] == 1
    @test track_state2.doppler_estimator.last_batch_prns == [2]

    track_state3 = merge_sats(
        track_state2,
        [
            TrackedSat(gpsl1, 3, 5.5, 80.0Hz; doppler_estimator = estimator),
            TrackedSat(gpsl1, 4, 6.5, 90.0Hz; doppler_estimator = estimator),
        ],
    )
    @test track_state3.doppler_estimator.sats_added[1] == 3
    @test sort(track_state3.doppler_estimator.last_batch_prns) == [3, 4]

    # Default fallback: estimators that don't override the hook keep the
    # estimator object unchanged through merge_sats.
    default_estimator = ConventionalAssistedPLLAndDLL()
    default_track_state = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = default_estimator)];
        doppler_estimator = default_estimator,
    )
    merged_default = merge_sats(
        default_track_state,
        1,
        TrackedSat(gpsl1, 2, 11.5, 20.0Hz; doppler_estimator = default_estimator),
    )
    @test merged_default.doppler_estimator === default_estimator
end

end
