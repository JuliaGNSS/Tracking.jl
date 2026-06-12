module TrackingStateTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz
using GNSSSignals:
    GNSSSignals, AbstractGNSSSignal, GPSL1CA, GPSL1C_D, GPSL1C_P, GPSL5I, GalileoE1B
using Dictionaries: Dictionary, dictionary
import Tracking
using Tracking:
    TrackedSat,
    TrackState,
    add_satellite!,
    get_signal,
    get_prn,
    get_num_ants,
    get_code_phase,
    get_code_doppler,
    get_carrier_phase,
    get_carrier_doppler,
    get_integrated_samples,
    get_signal_start_sample,
    get_correlator,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_filtered_prompt,
    get_post_corr_filter,
    get_cn0_estimator,
    get_bit_buffer,
    get_bits,
    get_num_bits,
    has_bit_or_secondary_code_been_found,
    estimate_cn0,
    get_sat_state,
    get_sat_states,
    merge_sats,
    remove_satellite!,
    remove_satellite,
    DefaultPostCorrFilter,
    MomentsCN0Estimator,
    BitBuffer,
    NumAnts,
    get_preferred_num_code_blocks_to_integrate,
    set_preferred_num_code_blocks_to_integrate!,
    reset_loop_filters!,
    get_doppler_estimator_state,
    ConventionalAssistedPLLAndDLL,
    ConventionalPLLAndDLL,
    SignalGroup

# No real GNSS signal pair mixes chip rates on a single band, so fake one to
# exercise the SignalGroup chip-rate invariant (issue #129): L1 band like
# GPS L1 C/A, but at double the chip rate. Only `get_band` and
# `get_code_frequency` are needed — validation throws before the group
# touches any other part of the signal API.
struct FakeDoubleRateL1Signal <: AbstractGNSSSignal{Matrix{Int16}} end
GNSSSignals.get_band(::FakeDoubleRateL1Signal) = GNSSSignals.L1()
GNSSSignals.get_code_frequency(::FakeDoubleRateL1Signal) = 2_046_000Hz

@testset "SignalGroup rejects mixed bands and mixed chip rates (issue #129)" begin
    # (a) Band homogeneity: every signal in a group is downconverted against
    # the single band measurement the group's `band` routes to, so an L5
    # signal in an L1 group would silently correlate against L1 samples.
    @test_throws ArgumentError SignalGroup((GPSL1CA(), GPSL5I()))
    @test_throws ArgumentError TrackState(; signals = (mix = (GPSL1CA(), GPSL5I()),))
    # An explicit `band` override that contradicts the signals is just as wrong.
    @test_throws ArgumentError SignalGroup((GPSL1CA(),); band = GNSSSignals.L5())

    # (b) Shared chip rate: the shared `code_phase` advances at `signals[1]`'s
    # code frequency, so a signal with a different chip rate silently mistracks.
    @test_throws ArgumentError SignalGroup((GPSL1CA(), FakeDoubleRateL1Signal()))
    @test_throws ArgumentError TrackState(;
        signals = (l1 = (GPSL1CA(), FakeDoubleRateL1Signal()),),
    )

    # Valid same-band, same-chip-rate combinations still construct.
    @test SignalGroup((GPSL1C_P(), GPSL1C_D(), GPSL1CA())) isa SignalGroup
    ts = TrackState(; signals = (l1 = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()), l5 = (GPSL5I(),)))
    @test ts isa TrackState
end

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

@testset "Positional TrackState(satellites::SatelliteDicts) infers default estimator" begin
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

@testset "Positional TrackState(dict) default estimator inference" begin
    # `_default_estimator_for_sats_dict` (non-empty path) sizes a default
    # estimator from the dict's first sat's driver signal.
    gpsl1 = GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL()
    sat = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator)
    ts = TrackState(dictionary([1 => sat]))
    @test length(get_sat_states(ts)) == 1
end

@testset "Positional TrackState constructor rejects sats with a different estimator" begin
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
    # The sat-state-level helper that the positional TrackState `merge_sats`
    # delegates to. Exercise it directly so the recursion is recorded.
    gpsl1 = GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL()
    a = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; doppler_estimator = estimator)
    b = TrackedSat(gpsl1, 2, 11.5, 20.0Hz; doppler_estimator = estimator)
    sats = (gps = dictionary([1 => a]),)
    new_sats = Tracking.merge_sats(sats, :gps, dictionary([2 => b]))
    @test length(new_sats.gps) == 2
    @test new_sats.gps[2].prn == 2

    # `_copy_slot_vectors` is the cheap per-iteration copy used inside
    # `track`'s loop (`downconvert_and_correlate` / `estimate`): it detaches
    # the slot *values* but deliberately *shares* the key set (`Indices`), so
    # the hash table is not copied on every loop iteration. The key set is
    # detached once at the `track` boundary instead (see below, #123).
    copies = Tracking._copy_slot_vectors(sats)
    @test length(copies.gps) == 1
    @test copies.gps.values !== sats.gps.values
    @test keys(copies.gps) === keys(sats.gps)

    # `_detach_slot_vector` / `_detach_groups_slot_vectors` are the boundary
    # copy: keys *and* values detached.
    detached = Tracking._detach_slot_vector(sats.gps)
    @test detached.values !== sats.gps.values
    @test keys(detached) !== keys(sats.gps)
end

@testset "Immutable boundary copies do not share Dictionary Indices (#123)" begin
    # `reset_start_sample_and_bit_buffer` is the boundary copy `track` makes
    # of the caller's live state. With a shared `Indices`, `add_satellite!`
    # on the copy would grow the original's key set without resizing its
    # values vector — leaving the original claiming keys it has no values for
    # (UndefRefError or silent garbage on access).
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1, carrier_doppler = 100.0Hz)

    ts2 = Tracking.reset_start_sample_and_bit_buffer(ts)
    ts2 = add_satellite!(ts2; prn = 2, carrier_doppler = 200.0Hz)

    # The original must be untouched: still exactly one satellite.
    @test length(get_sat_states(ts, :default)) == 1
    @test collect(keys(get_sat_states(ts, :default))) == [1]
    @test get_prn(ts, :default, 1) == 1
    # The copy got the new satellite.
    @test length(get_sat_states(ts2, :default)) == 2

    # Removal on the copy must not shrink the original's key set either.
    ts2 = remove_satellite!(ts2; prn = 1)
    @test length(get_sat_states(ts, :default)) == 1
    @test get_prn(ts, :default, 1) == 1
end

@testset "track output is structurally detached from input (#123)" begin
    # The public contract: `add_satellite!` on `track`'s output must not
    # corrupt the input state, even though the hot loop shares key sets among
    # its throwaway intermediates. The detach at the `track` boundary
    # (`reset_start_sample_and_bit_buffer`) makes the output's key set
    # independent of the input's.
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 10.5, 10.0Hz)])
    signal = rand(ComplexF64, 5000)

    ts2 = Tracking.track(signal, ts, sampling_frequency)
    ts2 = add_satellite!(ts2; prn = 2, carrier_doppler = 200.0Hz)

    # Input keeps exactly its one satellite; output has two.
    @test length(get_sat_states(ts, :default)) == 1
    @test collect(keys(get_sat_states(ts, :default))) == [1]
    @test get_prn(ts, :default, 1) == 1
    @test length(get_sat_states(ts2, :default)) == 2
end

@testset "remove_satellite! kwarg form mutates in place" begin
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 5, carrier_doppler = 100.0Hz)
    ts = add_satellite!(ts; prn = 6, carrier_doppler = 200.0Hz)
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

@testset "add_satellite! honors a rebuilt estimator on handoff" begin
    import Tracking
    using Tracking: AbstractDopplerEstimator

    # Spec-conforming estimator that *rebuilds itself* on handoff — same
    # concrete type, replaced field — the style the
    # `update_estimator_on_handoff` docstring explicitly sanctions.
    struct RebuildingEstimator <: AbstractDopplerEstimator
        num_registered::Int
    end

    struct SatRebuildingState end

    Tracking.init_estimator_state(::RebuildingEstimator, ::TrackedSat) =
        SatRebuildingState()

    Tracking.update_estimator_on_handoff(est::RebuildingEstimator, new_sats) =
        RebuildingEstimator(est.num_registered + length(new_sats))

    # Keyword form: the returned TrackState carries the rebuilt estimator.
    ts = TrackState(; signal = GPSL1CA(), doppler_estimator = RebuildingEstimator(0))
    ts2 = add_satellite!(ts; prn = 1, carrier_doppler = 100.0Hz)
    @test ts2.doppler_estimator.num_registered == 1
    # The dictionary mutation is still in place — shared with the input.
    @test length(get_sat_states(ts, :default)) == 1
    @test get_sat_states(ts2, :default) === get_sat_states(ts, :default)

    # Escape-hatch overload behaves the same.
    sat = get_sat_state(ts2, :default, 1)
    sat2 = TrackedSat(
        sat.prn + 1, sat.code_phase, sat.code_doppler,
        sat.carrier_phase, sat.carrier_doppler,
        sat.signal_start_sample, sat.signals, sat.doppler_estimator_state,
    )
    ts3 = add_satellite!(ts2, :default, sat2)
    @test ts3.doppler_estimator.num_registered == 2
    @test length(get_sat_states(ts3, :default)) == 2

    # Estimators that update in place (or don't override the hook) keep
    # the in-place contract: the identical TrackState comes back.
    conventional_ts = TrackState(; signal = GPSL1CA())
    ret = add_satellite!(conventional_ts; prn = 5, carrier_doppler = 100.0Hz)
    @test ret === conventional_ts
end

# Type stability matters because every accessor is on the path between
# user code and the hot tracking loop. A widened return type here can
# silently propagate into a dynamic dispatch in the caller; `@inferred`
# trips the moment the compiler's narrowest prediction stops matching
# the actual return.
@testset "Accessor type stability — single-group, single-signal" begin
    track_state = TrackState(; signal = GPSL1CA())
    track_state = add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz)

    @inferred get_prn(track_state, 1)
    @inferred get_num_ants(track_state, 1)
    @inferred get_code_phase(track_state, 1)
    @inferred get_code_doppler(track_state, 1)
    @inferred get_carrier_phase(track_state, 1)
    @inferred get_carrier_doppler(track_state, 1)
    @inferred get_signal_start_sample(track_state, 1)
    @inferred get_integrated_samples(track_state, 1)
    @inferred get_correlator(track_state, 1)
    @inferred get_last_fully_integrated_correlator(track_state, 1)
    @inferred get_last_fully_integrated_filtered_prompt(track_state, 1)
    @inferred get_post_corr_filter(track_state, 1)
    @inferred get_cn0_estimator(track_state, 1)
    @inferred get_bit_buffer(track_state, 1)
    @inferred get_bits(track_state, 1)
    @inferred get_num_bits(track_state, 1)
    @inferred has_bit_or_secondary_code_been_found(track_state, 1)
    @inferred estimate_cn0(track_state, 1)
    @inferred get_signal(track_state)
    @inferred get_sat_state(track_state, 1)
    @test true   # gives the testset a passing assertion if all `@inferred` succeed
end

@testset "preferred_num_code_blocks_to_integrate validation (issue #128)" begin
    gpsl1 = GPSL1CA()
    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 0.0Hz)])

    # GPS L1 C/A has 20 code blocks per bit; only divisors keep integrations
    # aligned to bit boundaries, so only divisors are accepted.
    for valid in (1, 2, 4, 5, 10, 20)
        set_preferred_num_code_blocks_to_integrate!(track_state, 1, valid)
        @test get_preferred_num_code_blocks_to_integrate(track_state, 1) == valid
    end
    for invalid in (0, -1, 3, 7, 21, 40)
        @test_throws ArgumentError set_preferred_num_code_blocks_to_integrate!(
            track_state,
            1,
            invalid,
        )
    end
    # A rejected value leaves the previous one untouched.
    @test get_preferred_num_code_blocks_to_integrate(track_state, 1) == 20

    # Validated at construction, too — not only via the setter.
    @test_throws ArgumentError Tracking.TrackedSignal(
        gpsl1;
        preferred_num_code_blocks_to_integrate = 3,
    )
    @test Tracking.get_preferred_num_code_blocks_to_integrate(
        Tracking.TrackedSignal(gpsl1; preferred_num_code_blocks_to_integrate = 10),
    ) == 10

    # GPS L5I has 10 blocks per bit — 4 doesn't divide it (issue #128), 10 does.
    gpsl5 = GPSL5I()
    l5_state = TrackState(gpsl5, [TrackedSat(gpsl5, 1, 0.0, 0.0Hz)])
    @test_throws ArgumentError set_preferred_num_code_blocks_to_integrate!(l5_state, 1, 4)
    set_preferred_num_code_blocks_to_integrate!(l5_state, 1, 10)
    @test get_preferred_num_code_blocks_to_integrate(l5_state, 1) == 10

    # Pilot signals carry no data bits, so there is no bit boundary to
    # straddle — any length of at least one block is accepted (the actual
    # integration is clamped to one block until secondary-code-bounded
    # integration is wired up, see `calc_num_code_blocks_to_integrate`).
    gpsl1c_p = GPSL1C_P()
    p_state = TrackState(gpsl1c_p, [TrackedSat(gpsl1c_p, 1, 0.0, 0.0Hz)])
    set_preferred_num_code_blocks_to_integrate!(p_state, 1, 7)
    @test get_preferred_num_code_blocks_to_integrate(p_state, 1) == 7
    @test_throws ArgumentError set_preferred_num_code_blocks_to_integrate!(p_state, 1, 0)
end

@testset "per-signal preferred_num_code_blocks_to_integrate + reset_loop_filters!" begin
    gpsl1 = GPSL1CA()
    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])

    # Defaults to one block; settable per signal; persists.
    @test get_preferred_num_code_blocks_to_integrate(track_state, 1) == 1
    set_preferred_num_code_blocks_to_integrate!(track_state, 1, 20)
    @test get_preferred_num_code_blocks_to_integrate(track_state, 1) == 20
    set_preferred_num_code_blocks_to_integrate!(track_state, :default, 1, GPSL1CA, 5)   # by signal type
    @test get_preferred_num_code_blocks_to_integrate(track_state, 1) == 5
    # Selector-less (group, prn) form mirrors `get_…(ts, group, prn)`.
    set_preferred_num_code_blocks_to_integrate!(track_state, :default, 1, 4)
    @test get_preferred_num_code_blocks_to_integrate(track_state, :default, 1) == 4
    # No-id form mirrors `get_…(ts)` on a 1-group/1-sat state.
    set_preferred_num_code_blocks_to_integrate!(track_state, 10)
    @test get_preferred_num_code_blocks_to_integrate(track_state) == 10

    # The estimator state starts seeded at the init Doppler (100 Hz).
    @test get_doppler_estimator_state(get_sat_state(track_state, 1)).init_carrier_doppler ==
          100.0Hz

    # Simulate the loop having converged to a different carrier/code Doppler.
    sats = get_sat_states(track_state)
    sats[1] = TrackedSat(sats[1]; carrier_doppler = 137.0Hz, code_doppler = 0.5Hz)
    @test get_carrier_doppler(track_state, 1) == 137.0Hz
    # The estimator state is still anchored at the original init (100 Hz).
    @test get_doppler_estimator_state(get_sat_state(track_state, 1)).init_carrier_doppler ==
          100.0Hz

    reset_loop_filters!(track_state, 1)
    # Carrier/code Doppler preserved across the reset.
    @test get_carrier_doppler(track_state, 1) == 137.0Hz
    @test get_code_doppler(track_state, 1) == 0.5Hz
    # …and the estimator is re-seeded from the converged Doppler.
    @test get_doppler_estimator_state(get_sat_state(track_state, 1)).init_carrier_doppler ==
          137.0Hz
    @test get_doppler_estimator_state(get_sat_state(track_state, 1)).init_code_doppler == 0.5Hz

    # Whole-TrackState form resets every satellite and preserves Doppler.
    sats[1] = TrackedSat(sats[1]; carrier_doppler = 201.0Hz)
    reset_loop_filters!(track_state)
    @test get_carrier_doppler(track_state, 1) == 201.0Hz
    @test get_doppler_estimator_state(get_sat_state(track_state, 1)).init_carrier_doppler ==
          201.0Hz
end

end
