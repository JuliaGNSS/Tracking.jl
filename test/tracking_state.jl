module TrackingStateTest

using Test: @test, @testset, @inferred, @test_throws
using GNSSSignals: CBOC, TMBOC
using TrackingLoopFilters: ThirdOrderAssistedBilinearLF
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
    get_code_phase,
    get_code_doppler,
    get_carrier_phase,
    get_carrier_doppler,
    get_num_integrated_samples,
    get_next_sample_index,
    get_correlator_bank,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_beamformed_prompt,
    get_beamformed_prompts,
    get_beamformer,
    get_cn0_estimator,
    get_bit_synchronizer,
    get_num_bits,
    get_sat_state,
    get_sat_states,
    merge_sats,
    remove_satellite!,
    remove_satellite,
    get_preferred_num_code_periods_to_integrate,
    set_preferred_num_code_periods_to_integrate!,
    reset_loop_filters!,
    get_loop_state,
    SignalGroup
import TrackingLoops
using TrackingLoops:
    get_num_antennas,
    VeryEarlyPromptLateCorrelatorBank,
    get_soft_bits,
    is_bit_or_secondary_code_synced,
    estimate_cn0,
    DefaultBeamformer,
    NoiseReferencedCN0Estimator,
    BitSynchronizer,
    NumAntennas,
    default_wide_carrier_loop_filter_bandwidth,
    ScalarTrackingLoop

# No real signal pair mixes chip rates on one band, so fake an L1 signal at
# double the chip rate (#129). Validation throws before any other signal API
# is needed.
struct FakeDoubleRateL1Signal <: AbstractGNSSSignal{Matrix{Int16}} end
GNSSSignals.get_band(::FakeDoubleRateL1Signal) = GNSSSignals.L1()
GNSSSignals.get_chip_rate(::FakeDoubleRateL1Signal) = 2_046_000Hz

@testset "SignalGroup rejects mixed bands and mixed chip rates (issue #129)" begin
    # A group correlates against one band's samples.
    @test_throws ArgumentError SignalGroup((GPSL1CA(), GPSL5I()))
    @test_throws ArgumentError TrackState(; signals = (mix = (GPSL1CA(), GPSL5I()),))
    @test_throws ArgumentError SignalGroup((GPSL1CA(),); band = GNSSSignals.L5())

    # The shared `code_phase` advances at `signals[1]`'s chip rate.
    @test_throws ArgumentError SignalGroup((GPSL1CA(), FakeDoubleRateL1Signal()))
    @test_throws ArgumentError TrackState(;
        signals = (l1 = (GPSL1CA(), FakeDoubleRateL1Signal()),),
    )

    @test SignalGroup((GPSL1C_P(TMBOC), GPSL1C_D(), GPSL1CA())) isa SignalGroup
    ts = TrackState(;
        signals = (l1 = (GPSL1C_P(TMBOC), GPSL1C_D(), GPSL1CA()), l5 = (GPSL5I(),)),
    )
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
    @test get_num_antennas(track_state, 1) == 1
    @test get_num_integrated_samples(track_state, 1) == 0
    @test get_next_sample_index(track_state, 1) == 1
    @test get_correlator_bank(track_state, 1).correlator_outputs == zeros(3)
    @test get_last_fully_integrated_correlator(track_state, 1).correlator_outputs ==
          zeros(3)
    @test get_beamformer(track_state, 1) isa DefaultBeamformer
    @test get_cn0_estimator(track_state, 1) isa NoiseReferencedCN0Estimator
    @test get_bit_synchronizer(track_state, 1) isa BitSynchronizer

    @test get_sat_states(track_state)[1].loop_state.reference_carrier_doppler == 10.0Hz
    @test get_sat_states(track_state)[2].loop_state.reference_carrier_doppler == 20.0Hz

    track_state2 = @inferred TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz), TrackedSat(gpsl1, 2, 11.5, 20.0Hz)],
    )

    @test @inferred(get_signal(track_state2)) isa GPSL1CA
    @test @inferred(get_sat_state(track_state2, 1)).prn == 1
    @test @inferred(get_sat_state(track_state2, 2)).prn == 2

    @test get_sat_states(track_state2)[1].loop_state.reference_carrier_doppler == 10.0Hz
    @test get_sat_states(track_state2)[2].loop_state.reference_carrier_doppler == 20.0Hz

    sat_states_num_antennas2 = [
        TrackedSat(gpsl1, 1, 10.5, 10.0Hz; num_antennas = NumAntennas(2)),
        TrackedSat(gpsl1, 2, 11.5, 20.0Hz; num_antennas = NumAntennas(2)),
    ]
end

@testset "Positional TrackState(satellites::SatelliteDicts) infers default tracking loop" begin
    # One shared auto-bandwidth tracking loop seeds each group's sats from that
    # group's own driver signal.
    gpsl1 = GPSL1CA()
    galileo = GalileoE1B(CBOC)
    tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    sats = (
        gps = dictionary([1 => TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop)]),
        gal = dictionary([2 => TrackedSat(galileo, 2, 11.5, 20.0Hz; tracking_loop)]),
    )
    ts = TrackState(sats)
    @test length(get_sat_states(ts, :gps)) == 1
    @test length(get_sat_states(ts, :gal)) == 1
    @test get_sat_states(ts, :gps)[1].loop_state.bandwidths.wide_carrier ==
          default_wide_carrier_loop_filter_bandwidth(gpsl1)
    @test get_sat_states(ts, :gal)[2].loop_state.bandwidths.wide_carrier ==
          default_wide_carrier_loop_filter_bandwidth(galileo)
end

@testset "Positional TrackState(dict) default tracking loop inference" begin
    gpsl1 = GPSL1CA()
    tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    sat = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop)
    ts = TrackState(dictionary([1 => sat]))
    @test length(get_sat_states(ts)) == 1
end

@testset "Positional TrackState constructor rejects sats with a different tracking loop" begin
    gpsl1 = GPSL1CA()
    sat_default = TrackedSat(gpsl1, 1, 10.5, 10.0Hz)  # default tracking loop
    different = ScalarTrackingLoop(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test_throws ArgumentError TrackState(gpsl1, [sat_default]; tracking_loop = different)
end

@testset "Slot-vector copy helpers share vs detach Dictionary Indices (#123)" begin
    gpsl1 = GPSL1CA()
    tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    a = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop)
    sats = dictionary([1 => a])

    # The per-iteration copy inside `track` shares the key set; the boundary
    # copy detaches it too.
    shared = Tracking._copy_slot_vector(sats)
    @test shared.values !== sats.values
    @test keys(shared) === keys(sats)

    detached = Tracking._detach_slot_vector(sats)
    @test detached.values !== sats.values
    @test keys(detached) !== keys(sats)
end

@testset "Immutable boundary copies do not share Dictionary Indices (#123)" begin
    # With a shared key set, `add_satellite!` on the copy would grow the
    # original's keys without its values.
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1, carrier_doppler = 100.0Hz)

    ts2 = Tracking._prepare_for_new_buffer(ts)
    ts2 = add_satellite!(ts2; prn = 2, carrier_doppler = 200.0Hz)

    @test length(get_sat_states(ts, :default)) == 1
    @test collect(keys(get_sat_states(ts, :default))) == [1]
    @test get_prn(ts, :default, 1) == 1
    @test length(get_sat_states(ts2, :default)) == 2

    ts2 = remove_satellite!(ts2; prn = 1)
    @test length(get_sat_states(ts, :default)) == 1
    @test get_prn(ts, :default, 1) == 1
end

@testset "track output is structurally detached from input (#123)" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 10.5, 10.0Hz)])
    samples = rand(ComplexF64, 5000)

    ts2 = Tracking.track(samples, ts, sampling_frequency)
    ts2 = add_satellite!(ts2; prn = 2, carrier_doppler = 200.0Hz)

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

@testset "remove_satellite(!) leaves no #undef hole → track! is safe (issue #182)" begin
    # `Dictionaries.delete!` compacts only past ~1/3 deletions, so one removal
    # out of ten leaves an `#undef` slot in `values`, which the hot loops
    # iterate directly.
    samples = rand(ComplexF32, 20000)
    build() = foldl(
        (ts, prn) -> add_satellite!(ts; prn = prn, carrier_doppler = (100.0 * prn)Hz),
        1:10;
        init = TrackState(; signal = GPSL1CA()),
    )

    for remove in (
        ts -> remove_satellite(ts; prn = 3),   # immutable variant
        ts -> remove_satellite!(ts; prn = 3),   # in-place variant
    )
        ts = remove(build())
        sats = get_sat_states(ts, :default)
        @test length(sats) == 9
        @test !haskey(sats, 3)
        vals = sats.values
        @test length(vals) == 9
        @test all(i -> isassigned(vals, i), eachindex(vals))
        ts2 = Tracking.track(samples, ts, 5e6Hz)
        @test length(get_sat_states(ts2, :default)) == 9
    end
end

@testset "remove_satellite! throws KeyError for a missing PRN" begin
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 5, carrier_doppler = 100.0Hz)
    @test_throws KeyError remove_satellite!(ts; prn = 99)
    @test_throws KeyError remove_satellite(ts; prn = 99)
    @test get_prn(ts, :default, 5) == 5
end

@testset "remove_satellite!/remove_satellite — multi-group omitted `group_key=` errors" begin
    ts = TrackState(; signals = (legacy = (GPSL1CA(),), galileo = (GalileoE1B(CBOC),)))
    ts = add_satellite!(ts; prn = 5, group_key = :legacy, carrier_doppler = 100.0Hz)
    @test_throws ArgumentError remove_satellite!(ts; prn = 5)
    @test_throws ArgumentError remove_satellite(ts; prn = 5)
    # A single named group still infers the group when omitted.
    ts2 = TrackState(; signals = (gps = (GPSL1CA(),),))
    ts2 = add_satellite!(ts2; prn = 5, carrier_doppler = 100.0Hz)
    @test remove_satellite!(ts2; prn = 5) === ts2
    @test isempty(get_sat_states(ts2, :gps))
end

@testset "Add and remove satellite state to and from track state" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    sat_states = [TrackedSat(gpsl1, 1, 10.5, 10.0Hz), TrackedSat(gpsl1, 2, 11.5, 20.0Hz)]

    track_state = @inferred TrackState(gpsl1, sat_states)

    new_track_state =
        @inferred merge_sats(track_state, 1, TrackedSat(gpsl1, 3, 5.5, 80.0Hz))

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

    filtered_new_new_track_state = @inferred remove_satellite(
        @inferred(remove_satellite(new_new_track_state; prn = 1));
        prn = 6,
    )

    @test length(@inferred(get_sat_states(filtered_new_new_track_state))) == 2
    @test @inferred(get_prn(filtered_new_new_track_state, 3)) == 3
    @test @inferred(get_prn(filtered_new_new_track_state, 8)) == 8
end

@testset "Add and remove satellite state to track state with multiple constellations" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    galileo_e1b = GalileoE1B(CBOC)

    tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    gps_sats = dictionary([
        1 => TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop),
        2 => TrackedSat(gpsl1, 2, 11.5, 20.0Hz; tracking_loop),
    ])
    gal_sats = dictionary([
        2 => TrackedSat(galileo_e1b, 2, 10.5, 10.0Hz; tracking_loop),
        3 => TrackedSat(galileo_e1b, 3, 11.5, 20.0Hz; tracking_loop),
    ])
    satellites = (gps = gps_sats, gal = gal_sats)

    track_state = @inferred TrackState(satellites; tracking_loop)

    new_track_state =
        @inferred merge_sats(track_state, 1, TrackedSat(gpsl1, 3, 5.5, 80.0Hz))

    @test length(@inferred(get_sat_states(new_track_state, Val(:gps)))) == 3
    @test length(@inferred(get_sat_states(new_track_state, Val(:gal)))) == 2
    @test length(@inferred(get_sat_states(new_track_state, Val(1)))) == 3
    @test length(@inferred(get_sat_states(new_track_state, Val(2)))) == 2
    @test @inferred(get_prn(new_track_state, Val(:gps), 3)) == 3
    @test @inferred(get_prn(new_track_state, Val(1), 3)) == 3

    filtered_new_track_state =
        @inferred remove_satellite(new_track_state; prn = 2, group_key = :gps)
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
        @inferred remove_satellite(new_new_track_state; prn = 2, group_key = :gps)
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

    filtered_new_new_new_track_state = @inferred remove_satellite(
        @inferred(remove_satellite(new_new_new_track_state; prn = 1, group_key = :gps));
        prn = 6,
        group_key = :gps,
    )

    @test length(get_sat_states(filtered_new_new_new_track_state, :gps)) == 2
    @test length(get_sat_states(filtered_new_new_new_track_state, :gal)) == 2
    @test get_prn(filtered_new_new_new_track_state, :gps, 3) == 3
    @test get_prn(filtered_new_new_new_track_state, :gps, 8) == 8
end

@testset "TrackedSat helpers" begin
    using Tracking: TrackedSat, get_loop_state, _prepare_for_new_buffer!
    using TrackingLoops: ScalarLoopState

    gpsl1 = GPSL1CA()
    tracking_loop = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)

    sat = TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop)
    sat2 = TrackedSat(gpsl1, 2, 11.5, 20.0Hz; tracking_loop)

    ts_from_vec = TrackState(gpsl1, [sat, sat2]; tracking_loop)
    @test length(get_sat_states(ts_from_vec)) == 2
    ts_from_one = TrackState(gpsl1, sat; tracking_loop)
    @test length(get_sat_states(ts_from_one)) == 1
    ts_from_dict = TrackState(dictionary([1 => sat, 2 => sat2]); tracking_loop)
    @test length(get_sat_states(ts_from_dict)) == 2

    @test sat.prn == 1
    @test get_loop_state(sat) isa ScalarLoopState

    track_state = TrackState(gpsl1, [sat, sat2]; tracking_loop)
    track_state2 = _prepare_for_new_buffer!(track_state)
    @test track_state2 === track_state
    for tsat in get_sat_states(track_state).values
        @test tsat.next_sample_index == 1
    end
end

@testset "merge_sats invokes update_tracking_loop_on_new_sats" begin
    import Tracking
    using TrackingLoops: AbstractTrackingLoop

    # Immutable, with its state in vectors the hook mutates in place.
    struct CountingTrackingLoop <: AbstractTrackingLoop
        sats_added::Vector{Int}
        last_batch_prns::Vector{Int}
    end
    CountingTrackingLoop() = CountingTrackingLoop(Int[0], Int[])

    struct CountingLoopState end

    TrackingLoops.init_loop_state(::CountingTrackingLoop, ::TrackedSat) =
        CountingLoopState()

    function Tracking.update_tracking_loop_on_new_sats(loop::CountingTrackingLoop, new_sats)
        loop.sats_added[1] += length(new_sats)
        empty!(loop.last_batch_prns)
        for s in new_sats
            push!(loop.last_batch_prns, s.prn)
        end
        return loop
    end

    gpsl1 = GPSL1CA()
    tracking_loop = CountingTrackingLoop()
    track_state = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop)];
        tracking_loop,
    )
    @test track_state.tracking_loop.sats_added[1] == 0

    track_state2 =
        merge_sats(track_state, 1, TrackedSat(gpsl1, 2, 11.5, 20.0Hz; tracking_loop))
    @test track_state2.tracking_loop === tracking_loop  # same object, mutated in place
    @test track_state2.tracking_loop.sats_added[1] == 1
    @test track_state2.tracking_loop.last_batch_prns == [2]

    track_state3 = merge_sats(
        track_state2,
        [
            TrackedSat(gpsl1, 3, 5.5, 80.0Hz; tracking_loop),
            TrackedSat(gpsl1, 4, 6.5, 90.0Hz; tracking_loop),
        ],
    )
    @test track_state3.tracking_loop.sats_added[1] == 3
    @test sort(track_state3.tracking_loop.last_batch_prns) == [3, 4]

    # Without an override the tracking loop is unchanged.
    default_estimator = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)
    default_track_state = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz; tracking_loop = default_estimator)];
        tracking_loop = default_estimator,
    )
    merged_default = merge_sats(
        default_track_state,
        1,
        TrackedSat(gpsl1, 2, 11.5, 20.0Hz; tracking_loop = default_estimator),
    )
    @test merged_default.tracking_loop === default_estimator
end

@testset "add_satellite! honors a tracking loop rebuilt for new sats" begin
    import Tracking
    using TrackingLoops: AbstractTrackingLoop

    # Rebuilds itself for new sats, same concrete type, as the
    # `update_tracking_loop_on_new_sats` docstring allows.
    struct RebuildingTrackingLoop <: AbstractTrackingLoop
        num_registered::Int
    end

    struct SatRebuildingState end

    TrackingLoops.init_loop_state(::RebuildingTrackingLoop, ::TrackedSat) =
        SatRebuildingState()

    Tracking.update_tracking_loop_on_new_sats(loop::RebuildingTrackingLoop, new_sats) =
        RebuildingTrackingLoop(loop.num_registered + length(new_sats))

    ts = TrackState(; signal = GPSL1CA(), tracking_loop = RebuildingTrackingLoop(0))
    ts2 = add_satellite!(ts; prn = 1, carrier_doppler = 100.0Hz)
    @test ts2.tracking_loop.num_registered == 1
    # The dictionary is still mutated in place.
    @test length(get_sat_states(ts, :default)) == 1
    @test get_sat_states(ts2, :default) === get_sat_states(ts, :default)

    sat = get_sat_state(ts2, :default, 1)
    sat2 = TrackedSat(
        sat.prn + 1,
        sat.code_phase,
        sat.code_doppler,
        sat.carrier_phase,
        sat.carrier_doppler,
        sat.next_sample_index,
        sat.signals,
        sat.loop_state,
    )
    ts3 = add_satellite!(ts2, :default, sat2)
    @test ts3.tracking_loop.num_registered == 2
    @test length(get_sat_states(ts3, :default)) == 2

    # Otherwise the identical TrackState comes back.
    conventional_ts = TrackState(; signal = GPSL1CA())
    ret = add_satellite!(conventional_ts; prn = 5, carrier_doppler = 100.0Hz)
    @test ret === conventional_ts
end

# A widened accessor return type turns into dynamic dispatch in the caller.
@testset "Accessor type stability — single-group, single-signal" begin
    track_state = TrackState(; signal = GPSL1CA())
    track_state =
        add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz)

    @inferred get_prn(track_state, 1)
    @inferred get_num_antennas(track_state, 1)
    @inferred get_code_phase(track_state, 1)
    @inferred get_code_doppler(track_state, 1)
    @inferred get_carrier_phase(track_state, 1)
    @inferred get_carrier_doppler(track_state, 1)
    @inferred get_next_sample_index(track_state, 1)
    @inferred get_num_integrated_samples(track_state, 1)
    @inferred get_correlator_bank(track_state, 1)
    @inferred get_last_fully_integrated_correlator(track_state, 1)
    @inferred get_last_fully_integrated_beamformed_prompt(track_state, 1)
    @inferred get_beamformer(track_state, 1)
    @inferred get_cn0_estimator(track_state, 1)
    @inferred get_bit_synchronizer(track_state, 1)
    @inferred get_soft_bits(track_state, 1)
    @inferred get_num_bits(track_state, 1)
    @inferred is_bit_or_secondary_code_synced(track_state, 1)
    @inferred estimate_cn0(track_state, 1)
    @inferred get_signal(track_state)
    @inferred get_sat_state(track_state, 1)
    @test true
end

@testset "preferred_num_code_periods_to_integrate validation (issue #128)" begin
    gpsl1 = GPSL1CA()
    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 0.0Hz)])

    # Only divisors of the 20 code periods per bit keep integrations on the bit grid.
    for valid in (1, 2, 4, 5, 10, 20)
        set_preferred_num_code_periods_to_integrate!(track_state, 1, valid)
        @test get_preferred_num_code_periods_to_integrate(track_state, 1) == valid
    end
    for invalid in (0, -1, 3, 7, 21, 40)
        @test_throws ArgumentError set_preferred_num_code_periods_to_integrate!(
            track_state,
            1,
            invalid,
        )
    end
    # A rejected value leaves the previous one untouched.
    @test get_preferred_num_code_periods_to_integrate(track_state, 1) == 20

    # Validated at construction, too.
    @test_throws ArgumentError Tracking.TrackedSignal(
        gpsl1;
        preferred_num_code_periods_to_integrate = 3,
    )
    @test Tracking.get_preferred_num_code_periods_to_integrate(
        Tracking.TrackedSignal(gpsl1; preferred_num_code_periods_to_integrate = 10),
    ) == 10

    # GPS L5I has 10 code periods per bit.
    gpsl5 = GPSL5I()
    l5_state = TrackState(gpsl5, [TrackedSat(gpsl5, 1, 0.0, 0.0Hz)])
    @test_throws ArgumentError set_preferred_num_code_periods_to_integrate!(l5_state, 1, 4)
    set_preferred_num_code_periods_to_integrate!(l5_state, 1, 10)
    @test get_preferred_num_code_periods_to_integrate(l5_state, 1) == 10

    # A pilot has no bit grid, so any length ≥ 1 is accepted; at runtime it is
    # clamped to a divisor of the secondary-code period (#134).
    gpsl1c_p = GPSL1C_P(TMBOC)
    p_state = TrackState(gpsl1c_p, [TrackedSat(gpsl1c_p, 1, 0.0, 0.0Hz)])
    set_preferred_num_code_periods_to_integrate!(p_state, 1, 7)
    @test get_preferred_num_code_periods_to_integrate(p_state, 1) == 7
    @test_throws ArgumentError set_preferred_num_code_periods_to_integrate!(p_state, 1, 0)
end

@testset "get_beamformed_prompts TrackState accessor (issue #134)" begin
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1)
    @test get_beamformed_prompts(ts, 1) isa Vector{ComplexF64}
    @test isempty(get_beamformed_prompts(ts, 1))
    @test get_beamformed_prompts(ts, :default, 1, 1) isa Vector{ComplexF64}
    @test get_beamformed_prompts(ts, :default, 1, GPSL1CA) isa Vector{ComplexF64}
end

@testset "get_signal works on a declared-but-unpopulated group (issue #134)" begin
    ts = TrackState(; signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(CBOC),)))
    @test get_signal(ts, :gps) isa GPSL1CA
    @test get_signal(ts, :gal) isa GalileoE1B
    @test get_signal(ts, 1) isa GPSL1CA
    @test @inferred(get_signal(ts, Val(:gal))) isa GalileoE1B
    single = TrackState(; signal = GPSL1CA())
    @test get_signal(single) isa GPSL1CA
end

@testset "remove_satellite throws KeyError carrying the prn (issue #134)" begin
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1)
    err = try
        remove_satellite(ts; prn = 99)
        nothing
    catch e
        e
    end
    @test err isa KeyError
    @test err.key == 99
end

@testset "TrackState from empty SatelliteDicts throws friendly error (issue #134)" begin
    gpsl1 = GPSL1CA()
    sat = TrackedSat(gpsl1, 1, 0.0, 0.0Hz)
    empty_dict = Dictionary{Int,typeof(sat)}(Int[], typeof(sat)[])
    # The ArgumentError from `_signal_group_from_dict`, not a BoundsError.
    @test_throws ArgumentError TrackState((gps = empty_dict,))
end

@testset "TrackState(satellites::SatelliteDicts) validates loop state types (issue #134)" begin
    gpsl1 = GPSL1CA()
    sat_default = TrackedSat(gpsl1, 1, 10.5, 10.0Hz)  # assisted default tracking loop
    different = ScalarTrackingLoop(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test_throws ArgumentError TrackState(
        (gps = dictionary([1 => sat_default]),);
        tracking_loop = different,
    )
end

@testset "TrackState(; signals = SignalGroup) validates pre-populated sats' tracking loop (issue #134)" begin
    ts0 = TrackState(; signal = GPSL1CA())
    ts0 = add_satellite!(ts0; prn = 1)
    populated_group = ts0.groups.default
    different = ScalarTrackingLoop(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test_throws ArgumentError TrackState(;
        signals = (gps = populated_group,),
        tracking_loop = different,
    )
end

@testset "merge_sats rejects sats whose full slot type mismatches (issue #134)" begin
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 10.5, 10.0Hz)])
    # Same tracking loop state, different correlator: an ArgumentError, not a
    # MethodError downstream.
    odd = TrackedSat(
        gpsl1,
        2,
        0.0,
        0.0Hz;
        correlator_bank = VeryEarlyPromptLateCorrelatorBank(),
    )
    @test_throws ArgumentError merge_sats(ts, 1, odd)
end

@testset "reset_loop_filters! keeps per-sat bandwidth override and clears FLL history (issue #134)" begin
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])
    sats = get_sat_states(ts)
    old = sats[1]
    overridden_state = TrackingLoops.ScalarLoopState(
        get_loop_state(old);
        bandwidths = TrackingLoops.LoopBandwidths(;
            wide_carrier = 5.0Hz,
            narrow_carrier = 4.0Hz,
            code = 0.25Hz,
        ),
    )
    fake_prompt_signals = map(
        s -> Tracking.TrackedSignal(
            s;
            last_fully_integrated_beamformed_prompt = 1.0 + 2.0im,
        ),
        old.signals,
    )
    sats[1] = TrackedSat(old; signals = fake_prompt_signals, loop_state = overridden_state)

    reset_loop_filters!(ts, 1)
    state = get_loop_state(get_sat_state(ts, 1))
    @test state.bandwidths.wide_carrier == 5.0Hz
    @test state.bandwidths.narrow_carrier == 4.0Hz
    @test state.bandwidths.code == 0.25Hz
    # The FLL's previous prompt is cleared, so the first discriminator after
    # the reset does not span the old integration.
    @test get_last_fully_integrated_beamformed_prompt(ts, 1) == 0.0 + 0.0im
end

@testset "per-signal preferred_num_code_periods_to_integrate + reset_loop_filters!" begin
    gpsl1 = GPSL1CA()
    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])

    @test get_preferred_num_code_periods_to_integrate(track_state, 1) == 1
    set_preferred_num_code_periods_to_integrate!(track_state, 1, 20)
    @test get_preferred_num_code_periods_to_integrate(track_state, 1) == 20
    set_preferred_num_code_periods_to_integrate!(track_state, :default, 1, GPSL1CA, 5)   # by signal type
    @test get_preferred_num_code_periods_to_integrate(track_state, 1) == 5
    set_preferred_num_code_periods_to_integrate!(track_state, :default, 1, 4)
    @test get_preferred_num_code_periods_to_integrate(track_state, :default, 1) == 4
    set_preferred_num_code_periods_to_integrate!(track_state, 10)
    @test get_preferred_num_code_periods_to_integrate(track_state) == 10

    @test get_loop_state(get_sat_state(track_state, 1)).reference_carrier_doppler == 100.0Hz

    # As if the loop had converged elsewhere.
    sats = get_sat_states(track_state)
    sats[1] = TrackedSat(sats[1]; carrier_doppler = 137.0Hz, code_doppler = 0.5Hz)
    @test get_carrier_doppler(track_state, 1) == 137.0Hz
    @test get_loop_state(get_sat_state(track_state, 1)).reference_carrier_doppler == 100.0Hz

    reset_loop_filters!(track_state, 1)
    # The reset keeps the Doppler and re-seeds the tracking loop from it.
    @test get_carrier_doppler(track_state, 1) == 137.0Hz
    @test get_code_doppler(track_state, 1) == 0.5Hz
    @test get_loop_state(get_sat_state(track_state, 1)).reference_carrier_doppler == 137.0Hz
    @test get_loop_state(get_sat_state(track_state, 1)).reference_code_doppler == 0.5Hz

    sats[1] = TrackedSat(sats[1]; carrier_doppler = 201.0Hz)
    reset_loop_filters!(track_state)
    @test get_carrier_doppler(track_state, 1) == 201.0Hz
    @test get_loop_state(get_sat_state(track_state, 1)).reference_carrier_doppler == 201.0Hz
end

@testset "non-mutating constructors give their result its own sample offsets" begin
    gpsl1 = GPSL1CA()
    band = GNSSSignals.get_band_id(GNSSSignals.get_band(gpsl1))
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 10.5, 10.0Hz)])
    ts.sample_index_offsets[band] = 4000
    derived_states = (
        merge_sats(ts, TrackedSat(gpsl1, 2, 11.5, 20.0Hz)),
        Tracking.add_satellite(ts, :default, TrackedSat(gpsl1, 3, 12.5, 30.0Hz)),
        remove_satellite(ts; prn = 1),
        Tracking._prepare_for_new_buffer(ts),
    )
    for derived in derived_states
        @test derived.sample_index_offsets == ts.sample_index_offsets
        derived.sample_index_offsets[band] += 4000
        @test ts.sample_index_offsets[band] == 4000
    end
    # The in-place mutators keep sharing them with the state they return.
    @test add_satellite!(ts; prn = 4).sample_index_offsets === ts.sample_index_offsets
end

@testset "constructors carry sample offsets over" begin
    gpsl1 = GPSL1CA()
    l1 = GNSSSignals.get_band_id(GNSSSignals.get_band(gpsl1))
    l5 = GNSSSignals.get_band_id(GNSSSignals.get_band(GPSL5I()))
    @test TrackState(; signal = gpsl1).sample_index_offsets == Dict(l1 => 0)
    @test TrackState(; signal = gpsl1, sample_index_offsets = Dict(l1 => 8000)).sample_index_offsets ==
          Dict(l1 => 8000)
    # A rebuilt state with a new band keeps the old offset, starts the new band
    # at 0, and does not alias the input.
    old = TrackState(; signal = gpsl1, sample_index_offsets = (; l1 => 8000))
    rebuilt = TrackState(;
        signals = (l1 = (gpsl1,), l5 = (GPSL5I(),)),
        tracking_loop = old.tracking_loop,
        sample_index_offsets = old.sample_index_offsets,
    )
    @test rebuilt.sample_index_offsets == Dict(l1 => 8000, l5 => 0)
    @test rebuilt.sample_index_offsets !== old.sample_index_offsets
    @test TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 10.5, 10.0Hz)];
        sample_index_offsets = (; l1 => 4000),
    ).sample_index_offsets == Dict(l1 => 4000)
    @test TrackState(
        dictionary((1 => TrackedSat(gpsl1, 1, 10.5, 10.0Hz),));
        sample_index_offsets = (; l1 => 4000),
    ).sample_index_offsets == Dict(l1 => 4000)
    @test TrackState(
        (default = dictionary((1 => TrackedSat(gpsl1, 1, 10.5, 10.0Hz),)),);
        sample_index_offsets = (; l1 => 4000),
    ).sample_index_offsets == Dict(l1 => 4000)
    @test_throws ArgumentError TrackState(;
        signal = gpsl1,
        sample_index_offsets = (; l5 => 0),
    )
    @test_throws ArgumentError TrackState(;
        signal = gpsl1,
        sample_index_offsets = (; l1 => -1),
    )
end

end
