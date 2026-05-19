module AddSatelliteTest

# Tests for the v2 user-facing construction API:
#   * `TrackState(; signals = (...))` with NamedTuple and bare-tuple forms
#   * `add_satellite!` / `add_satellite` keyword API + escape-hatch overload
#   * Overwrite semantics on duplicate PRN
#   * Doppler estimator re-seeding when sats arrive with a mismatched type

using Test: @test, @testset, @test_throws, @inferred
using Unitful: Hz
using GNSSSignals:
    GPSL1CA, GalileoE1B, get_code_center_frequency_ratio
using Tracking:
    TrackState,
    TrackedSat,
    TrackedSignal,
    SignalGroup,
    NumAnts,
    add_satellite!,
    add_satellite,
    get_sat_state,
    get_sat_states,
    get_prn,
    get_carrier_doppler,
    get_code_doppler,
    ConventionalAssistedPLLAndDLL,
    ConventionalPLLAndDLL,
    SatConventionalPLLAndDLL,
    init_estimator_state

@testset "TrackState(; signals = ...) — single-group shortcut" begin
    track_state = TrackState(; signal = GPSL1CA())
    # Implicit `:default` group.
    @test isempty(get_sat_states(track_state, :default))
end

@testset "TrackState constructor argument errors" begin
    # Neither `signal` nor `signals` supplied.
    @test_throws ArgumentError TrackState()
    # Both supplied — must reject.
    @test_throws ArgumentError TrackState(;
        signal = GPSL1CA(),
        signals = (default = (GPSL1CA(),),),
    )
end

@testset "TrackState(; signals = (bare tuple,)) wraps to :default group" begin
    # When `signals` is a bare `Tuple{Vararg{AbstractGNSSSignal}}`,
    # `_normalize_signal_groups` wraps it as `(default = signals,)`.
    ts = TrackState(; signals = (GPSL1CA(),))
    @test isempty(get_sat_states(ts, :default))
end

@testset "TrackState rebuilds SignalGroup slot when estimator types differ" begin
    # `_normalize_group_entry` for a pre-built `SignalGroup` rebuilds the
    # empty dict when the template's estimator-state type doesn't match
    # the TrackState's `doppler_estimator` kwarg.
    custom = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sg = SignalGroup((GPSL1CA(),); num_ants = NumAnts(1))  # default estimator
    ts = TrackState(;
        signals = (legacy = sg,),
        doppler_estimator = custom,
    )
    add_satellite!(ts; prn = 1, group = :legacy, carrier_doppler = 0.0Hz)
    de_state = get_sat_state(ts, :legacy, 1).doppler_estimator_state
    @test de_state.carrier_loop_filter_bandwidth == 22.0Hz
end

@testset "TrackState passes through pre-populated SignalGroup unchanged" begin
    # A SignalGroup whose dict is already populated is trusted as-is — the
    # `!isempty(sats)` branch in `_normalize_group_entry`.
    using Dictionaries: insert!
    estimator = ConventionalAssistedPLLAndDLL()
    sat = TrackedSat(GPSL1CA(), 1, 10.5, 10.0Hz; doppler_estimator = estimator)
    sg_template = SignalGroup((GPSL1CA(),); num_ants = NumAnts(1),
                              doppler_estimator = estimator)
    insert!(sg_template.satellites, 1, sat)
    ts = TrackState(; signals = (legacy = sg_template,), doppler_estimator = estimator)
    @test length(get_sat_states(ts, :legacy)) == 1
end

@testset "TrackState(; signals = ...) — multi-group NamedTuple" begin
    track_state = TrackState(;
        signals = (
            legacy_gps = (GPSL1CA(),),
            galileo = (GalileoE1B(),),
        ),
    )
    @test isempty(get_sat_states(track_state, :legacy_gps))
    @test isempty(get_sat_states(track_state, :galileo))
end

@testset "add_satellite! — single-group shortcut (no `group=`)" begin
    track_state = TrackState(; signal = GPSL1CA())
    add_satellite!(track_state;
        prn = 11,
        code_phase = 10.5,
        carrier_doppler = 1234.0Hz,
    )
    @test length(get_sat_states(track_state, :default)) == 1
    @test get_prn(track_state, :default, 11) == 11
    @test get_carrier_doppler(track_state, :default, 11) == 1234.0Hz
    # code_doppler defaults to carrier_doppler × code-center-freq ratio.
    @test get_code_doppler(track_state, :default, 11) ≈
        1234.0Hz * get_code_center_frequency_ratio(GPSL1CA())
end

@testset "add_satellite! — multi-group with explicit `group=`" begin
    track_state = TrackState(;
        signals = (
            legacy = (GPSL1CA(),),
            galileo = (GalileoE1B(),),
        ),
    )
    add_satellite!(track_state;
        prn = 5, group = :legacy, carrier_doppler = 500.0Hz,
    )
    add_satellite!(track_state;
        prn = 11, group = :galileo, carrier_doppler = 2000.0Hz,
    )
    @test length(get_sat_states(track_state, :legacy)) == 1
    @test length(get_sat_states(track_state, :galileo)) == 1
    @test get_carrier_doppler(track_state, :legacy, 5) == 500.0Hz
    @test get_carrier_doppler(track_state, :galileo, 11) == 2000.0Hz
end

@testset "add_satellite! — overwrites on duplicate PRN" begin
    track_state = TrackState(; signal = GPSL1CA())
    add_satellite!(track_state;
        prn = 7, carrier_doppler = 100.0Hz,
    )
    add_satellite!(track_state;
        prn = 7, carrier_doppler = 200.0Hz,
    )
    @test length(get_sat_states(track_state, :default)) == 1
    @test get_carrier_doppler(track_state, :default, 7) == 200.0Hz
end

@testset "add_satellite — immutable variant" begin
    track_state = TrackState(; signal = GPSL1CA())
    new_track_state = add_satellite(track_state;
        prn = 3, carrier_doppler = 50.0Hz,
    )
    # Original is unchanged.
    @test isempty(get_sat_states(track_state, :default))
    # New one has the sat.
    @test length(get_sat_states(new_track_state, :default)) == 1
    @test get_carrier_doppler(new_track_state, :default, 3) == 50.0Hz
end

@testset "add_satellite!(track_state, group, sat) — escape hatch" begin
    track_state = TrackState(; signal = GPSL1CA())
    sat = TrackedSat(GPSL1CA(), 9, 5.0, 300.0Hz;
        doppler_estimator = ConventionalAssistedPLLAndDLL())
    add_satellite!(track_state, :default, sat)
    @test get_prn(track_state, :default, 9) == 9
    @test get_carrier_doppler(track_state, :default, 9) == 300.0Hz
end

@testset "Escape-hatch add_satellite! rejects sat of wrong slot type" begin
    # `_assert_sat_matches_slot_type` should raise an ArgumentError when
    # the sat's concrete type doesn't match the group's fixed slot type.
    track_state = TrackState(; signal = GPSL1CA())
    galileo_sat = TrackedSat(GalileoE1B(), 1, 10.5, 100.0Hz)
    @test_throws ArgumentError add_satellite!(track_state, :default, galileo_sat)
end

@testset "Custom doppler estimator config seeds per-sat state correctly" begin
    custom = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    track_state = TrackState(; signal = GPSL1CA(), doppler_estimator = custom)
    add_satellite!(track_state;
        prn = 4, carrier_doppler = 0.0Hz,
    )
    de_state = get_sat_state(track_state, :default, 4).doppler_estimator_state
    @test de_state isa SatConventionalPLLAndDLL
    @test de_state.carrier_loop_filter_bandwidth == 22.0Hz
    @test de_state.code_loop_filter_bandwidth == 1.5Hz
end

end
