module AddSatelliteTest

# Tests for the v2 user-facing construction API:
#   * `TrackState(; signals = (...))` with NamedTuple and bare-tuple forms
#   * `add_satellite!` / `add_satellite` keyword API + escape-hatch overload
#   * Overwrite semantics on duplicate PRN
#   * Doppler estimator re-seeding when sats arrive with a mismatched type

using Test: @test, @testset, @test_throws, @inferred
using Unitful: Hz
using GNSSSignals:
    GPSL1CA, GPSL1C_D, GPSL1C_P, GalileoE1B, get_code_center_frequency_ratio
using Acquisition: Acquisition, AcquisitionResults
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

# Acquisition v2 added num_blocks/block_size fields (FM-DBZP backend).
# Keep both constructors supported as long as Project.toml's Acquisition
# compat range does.
function _make_acq(signal, prn, code_phase, carrier_doppler)
    args = (
        signal, prn, 5e6Hz, carrier_doppler, code_phase,
        45.0, 1.0, 10.0, 1, randn(10, 10), -500:100.0:500,
    )
    pkgversion(Acquisition) >= v"2" ?
        AcquisitionResults(args..., 1, length(-500:100.0:500)) :
        AcquisitionResults(args...)
end

@testset "add_satellite!(ts, acq) — single-group shortcut" begin
    ts = TrackState(; signal = GPSL1CA())
    acq = _make_acq(GPSL1CA(), 7, 524.6, 100.0Hz)
    ret = add_satellite!(ts, acq)
    @test ret === ts
    @test get_prn(ts, :default, 7) == 7
    @test get_carrier_doppler(ts, :default, 7) == 100.0Hz
    sat = get_sat_state(ts, :default, 7)
    @test sat.code_phase == 524.6
end

@testset "add_satellite(ts, acq) immutable returns a new TrackState" begin
    ts = TrackState(; signal = GPSL1CA())
    acq = _make_acq(GPSL1CA(), 9, 12.5, -250.0Hz)
    new_ts = add_satellite(ts, acq)
    @test new_ts !== ts
    @test isempty(get_sat_states(ts, :default))
    @test get_prn(new_ts, :default, 9) == 9
    @test get_carrier_doppler(new_ts, :default, 9) == -250.0Hz
end

@testset "add_satellite!(ts, acqs::Vector) batch form" begin
    ts = TrackState(; signal = GPSL1CA())
    acqs = [
        _make_acq(GPSL1CA(), 1, 0.0, 100.0Hz),
        _make_acq(GPSL1CA(), 2, 10.0, 200.0Hz),
        _make_acq(GPSL1CA(), 3, 20.0, 300.0Hz),
    ]
    add_satellite!(ts, acqs)
    @test length(get_sat_states(ts, :default)) == 3
    @test get_carrier_doppler(ts, :default, 2) == 200.0Hz
end

@testset "add_satellite!(ts, acq; group) — multi-group with explicit `group=`" begin
    ts = TrackState(; signals = (
        gps = (GPSL1CA(),),
        gal = (GalileoE1B(),),
    ))
    gps_acq = _make_acq(GPSL1CA(), 11, 1.5, 50.0Hz)
    gal_acq = _make_acq(GalileoE1B(), 12, 2.5, 60.0Hz)
    add_satellite!(ts, gps_acq; group = :gps)
    add_satellite!(ts, gal_acq; group = :gal)
    @test get_prn(ts, :gps, 11) == 11
    @test get_prn(ts, :gal, 12) == 12
end

@testset "add_satellite!(ts, acq) — multi-group auto-routes by signal type" begin
    ts = TrackState(; signals = (
        gps = (GPSL1CA(),),
        gal = (GalileoE1B(),),
    ))
    gps_acq = _make_acq(GPSL1CA(), 11, 1.5, 50.0Hz)
    gal_acq = _make_acq(GalileoE1B(), 12, 2.5, 60.0Hz)
    # No `group=`: each acq lands in the matching group automatically.
    add_satellite!(ts, gps_acq)
    add_satellite!(ts, gal_acq)
    @test get_prn(ts, :gps, 11) == 11
    @test get_prn(ts, :gal, 12) == 12
end

@testset "add_satellite!(ts, acqs::Vector) — mixed-system batch auto-routes" begin
    ts = TrackState(; signals = (
        gps = (GPSL1CA(),),
        gal = (GalileoE1B(),),
    ))
    acqs = [
        _make_acq(GPSL1CA(),    11, 1.5, 50.0Hz),
        _make_acq(GalileoE1B(), 12, 2.5, 60.0Hz),
        _make_acq(GPSL1CA(),    13, 3.5, 70.0Hz),
    ]
    add_satellite!(ts, acqs)  # no group= — routes per-entry
    @test get_prn(ts, :gps, 11) == 11
    @test get_prn(ts, :gps, 13) == 13
    @test get_prn(ts, :gal, 12) == 12
end

@testset "add_satellite!(ts, acq) — unknown signal type errors with auto-routing" begin
    ts = TrackState(; signals = (gps = (GPSL1CA(),),))
    bad = _make_acq(GalileoE1B(), 5, 10.0, 100.0Hz)
    @test_throws ArgumentError add_satellite!(ts, bad)
end

@testset "add_satellite!(ts, acq) rejects wrong-signal acq" begin
    ts = TrackState(; signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(),)))
    # Galileo acq into the GPS group is the canonical footgun.
    bad = _make_acq(GalileoE1B(), 5, 10.0, 100.0Hz)
    @test_throws ArgumentError add_satellite!(ts, bad; group = :gps)
end

@testset "add_satellite!(ts, acq) requires longest-code signal in mixed group" begin
    # Group containing GPS L1 C/A (1023 chips) + GPS L1C-P (10230 chips):
    # the longest primary is L1C-P, so handing over an L1CA acquisition
    # would alias its 1023-chip code-phase inside L1C-P's 10230-chip
    # primary period. Reject.
    ts = TrackState(; signal = GPSL1C_P())  # placeholder single-group reset
    ts = TrackState(; signals = (mix = (GPSL1C_P(), GPSL1CA()),))
    short_acq = _make_acq(GPSL1CA(), 7, 100.0, 50.0Hz)
    @test_throws ArgumentError add_satellite!(ts, short_acq; group = :mix)
    # An L1C-P acquisition is accepted.
    long_acq = _make_acq(GPSL1C_P(), 7, 100.0, 50.0Hz)
    add_satellite!(ts, long_acq; group = :mix)
    @test get_prn(ts, :mix, 7) == 7
end

@testset "TrackState(acq) — single-acq convenience" begin
    acq = _make_acq(GPSL1CA(), 7, 524.6, 100.0Hz)
    ts = TrackState(acq)
    @test length(get_sat_states(ts, :default)) == 1
    @test get_prn(ts, :default, 7) == 7
    @test get_carrier_doppler(ts, :default, 7) == 100.0Hz
end

@testset "TrackState(acqs::Vector) — batch single-signal convenience" begin
    acqs = [
        _make_acq(GPSL1CA(), 1, 0.0, 100.0Hz),
        _make_acq(GPSL1CA(), 2, 10.0, 200.0Hz),
        _make_acq(GPSL1CA(), 3, 20.0, 300.0Hz),
    ]
    ts = TrackState(acqs)
    @test length(get_sat_states(ts, :default)) == 3
    @test get_carrier_doppler(ts, :default, 2) == 200.0Hz
end

@testset "TrackState(empty acq vector) errors" begin
    @test_throws ArgumentError TrackState(AcquisitionResults[])
end

@testset "TrackState(mixed acq vector) without `signals=` errors" begin
    acqs = [
        _make_acq(GPSL1CA(), 1, 0.0, 100.0Hz),
        _make_acq(GalileoE1B(), 2, 10.0, 200.0Hz),
    ]
    @test_throws ArgumentError TrackState(acqs)
end

@testset "TrackState(acqs; signals = ...) — multi-group routing" begin
    acqs = [
        _make_acq(GPSL1CA(), 11, 1.5, 50.0Hz),
        _make_acq(GalileoE1B(), 12, 2.5, 60.0Hz),
        _make_acq(GPSL1CA(), 13, 3.5, 70.0Hz),
    ]
    ts = TrackState(acqs;
        signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(),)),
    )
    @test length(get_sat_states(ts, :gps)) == 2
    @test length(get_sat_states(ts, :gal)) == 1
    @test get_carrier_doppler(ts, :gps, 11) == 50.0Hz
    @test get_carrier_doppler(ts, :gal, 12) == 60.0Hz
end

@testset "TrackState(acqs; signals = ...) errors on unmatched system" begin
    acqs = [_make_acq(GalileoE1B(), 5, 10.0, 100.0Hz)]
    @test_throws ArgumentError TrackState(acqs;
        signals = (gps = (GPSL1CA(),),),
    )
end

@testset "TrackState(acq; doppler_estimator) forwards kwargs" begin
    custom = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    acq = _make_acq(GPSL1CA(), 9, 5.0, 80.0Hz)
    ts = TrackState(acq; doppler_estimator = custom)
    de_state = get_sat_state(ts, :default, 9).doppler_estimator_state
    @test de_state isa SatConventionalPLLAndDLL
    @test de_state.carrier_loop_filter_bandwidth == 22.0Hz
end

end
