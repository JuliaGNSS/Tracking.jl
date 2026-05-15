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

@testset "TrackState(; signals = ...) — single-capability shortcut" begin
    track_state = TrackState(; signal = GPSL1CA())
    # Implicit `:default` capability.
    @test isempty(get_sat_states(track_state, :default))
end

@testset "TrackState(; signals = ...) — multi-capability NamedTuple" begin
    track_state = TrackState(;
        signals = (
            legacy_gps = (GPSL1CA(),),
            galileo = (GalileoE1B(),),
        ),
    )
    @test isempty(get_sat_states(track_state, :legacy_gps))
    @test isempty(get_sat_states(track_state, :galileo))
end

@testset "add_satellite! — single-capability shortcut (no `capability=`)" begin
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

@testset "add_satellite! — multi-capability with explicit `capability=`" begin
    track_state = TrackState(;
        signals = (
            legacy = (GPSL1CA(),),
            galileo = (GalileoE1B(),),
        ),
    )
    add_satellite!(track_state;
        prn = 5, capability = :legacy, carrier_doppler = 500.0Hz,
    )
    add_satellite!(track_state;
        prn = 11, capability = :galileo, carrier_doppler = 2000.0Hz,
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

@testset "add_satellite!(track_state, capability, sat) — escape hatch" begin
    track_state = TrackState(; signal = GPSL1CA())
    sat = TrackedSat(GPSL1CA(), 9, 5.0, 300.0Hz;
        doppler_estimator = ConventionalAssistedPLLAndDLL())
    add_satellite!(track_state, :default, sat)
    @test get_prn(track_state, :default, 9) == 9
    @test get_carrier_doppler(track_state, :default, 9) == 300.0Hz
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
