module PositionVelocityTimeTest

using Test: @test, @testset
using Unitful: Hz
using GNSSSignals: GPSL1CA, GalileoE1B
using GNSSDecoder: GNSSDecoderState
using Dictionaries: Dictionary
import PositionVelocityTime
import Tracking
using Tracking:
    TrackedSat, TrackState, get_code_phase, get_carrier_doppler, get_carrier_phase

@testset "PositionVelocityTime: SatelliteState from a tracked satellite" begin
    gpsl1 = GPSL1CA()
    decoder = GNSSDecoderState(gpsl1, 1)
    sat_state = TrackedSat(gpsl1, 1, 123.0, 1234.0Hz; carrier_phase = 0.25)

    state = PositionVelocityTime.SatelliteState(decoder, gpsl1, sat_state)

    @test state isa PositionVelocityTime.SatelliteState
    @test state.decoder === decoder
    @test state.signal === gpsl1
    @test state.code_phase == get_code_phase(sat_state)
    @test state.carrier_doppler == get_carrier_doppler(sat_state)
    @test state.carrier_phase == get_carrier_phase(sat_state)
end

@testset "PositionVelocityTime: signal groups from a TrackState" begin
    track_state = TrackState(; signals = (gps = (GPSL1CA(),), galileo = (GalileoE1B(),)))
    for (group_key, prn, code_phase) in
        ((:gps, 7, 12.0), (:gps, 8, 13.0), (:galileo, 2, 5.0))
        track_state = Tracking.add_satellite!(
            track_state;
            group_key,
            prn,
            carrier_doppler = 100.0Hz,
            code_phase,
        )
    end
    # PRN 8 has no decoder yet, as in a satellite's first minute: skipped.
    decoders = (
        gps = Dictionary([7], [GNSSDecoderState(GPSL1CA(), 7)]),
        galileo = Dictionary([2], [GNSSDecoderState(GalileoE1B(), 2)]),
    )

    groups = PositionVelocityTime.signal_groups(track_state, decoders)
    @test keys(groups) == keys(track_state.groups)
    @test collect(keys(groups.gps.satellites)) == [7]
    @test collect(keys(groups.galileo.satellites)) == [2]
    # The ranging signal defaults to the group's first.
    @test groups.gps.signal === first(track_state.groups.gps.signals)
    @test groups.gps.satellites[7].signal === groups.gps.signal
    @test groups.gps.satellites[7].decoder === decoders.gps[7]
    tracked = Tracking.get_sat_state(track_state, :gps, 7)
    @test groups.gps.satellites[7].code_phase == get_code_phase(tracked)
    @test groups.gps.satellites[7].carrier_doppler == get_carrier_doppler(tracked)

    # Concretely typed per group and inferable from `track_state`'s type
    # (`Base.infer_return_type` needs Julia 1.11).
    @test isconcretetype(eltype(groups.gps.satellites))
    @test only(
        Base.return_types(
            PositionVelocityTime.signal_groups,
            Tuple{typeof(track_state),typeof(decoders)},
        ),
    ) === typeof(groups)

    # The `signal` keyword ranges on a signal other than the driver.
    one_group = PositionVelocityTime.SignalGroup(track_state.groups.gps, decoders.gps)
    @test collect(keys(one_group.satellites)) == [7]
    @test one_group.signal === groups.gps.signal
    ranging_on = GPSL1CA()
    @test PositionVelocityTime.SignalGroup(
        track_state.groups.gps,
        decoders.gps;
        signal = ranging_on,
    ).signal === ranging_on
    # No decoders at all: an empty group.
    @test isempty(
        PositionVelocityTime.SignalGroup(
            track_state.groups.gps,
            Dictionary{Int,GNSSDecoderState}(),
        ).satellites,
    )
end

end
