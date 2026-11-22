module PositionVelocityTimeTest

using Test: @test, @testset
using GNSSSignals: CBOC
using Unitful: Hz
using GNSSSignals: GPSL1CA, GalileoE1B
using GNSSDecoder: GNSSDecoderState
using Dictionaries: Dictionary
import PositionVelocityTime
import Tracking
using Tracking:
    TrackedSat, TrackState, get_code_phase, get_carrier_doppler, get_signal_observations

@testset "PositionVelocityTime: SignalObservation from a tracked satellite" begin
    gpsl1 = GPSL1CA()
    decoder = GNSSDecoderState(gpsl1, 1)
    tracked_sat = TrackedSat(gpsl1, 1, 123.0, 1234.0Hz)

    observation = PositionVelocityTime.SignalObservation(decoder, gpsl1, tracked_sat)

    @test observation isa PositionVelocityTime.SignalObservation
    @test observation.decoder === decoder
    @test observation.signal === gpsl1
    @test observation.code_phase == get_code_phase(tracked_sat)
    @test observation.carrier_doppler == get_carrier_doppler(tracked_sat)
end

@testset "get_signal_observations from a TrackState" begin
    track_state =
        TrackState(; signals = (gps = (GPSL1CA(),), galileo = (GalileoE1B(CBOC),)))
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
        galileo = Dictionary([2], [GNSSDecoderState(GalileoE1B(CBOC), 2)]),
    )

    observations = get_signal_observations(track_state, decoders)
    @test keys(observations) == keys(track_state.groups)
    @test collect(keys(observations.gps)) == [7]
    @test collect(keys(observations.galileo)) == [2]
    # The ranging signal defaults to the group's first.
    @test observations.gps[7].signal === first(track_state.groups.gps.signals)
    @test observations.gps[7].decoder === decoders.gps[7]
    tracked = Tracking.get_tracked_sat(track_state, :gps, 7)
    @test observations.gps[7].code_phase == get_code_phase(tracked)
    @test observations.gps[7].carrier_doppler == get_carrier_doppler(tracked)

    # Concretely typed per group, as `calc_pvt` requires, and inferable from
    # `track_state`'s type (`Base.infer_return_type` needs Julia 1.11).
    @test isconcretetype(eltype(observations.gps))
    @test only(
        Base.return_types(
            get_signal_observations,
            Tuple{typeof(track_state),typeof(decoders)},
        ),
    ) === typeof(observations)

    # The `signal` keyword ranges on a signal other than the driver.
    one_group = get_signal_observations(track_state.groups.gps, decoders.gps)
    @test collect(keys(one_group)) == [7]
    ranging_on = GPSL1CA()
    @test get_signal_observations(
        track_state.groups.gps,
        decoders.gps;
        signal = ranging_on,
    )[7].signal === ranging_on
    # No decoders at all: no observations.
    @test isempty(
        get_signal_observations(track_state.groups.gps, Dictionary{Int,GNSSDecoderState}()),
    )
end

end
