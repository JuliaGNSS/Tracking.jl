# The hand-over to PositionVelocityTime.jl: tracked satellites paired with their
# decoded navigation data, as the measurements `PositionVelocityTime.calc_pvt`
# solves. Not a PositionVelocityTime extension: this package already depends on
# it (through TrackingLoops), so that would close a dependency cycle. Both
# packages define a `SignalGroup`, so it is always qualified below.

"""
    PositionVelocityTime.SatelliteState(decoder::GNSSDecoderState, signal, sat_state::TrackedSat)

One tracked satellite's code phase, carrier Doppler and carrier phase (from
[`get_sat_state`](@ref)) paired with its decoded navigation data, ranging on `signal`.
"""
function PositionVelocityTime.SatelliteState(
    decoder::GNSSDecoderState,
    signal::AbstractGNSSSignal,
    sat_state::TrackedSat,
)
    PositionVelocityTime.SatelliteState(
        decoder,
        signal,
        get_code_phase(sat_state),
        get_carrier_doppler(sat_state),
        get_carrier_phase(sat_state),
    )
end

"""
    PositionVelocityTime.SignalGroup(group::Tracking.SignalGroup, decoders;
                                     signal = first(group.signals))

One `Tracking.SignalGroup`'s measurements for an epoch, each satellite paired with
its decoded navigation data.

`decoders` maps the group's PRNs to `GNSSDecoderState`s (a `Dictionary`, or any
`getindex`/`haskey` mapping); the result keeps those keys. Satellites without a
decoder are skipped, as tracking starts long before the message is complete.

`signal` is the ranging signal, by default the group's first (the estimator
driver, and GNSSReceiver's `RANGING_SIGNAL_INDEX = 1`).
"""
function PositionVelocityTime.SignalGroup(
    group::SignalGroup,
    decoders;
    signal::AbstractGNSSSignal = first(group.signals),
)
    # `map` over the kept keys gives a concretely typed vector of states.
    kept = [key for key in keys(group.satellites) if haskey(decoders, key)]
    states = map(
        key -> PositionVelocityTime.SatelliteState(
            decoders[key],
            signal,
            group.satellites[key],
        ),
        kept,
    )
    PositionVelocityTime.SignalGroup(signal, Dictionary(kept, states))
end

"""
    PositionVelocityTime.signal_groups(track_state::TrackState, decoders)
        -> NamedTuple of PositionVelocityTime.SignalGroup

One epoch's measurements from a whole [`TrackState`](@ref) for
`PositionVelocityTime.calc_pvt`. `decoders` has the group keys of
`track_state.groups`, each a PRN → `GNSSDecoderState` mapping; both levels of key
are kept. Each group ranges on its first signal (see
`PositionVelocityTime.SignalGroup(::Tracking.SignalGroup, decoders)`). Inferable.
"""
function PositionVelocityTime.signal_groups(track_state::TrackState, decoders::NamedTuple)
    map(track_state.groups, decoders) do group, group_decoders
        PositionVelocityTime.SignalGroup(group, group_decoders)
    end
end
