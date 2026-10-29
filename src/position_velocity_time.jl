# The hand-over to PositionVelocityTime.jl: tracked satellites paired with their
# decoded navigation data, as the observations `PositionVelocityTime.calc_pvt`
# solves. Not a PositionVelocityTime extension: this package already depends on
# it (through TrackingLoops), so that would close a dependency cycle.

"""
    PositionVelocityTime.SignalObservation(decoder::GNSSDecoderState, signal, sat_state::TrackedSat)

One tracked satellite's code phase and carrier Doppler (from [`get_sat_state`](@ref))
paired with its decoded navigation data, ranging on `signal`.
"""
function PositionVelocityTime.SignalObservation(
    decoder::GNSSDecoderState,
    signal::AbstractGNSSSignal,
    sat_state::TrackedSat,
)
    PositionVelocityTime.SignalObservation(
        decoder,
        signal,
        get_code_phase(sat_state),
        get_carrier_doppler(sat_state),
    )
end

"""
    get_signal_observations(group::SignalGroup, decoders; signal = first(group.signals))
        -> Dictionary of PositionVelocityTime.SignalObservation

One [`SignalGroup`](@ref)'s observations for an epoch, each satellite paired with
its decoded navigation data.

`decoders` maps the group's PRNs to `GNSSDecoderState`s (a `Dictionary`, or any
`getindex`/`haskey` mapping); the result keeps those keys. Satellites without a
decoder are skipped, as tracking starts long before the message is complete.

`signal` is the ranging signal, by default the group's first (the loop's driver
signal, and GNSSReceiver's `RANGING_SIGNAL_INDEX = 1`).
"""
function get_signal_observations(
    group::SignalGroup,
    decoders;
    signal::AbstractGNSSSignal = first(group.signals),
)
    # `map` over the kept keys gives a concretely typed vector of observations.
    kept = [key for key in keys(group.satellites) if haskey(decoders, key)]
    observations = map(
        key -> PositionVelocityTime.SignalObservation(
            decoders[key],
            signal,
            group.satellites[key],
        ),
        kept,
    )
    Dictionary(kept, observations)
end

"""
    get_signal_observations(track_state::TrackState, decoders::NamedTuple)
        -> NamedTuple of Dictionary of PositionVelocityTime.SignalObservation

One epoch's observations from a whole [`TrackState`](@ref) for
`PositionVelocityTime.calc_pvt`. `decoders` has the group keys of
`track_state.groups`, each a PRN → `GNSSDecoderState` mapping; both levels of key
are kept. Each group ranges on its first signal (see
[`get_signal_observations(::SignalGroup, ::Any)`](@ref)). Inferable.
"""
function get_signal_observations(track_state::TrackState, decoders::NamedTuple)
    map(track_state.groups, decoders) do group, group_decoders
        get_signal_observations(group, group_decoders)
    end
end
