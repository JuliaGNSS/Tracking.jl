# Shared test helper. The including module must have `Tracking` in scope.

# A track state around one satellite whose signals repeat. A signal group rejects a
# repeated signal, so the group is built through its unchecked constructor: only
# the tile-share kernel tests, which check that every copy gets the single-signal
# result, may build one.
function _repeated_signal_track_state(sat::Tracking.TrackedSat, doppler_estimator)
    sats = Tracking.to_dictionary(sat)
    signals = map(s -> s.signal, sat.signals)
    band = Tracking.get_band(first(signals))
    num_ants = Tracking.NumAnts(Tracking.get_num_ants(sat))
    group =
        Tracking.SignalGroup{typeof(band),typeof(sats),typeof(signals),typeof(num_ants)}(
            band,
            sats,
            signals,
            num_ants,
        )
    groups = (default = group,)
    Tracking.TrackState(
        groups,
        doppler_estimator,
        Tracking._resolve_noise_estimators(nothing, groups),
    )
end
