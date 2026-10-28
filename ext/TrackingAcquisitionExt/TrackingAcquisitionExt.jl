module TrackingAcquisitionExt

using Acquisition: AcquisitionResults
using GNSSSignals: AbstractGNSSSignal, get_code_length, get_signal_id, get_signal_name
using Tracking: Tracking, TrackState, TrackedSat

function Tracking.TrackedSat(acq::AcquisitionResults; args...)
    TrackedSat(acq.system, acq.prn, acq.code_phase, acq.carrier_doppler; args...)
end

"""
    TrackState(acq::AcquisitionResults; doppler_estimator, num_ants)

Build a single-group `:default` [`TrackState`](@ref) on `acq.system`, holding the
satellite of `acq` (its `prn`, `code_phase` and `carrier_doppler`).

```julia
acq = acquire(GPSL1CA(), data, fs, 7)
ts = TrackState(acq)
```
"""
function Tracking.TrackState(acq::AcquisitionResults; kwargs...)
    ts = Tracking.TrackState(; signal = acq.system, kwargs...)
    Tracking.add_satellite!(ts, acq)
end

"""
    TrackState(acqs::AbstractVector{<:AcquisitionResults}; signals = nothing, kwargs...)

Build a [`TrackState`](@ref) holding the satellites of a batch of acquisition results.

  - **Single group**: omit `signals`. All entries must share one `acq.system`, which
    becomes the `:default` group's signal; an empty vector errors.
  - **Multi-group**: pass `signals` as to the regular [`TrackState`](@ref)
    constructor. Each acquisition goes to the group whose longest-primary-code
    signal (ties included) matches `acq.system`; it errors if no group or more than
    one matches (then use `add_satellite!` with `group =`).

Other `kwargs` (`doppler_estimator`, `num_ants`, ...) go to the regular constructor.

```julia
# Single-signal
acqs = acquire(GPSL1CA(), data, fs, 1:32)
ts = TrackState(filter(is_detected, acqs))

# Multi-signal
acqs = vcat(acquire(GPSL1CA(), data, fs, 1:32), acquire(GalileoE1B(), data, fs, 1:36))
ts = TrackState(
    filter(is_detected, acqs);
    signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(),)),
)
```
"""
function Tracking.TrackState(
    acqs::AbstractVector{<:AcquisitionResults};
    signals = nothing,
    kwargs...,
)
    if signals === nothing
        isempty(acqs) && throw(
            ArgumentError(
                "Cannot infer the group's signal from an empty acquisition-results " *
                "vector. Pass `signals = (...)` to declare signals explicitly, " *
                "or pass at least one acquisition.",
            ),
        )
        sig = first(acqs).system
        # By signal id, not `typeof` (see `_acq_signal_matches`).
        sig_id = get_signal_id(sig)
        for a in acqs
            get_signal_id(a.system) === sig_id || throw(
                ArgumentError(
                    string(
                        "All acquisition results must share the same `system`. ",
                        "Got ",
                        get_signal_name(sig),
                        " and ",
                        get_signal_name(a.system),
                        ". ",
                        "Pass `signals = (...)` to track multiple signals.",
                    ),
                ),
            )
        end
        ts = Tracking.TrackState(; signal = sig, kwargs...)
        return Tracking.add_satellite!(ts, acqs)
    end
    ts = Tracking.TrackState(; signals, kwargs...)
    for acq in acqs
        group = _find_group_for_acq(ts, acq)
        ts = Tracking.add_satellite!(ts, acq; group)
    end
    ts
end

# Whether `system` is in `sig_tuple` with the group's longest primary code length
# (ties qualify: GPS L1C-D and L1C-P are both 10230 chips). Compared by signal id,
# not `typeof`, since the code-matrix type parameter is not part of the identity.
@inline function _acq_signal_matches(sig_tuple::Tuple, system::AbstractGNSSSignal)
    longest_len = get_code_length(_longest_code_signal(sig_tuple))
    sys_id = get_signal_id(system)
    any(s -> get_signal_id(s) === sys_id && get_code_length(s) == longest_len, sig_tuple)
end

# The key of the one group matching `acq.system`; errors on none or several.
@inline function _find_group_for_acq(track_state::TrackState, acq::AcquisitionResults)
    keys_tuple = keys(track_state.groups)
    acq_name = get_signal_name(acq.system)
    matches = filter(
        k -> _acq_signal_matches(track_state.groups[k].signals, acq.system),
        keys_tuple,
    )
    isempty(matches) && throw(
        ArgumentError(
            string(
                "No group's longest-primary-code signal matches `acq.system` (",
                acq_name,
                ", PRN ",
                acq.prn,
                "). Declared groups: ",
                keys_tuple,
                ".",
            ),
        ),
    )
    length(matches) > 1 && throw(
        ArgumentError(
            string(
                "Acquisition routing is ambiguous: `acq.system` (",
                acq_name,
                ", PRN ",
                acq.prn,
                ") matches groups ",
                matches,
                ". Pass `group = <key>` explicitly to pick one.",
            ),
        ),
    )
    only(matches)
end

"""
    add_satellite!(track_state, acq::AcquisitionResults; group = nothing)

Add (or replace) a satellite in `track_state` from an acquisition result: `prn`,
`code_phase` and `carrier_doppler` from `acq`, the rest at the group's defaults.

`acq.system` must have the **longest code** in the group, the only unambiguous code
phase: for `(GPS L1C-P, GPS L1C-D, GPS L1 C/A)` hand over L1C-P or L1C-D (10230
chips); L1 C/A's 1023 chips would alias inside the 10230-chip primary.

With `group = nothing` the group is inferred as in
[`TrackState(acqs; signals = ...)`](@ref); an explicit `group` is checked instead.
"""
function Tracking.add_satellite!(
    track_state::TrackState,
    acq::AcquisitionResults;
    group::Union{Symbol,Nothing} = nothing,
)
    resolved, sat = _group_and_sat_for_acq(track_state, acq, group)
    Tracking.add_satellite!(track_state, resolved, sat)
end

"""
    add_satellite(track_state, acq::AcquisitionResults; group = nothing)

Immutable variant of [`add_satellite!`](@ref). Returns a new
[`TrackState`](@ref); the input is left unchanged.
"""
function Tracking.add_satellite(
    track_state::TrackState,
    acq::AcquisitionResults;
    group::Union{Symbol,Nothing} = nothing,
)
    resolved, sat = _group_and_sat_for_acq(track_state, acq, group)
    Tracking.add_satellite(track_state, resolved, sat)
end

# Shared by the two forms above: resolve or check the group, build the sat.
@inline function _group_and_sat_for_acq(
    track_state::TrackState,
    acq::AcquisitionResults,
    group::Union{Symbol,Nothing},
)
    resolved = isnothing(group) ? _find_group_for_acq(track_state, acq) : group
    _assert_acq_matches_group(track_state, resolved, acq)
    sat = Tracking._make_default_tracked_sat_for_group(
        track_state,
        resolved;
        prn = acq.prn,
        code_phase = acq.code_phase,
        carrier_doppler = acq.carrier_doppler,
    )
    resolved, sat
end

"""
    add_satellite!(track_state, acqs::AbstractVector{<:AcquisitionResults}; group = nothing)

Batch variant: add each entry of `acqs`, validated as in the single form. With
`group = nothing` each is routed to its group, so a mixed-constellation vector works;
an explicit `group` applies to every entry. Returns the track state carrying the
estimator after every [`Tracking.update_estimator_on_handoff`](@ref); keep using the
return value.
"""
function Tracking.add_satellite!(
    track_state::TrackState,
    acqs::AbstractVector{<:AcquisitionResults};
    group::Union{Symbol,Nothing} = nothing,
)
    for acq in acqs
        track_state = Tracking.add_satellite!(track_state, acq; group)
    end
    track_state
end

"""
    add_satellite(track_state, acqs::AbstractVector{<:AcquisitionResults}; group = nothing)

Immutable batch variant of [`add_satellite!`](@ref): returns a new
[`TrackState`](@ref) with every entry of `acqs` added, the input unchanged.
"""
function Tracking.add_satellite(
    track_state::TrackState,
    acqs::AbstractVector{<:AcquisitionResults};
    group::Union{Symbol,Nothing} = nothing,
)
    foldl(acqs; init = track_state) do ts, acq
        Tracking.add_satellite(ts, acq; group)
    end
end

# A shorter signal's code phase would alias inside the longest primary code, which
# sets the group's `max_code_wrap` wrap.
@inline function _assert_acq_matches_group(
    track_state::TrackState,
    group::Symbol,
    acq::AcquisitionResults,
)
    sig_tuple = track_state.groups[group].signals
    _acq_signal_matches(sig_tuple, acq.system) && return nothing
    longest = _longest_code_signal(sig_tuple)
    throw(
        ArgumentError(
            string(
                "Acquisition signal does not match a longest-code signal ",
                "in group `:",
                group,
                "`. ",
                "Got `acq.system` (",
                get_signal_name(acq.system),
                ") but the group's longest-code signal is ",
                get_signal_name(longest),
                ". ",
                "Hand over an acquisition for a signal with the longest code ",
                "length — its code phase is the only one that is unambiguous ",
                "when the group tracks multiple signals.",
            ),
        ),
    )
end

# The signal with the longest primary code (tracking finds the secondary itself);
# folds at compile time. What matters is the code *period*; chips stand in for it
# because a group shares one chip rate (`_validate_signal_group`). If that is
# relaxed (#151), compare `get_code_length(s) / get_chip_rate(s)`.
@inline _longest_code_signal(t::Tuple{AbstractGNSSSignal}) = only(t)
@inline function _longest_code_signal(
    t::Tuple{AbstractGNSSSignal,AbstractGNSSSignal,Vararg{AbstractGNSSSignal}},
)
    head = first(t)
    rest_longest = _longest_code_signal(Base.tail(t))
    get_code_length(head) >= get_code_length(rest_longest) ? head : rest_longest
end

end
