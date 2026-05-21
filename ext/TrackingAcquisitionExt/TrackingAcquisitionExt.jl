module TrackingAcquisitionExt

using Acquisition: AcquisitionResults
using GNSSSignals: AbstractGNSSSignal, get_code_length
using Tracking: Tracking, TrackState, TrackedSat

function Tracking.TrackedSat(acq::AcquisitionResults; args...)
    TrackedSat(acq.system, acq.prn, acq.code_phase, acq.carrier_doppler; args...)
end

"""
    TrackState(acq::AcquisitionResults; doppler_estimator, num_ants)

Convenience constructor: build a single-group `:default` [`TrackState`](@ref)
from one acquisition result. The group's signal is `acq.system`; the
satellite is pre-populated with `prn`, `code_phase`, and `carrier_doppler`
from `acq`. Loop-filter bandwidths default to the recommended values for
`acq.system` and can be overridden via `doppler_estimator`.

```julia
acq = acquire(GPSL1CA(), data, fs, 7)
ts  = TrackState(acq)
```
"""
function Tracking.TrackState(acq::AcquisitionResults; kwargs...)
    ts = Tracking.TrackState(; signal = acq.system, kwargs...)
    Tracking.add_satellite!(ts, acq)
    ts
end

"""
    TrackState(acqs::AbstractVector{<:AcquisitionResults}; signals = nothing, kwargs...)

Convenience constructor: build a [`TrackState`](@ref) pre-populated from
a batch of acquisition results.

* **Default (single-group)**: omit `signals`. All entries must share the
  same `acq.system`; that shared signal becomes the `:default` group's
  signal. Errors on an empty vector (no signal to infer).
* **Multi-group**: pass `signals = (group_a = (...), group_b = (...))` as
  in the regular [`TrackState`](@ref) constructor. Each acquisition is
  routed to the group whose longest-primary-code signal matches
  `acq.system`. Errors on an acq whose system doesn't match any group.

Other `kwargs` (`doppler_estimator`, `num_ants`, ...) forward to the
regular `TrackState` constructor.

```julia
# Single-signal
acqs = acquire(GPSL1CA(), data, fs, 1:32)
ts   = TrackState(filter(is_detected, acqs))

# Multi-signal
acqs = vcat(
    acquire(GPSL1CA(), data, fs, 1:32),
    acquire(GalileoE1B(), data, fs, 1:36),
)
ts = TrackState(filter(is_detected, acqs);
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
        isempty(acqs) && throw(ArgumentError(
            "Cannot infer the group's signal from an empty acquisition-results " *
            "vector. Pass `signals = (...)` to declare signals explicitly, " *
            "or pass at least one acquisition.",
        ))
        sig = first(acqs).system
        for a in acqs
            typeof(a.system) === typeof(sig) || throw(ArgumentError(string(
                "All acquisition results must share the same `system`. ",
                "Got `", typeof(sig), "` and `", typeof(a.system), "`. ",
                "Pass `signals = (...)` to track multiple signals.",
            )))
        end
        ts = Tracking.TrackState(; signal = sig, kwargs...)
        Tracking.add_satellite!(ts, acqs)
        return ts
    end
    ts = Tracking.TrackState(; signals, kwargs...)
    for acq in acqs
        group = _find_group_for_acq(ts, acq)
        Tracking.add_satellite!(ts, acq; group)
    end
    ts
end

# Walk the TrackState's groups and return the key of the one whose
# longest-primary-code signal matches `acq.system`. Errors if no group
# matches.
@inline function _find_group_for_acq(track_state::TrackState, acq::AcquisitionResults)
    keys_tuple = keys(track_state.groups)
    acq_type = typeof(acq.system)
    for k in keys_tuple
        longest = _longest_code_signal(track_state.groups[k].signals)
        typeof(longest) === acq_type && return k
    end
    throw(ArgumentError(string(
        "No group's longest-primary-code signal matches `acq.system::",
        acq_type, "` (PRN ", acq.prn, "). Declared groups: ", keys_tuple, ".",
    )))
end

"""
    add_satellite!(track_state, acq::AcquisitionResults; group = nothing)

Add (or replace) a satellite in `track_state` from an acquisition result.
The `prn`, `code_phase`, and `carrier_doppler` are read off `acq`; the
remaining tracking state (correlator, post-corr filter, doppler-estimator
state) is initialized to the group's defaults.

`acq.system` must be the signal with the **longest code** in the
group's signal tuple — its code-phase scaling is the only one that's
unambiguous when the group tracks multiple signals on shared chips.
For a group tracking `(GPS L1C-P, GPS L1C-D, GPS L1 C/A)`, hand over a
GPS L1C-P acquisition; L1CA's 1023-chip period would alias inside
L1C-P's 10230-chip primary.

If `group` is left as `nothing` (the default), the routing is inferred
by matching `acq.system` against each group's longest-primary-code
signal — the same rule the [`TrackState(acqs; signals = ...)`](@ref)
constructor uses. Passing an explicit `group =` keyword bypasses the
inference and asserts the match against that specific group.
"""
function Tracking.add_satellite!(
    track_state::TrackState,
    acq::AcquisitionResults;
    group::Union{Symbol,Nothing} = nothing,
)
    resolved = isnothing(group) ? _find_group_for_acq(track_state, acq) : group
    _assert_acq_matches_group(track_state, resolved, acq)
    Tracking.add_satellite!(
        track_state;
        prn = acq.prn,
        group = resolved,
        code_phase = acq.code_phase,
        carrier_doppler = acq.carrier_doppler,
    )
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
    resolved = isnothing(group) ? _find_group_for_acq(track_state, acq) : group
    _assert_acq_matches_group(track_state, resolved, acq)
    Tracking.add_satellite(
        track_state;
        prn = acq.prn,
        group = resolved,
        code_phase = acq.code_phase,
        carrier_doppler = acq.carrier_doppler,
    )
end

"""
    add_satellite!(track_state, acqs::AbstractVector{<:AcquisitionResults}; group = nothing)

Batch variant: add each entry of `acqs` to `track_state`. Same validation
as the single-acq form runs per entry. With `group = nothing` (the
default) each acq is routed to the matching group by signal type, so a
mixed `Vector{AcquisitionResults}` from multiple constellations lands in
the right place. Passing an explicit `group =` keyword applies that
group to every entry (use the single-acq form per entry if your vector
mixes groups). Returns `track_state`.
"""
function Tracking.add_satellite!(
    track_state::TrackState,
    acqs::AbstractVector{<:AcquisitionResults};
    group::Union{Symbol,Nothing} = nothing,
)
    for acq in acqs
        Tracking.add_satellite!(track_state, acq; group)
    end
    track_state
end

# Check that `acq.system` is the longest-code signal in the group's
# signal tuple. The longest signal is what defines the group's
# `max_code_length` wrap point — feeding code_phase from a shorter
# signal would alias into the wrong place inside the longer signal's
# primary code period.
@inline function _assert_acq_matches_group(
    track_state::TrackState, group::Symbol, acq::AcquisitionResults,
)
    sig_tuple = track_state.groups[group].signals
    longest = _longest_code_signal(sig_tuple)
    typeof(acq.system) === typeof(longest) && return nothing
    throw(ArgumentError(string(
        "Acquisition signal type does not match the longest-code signal ",
        "in group `:", group, "`. ",
        "Got `acq.system::", typeof(acq.system), "` ",
        "but the group's longest-code signal is `", typeof(longest), "`. ",
        "Hand over an acquisition for the longest signal — its code phase ",
        "is the only one that is unambiguous when the group tracks ",
        "multiple signals.",
    )))
end

# Recursive walk to find the signal with the largest primary code
# length in a signal tuple. Secondary codes are negotiated by tracking
# itself, so the wrap point that matters for handoff is the primary.
# Type-stable on concrete tuple types, folds at compile time.
@inline _longest_code_signal(t::Tuple{AbstractGNSSSignal}) = only(t)
@inline function _longest_code_signal(t::Tuple{AbstractGNSSSignal,AbstractGNSSSignal,Vararg{AbstractGNSSSignal}})
    head = first(t)
    rest_longest = _longest_code_signal(Base.tail(t))
    get_code_length(head) >= get_code_length(rest_longest) ? head : rest_longest
end

end
