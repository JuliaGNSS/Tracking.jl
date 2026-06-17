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
ts = TrackState(acq)
```
"""
function Tracking.TrackState(acq::AcquisitionResults; kwargs...)
    ts = Tracking.TrackState(; signal = acq.system, kwargs...)
    Tracking.add_satellite!(ts, acq)
end

"""
    TrackState(acqs::AbstractVector{<:AcquisitionResults}; signals = nothing, kwargs...)

Convenience constructor: build a [`TrackState`](@ref) pre-populated from
a batch of acquisition results.

  - **Default (single-group)**: omit `signals`. All entries must share the
    same `acq.system`; that shared signal becomes the `:default` group's
    signal. Errors on an empty vector (no signal to infer).
  - **Multi-group**: pass `signals = (group_a = (...), group_b = (...))` as
    in the regular [`TrackState`](@ref) constructor. Each acquisition is
    routed to the group whose longest-primary-code signal matches
    `acq.system` (signals tied at the maximum code length all qualify).
    Errors on an acq whose system doesn't match any group, and on an acq
    that matches more than one group (use `add_satellite!` with an
    explicit `group =` instead).

Other `kwargs` (`doppler_estimator`, `num_ants`, ...) forward to the
regular `TrackState` constructor.

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
        for a in acqs
            typeof(a.system) === typeof(sig) || throw(
                ArgumentError(
                    string(
                        "All acquisition results must share the same `system`. ",
                        "Got `",
                        typeof(sig),
                        "` and `",
                        typeof(a.system),
                        "`. ",
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

# Does `system` qualify as a handoff signal for a group declaring
# `sig_tuple`? It must appear in the tuple with a primary code length
# equal to the group's longest — the code-phase scaling is unambiguous
# for any signal tied at the maximum length (e.g. GPS L1C-D and L1C-P
# are both 10230 chips, so either acquisition is a valid handoff).
@inline function _acq_signal_matches(sig_tuple::Tuple, system::AbstractGNSSSignal)
    longest_len = get_code_length(_longest_code_signal(sig_tuple))
    any(s -> typeof(s) === typeof(system) && get_code_length(s) == longest_len, sig_tuple)
end

# Walk the TrackState's groups and return the key of the one whose
# longest-code signal matches `acq.system` (code-length ties all
# qualify). Errors if no group matches, and also if more than one group
# matches — silently routing to the first declared group would be a
# footgun when two groups share a signature.
@inline function _find_group_for_acq(track_state::TrackState, acq::AcquisitionResults)
    keys_tuple = keys(track_state.groups)
    acq_type = typeof(acq.system)
    matches = filter(
        k -> _acq_signal_matches(track_state.groups[k].signals, acq.system),
        keys_tuple,
    )
    isempty(matches) && throw(
        ArgumentError(
            string(
                "No group's longest-primary-code signal matches `acq.system::",
                acq_type,
                "` (PRN ",
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
                "Acquisition routing is ambiguous: `acq.system::",
                acq_type,
                "` (PRN ",
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

Add (or replace) a satellite in `track_state` from an acquisition result.
The `prn`, `code_phase`, and `carrier_doppler` are read off `acq`; the
remaining tracking state (correlator, post-corr filter, doppler-estimator
state) is initialized to the group's defaults.

`acq.system` must be a signal with the **longest code** in the
group's signal tuple — its code-phase scaling is the only one that's
unambiguous when the group tracks multiple signals on shared chips.
For a group tracking `(GPS L1C-P, GPS L1C-D, GPS L1 C/A)`, hand over a
GPS L1C-P or L1C-D acquisition (both 10230 chips); L1CA's 1023-chip
period would alias inside the 10230-chip primary.

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

# Shared body of the mutable/immutable acq pair above: resolve (or assert)
# the target group for `acq.system`, then build the group's
# default-correlator TrackedSat from the acquisition handoff values.
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

Batch variant: add each entry of `acqs` to `track_state`. Same validation
as the single-acq form runs per entry. With `group = nothing` (the
default) each acq is routed to the matching group by signal type, so a
mixed `Vector{AcquisitionResults}` from multiple constellations lands in
the right place. Passing an explicit `group =` keyword applies that
group to every entry (use the single-acq form per entry if your vector
mixes groups). Returns the track state carrying the estimator after all
per-entry [`Tracking.update_estimator_on_handoff`](@ref) updates — keep
using the return value.
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

Immutable batch variant of [`add_satellite!`](@ref). Returns a new
[`TrackState`](@ref) with every entry of `acqs` added; the input is left
unchanged. Routing and validation match the mutable batch form.
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

# Check that `acq.system` is a longest-code signal of the group's
# signal tuple (ties at the maximum code length all qualify). The
# longest signals are what define the group's `max_code_length` wrap
# point — feeding code_phase from a shorter signal would alias into the
# wrong place inside the longer signal's primary code period.
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
                "Acquisition signal type does not match a longest-code signal ",
                "in group `:",
                group,
                "`. ",
                "Got `acq.system::",
                typeof(acq.system),
                "` ",
                "but the group's longest-code signal is `",
                typeof(longest),
                "`. ",
                "Hand over an acquisition for a signal with the longest code ",
                "length — its code phase is the only one that is unambiguous ",
                "when the group tracks multiple signals.",
            ),
        ),
    )
end

# Recursive walk to find the signal with the largest primary code
# length in a signal tuple. Secondary codes are negotiated by tracking
# itself, so the wrap point that matters for handoff is the primary.
# Type-stable on concrete tuple types, folds at compile time.
#
# The physically meaningful criterion is the longest code *period* (in
# time): a shorter-period code's phase is ambiguous about which
# repetition it sits in within a longer-period code. Comparing raw code
# *length* (chips) is equivalent only because every signal in a group
# shares one chip rate (enforced by the SignalGroup constructor — see
# `_validate_signal_group` in sat_state.jl), so period ∝ length. If that
# equal-chip-rate invariant is ever relaxed (see issue #151), switch this
# to compare `get_code_length(s) / get_code_frequency(s)`.
@inline _longest_code_signal(t::Tuple{AbstractGNSSSignal}) = only(t)
@inline function _longest_code_signal(
    t::Tuple{AbstractGNSSSignal,AbstractGNSSSignal,Vararg{AbstractGNSSSignal}},
)
    head = first(t)
    rest_longest = _longest_code_signal(Base.tail(t))
    get_code_length(head) >= get_code_length(rest_longest) ? head : rest_longest
end

end
