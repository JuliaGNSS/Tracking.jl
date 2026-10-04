module TrackingAcquisitionExt

using Acquisition: AcquisitionResults
using GNSSSignals: AbstractGNSSSignal, get_code_length, get_signal_id, get_signal_name
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
    routed as by [`add_satellite!`](@ref); an acq matching no group or more
    than one group errors.

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
        # Compare by signal id, not `typeof`: the code-matrix type parameter
        # is not part of a signal's identity (see `_acq_signal_matches`).
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

# Does `system` qualify as a handoff signal for `sig_tuple`, i.e. is it one of
# the signals tied at the longest primary code length? Compared by signal id, not
# `typeof`: the code-matrix type parameter is not part of a signal's identity.
@inline function _acq_signal_matches(sig_tuple::Tuple, system::AbstractGNSSSignal)
    longest_len = get_code_length(_longest_code_signal(sig_tuple))
    sys_id = get_signal_id(system)
    any(s -> get_signal_id(s) === sys_id && get_code_length(s) == longest_len, sig_tuple)
end

# Key of the unique group `acq.system` qualifies for; errors on no match and on
# an ambiguous match rather than silently picking the first group.
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

Add (or replace) a satellite in `track_state` from an acquisition result.
The `prn`, `code_phase`, and `carrier_doppler` are read off `acq`; the
remaining tracking state (correlator, post-corr filter, doppler-estimator
state) is initialized to the group's defaults.

`acq.system` must be a signal with the **longest code** in the
group's signal tuple — its code phase is the only unambiguous one.
For a group tracking `(GPS L1C-P, GPS L1C-D, GPS L1 C/A)`, hand over a
GPS L1C-P or L1C-D acquisition (both 10230 chips); L1CA's 1023-chip
period would alias inside the 10230-chip primary.

With `group = nothing` (the default) the group is inferred by this rule;
an explicit `group =` bypasses the inference and asserts the match.
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

# Shared body of the acq pair above: resolve the group and build its default sat.
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

Batch variant: add each entry of `acqs` to `track_state`, validated and (with
`group = nothing`) routed per entry like the single-acq form, so a vector mixing
constellations lands in the right groups. An explicit `group =` applies to every
entry. Keep using the returned track state; it carries the estimator after
all [`Tracking.update_estimator_on_handoff`](@ref) updates.
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

# Check that `acq.system` is a longest-code signal of the group (see the
# `add_satellite!` docstring above).
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

# The signal with the longest primary code (secondary codes are resolved by
# tracking itself); folds at compile time. What matters is the code *period*;
# comparing chips is equivalent only because a group shares one chip rate (see
# `_validate_signal_group` in sat_state.jl). If that is relaxed (issue #151),
# compare `get_code_length(s) / get_code_frequency(s)` instead.
@inline _longest_code_signal(t::Tuple{AbstractGNSSSignal}) = only(t)
@inline function _longest_code_signal(
    t::Tuple{AbstractGNSSSignal,AbstractGNSSSignal,Vararg{AbstractGNSSSignal}},
)
    head = first(t)
    rest_longest = _longest_code_signal(Base.tail(t))
    get_code_length(head) >= get_code_length(rest_longest) ? head : rest_longest
end

end
