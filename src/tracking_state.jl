"""
$(SIGNATURES)

Construct a fresh `TrackState` from a declaration of which signals each group
tracks. A *group* is a set of satellites with the same signal tuple, hence the
same concrete `TrackedSat` type (type stability of the hot loop). The first signal
of each tuple is the driver signal: its integrations close the tracking loop,
which updates the satellite-shared carrier and code Doppler. `add_satellite!`
refers to a group by its key (`:modern_gps`, …).

```julia
track_state = TrackState(;
    signals = (
        legacy_gps = (GPSL1CA(),),
        modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
        galileo = (GalileoE1B(),),
    ),
)
```

For one group tracking one signal, use the singular `signal` keyword:

```julia
track_state = TrackState(; signal = GPSL1CA())
```

This is `TrackState(; signals = (default = (GPSL1CA(),),))`; `add_satellite!`
may then omit `group_key=`.

The per-group `TrackedSat` type (correlator bank, post-correlation filter and loop-state
types) is frozen at construction. For non-default correlator bank or post-correlation filter *types*,
build the `TrackedSat`s yourself and construct from them (e.g.
`TrackState(sats)`).

`noise_estimators` declares the per-signal noise sources
([`AbstractNoiseDensityEstimator`](@extref TrackingLoops.AbstractNoiseDensityEstimator)s keyed by
`GNSSSignals.get_signal_id`, e.g. `:GPSL1CA`). Left at `nothing`, every signal whose
C/N₀ estimator reads a noise density (see
[`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:TrackingLoops.AbstractCN0Estimator}}))
gets a [`CorrelatorNoiseDensityEstimator`](@extref TrackingLoops.CorrelatorNoiseDensityEstimator)
and every other signal none, so a state on
[`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator) despreads nothing. Pass a
NamedTuple to configure the window, or to declare a source you fill with
[`append_noise_observation!`](@ref) on a correlator-ingest path.

`sample_index_offsets` sets where each band's sample grid starts, keyed by `get_band_id`
(`NamedTuple` or `AbstractDict`; missing bands start at `0`). Pass the old state's
when rebuilding around a tracking loop that keeps state of its own; see the
[`TrackState`](@ref) struct. The other constructors take the same keyword.
"""
function TrackState(;
    signal::Maybe{AbstractGNSSSignal} = nothing,
    signals = nothing,
    tracking_loop::Maybe{AbstractTrackingLoop} = nothing,
    num_ants::NumAnts = NumAnts(1),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_index_offsets = nothing,
)
    if isnothing(signal) && isnothing(signals)
        throw(
            ArgumentError(
                "TrackState requires either `signal = <AbstractGNSSSignal>` or " *
                "`signals = (...)`. See the docstring for examples.",
            ),
        )
    end
    if !isnothing(signal) && !isnothing(signals)
        throw(
            ArgumentError(
                "Pass either `signal` (singular, one AbstractGNSSSignal) or " *
                "`signals` (plural, a NamedTuple of signal tuples) — not both.",
            ),
        )
    end
    sig_groups_nt =
        isnothing(signal) ? _normalize_signal_groups(signals) : (default = (signal,),)
    # Each satellite's loop bandwidth is seeded from its own group's driver
    # signal (`init_loop_state`), so groups need no shared compromise.
    tracking_loop = something(tracking_loop, _default_tracking_loop())
    # An entry is a bare signal tuple (takes `num_ants`) or a pre-built
    # `SignalGroup` (carries its own band and `num_ants`).
    groups = map(sig_groups_nt) do entry
        _normalize_group_entry(entry, tracking_loop, num_ants)
    end
    _validate_same_band_num_ants(groups)
    TrackState(
        groups,
        tracking_loop,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_index_offsets,
    )
end

# `nothing` provisions a `CorrelatorNoiseDensityEstimator` only for the signals that
# read a density, so nobody else pays for the despread.
@inline function _resolve_noise_estimators(
    noise_estimators::NamedTuple,
    groups::SignalGroups,
)
    _validate_noise_estimator_num_ants(noise_estimators, groups)
    noise_estimators
end
@inline function _resolve_noise_estimators(::Nothing, groups::SignalGroups)
    entries = _noise_reference_signal_entries(Tuple(groups), ())
    signal_ids = map(first, entries)
    NamedTuple{signal_ids}(
        map(entry -> CorrelatorNoiseDensityEstimator(; num_ants = last(entry)), entries),
    )
end

# An explicit estimator must despread as many antennas as its signal's group, or
# the floor and the prompt describe different arrays. Checked here, where both are
# known, rather than as a shape error deep in the loop update.
@inline function _validate_noise_estimator_num_ants(
    noise_estimators::NamedTuple,
    groups::SignalGroups,
)
    entries = _noise_reference_signal_entries(Tuple(groups), ())
    foreach(entries) do (signal_id, num_ants)
        haskey(noise_estimators, signal_id) || return nothing
        estimator = noise_estimators[signal_id]
        measured = NumAnts(estimator)
        measured === num_ants && return nothing
        throw(
            ArgumentError(
                string(
                    "The noise estimator for signal `:",
                    signal_id,
                    "` measures ",
                    _num_ants_count(measured),
                    " antenna(s), but its signal group declares NumAnts(",
                    _num_ants_count(num_ants),
                    "). Construct it for the group's antenna count, e.g. ",
                    "`CorrelatorNoiseDensityEstimator(; num_ants = NumAnts(",
                    _num_ants_count(num_ants),
                    "))` — the measured noise floor and the prompt it is divided ",
                    "into have to come from the same array.",
                ),
            ),
        )
    end
    nothing
end

@inline _num_ants_count(::NumAnts{M}) where {M} = M

# `(signal id, num_ants)` for the signals that need a noise density, deduplicated
# (one despread per signal per NCO update span) in first-encounter order. A signal fixes its
# band and `_validate_same_band_num_ants` makes groups on one band agree on
# `num_ants`, so the pairing is unambiguous. Folds at compile time.
@inline _noise_reference_signal_entries(::Tuple{}, acc::Tuple) = acc
@inline _noise_reference_signal_entries(t::Tuple, acc::Tuple) =
    _noise_reference_signal_entries(
        Base.tail(t),
        _group_noise_reference_keys(eltype(first(t).satellites), first(t).num_ants, acc),
    )

# Which of a group's signals read a noise density? Asked of the group's slot
# *type* (fixed at construction, by the user's `TrackedSat`s or the template), a
# compile-time constant; branching on `isempty(satellites)` would make the
# `noise_estimators` NamedTuple, and so `TrackState`'s type, a union.
@inline _group_noise_reference_keys(
    ::Type{<:TrackedSat{Signals}},
    num_ants::NumAnts,
    acc::Tuple,
) where {Signals} = _signal_type_noise_reference_keys(Signals, num_ants, acc)

# Signal id and C/N₀ estimator are type parameters of `TrackedSignal`, so this
# walk collapses to a literal tuple.
@inline _signal_type_noise_reference_keys(::Type{Tuple{}}, ::NumAnts, acc::Tuple) = acc
@inline function _signal_type_noise_reference_keys(
    ::Type{T},
    num_ants::NumAnts,
    acc::Tuple,
) where {T<:Tuple}
    head = Base.tuple_type_head(T)
    k = get_signal_id(_signal_type(head))
    seen = any(entry -> first(entry) === k, acc)
    new_acc =
        (seen || !requires_noise_density(_cn0_estimator_type(head))) ? acc :
        (acc..., (k, num_ants))
    _signal_type_noise_reference_keys(Base.tuple_type_tail(T), num_ants, new_acc)
end

@inline _signal_type(::Type{<:TrackedSignal{S}}) where {S} = S

@inline _cn0_estimator_type(
    ::Type{
        <:TrackedSignal{<:Any,<:Any,<:PostCorrelationState{<:Any,<:Any,<:Any,<:Any,CN0}},
    },
) where {CN0} = CN0

# Bare tuple of AbstractGNSSSignal → single :default group NamedTuple.
@inline _normalize_signal_groups(signals::Tuple{Vararg{AbstractGNSSSignal}}) =
    (default = signals,)
# NamedTuple of signal tuples / SignalGroups → pass through unchanged.
@inline _normalize_signal_groups(signals::NamedTuple) = signals

# Build a freshly-templated SignalGroup from a bare signal tuple.
@inline function _normalize_group_entry(
    sig_tuple::Tuple{Vararg{AbstractGNSSSignal}},
    tracking_loop::AbstractTrackingLoop,
    num_ants::NumAnts,
)
    band = get_band(first(sig_tuple))
    _validate_signal_group(sig_tuple, band)
    template = _make_template_tracked_sat(sig_tuple, tracking_loop, num_ants)
    sats = Dictionary{Int,typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, sig_tuple, num_ants)
end

# Pre-built SignalGroup: passes through, except that an empty one whose slot
# type does not fit `tracking_loop` is rebuilt from a fresh template. Its own
# `num_ants` wins over the `TrackState` keyword.
@inline function _normalize_group_entry(
    g::SignalGroup,
    tracking_loop::AbstractTrackingLoop,
    _num_ants::NumAnts,
)
    sats = g.satellites
    if !isempty(sats)
        # The existing sats fix the slot type; it must already match.
        _assert_tracking_loop_types_match(sats, tracking_loop)
        return g
    end
    template = _make_template_tracked_sat(g.signals, tracking_loop, g.num_ants)
    if typeof(template) === eltype(sats)
        return g
    end
    new_sats = Dictionary{Int,typeof(template)}(Int[], typeof(template)[])
    SignalGroup(g.band, new_sats, g.signals, g.num_ants)
end

# Groups on one band share a front end, so they must declare the same `num_ants`.
@inline function _validate_same_band_num_ants(groups::NamedTuple)
    _check_same_band_num_ants(Tuple(groups), ())
end

@inline _check_same_band_num_ants(::Tuple{}, ::Tuple) = nothing
@inline function _check_same_band_num_ants(t::Tuple, seen::Tuple)
    g = first(t)
    bk = get_band_id(g.band)
    for (sk, sna) in seen
        if sk === bk && sna !== g.num_ants
            throw(
                ArgumentError(
                    string(
                        "Two groups on band `:",
                        bk,
                        "` declare different antenna counts (",
                        sna,
                        " vs ",
                        g.num_ants,
                        "). Groups sharing a band must share NumAnts ",
                        "— they're sampled by the same front-end.",
                    ),
                ),
            )
        end
    end
    _check_same_band_num_ants(Base.tail(t), (seen..., (bk, g.num_ants)))
end

"""
$(SIGNATURES)

Build a template `TrackedSat` for `signal_tuple` whose concrete type fixes the
group's dictionary value type; its data (PRN 0, zero Doppler) is meaningless.
"""
function _make_template_tracked_sat(
    signal_tuple::Tuple{Vararg{AbstractGNSSSignal}},
    tracking_loop::AbstractTrackingLoop,
    num_ants::NumAnts,
)
    TrackedSat(signal_tuple, 0, 0.0, 0.0Hz; tracking_loop, num_ants)
end

function TrackState(
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}};
    tracking_loop::AbstractTrackingLoop = _default_tracking_loop(),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_index_offsets = nothing,
)
    sats_dict = to_dictionary(tracked_sats)
    _assert_tracking_loop_types_match(sats_dict, tracking_loop)
    groups = (default = _signal_group_from_dict(sats_dict),)
    TrackState(
        groups,
        tracking_loop,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_index_offsets,
    )
end

# Type-only check that the sats' `loop_state` is what `tracking_loop`
# would produce, with an error naming both types.
@inline function _assert_tracking_loop_types_match(
    dict::Dictionary{<:Any,<:TrackedSat},
    tracking_loop::AbstractTrackingLoop,
)
    isempty(dict) && return nothing
    sat = first(dict.values)
    expected_state_type = typeof(_init_loop_state(tracking_loop, sat))
    actual_state_type = typeof(sat.loop_state)
    expected_state_type === actual_state_type && return nothing
    throw(
        ArgumentError(
            string(
                "TrackedSat has loop_state of type ",
                actual_state_type,
                ", but the configured tracking_loop (",
                typeof(tracking_loop),
                ") would produce ",
                expected_state_type,
                ". Construct each TrackedSat with `tracking_loop = <same instance>` ",
                "as the one passed here.",
            ),
        ),
    )
end

function TrackState(
    tracked_sats::TrackedSatDicts;
    tracking_loop::AbstractTrackingLoop = _default_tracking_loop(),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_index_offsets = nothing,
)
    foreach(d -> _assert_tracking_loop_types_match(d, tracking_loop), Tuple(tracked_sats))
    groups = map(_signal_group_from_dict, tracked_sats)
    TrackState(
        groups,
        tracking_loop,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_index_offsets,
    )
end

# Copy with overrides; the overrides keep the input's concrete type parameters.
function TrackState(
    track_state::TrackState{G,TL,NE};
    groups::Maybe{G} = nothing,
    tracking_loop::Maybe{TL} = nothing,
    noise_estimators::Maybe{NE} = nothing,
    sample_index_offsets::Maybe{Dict{Symbol,Int}} = nothing,
) where {G<:SignalGroups,TL<:AbstractTrackingLoop,NE<:NoiseDensityEstimators}
    TrackState{G,TL,NE}(
        isnothing(groups) ? track_state.groups : groups,
        isnothing(tracking_loop) ? track_state.tracking_loop : tracking_loop,
        isnothing(noise_estimators) ? track_state.noise_estimators : noise_estimators,
        track_state.noise_descriptor,
        isnothing(sample_index_offsets) ? track_state.sample_index_offsets :
        sample_index_offsets,
    )
end

# Recover a SignalGroup (signals, band, antenna count) from a non-empty
# satellites dictionary, for the positional constructors.
@inline function _signal_group_from_dict(dict::Dictionary{<:Any,<:TrackedSat})
    isempty(dict) && throw(
        ArgumentError(
            "Cannot recover the signal-instance tuple from an empty " *
            "satellites dictionary. Use the `TrackState(; signals = ...)` " *
            "constructor to declare signal groups before populating sats.",
        ),
    )
    sat = first(dict.values)
    sig_tuple = map(s -> get_signal(s), sat.signals)
    band = get_band(first(sig_tuple))
    num_ants = NumAnts(get_num_ants(sat))
    SignalGroup(band, dict, sig_tuple, num_ants)
end

# The first copy `track` makes of the caller's state. Detaches the key set as
# well as the values (#123), so adding or removing satellites on the result
# cannot corrupt the input; `track`'s later steps share these keys
# (`_copy_groups_slot_vectors`).
function _prepare_for_new_buffer(track_state::TrackState)
    new_groups = _detach_groups_slot_vectors(track_state.groups)
    _prepare_for_new_buffer!(new_groups)
    TrackState(
        track_state;
        groups = new_groups,
        sample_index_offsets = copy(track_state.sample_index_offsets),
    )
end

function _prepare_for_new_buffer!(track_state::TrackState)
    _prepare_for_new_buffer!(track_state.groups)
    return track_state
end

"""
$(SIGNATURES)

Get the group `group_key`'s dictionary of [`TrackedSat`](@ref)s, keyed by PRN.
"""
function get_tracked_sats(
    tracked_sats::TrackedSatDicts{N},
    group_key::Union{Symbol,Integer,Val},
) where {N}
    _index_group(tracked_sats, group_key)
end

# Works for any group-shaped (named) tuple — per-group satellite dicts
# as well as the `SignalGroup`s NamedTuple itself.
@inline _index_group(t, i::Union{Symbol,Integer}) = t[i]
@inline _index_group(t, ::Val{M}) where {M} = t[M]

function get_tracked_sats(tracked_sats::TrackedSatDicts{1})
    get_tracked_sats(tracked_sats, 1)
end

function get_tracked_sats(
    track_state::TrackState{<:SignalGroups{N}},
    group_key::Union{Symbol,Integer,Val},
) where {N}
    get_tracked_sats(map(g -> g.satellites, track_state.groups), group_key)
end

function get_tracked_sats(track_state::TrackState{<:SignalGroups{1}})
    get_tracked_sats(map(g -> g.satellites, track_state.groups))
end

"""
$(SIGNATURES)

Return the first (driver) signal of the given group, also for a group without
satellites. Resolves at compile time for a literal `group_key`.

For a single-group `TrackState` the group can be omitted.
"""
function get_signal(
    track_state::TrackState{<:SignalGroups{N}},
    group_key::Union{Symbol,Integer,Val},
) where {N}
    first(_index_group(track_state.groups, group_key).signals)
end

function get_signal(track_state::TrackState{<:SignalGroups{1}})
    get_signal(track_state, 1)
end

function get_tracked_sat(
    track_state::TrackState{<:SignalGroups{N}},
    group_key::Union{Symbol,Integer,Val},
    prn,
) where {N}
    get_tracked_sat(get_tracked_sats(track_state, group_key), prn)
end

function get_tracked_sat(track_state::TrackState{<:SignalGroups{1}}, prn)
    get_tracked_sat(track_state, 1, prn)
end

function get_tracked_sat(track_state::TrackState{<:SignalGroups{1}})
    only(get_tracked_sats(track_state, 1))
end

"""
$(SIGNATURES)

Merge already-built [`TrackedSat`](@ref)s into the group `group_key` of
`track_state`, returning a new [`TrackState`](@ref) (the input is left
unchanged). `tracked_sats` may be a single `TrackedSat`, a `Vector`, or a
`Dictionary` keyed by PRN; existing PRNs in the group are overwritten.

Each sat's type must match the group's slot type, including the
`loop_state` the configured tracking loop produces. The tracking loop's
[`update_tracking_loop_on_new_sats`](@ref Tracking.update_tracking_loop_on_new_sats) hook is called once with the incoming sats.

For a single-group `TrackState` the `group_key` may be omitted.
"""
function merge_sats(
    track_state::TrackState{G,TL,NE},
    group_key::Union{Symbol,Integer},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
) where {G<:SignalGroups,TL<:AbstractTrackingLoop,NE<:NoiseDensityEstimators}
    new_sats_dict = to_dictionary(tracked_sats)
    _assert_tracking_loop_types_match(new_sats_dict, track_state.tracking_loop)
    new_tracking_loop =
        update_tracking_loop_on_new_sats(track_state.tracking_loop, new_sats_dict)
    groups = track_state.groups
    g = groups[group_key]
    _assert_sats_match_slot_type(g, new_sats_dict, group_key)
    new_group = SignalGroup(g; satellites = merge(g.satellites, new_sats_dict))
    new_groups = @set groups[group_key] = new_group
    TrackState{G,TL,NE}(
        new_groups,
        new_tracking_loop,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        copy(track_state.sample_index_offsets),
    )
end

function merge_sats(
    track_state::TrackState{<:SignalGroups{1}},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
)
    merge_sats(track_state, 1, tracked_sats)
end

"""
$(SIGNATURES)

Add (or replace, by `prn`) a satellite in `track_state` in place. Builds a
[`TrackedSat`](@ref) for the group with the default correlator bank and post-correlation filter from the acquisition handoff (`prn`, `code_phase`, `code_doppler`,
`carrier_phase` in cycles, `carrier_doppler`), its loop state from TrackingLoops'
`init_loop_state`.

`group_key` may be omitted on a single-group `TrackState`; otherwise an
`ArgumentError` names the available groups.

Keep using the *returned* `TrackState`: a tracking loop whose
[`update_tracking_loop_on_new_sats`](@ref Tracking.update_tracking_loop_on_new_sats) returns a rebuilt tracking loop can only hand
it back there (`TrackState` is immutable). Tracking loops that update in place,
the scalar ones included, return the same `track_state`.

```julia
track_state = TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
add_satellite!(
    track_state;
    prn = 11,
    group_key = :modern_gps,
    code_phase = 0.0,
    carrier_doppler = 1234.0Hz,
)
```

For non-default correlator bank or post-corr-filter *types*, build the `TrackState`
from such `TrackedSat`s and add more with the
`add_satellite!(track_state, group, sat)` overload.
"""
function add_satellite!(
    track_state::TrackState;
    prn::Int,
    group_key::Union{Symbol,Nothing} = nothing,
    kwargs...,
)
    resolved = _resolve_group(track_state, group_key)
    sat = _make_default_tracked_sat_for_group(track_state, resolved; prn, kwargs...)
    add_satellite!(track_state, resolved, sat)
end

"""
$(SIGNATURES)

In-place add (or replace) with a pre-built [`TrackedSat`](@ref). Its type must
match the group's slot type fixed at [`TrackState`](@ref) construction, so
non-default correlator bank or post-corr-filter types need a state built from such
sats. Keep using the return value, as for the keyword form.
"""
function add_satellite!(
    track_state::TrackState{G,TL,NE},
    group_key::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,TL<:AbstractTrackingLoop,NE<:NoiseDensityEstimators}
    _assert_sat_matches_slot_type(track_state, group_key, sat)
    insert_or_set!(_dict_for_group(track_state, group_key), sat.prn, sat)
    new_tracking_loop = update_tracking_loop_on_new_sats(
        track_state.tracking_loop,
        dictionary((sat.prn => sat,)),
    )
    new_tracking_loop === track_state.tracking_loop && return track_state
    TrackState{G,TL,NE}(
        track_state.groups,
        new_tracking_loop,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        track_state.sample_index_offsets,
    )
end

# Name the mismatching types before Dictionaries.jl's `set!` raises a deep
# `convert` MethodError.
@inline function _assert_sat_matches_slot_type(
    track_state::TrackState,
    group_key::Symbol,
    sat::TrackedSat,
)
    SlotT = eltype(track_state.groups[group_key].satellites)
    typeof(sat) === SlotT && return nothing
    throw(
        ArgumentError(
            string(
                "TrackedSat type does not match the `:",
                group_key,
                "` group's slot. ",
                "Got: ",
                typeof(sat),
                ". Expected: ",
                SlotT,
                ". The slot type is fixed at TrackState construction; rebuild the ",
                "sat with the matching correlator bank, post_corr_filter, and ",
                "tracking_loop types.",
            ),
        ),
    )
end

# Dictionary-level variant for `merge_sats`.
@inline function _assert_sats_match_slot_type(
    g::SignalGroup,
    new_sats_dict::Dictionary{<:Any,<:TrackedSat},
    group_key,
)
    SlotT = eltype(g.satellites)
    T = eltype(new_sats_dict)
    T === SlotT && return nothing
    throw(
        ArgumentError(
            string(
                "TrackedSat type does not match the `",
                group_key,
                "` group's slot. ",
                "Got: ",
                T,
                ". Expected: ",
                SlotT,
                ". The slot type is fixed at TrackState construction; rebuild the ",
                "sats with the matching signals, correlator bank, post_corr_filter, and ",
                "tracking_loop types.",
            ),
        ),
    )
end

# Overwrite semantics, like `merge_sats` (`insert!` errors on an existing key).
@inline insert_or_set!(d::Dictionary, k, v) = set!(d, k, v)

"""
$(SIGNATURES)

Immutable variant of [`add_satellite!`](@ref): returns a new
[`TrackState`](@ref) with the satellite added, the input unchanged.
"""
function add_satellite(
    track_state::TrackState;
    prn::Int,
    group_key::Union{Symbol,Nothing} = nothing,
    kwargs...,
)
    resolved = _resolve_group(track_state, group_key)
    sat = _make_default_tracked_sat_for_group(track_state, resolved; prn, kwargs...)
    add_satellite(track_state, resolved, sat)
end

function add_satellite(
    track_state::TrackState{G,TL,NE},
    group_key::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,TL<:AbstractTrackingLoop,NE<:NoiseDensityEstimators}
    _assert_sat_matches_slot_type(track_state, group_key, sat)
    new_tracking_loop = update_tracking_loop_on_new_sats(
        track_state.tracking_loop,
        dictionary((sat.prn => sat,)),
    )
    groups = track_state.groups
    g = groups[group_key]
    new_dict = merge(g.satellites, dictionary((sat.prn => sat,)))
    new_group = SignalGroup(g; satellites = new_dict)
    new_groups = @set groups[group_key] = new_group
    TrackState{G,TL,NE}(
        new_groups,
        new_tracking_loop,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        copy(track_state.sample_index_offsets),
    )
end

"""
$(SIGNATURES)

Remove a satellite from `track_state` in place; a `KeyError` if the group has no
such `prn`.

```julia
remove_satellite!(track_state; prn = 11, group_key = :modern_gps)
```

`group_key` may be omitted as for [`add_satellite!`](@ref). Returns `track_state`.
"""
function remove_satellite!(
    track_state::TrackState;
    prn::Int,
    group_key::Union{Symbol,Nothing} = nothing,
)
    dict = _dict_for_group(track_state, _resolve_group(track_state, group_key))
    haskey(dict, prn) || throw(KeyError(prn))
    delete!(dict, prn)
    # `delete!` may leave a `#undef` hole (see `_satellites_without`); only then
    # rebuild in place (the `SignalGroup.satellites` field is immutable).
    if length(dict) != length(dict.values)
        keys_kept, sats_kept = _satellites_without(dict, prn)
        empty!(dict)
        for (key, sat) in zip(keys_kept, sats_kept)
            insert!(dict, key, sat)
        end
    end
    track_state
end

"""
$(SIGNATURES)

Immutable variant of [`remove_satellite!`](@ref): returns a new
[`TrackState`](@ref) with the satellite removed, the input unchanged. Errors if no
satellite with the given `prn` exists.
"""
function remove_satellite(
    track_state::TrackState{G,TL,NE};
    prn::Int,
    group_key::Union{Symbol,Nothing} = nothing,
) where {G<:SignalGroups,TL<:AbstractTrackingLoop,NE<:NoiseDensityEstimators}
    resolved = _resolve_group(track_state, group_key)
    dict = _dict_for_group(track_state, resolved)
    haskey(dict, prn) || throw(KeyError(prn))
    new_dict = Dictionary(_satellites_without(dict, prn)...)
    groups = track_state.groups
    g = groups[resolved]
    new_group = SignalGroup(g; satellites = new_dict)
    new_groups = @set groups[resolved] = new_group
    TrackState{G,TL,NE}(
        new_groups,
        track_state.tracking_loop,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        copy(track_state.sample_index_offsets),
    )
end

# `(keys, values)` of every satellite except `prn`, in concretely typed vectors.
# Removal must be hole-free (#182): `Dictionaries.delete!` can leave a `#undef` slot
# in `values`, which the hot paths iterate and hit as an `UndefRefError`.
function _satellites_without(dict::Dictionary{I,T}, prn) where {I,T}
    keys_kept = Vector{I}(undef, 0)
    sats_kept = Vector{T}(undef, 0)
    sizehint!(keys_kept, length(dict) - 1)
    sizehint!(sats_kept, length(dict) - 1)
    for (key, sat) in pairs(dict)
        key == prn && continue
        push!(keys_kept, key)
        push!(sats_kept, sat)
    end
    return keys_kept, sats_kept
end

@inline _dict_for_group(track_state::TrackState, group_key::Symbol) =
    track_state.groups[group_key].satellites

# Resolve `group_key=`: `nothing` is the only group of a single-group state (folds to
# a constant), else an error naming the candidates.
@inline _resolve_group(::TrackState, group_key::Symbol) = group_key
@inline function _resolve_group(track_state::TrackState, ::Nothing)
    g = keys(track_state.groups)
    length(g) == 1 && return only(g)
    throw(
        ArgumentError(
            string(
                "track_state has ",
                length(g),
                " groups ",
                g,
                "; pass `group_key = :one_of_them` to select one.",
            ),
        ),
    )
end

# Default-correlator, default-filter TrackedSat for the group; holds the handoff
# keyword defaults of `add_satellite!`/`add_satellite`.
function _make_default_tracked_sat_for_group(
    track_state::TrackState,
    group_key::Symbol;
    prn::Int,
    code_phase = 0.0,
    code_doppler = nothing,
    carrier_phase = 0.0,
    carrier_doppler = 0.0Hz,
)
    g = track_state.groups[group_key]
    TrackedSat(
        g.signals,
        prn,
        code_phase,
        carrier_doppler;
        tracking_loop = track_state.tracking_loop,
        num_ants = g.num_ants,
        carrier_phase,
        code_doppler,
    )
end

# Sat-level accessors: shared by all signals on a sat.
get_prn(track_state::TrackState, id...) = get_prn(get_tracked_sat(track_state, id...))
get_num_ants(track_state::TrackState, id...) =
    get_num_ants(get_tracked_sat(track_state, id...))
get_code_phase(track_state::TrackState, id...) =
    get_code_phase(get_tracked_sat(track_state, id...))
get_code_doppler(track_state::TrackState, id...) =
    get_code_doppler(get_tracked_sat(track_state, id...))
get_carrier_phase(track_state::TrackState, id...) =
    get_carrier_phase(get_tracked_sat(track_state, id...))
get_carrier_doppler(track_state::TrackState, id...) =
    get_carrier_doppler(get_tracked_sat(track_state, id...))
get_next_sample_index(track_state::TrackState, id...) =
    get_next_sample_index(get_tracked_sat(track_state, id...))
get_carrier_loop_stage(track_state::TrackState, id...) =
    get_carrier_loop_stage(get_tracked_sat(track_state, id...))
is_phase_locked(track_state::TrackState, id...) =
    is_phase_locked(get_tracked_sat(track_state, id...))
get_phase_lock_indicator(track_state::TrackState, id...) =
    get_phase_lock_indicator(get_tracked_sat(track_state, id...))
get_loop_bandwidths(track_state::TrackState, id...) =
    get_loop_bandwidths(get_tracked_sat(track_state, id...))

"""
$(SIGNATURES)

The carrier phase (cycles, not wrapped) the tracking loop adds to the replica's
[`get_carrier_phase`](@ref) for its estimate of the incoming carrier phase, the
carrier-phase measurement (RINEX `L`); `0.0` for the scalar loops and for a
satellite out of the vector loop.

In TrackingLoops' `VectorTrackingLoop` the navigation filter sets the carrier NCO
and the residual PLL follows the carrier without steering it, so the replica drifts
off the carrier by whole cycles. This is the residual PLL's prediction at the end of
the latest [`track`](@ref)/[`track!`](@ref) call, every cycle counted (see
`TrackingLoops.calc_carrier_phase_correction`, which also names where it steps), at the
band's sampling frequency the latest integration stored. The group can be omitted on a
single-group `TrackState`.

The whole cycles only count against a replica phase that counts them too:
`get_carrier_phase` wraps, so a RINEX host builds `L` from its own accumulated
replica phase (e.g. the integrated carrier Doppler) plus this.

Only valid right after `track`/`track!`: a host feeding correlator outputs through
[`append_correlator_bank_output!`](@ref) calls `TrackingLoops.calc_carrier_phase_correction`
with its own NCO history and sample index instead.
"""
function calc_carrier_phase_correction(
    track_state::TrackState,
    group_key::Union{Symbol,Integer,Val},
    prn,
)
    group = _index_group(track_state.groups, group_key)
    sat = get_tracked_sat(group.satellites, prn)
    # `track!` leaves the replica at the band's sample offset (the end of the
    # samples), having run on the satellite's Doppler since its latest integration.
    band_grid_index = get(track_state.sample_index_offsets, get_band_id(group.band), 0)
    # Every integration stores its band's sampling frequency in the post-correlation
    # state; before the first it is still NaN, but TrackingLoops returns a correction of
    # 0 then without reading it.
    sampling_frequency = first(sat.signals).post_correlation_state.sampling_frequency
    calc_carrier_phase_correction(
        sat.loop_state,
        _software_nco_command(sat),
        band_grid_index,
        sampling_frequency,
    )
end

calc_carrier_phase_correction(track_state::TrackState{<:SignalGroups{1}}, prn) =
    calc_carrier_phase_correction(track_state, 1, prn)

# Per-signal accessors, addressed as (using `get_correlator_bank`):
#   * `get_correlator_bank(track_state)` — 1 group, 1 sat, 1 signal.
#   * `get_correlator_bank(track_state, prn)` — 1 group, 1 signal.
#   * `get_correlator_bank(track_state, group_key, prn)` — multi-group, 1 signal.
#   * `get_correlator_bank(track_state, group_key, prn, sig)` — per-signal, `sig` an index
#     or signal type; always names the group (`:default` or `1` for one group).
# Without `sig` they are only valid on single-signal sats.
const _SignalSelector = Union{Integer,Type{<:AbstractGNSSSignal}}

for fn in (
    :get_num_integrated_samples,
    :get_correlator_bank,
    :get_correlator_bank_outputs,
    :get_last_correlator_bank,
    :get_last_filtered_prompt,
    :get_last_num_code_periods,
    :get_last_coherent_integration_time,
    :get_filtered_prompts,
    :get_post_corr_filter,
    :get_cn0_estimator,
    :get_bit_synchronizer,
    :get_soft_bits,
    :get_num_bits,
    :is_bit_or_secondary_code_synced,
    :get_bit_polarity,
    :get_num_code_periods_into_bit,
    :estimate_cn0,
    :get_preferred_num_code_periods_to_integrate,
)
    @eval begin
        $fn(track_state::TrackState, id...) = $fn(get_tracked_sat(track_state, id...))
        $fn(
            track_state::TrackState{<:SignalGroups},
            group_key::Union{Symbol,Integer,Val},
            prn,
            sig::_SignalSelector,
        ) = $fn(get_tracked_sat(track_state, group_key, prn), sig)
    end
end

"""
$(SIGNATURES)

Append an externally built
[`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput) to the addressed
signal's buffer and return `track_state`. Addressed like the per-signal accessors
(e.g. [`get_correlator_bank_outputs`](@ref)):

  - `append_correlator_bank_output!(track_state, output)` — one group, one sat, one signal.
  - `append_correlator_bank_output!(track_state, output, prn)` — one group, one signal.
  - `append_correlator_bank_output!(track_state, output, group_key, prn)` — multi-group, one signal.
  - `append_correlator_bank_output!(track_state, output, group_key, prn, sig)` — per-signal.

The contract is on the signal-level method
(`append_correlator_bank_output!(::TrackedSignal, output)`); see it and
[External correlator producers](@ref).
"""
append_correlator_bank_output!(
    track_state::TrackState,
    output::CorrelatorBankOutput,
    id...,
) = (
    append_correlator_bank_output!(get_tracked_sat(track_state, id...), output);
    track_state
)
append_correlator_bank_output!(
    track_state::TrackState{<:SignalGroups},
    output::CorrelatorBankOutput,
    group_key::Union{Symbol,Integer,Val},
    prn,
    sig::_SignalSelector,
) = (
    append_correlator_bank_output!(
        get_tracked_sat(track_state, group_key, prn),
        output,
        sig,
    );
    track_state
)

"""
$(SIGNATURES)

Append an externally built [`NoiseObservation`](@extref TrackingLoops.NoiseObservation) to the addressed
**signal**'s noise estimator and return `track_state`:

  - `append_noise_observation!(track_state, obs)` — single-signal `TrackState`.
  - `append_noise_observation!(track_state, obs, :GPSL1CA)` — explicit signal id.
  - `append_noise_observation!(track_state, obs, GPSL1CA)` — or the signal type /
    an instance of it, which is what a caller usually has to hand.

The hardware/FPGA fill path: [`append_correlator_bank_output!`](@ref) and
`append_noise_observation!` per signal, then
[`step_tracking_loops!`](@ref) to step the loops (separate mechanisms; see the
estimator-level `append_noise_observation!` for how they differ). Per signal, not
per band, because the post-correlation floor depends on the despreading modulation:
a noise channel is a tracking channel with a wrong PRN (see [`AbstractNoiseDensityEstimator`](@extref TrackingLoops.AbstractNoiseDensityEstimator)).

The signal must have a noise estimator: provisioned when its C/N₀ estimator reads
a density (see
[`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:TrackingLoops.AbstractCN0Estimator}})),
or declared through `TrackState`'s `noise_estimators` keyword.
"""
append_noise_observation!(
    track_state::TrackState,
    observation::NoiseObservation,
    signal_id::Symbol,
) = (
    append_noise_observation!(
        _noise_estimator_for_signal(track_state, signal_id),
        observation,
    );
    track_state
)

append_noise_observation!(
    track_state::TrackState,
    observation::NoiseObservation,
    signal::Union{AbstractGNSSSignal,Type{<:AbstractGNSSSignal}},
) = append_noise_observation!(track_state, observation, get_signal_id(signal))

function append_noise_observation!(track_state::TrackState, observation::NoiseObservation)
    noise_estimators = track_state.noise_estimators
    length(noise_estimators) == 1 || throw(
        ArgumentError(
            string(
                "track_state has ",
                length(noise_estimators),
                " noise estimators ",
                keys(noise_estimators),
                "; pass the signal id, e.g. ",
                "`append_noise_observation!(track_state, obs, :GPSL1CA)`.",
            ),
        ),
    )
    append_noise_observation!(only(Tuple(noise_estimators)), observation)
    track_state
end

# Error naming what is configured, not a bare `KeyError`.
@inline function _noise_estimator_for_signal(track_state::TrackState, signal_id::Symbol)
    noise_estimators = track_state.noise_estimators
    haskey(noise_estimators, signal_id) ||
        _throw_no_noise_estimator(signal_id, keys(noise_estimators))
    noise_estimators[signal_id]
end

@noinline _throw_no_noise_estimator(signal_id::Symbol, configured) = throw(
    ArgumentError(
        string(
            "no noise estimator configured for signal `:",
            signal_id,
            "`; this TrackState has ",
            isempty(configured) ? "none at all" : string(configured),
            ". A signal is provisioned automatically only where its ",
            "C/N₀ estimator reads a noise density (see `requires_noise_density`) ",
            "— pass `noise_estimators = (",
            signal_id,
            " = CorrelatorNoiseDensityEstimator(),)` to `TrackState` to declare one.",
        ),
    ),
)

# Index of the addressed signal within a sat's signals (config time only).
_signal_index(signals::Tuple) = length(signals) == 1 ? 1 : _throw_needs_signal_selector()
_signal_index(::Tuple, i::Integer) = Int(i)
function _signal_index(signals::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    idx = findfirst(s -> get_signal(s) isa T, signals)
    isnothing(idx) && throw(ArgumentError("no signal of type $T on this satellite"))
    # A type selector must be unambiguous, as in `_find_signal_by_type`.
    isnothing(findnext(s -> get_signal(s) isa T, signals, idx + 1)) || throw(
        ArgumentError(
            "more than one signal of type $T on this satellite — " *
            "pass an integer index to address one of them.",
        ),
    )
    idx
end

# Rebuild `sat` with the addressed signal's integration length set to `N`.
function _set_sat_signal_preferred_code_periods(sat::TrackedSat, N::Int, sel...)
    idx = _signal_index(sat.signals, sel...)
    idx_tuple = ntuple(identity, length(sat.signals))
    new_signals = map(sat.signals, idx_tuple) do s, i
        i == idx ? TrackedSignal(s; preferred_num_code_periods_to_integrate = N) : s
    end
    TrackedSat(sat; signals = new_signals)
end

"""
$(SIGNATURES)

Set the preferred coherent-integration length, in primary code periods, of one
signal on one satellite (the addressed [`TrackedSignal`](@ref)'s
`preferred_num_code_periods_to_integrate`). The actual length is capped by the
bit/secondary-code period and held at 1 until sync (see
`calc_num_code_periods_to_integrate`); the scalar tracking loops cap their
bandwidths against the integration time to stay stable.

For data-bearing signals the length must divide the code periods per bit (20 for
GPS L1 C/A, 10 for GPS L5I), or an `ArgumentError` is thrown (#128). Pilots accept
any length of at least one code period.

Addressed like the per-signal accessors; the per-signal form names the group even
on a single-group `TrackState`:

```julia
set_preferred_num_code_periods_to_integrate!(ts, :gps_l5, 1, GPSL5I, 10)  # (group_key, prn, signal)
set_preferred_num_code_periods_to_integrate!(ts, :gps_l5, 1, 10)          # single-signal sat
set_preferred_num_code_periods_to_integrate!(ts, 1, 10)                   # single-group state
set_preferred_num_code_periods_to_integrate!(ts, 10)                      # 1 group, 1 sat, 1 signal
```

Mutates `track_state` in place and returns it.
"""
function set_preferred_num_code_periods_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group_key::Union{Symbol,Integer,Val},
    prn::Integer,
    sig::_SignalSelector,
    num_code_periods::Integer,
)
    _set_preferred_code_periods!(
        track_state,
        get_tracked_sats(track_state, group_key),
        prn,
        num_code_periods,
        sig,
    )
end

function set_preferred_num_code_periods_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group_key::Union{Symbol,Integer,Val},
    prn::Integer,
    num_code_periods::Integer,
)
    _set_preferred_code_periods!(
        track_state,
        get_tracked_sats(track_state, group_key),
        prn,
        num_code_periods,
    )
end

function set_preferred_num_code_periods_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    prn::Integer,
    num_code_periods::Integer,
)
    _set_preferred_code_periods!(
        track_state,
        get_tracked_sats(track_state),
        prn,
        num_code_periods,
    )
end

# Shared body for the three addressing overloads above.
@inline function _set_preferred_code_periods!(
    track_state::TrackState,
    sats,
    prn,
    num_code_periods::Integer,
    sel...,
)
    sats[prn] =
        _set_sat_signal_preferred_code_periods(sats[prn], Int(num_code_periods), sel...)
    track_state
end

function set_preferred_num_code_periods_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    num_code_periods::Integer,
)
    sats = get_tracked_sats(track_state)
    prn = only(keys(sats))
    sats[prn] = _set_sat_signal_preferred_code_periods(sats[prn], Int(num_code_periods))
    track_state
end

# Re-initialise one satellite's loop state (`reinit_loop_state`) and zero each signal's
# latest filtered correlator bank, so the next integration's previous prompt is 0, for
# which `fll_disc` returns 0 (see `reinit_loop_states!`).
@inline function _reinit_sat_loop_states(track_state::TrackState, sat::TrackedSat)
    new_signals = map(sat.signals) do s
        state = s.post_correlation_state
        TrackedSignal(
            s;
            post_correlation_state = PostCorrelationState(
                state;
                last_filtered_correlator_bank = zero(state.last_filtered_correlator_bank),
            ),
        )
    end
    TrackedSat(
        sat;
        signals = new_signals,
        loop_state = reinit_loop_state(
            track_state.tracking_loop,
            sat.loop_state,
            sat.carrier_doppler,
            sat.code_doppler,
        ),
    )
end

"""
$(SIGNATURES)

Re-initialise the tracking loop's state of every satellite (or one addressed
satellite) from its current Doppler. For the scalar tracking loops this zeroes the
loop-filter integrators and keeps `carrier_doppler`, `code_doppler`, any
per-satellite bandwidth override and the carrier-loop stage. Each signal's
latest filtered prompt is cleared, so the first FLL update does not measure a rotation
spanning the old integration interval.

Use it when a signal's integration length changes mid-track, e.g. GPS L5I from
1 ms to 10 ms via [`set_preferred_num_code_periods_to_integrate!`](@ref): the
integrator state does not carry over to the new update interval, and the
transient can drag the loop out of lock.

No satellite id re-initialises every satellite; `(group_key, prn)` or (single-group) `prn`
re-initialises one. Mutates `track_state` and returns it. Works for any
[`AbstractTrackingLoop`](@extref TrackingLoops.AbstractTrackingLoop)
through its `reinit_loop_state` method.

```julia
set_preferred_num_code_periods_to_integrate!(track_state, 1, GPSL5I, 10)
reinit_loop_states!(track_state, 1)          # clean handoff for PRN 1
reinit_loop_states!(track_state)             # …or every satellite
```
"""
function reinit_loop_states!(track_state::TrackState)
    for g in Tuple(track_state.groups)
        vals = g.satellites.values
        @inbounds for i in eachindex(vals)
            vals[i] = _reinit_sat_loop_states(track_state, vals[i])
        end
    end
    track_state
end

function reinit_loop_states!(
    track_state::TrackState{<:SignalGroups},
    group_key::Union{Symbol,Integer,Val},
    prn,
)
    sats = get_tracked_sats(track_state, group_key)
    sats[prn] = _reinit_sat_loop_states(track_state, sats[prn])
    track_state
end

reinit_loop_states!(track_state::TrackState{<:SignalGroups{1}}, prn) =
    reinit_loop_states!(track_state, 1, prn)

# Distinct band ids across `groups`, in first-encounter order; folds at compile time.
@inline _band_ids_in_groups(::Tuple{}, acc::Tuple{Vararg{Symbol}}) = acc
@inline function _band_ids_in_groups(t::Tuple, acc::Tuple{Vararg{Symbol}})
    k = get_band_id(first(t).band)
    new_acc = k in acc ? acc : (acc..., k)
    _band_ids_in_groups(Base.tail(t), new_acc)
end

"""
$(SIGNATURES)

The distinct band ids (`GNSSSignals.get_band_id`) of `track_state`'s groups, in
first-encounter order: the keys a [`SamplesByBand`](@ref Tracking.SamplesByBand) NamedTuple must use.

```julia
ts = TrackState(; signals = (legacy_gps = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
get_band_ids(ts) == (:L1, :L5)
```
"""
@inline get_band_ids(track_state::TrackState) =
    _band_ids_in_groups(Tuple(track_state.groups), ())

"""
$(SIGNATURES)

The tracking loop `track_state` runs. TrackingLoops' readers of a tracking loop's
configuration read it, e.g. `get_loop_bandwidths(get_tracking_loop(track_state), signal)` or, of a `VectorTrackingLoop`, `get_signal_pairs`, `get_scalar_loop`,
`get_navigation_update_interval`, `get_min_tracking_cn0` and `get_nav_filter_config`.
"""
get_tracking_loop(track_state::TrackState) = track_state.tracking_loop

# TrackingLoops' navigation-solution interface, forwarded to the tracking loop (the
# scalar loops answer `nothing`) and re-exported, with the vector loop's members.
TrackingLoops.get_navigation_solution(track_state::TrackState) =
    TrackingLoops.get_navigation_solution(track_state.tracking_loop)
TrackingLoops.get_navigation_status(track_state::TrackState) =
    TrackingLoops.get_navigation_status(track_state.tracking_loop)
TrackingLoops.get_navigation_update_index(track_state::TrackState) =
    TrackingLoops.get_navigation_update_index(track_state.tracking_loop)
TrackingLoops.get_navigation_elapsed_time(track_state::TrackState) =
    TrackingLoops.get_navigation_elapsed_time(track_state.tracking_loop)
TrackingLoops.get_satellite_report(
    track_state::TrackState,
    signal::AbstractGNSSSignal,
    prn::Integer,
) = TrackingLoops.get_satellite_report(track_state.tracking_loop, signal, prn)
TrackingLoops.get_decoder(
    track_state::TrackState,
    signal::AbstractGNSSSignal,
    prn::Integer,
) = TrackingLoops.get_decoder(track_state.tracking_loop, signal, prn)
TrackingLoops.get_vector_tracking_members(track_state::TrackState) =
    TrackingLoops.get_vector_tracking_members(track_state.tracking_loop)
