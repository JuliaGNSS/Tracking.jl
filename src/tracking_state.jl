"""
$(SIGNATURES)

Construct a fresh `TrackState` from a declaration of which signals each
group tracks. A *group* is a set of satellites with the same signal tuple, and
so the same concrete `TrackedSat` type (type stability of the hot loop). The
first signal of each tuple drives the Doppler estimator, which updates the
satellite-shared carrier and code Doppler. The group key (`:modern_gps`, …) is
what `add_satellite!` refers to.

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

The per-group `TrackedSat` type (correlator, post-corr filter and
doppler-estimator-state types) is frozen at construction, so `add_satellite!`
cannot change it. For non-default correlator or PCF *types*, build the
`TrackedSat`s yourself and construct from them (e.g. `TrackState(signal, sats)`).

`noise_estimators` declares the per-signal noise sources
([`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)s keyed by
`GNSSSignals.get_signal_id`, e.g. `:GPSL1CA`). Left at `nothing`, every signal whose
C/N₀ estimator reads a noise density (see
[`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:TrackingLoops.AbstractCN0Estimator}}))
gets a [`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator)
and every other signal none, so a state on
[`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator) despreads nothing. Pass a
NamedTuple to configure the window, or to declare a source you fill with
[`append_noise_observation!`](@ref) on a correlator-ingest path.

`sample_offsets` sets where each band's sample grid starts, keyed by `get_band_id`
(`NamedTuple` or `AbstractDict`; missing bands start at `0`). Pass the old state's
when rebuilding around a Doppler estimator that keeps state of its own; see the
[`TrackState`](@ref) struct. The other constructors take the same keyword.
"""
function TrackState(;
    signal::Maybe{AbstractGNSSSignal} = nothing,
    signals = nothing,
    doppler_estimator::Maybe{AbstractDopplerEstimator} = nothing,
    num_ants::NumAnts = NumAnts(1),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_offsets = nothing,
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
    # signal (`init_estimator_state`), so groups need no shared compromise.
    estimator =
        isnothing(doppler_estimator) ? ConventionalAssistedPLLAndDLL() : doppler_estimator
    # An entry is a bare signal tuple (takes `num_ants`) or a pre-built
    # `SignalGroup` (carries its own band and `num_ants`).
    groups = map(sig_groups_nt) do entry
        _normalize_group_entry(entry, estimator, num_ants)
    end
    _validate_same_band_num_ants(groups)
    TrackState(
        groups,
        estimator,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_offsets,
    )
end

# `nothing` provisions a `CorrelatorNoiseEstimator` only for the signals that
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
        map(entry -> CorrelatorNoiseEstimator(; num_ants = last(entry)), entries),
    )
end

# An explicitly passed estimator must despread as many antennas as its signal's
# group, or the floor and the prompt describe different arrays. Checked here,
# where both are known, rather than as a shape error deep in the fold;
# `NumAnts(estimator)` reads the count of any `AbstractNoiseEstimator`.
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
                    "`CorrelatorNoiseEstimator(; num_ants = NumAnts(",
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
# (one despread per signal per chunk) in first-encounter order. The pairing is
# unambiguous: a signal fixes its band, and `_validate_same_band_num_ants` makes
# groups on one band agree on `num_ants`. Folds at compile time.
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

@inline _signal_type(::Type{<:TrackedSignal{Sig}}) where {Sig} = Sig

@inline _cn0_estimator_type(
    ::Type{<:TrackedSignal{Sig,B,C,PCF,CN0}},
) where {Sig,B,C,PCF,CN0} = CN0

# A signal type's band id, folded to a literal symbol (bands are singletons).
@inline _signal_band_id(::Type{Sig}) where {Sig<:AbstractGNSSSignal} =
    get_band_id(get_band(Sig))

# Bare tuple of AbstractGNSSSignal → single :default group NamedTuple.
@inline _normalize_signal_groups(signals::Tuple{Vararg{AbstractGNSSSignal}}) =
    (default = signals,)
# NamedTuple of signal tuples / SignalGroups → pass through unchanged.
@inline _normalize_signal_groups(signals::NamedTuple) = signals

# Build a freshly-templated SignalGroup from a bare signal tuple.
@inline function _normalize_group_entry(
    sig_tuple::Tuple{Vararg{AbstractGNSSSignal}},
    doppler_estimator::AbstractDopplerEstimator,
    num_ants::NumAnts,
)
    band = get_band(first(sig_tuple))
    _validate_signal_group(sig_tuple, band)
    template = _make_template_tracked_sat(sig_tuple, doppler_estimator, num_ants)
    sats = Dictionary{Int,typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, sig_tuple, num_ants)
end

# Pre-built SignalGroup: passes through, except that an empty one whose slot
# type does not fit `doppler_estimator` is rebuilt from a fresh template. Its own
# `num_ants` wins over the `TrackState` keyword.
@inline function _normalize_group_entry(
    g::SignalGroup,
    doppler_estimator::AbstractDopplerEstimator,
    _num_ants::NumAnts,
)
    sats = g.satellites
    if !isempty(sats)
        # The existing sats fix the slot type; it must already match.
        _assert_doppler_estimator_types_match(sats, doppler_estimator)
        return g
    end
    template = _make_template_tracked_sat(g.signals, doppler_estimator, g.num_ants)
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
    doppler_estimator::AbstractDopplerEstimator,
    num_ants::NumAnts,
)
    TrackedSat(signal_tuple, 0, 0.0, 0.0Hz; doppler_estimator, num_ants)
end

function TrackState(
    signal::AbstractGNSSSignal,
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}};
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_offsets = nothing,
)
    # `signal` is unused (each sat carries its signals); kept for compatibility.
    sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(sats_dict, doppler_estimator)
    groups = (default = _signal_group_from_dict(sats_dict),)
    TrackState(
        groups,
        doppler_estimator,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_offsets,
    )
end

function TrackState(
    tracked_sats::Dictionary{<:Any,<:TrackedSat};
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_offsets = nothing,
)
    _assert_doppler_estimator_types_match(tracked_sats, doppler_estimator)
    groups = (default = _signal_group_from_dict(tracked_sats),)
    TrackState(
        groups,
        doppler_estimator,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_offsets,
    )
end

# Type-only check that the sats' `doppler_estimator_state` is what `estimator`
# would produce, with an error naming both types.
@inline function _assert_doppler_estimator_types_match(
    dict::Dictionary{<:Any,<:TrackedSat},
    estimator::AbstractDopplerEstimator,
)
    isempty(dict) && return nothing
    sat = first(dict.values)
    expected_state_type = typeof(init_estimator_state(estimator, sat))
    actual_state_type = typeof(sat.doppler_estimator_state)
    expected_state_type === actual_state_type && return nothing
    throw(
        ArgumentError(
            string(
                "TrackedSat has doppler_estimator_state of type ",
                actual_state_type,
                ", but the configured doppler_estimator (",
                typeof(estimator),
                ") would produce ",
                expected_state_type,
                ". Construct each TrackedSat with `doppler_estimator = <same instance>` ",
                "as the one passed here.",
            ),
        ),
    )
end

function TrackState(
    satellites::SatelliteDicts;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    noise_estimators::Maybe{NamedTuple} = nothing,
    sample_offsets = nothing,
)
    foreach(
        d -> _assert_doppler_estimator_types_match(d, doppler_estimator),
        Tuple(satellites),
    )
    groups = map(_signal_group_from_dict, satellites)
    TrackState(
        groups,
        doppler_estimator,
        _resolve_noise_estimators(noise_estimators, groups);
        sample_offsets,
    )
end

# Copy with overrides; the overrides keep the input's concrete type parameters.
function TrackState(
    track_state::TrackState{G,DE,NE};
    groups::Maybe{G} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
    noise_estimators::Maybe{NE} = nothing,
    sample_offsets::Maybe{Dict{Symbol,Int}} = nothing,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    TrackState{G,DE,NE}(
        isnothing(groups) ? track_state.groups : groups,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
        isnothing(noise_estimators) ? track_state.noise_estimators : noise_estimators,
        track_state.noise_descriptor,
        isnothing(sample_offsets) ? track_state.sample_offsets : sample_offsets,
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
    sig_tuple = map(s -> s.signal, sat.signals)
    band = get_band(first(sig_tuple))
    num_ants = NumAnts(get_num_ants(sat))
    SignalGroup(band, dict, sig_tuple, num_ants)
end

# The first copy `track` makes of the caller's state. Detaches the key set as
# well as the values (#123), so adding or removing satellites on the result
# cannot corrupt the input; `track`'s later steps share these keys
# (`_copy_groups_slot_vectors`).
function reset_start_sample_and_bit_buffer(track_state::TrackState)
    new_groups = _detach_groups_slot_vectors(track_state.groups)
    reset_start_sample_and_bit_buffer!(new_groups)
    TrackState(
        track_state;
        groups = new_groups,
        sample_offsets = copy(track_state.sample_offsets),
    )
end

function reset_start_sample_and_bit_buffer!(track_state::TrackState)
    reset_start_sample_and_bit_buffer!(track_state.groups)
    return track_state
end

"""
$(SIGNATURES)

Get the per-group satellite dictionary for the given `group_key`.
"""
function get_sat_states(
    satellites::SatelliteDicts{N},
    group_key::Union{Symbol,Integer,Val},
) where {N}
    _index_group(satellites, group_key)
end

# Works for any group-shaped (named) tuple — per-group satellite dicts
# as well as the `SignalGroup`s NamedTuple itself.
@inline _index_group(t, i::Union{Symbol,Integer}) = t[i]
@inline _index_group(t, ::Val{M}) where {M} = t[M]

function get_sat_states(satellites::SatelliteDicts{1})
    get_sat_states(satellites, 1)
end

function get_sat_states(
    track_state::TrackState{<:SignalGroups{N}},
    group_key::Union{Symbol,Integer,Val},
) where {N}
    get_sat_states(map(g -> g.satellites, track_state.groups), group_key)
end

function get_sat_states(track_state::TrackState{<:SignalGroups{1}})
    get_sat_states(map(g -> g.satellites, track_state.groups))
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

function get_sat_state(
    track_state::TrackState{<:SignalGroups{N}},
    group_key::Union{Symbol,Integer,Val},
    prn,
) where {N}
    get_sat_state(get_sat_states(track_state, group_key), prn)
end

function get_sat_state(track_state::TrackState{<:SignalGroups{1}}, prn)
    get_sat_state(track_state, 1, prn)
end

function get_sat_state(track_state::TrackState{<:SignalGroups{1}})
    only(get_sat_states(track_state, 1))
end

"""
$(SIGNATURES)

Merge already-built [`TrackedSat`](@ref)s into the group `group_key` of
`track_state`, returning a new [`TrackState`](@ref) (the input is left
unchanged). `tracked_sats` may be a single `TrackedSat`, a `Vector`, or a
`Dictionary` keyed by PRN; existing PRNs in the group are overwritten.

Each sat's type must match the group's slot type, including the
`doppler_estimator_state` the configured estimator produces. The estimator's
[`update_estimator_on_handoff`](@ref Tracking.update_estimator_on_handoff) hook is called once with the incoming sats.

For a single-group `TrackState` the `group_key` may be omitted.
"""
function merge_sats(
    track_state::TrackState{G,DE,NE},
    group_key::Union{Symbol,Integer},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    new_sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(new_sats_dict, track_state.doppler_estimator)
    new_estimator =
        update_estimator_on_handoff(track_state.doppler_estimator, new_sats_dict)
    groups = track_state.groups
    g = groups[group_key]
    _assert_sats_match_slot_type(g, new_sats_dict, group_key)
    new_group = SignalGroup(g; satellites = merge(g.satellites, new_sats_dict))
    new_groups = @set groups[group_key] = new_group
    TrackState{G,DE,NE}(
        new_groups,
        new_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        copy(track_state.sample_offsets),
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
[`TrackedSat`](@ref) for the group with the default correlator and post-corr
filter from the acquisition handoff (`prn`, `code_phase`, `code_doppler`,
`carrier_phase`, `carrier_doppler`), its estimator state from
[`init_estimator_state`](@ref).

`group_key` may be omitted on a single-group `TrackState`; otherwise an
`ArgumentError` names the available groups.

Keep using the *returned* `TrackState`: an estimator whose
[`update_estimator_on_handoff`](@ref Tracking.update_estimator_on_handoff) returns a rebuilt estimator can only hand
it back there (`TrackState` is immutable). Estimators that update in place,
the conventional ones included, return the same `track_state`.

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

For non-default correlator or post-corr-filter *types*, build the `TrackState`
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
non-default correlator or post-corr-filter types need a state built from such
sats. Keep using the
return value, as for the keyword form.
"""
function add_satellite!(
    track_state::TrackState{G,DE,NE},
    group_key::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    _assert_sat_matches_slot_type(track_state, group_key, sat)
    insert_or_set!(_dict_for_group(track_state, group_key), sat.prn, sat)
    new_estimator = update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((sat.prn => sat,)),
    )
    new_estimator === track_state.doppler_estimator && return track_state
    TrackState{G,DE,NE}(
        track_state.groups,
        new_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        track_state.sample_offsets,
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
                "sat with the matching correlator, post_corr_filter, and ",
                "doppler_estimator types.",
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
                "sats with the matching signals, correlator, post_corr_filter, and ",
                "doppler_estimator types.",
            ),
        ),
    )
end

# Overwrite semantics, like `merge_sats` (`insert!` errors on an existing key).
@inline insert_or_set!(d::Dictionary, k, v) = set!(d, k, v)

"""
$(SIGNATURES)

Immutable variant of [`add_satellite!`](@ref). Returns a new
[`TrackState`](@ref) with the satellite added; the input is left
unchanged.
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
    track_state::TrackState{G,DE,NE},
    group_key::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    _assert_sat_matches_slot_type(track_state, group_key, sat)
    new_estimator = update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((sat.prn => sat,)),
    )
    groups = track_state.groups
    g = groups[group_key]
    new_dict = merge(g.satellites, dictionary((sat.prn => sat,)))
    new_group = SignalGroup(g; satellites = new_dict)
    new_groups = @set groups[group_key] = new_group
    TrackState{G,DE,NE}(
        new_groups,
        new_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        copy(track_state.sample_offsets),
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
    # `delete!` may leave a `#undef` hole (see `_satellites_without`); rebuild in
    # place (the `SignalGroup.satellites` field is immutable) only when it did.
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

Immutable variant of [`remove_satellite!`](@ref). Returns a new
[`TrackState`](@ref) with the satellite removed; the input is left
unchanged. Errors if no satellite with the given `prn` exists.
"""
function remove_satellite(
    track_state::TrackState{G,DE,NE};
    prn::Int,
    group_key::Union{Symbol,Nothing} = nothing,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    resolved = _resolve_group(track_state, group_key)
    dict = _dict_for_group(track_state, resolved)
    haskey(dict, prn) || throw(KeyError(prn))
    new_dict = Dictionary(_satellites_without(dict, prn)...)
    groups = track_state.groups
    g = groups[resolved]
    new_group = SignalGroup(g; satellites = new_dict)
    new_groups = @set groups[resolved] = new_group
    TrackState{G,DE,NE}(
        new_groups,
        track_state.doppler_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
        copy(track_state.sample_offsets),
    )
end

# `(keys, values)` of every satellite except `prn`, in concretely typed vectors.
# Removal must be hole-free (#182): `Dictionaries.delete!` can leave a `#undef`
# slot in `values` that the hot paths, which iterate `values`, hit as an
# `UndefRefError`.
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

# Resolve `group_key=`: `nothing` means the only group of a single-group state
# (folds to a constant), and errors with the candidates otherwise.
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

# Default-correlator, default-PCF TrackedSat for the group; holds the handoff
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
        doppler_estimator = track_state.doppler_estimator,
        num_ants = g.num_ants,
        carrier_phase,
        code_doppler,
    )
end

# Sat-level accessors: shared by all signals on a sat.
get_prn(s::TrackState, id...) = get_prn(get_sat_state(s, id...))
get_num_ants(s::TrackState, id...) = get_num_ants(get_sat_state(s, id...))
get_code_phase(s::TrackState, id...) = get_code_phase(get_sat_state(s, id...))
get_code_doppler(s::TrackState, id...) = get_code_doppler(get_sat_state(s, id...))
get_carrier_phase(s::TrackState, id...) = get_carrier_phase(get_sat_state(s, id...))
get_carrier_doppler(s::TrackState, id...) = get_carrier_doppler(get_sat_state(s, id...))
get_signal_start_sample(s::TrackState, id...) =
    get_signal_start_sample(get_sat_state(s, id...))

"""
$(SIGNATURES)

The carrier phase (radians, not wrapped) the Doppler estimator adds to the
replica's [`get_carrier_phase`](@ref) for its estimate of the incoming carrier
phase, the carrier-phase measurement (RINEX `L`); `0.0` for the scalar loops and
for a satellite out of the vector loop.

In TrackingLoops' `VectorPLLAndDLL` loop the navigation filter sets the carrier
NCO, and the residual PLL follows the carrier without steering it, so the
replica drifts off the carrier by whole cycles; this is the residual PLL's
prediction at the end of the latest [`track`](@ref)/[`track!`](@ref) call, every
cycle counted (see `TrackingLoops.carrier_phase_correction`, which also names
where it steps). `sampling_frequency` is the band's. The group can be omitted on
a single-group `TrackState`.

The whole cycles only count against a replica phase that counts them too:
`get_carrier_phase` wraps, so a RINEX host builds `L` from its own accumulated
replica phase (e.g. the integrated carrier Doppler) plus this.

Only valid right after `track`/`track!`: a host feeding records through
[`append_correlator_output!`](@ref) calls `TrackingLoops.carrier_phase_correction`
with its own NCO words and sample instead.
"""
function get_carrier_phase_correction(
    track_state::TrackState{<:SignalGroups{N}},
    group_key::Union{Symbol,Integer,Val},
    prn,
    sampling_frequency,
) where {N}
    signal_group = _index_group(track_state.groups, group_key)
    sat = get_sat_state(signal_group.satellites, prn)
    # `track!` leaves the replica at the end of the measurement, the band's sample
    # offset, having run on the satellite's Doppler since its latest record.
    sample = get(track_state.sample_offsets, get_band_id(signal_group.band), 0)
    2π * carrier_phase_correction(
        sat.doppler_estimator_state,
        _software_words(sat),
        sample,
        sampling_frequency,
    )
end

get_carrier_phase_correction(
    track_state::TrackState{<:SignalGroups{1}},
    prn,
    sampling_frequency,
) = get_carrier_phase_correction(track_state, 1, prn, sampling_frequency)

# Per-signal accessors, addressed as (using `get_correlator`):
#   * `get_correlator(track_state)` — 1 group, 1 sat, 1 signal.
#   * `get_correlator(track_state, prn)` — 1 group, 1 signal.
#   * `get_correlator(track_state, group_key, prn)` — multi-group, 1 signal.
#   * `get_correlator(track_state, group_key, prn, sig)` — per-signal, `sig` an index
#     or signal type; always names the group (`:default` or `1` for one group).
# Without `sig` they are only valid on single-signal sats.
const _SignalSelector = Union{Integer,Type{<:AbstractGNSSSignal}}

for fn in (
    :get_integrated_samples,
    :get_correlator,
    :get_correlator_outputs,
    :get_last_fully_integrated_correlator,
    :get_last_fully_integrated_filtered_prompt,
    :get_last_fully_integrated_num_code_blocks,
    :get_last_fully_integrated_integration_time,
    :get_filtered_prompts,
    :get_post_corr_filter,
    :get_cn0_estimator,
    :get_bit_buffer,
    :get_soft_bits,
    :get_num_bits,
    :has_bit_or_secondary_code_been_found,
    :estimate_cn0,
    :get_preferred_num_code_blocks_to_integrate,
)
    @eval begin
        $fn(s::TrackState, id...) = $fn(get_sat_state(s, id...))
        $fn(
            s::TrackState{<:SignalGroups},
            group_key::Union{Symbol,Integer,Val},
            prn,
            sig::_SignalSelector,
        ) = $fn(get_sat_state(s, group_key, prn), sig)
    end
end

"""
$(SIGNATURES)

Append an externally built
[`CorrelatorOutput`](@extref TrackingLoops.CorrelatorOutput) to the addressed
signal's buffer and return `track_state`. Addressed like the per-signal accessors
(e.g. [`get_correlator_outputs`](@ref)):

  - `append_correlator_output!(track_state, output)` — one group, one sat, one signal.
  - `append_correlator_output!(track_state, output, prn)` — one group, one signal.
  - `append_correlator_output!(track_state, output, group_key, prn)` — multi-group, one signal.
  - `append_correlator_output!(track_state, output, group_key, prn, sig)` — per-signal.

The signal-level method (`append_correlator_output!(::TrackedSignal, output)`)
carries the contract; see it and [External correlator producers](@ref).
"""
append_correlator_output!(s::TrackState, output::CorrelatorOutput, id...) =
    (append_correlator_output!(get_sat_state(s, id...), output); s)
append_correlator_output!(
    s::TrackState{<:SignalGroups},
    output::CorrelatorOutput,
    group_key::Union{Symbol,Integer,Val},
    prn,
    sig::_SignalSelector,
) = (append_correlator_output!(get_sat_state(s, group_key, prn), output, sig); s)

"""
$(SIGNATURES)

Append an externally built [`NoiseObservation`](@extref TrackingLoops.NoiseObservation) to the addressed
**signal**'s noise estimator and return `track_state`:

  - `append_noise_observation!(track_state, obs)` — single-signal `TrackState`.
  - `append_noise_observation!(track_state, obs, :GPSL1CA)` — explicit signal id.
  - `append_noise_observation!(track_state, obs, GPSL1CA)` — or the signal type /
    an instance of it, which is what a caller usually has to hand.

The hardware/FPGA fill path: [`append_correlator_output!`](@ref) and
`append_noise_observation!` per signal, then
[`estimate_dopplers_and_filter_prompt!`](@ref) to fold (deliberately separate
mechanisms; see the estimator-level `append_noise_observation!` for how they
differ). Per signal, not per band, because the post-correlation floor depends on
the despreading modulation: a noise channel is a tracking channel with a wrong
PRN (see [`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)).

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
            " = CorrelatorNoiseEstimator(),)` to `TrackState` to declare one.",
        ),
    ),
)

# Index of the addressed signal within a sat's signals (config time only).
_signal_index(signals::Tuple) = length(signals) == 1 ? 1 : _throw_needs_signal_selector()
_signal_index(::Tuple, i::Integer) = Int(i)
function _signal_index(signals::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    idx = findfirst(s -> s.signal isa T, signals)
    isnothing(idx) && throw(ArgumentError("no signal of type $T on this satellite"))
    # A type selector must be unambiguous, as in `_find_signal_by_type`.
    isnothing(findnext(s -> s.signal isa T, signals, idx + 1)) || throw(
        ArgumentError(
            "more than one signal of type $T on this satellite — " *
            "pass an integer index to address one of them.",
        ),
    )
    idx
end

# Rebuild `sat` with the addressed signal's integration length set to `N`.
function _set_sat_signal_preferred_blocks(sat::TrackedSat, N::Int, sel...)
    idx = _signal_index(sat.signals, sel...)
    idx_tuple = ntuple(identity, length(sat.signals))
    new_signals = map(sat.signals, idx_tuple) do s, i
        i == idx ? TrackedSignal(s; preferred_num_code_blocks_to_integrate = N) : s
    end
    TrackedSat(sat; signals = new_signals)
end

"""
$(SIGNATURES)

Set the preferred coherent-integration length, in primary code blocks, of one
signal on one satellite (the addressed [`TrackedSignal`](@ref)'s
`preferred_num_code_blocks_to_integrate`). The actual length is capped by the
bit/secondary-code period and held at 1 until sync (see
`calc_num_code_blocks_to_integrate`); the conventional estimator caps its loop
bandwidths against the integration time, so the loop stays stable.

For data-bearing signals the length must divide the code blocks per bit (20 for
GPS L1 C/A, 10 for GPS L5I), or an `ArgumentError` is thrown (#128). Pilots accept
any length of at least one block.

Addressed like the per-signal accessors; the per-signal form names the group even
on a single-group `TrackState`:

```julia
set_preferred_num_code_blocks_to_integrate!(ts, :gps_l5, 1, GPSL5I, 10)  # (group_key, prn, signal)
set_preferred_num_code_blocks_to_integrate!(ts, :gps_l5, 1, 10)          # single-signal sat
set_preferred_num_code_blocks_to_integrate!(ts, 1, 10)                   # single-group state
set_preferred_num_code_blocks_to_integrate!(ts, 10)                      # 1 group, 1 sat, 1 signal
```

Mutates `track_state` in place and returns it.
"""
function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group_key::Union{Symbol,Integer,Val},
    prn::Integer,
    sig::_SignalSelector,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state,
        get_sat_states(track_state, group_key),
        prn,
        num_code_blocks,
        sig,
    )
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group_key::Union{Symbol,Integer,Val},
    prn::Integer,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state,
        get_sat_states(track_state, group_key),
        prn,
        num_code_blocks,
    )
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    prn::Integer,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(track_state, get_sat_states(track_state), prn, num_code_blocks)
end

# Shared body for the three addressing overloads above.
@inline function _set_preferred_blocks!(
    track_state::TrackState,
    sats,
    prn,
    num_code_blocks::Integer,
    sel...,
)
    sats[prn] = _set_sat_signal_preferred_blocks(sats[prn], Int(num_code_blocks), sel...)
    track_state
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    num_code_blocks::Integer,
)
    sats = get_sat_states(track_state)
    prn = only(keys(sats))
    sats[prn] = _set_sat_signal_preferred_blocks(sats[prn], Int(num_code_blocks))
    track_state
end

# Re-seed one satellite's estimator state (`_reset_estimator_state`) and zero
# each signal's previous prompt, for which `fll_disc` returns 0 (see
# `reset_loop_filters!`).
@inline function _reset_sat_loop_filters(track_state::TrackState, sat::TrackedSat)
    new_signals = map(
        s ->
            TrackedSignal(s; last_fully_integrated_filtered_prompt = complex(0.0, 0.0)),
        sat.signals,
    )
    TrackedSat(
        sat;
        signals = new_signals,
        doppler_estimator_state = _reset_estimator_state(
            track_state.doppler_estimator,
            sat,
        ),
    )
end

"""
$(SIGNATURES)

Re-seed the Doppler-estimator state of every satellite (or one addressed
satellite) from its current Doppler. For the conventional PLL/DLL this zeroes the
loop-filter integrators and keeps `carrier_doppler`, `code_doppler`, any
per-satellite bandwidth override and the carrier-loop stage. Each signal's
`last_fully_integrated_filtered_prompt` is cleared, so the first FLL update does
not measure a rotation spanning the old integration interval.

Use it when a signal's integration length changes mid-track, e.g. GPS L5I from
1 ms to 10 ms via [`set_preferred_num_code_blocks_to_integrate!`](@ref): the
integrator state does not carry over to the new update interval, and the
transient can drag the loop out of lock.

No satellite id resets every satellite; `(group_key, prn)` or (single-group) `prn`
resets one. Mutates `track_state` and returns it. Works for any
[`AbstractDopplerEstimator`](@extref TrackingLoops.AbstractDopplerEstimator)
through its [`init_estimator_state`](@ref) hook.

```julia
set_preferred_num_code_blocks_to_integrate!(track_state, 1, GPSL5I, 10)
reset_loop_filters!(track_state, 1)          # clean handoff for PRN 1
reset_loop_filters!(track_state)             # …or reset every satellite
```
"""
function reset_loop_filters!(track_state::TrackState)
    for g in Tuple(track_state.groups)
        vals = g.satellites.values
        @inbounds for i in eachindex(vals)
            vals[i] = _reset_sat_loop_filters(track_state, vals[i])
        end
    end
    track_state
end

function reset_loop_filters!(
    track_state::TrackState{<:SignalGroups},
    group_key::Union{Symbol,Integer,Val},
    prn,
)
    sats = get_sat_states(track_state, group_key)
    sats[prn] = _reset_sat_loop_filters(track_state, sats[prn])
    track_state
end

reset_loop_filters!(track_state::TrackState{<:SignalGroups{1}}, prn) =
    reset_loop_filters!(track_state, 1, prn)

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
first-encounter order: the keys a multi-band measurements NamedTuple must use.

```julia
ts = TrackState(; signals = (legacy_gps = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
get_band_ids(ts) == (:L1, :L5)
```
"""
@inline get_band_ids(track_state::TrackState) =
    _band_ids_in_groups(Tuple(track_state.groups), ())

# The one band of a single-band state, for the bare-buffer `track` entry points.
@inline function _single_band(track_state::TrackState)
    keys_tuple = get_band_ids(track_state)
    if length(keys_tuple) != 1
        throw(
            ArgumentError(
                string(
                    "Bare-buffer `track`/`track!` requires a single-band TrackState, ",
                    "but this TrackState spans bands ",
                    keys_tuple,
                    ". Pass a NamedTuple of `BandMeasurement`s instead, ",
                    "one per band.",
                ),
            ),
        )
    end
    first(track_state.groups).band
end

# `Set(a) == Set(b)` without allocating, for short tuples.
@inline function _tuple_sets_equal(a::Tuple, b::Tuple)
    length(a) == length(b) || return false
    all(x -> x in b, a) && all(y -> y in a, b)
end

# Measurements must match the band keys and antenna shapes and have exactly equal
# durations.
@inline function _validate_measurements(
    track_state::TrackState,
    measurements::BandMeasurements,
)
    expected_keys = get_band_ids(track_state)
    got_keys = keys(measurements)
    if !_tuple_sets_equal(expected_keys, got_keys)
        throw(
            ArgumentError(
                string(
                    "BandMeasurement keys do not match the TrackState's band set. ",
                    "Expected: ",
                    expected_keys,
                    ". Got: ",
                    got_keys,
                    ".",
                ),
            ),
        )
    end
    _validate_antenna_shapes(track_state, measurements)
    _validate_equal_durations(measurements)
    return nothing
end

# Vector → 1 antenna; Matrix → cols = num_ants. Tuple recursion, since a `for`
# loop over heterogeneous groups boxes.
@inline function _validate_antenna_shapes(
    track_state::TrackState,
    measurements::BandMeasurements,
)
    _validate_antenna_shapes_walk(Tuple(track_state.groups), measurements)
end

@inline _validate_antenna_shapes_walk(::Tuple{}, ::BandMeasurements) = nothing
@inline function _validate_antenna_shapes_walk(
    groups::Tuple,
    measurements::BandMeasurements,
)
    g = first(groups)
    m = measurements[get_band_id(g.band)]
    _assert_antenna_shape(g, m)
    _validate_antenna_shapes_walk(Base.tail(groups), measurements)
end

@inline function _assert_antenna_shape(
    g::SignalGroup{B,S,Sigs,NumAnts{M}},
    m::BandMeasurement,
) where {B,S,Sigs,M}
    cols = m.samples isa AbstractMatrix ? size(m.samples, 2) : 1
    cols == M && return nothing
    throw(
        ArgumentError(
            string(
                "Antenna shape mismatch for band `:",
                get_band_id(g.band),
                "`. Group declares NumAnts(",
                M,
                ") but the measurement's ",
                "`samples` has ",
                cols,
                " column(s).",
            ),
        ),
    )
end

# Cross-multiplied with promoted rates, so mixed-precision sampling rates with
# identical durations compare equal.
@inline function _validate_equal_durations(measurements::BandMeasurements)
    ms = Tuple(measurements)
    isempty(ms) && return nothing
    first_m = first(ms)
    ref_num_samples = get_num_samples(first_m)
    for m in Base.tail(ms)
        ref_fs, fs = promote(first_m.sampling_frequency, m.sampling_frequency)
        get_num_samples(m) * ref_fs == ref_num_samples * fs && continue
        throw(
            ArgumentError(
                string(
                    "BandMeasurement durations must be exactly equal across bands. ",
                    "Got ",
                    uconvert(s, get_num_samples(m) / m.sampling_frequency),
                    " and ",
                    uconvert(s, ref_num_samples / first_m.sampling_frequency),
                    ". ",
                    "Check `num_samples / sampling_frequency` for each band.",
                ),
            ),
        )
    end
    return nothing
end

# TrackingLoops' navigation-solution interface, forwarded to the estimator (the
# scalar loops answer `nothing`).
TrackingLoops.navigation_solution(track_state::TrackState) =
    TrackingLoops.navigation_solution(track_state.doppler_estimator)
TrackingLoops.navigation_status(track_state::TrackState) =
    TrackingLoops.navigation_status(track_state.doppler_estimator)
TrackingLoops.navigation_cycle(track_state::TrackState) =
    TrackingLoops.navigation_cycle(track_state.doppler_estimator)
TrackingLoops.navigation_epoch(track_state::TrackState) =
    TrackingLoops.navigation_epoch(track_state.doppler_estimator)
TrackingLoops.satellite_report(
    track_state::TrackState,
    signal::AbstractGNSSSignal,
    prn::Integer,
) = TrackingLoops.satellite_report(track_state.doppler_estimator, signal, prn)
