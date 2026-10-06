"""
$(SIGNATURES)

Construct a fresh `TrackState` from a declaration of which signals each
group tracks. Each entry in `signals` is a tuple of `AbstractGNSSSignal`s
whose first signal is the [Estimator-driver signal](@ref), whose correlator feeds
the Doppler estimator;
the group key (`:modern_gps`, …) is what `add_satellite!` later refers to.

```julia
track_state = TrackState(;
    signals = (
        legacy_gps = (GPSL1CA(),),
        modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
        galileo = (GalileoE1B(),),
    ),
)
```

For the common case of one group tracking one signal, use the singular
`signal` keyword instead:

```julia
track_state = TrackState(; signal = GPSL1CA())
```

This is equivalent to `TrackState(; signals = (default = (GPSL1CA(),),))`.
`add_satellite!` may then omit the `group=` keyword.

Each group's `TrackedSat` type (correlator, post-corr filter, estimator state)
is frozen at construction; for non-default correlator or PCF *types*, build the
`TrackedSat`s yourself and use `add_satellite!(track_state, group, sat)`.

`noise_estimators` declares the per-signal noise sources
([`AbstractNoiseEstimator`](@ref)s keyed by signal id, `GNSSSignals.get_signal_id`
— `:GPSL1CA`, `:GalileoE1B`, …). Left at `nothing` it is derived: a
[`CorrelatorNoiseEstimator`](@ref) for every signal whose C/N₀ estimator reads a
noise density (see [`requires_noise_density`](@ref)), and **no entry** for any
other signal, so a state that stays on [`NWPRCN0Estimator`](@ref) runs no
despread at all. Pass an explicit NamedTuple to configure the window, or to
declare a signal's source on a correlator-ingest path where you fill it with
[`append_noise_observation!`](@ref) rather than from samples.
"""
function TrackState(;
    signal::Maybe{AbstractGNSSSignal} = nothing,
    signals = nothing,
    doppler_estimator::Maybe{AbstractDopplerEstimator} = nothing,
    num_ants::NumAnts = NumAnts(1),
    noise_estimators::Maybe{NamedTuple} = nothing,
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
    # Auto-bandwidth default: `init_estimator_state` sizes each sat's loop from
    # its own group's driver signal, so groups need no cross-group compromise.
    estimator =
        isnothing(doppler_estimator) ? ConventionalAssistedPLLAndDLL() : doppler_estimator
    # Entries are bare signal tuples (use the `num_ants` kwarg) or pre-built
    # `SignalGroup`s (carry their own band / num_ants).
    groups = map(sig_groups_nt) do entry
        _normalize_group_entry(entry, estimator, num_ants)
    end
    _validate_same_band_num_ants(groups)
    TrackState(groups, estimator, _resolve_noise_estimators(noise_estimators, groups))
end

# Resolve the `noise_estimators` kwarg; `nothing` provisions a
# `CorrelatorNoiseEstimator` only for signals that need a density.
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

# An explicitly passed noise estimator must cover as many antennas as its signal's
# group, or the measured floor and the prompt describe different arrays. Checked at
# construction (not in the per-chunk fold) for an actionable error; reads the count
# off `noise_density_type`, so it works for any `AbstractNoiseEstimator`.
@inline function _validate_noise_estimator_num_ants(
    noise_estimators::NamedTuple,
    groups::SignalGroups,
)
    entries = _noise_reference_signal_entries(Tuple(groups), ())
    foreach(entries) do (signal_id, num_ants)
        haskey(noise_estimators, signal_id) || return nothing
        estimator = noise_estimators[signal_id]
        measured = _num_ants_of_density_type(noise_density_type(estimator))
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
# same-band groups agree on `num_ants`. Tuple recursion folds at compile time.
@inline _noise_reference_signal_entries(::Tuple{}, acc::Tuple) = acc
@inline _noise_reference_signal_entries(t::Tuple, acc::Tuple) =
    _noise_reference_signal_entries(
        Base.tail(t),
        _group_noise_reference_keys(eltype(first(t).satellites), first(t).num_ants, acc),
    )

# Which of a group's signals use a C/N₀ estimator that reads a noise density?
# Asked of the group's slot type, never of a satellite value: the slot type is
# right for empty and populated groups alike and is a compile-time constant, so
# the `noise_estimators` NamedTuple (and thus `TrackState`'s type) stays inferable.
@inline _group_noise_reference_keys(
    ::Type{<:TrackedSat{Signals}},
    num_ants::NumAnts,
    acc::Tuple,
) where {Signals} = _signal_type_noise_reference_keys(Signals, num_ants, acc)

# Recurse down the signals' tuple type; folds to a literal `(Symbol, NumAnts)` tuple.
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

# The band id of a signal type, folded to a literal symbol at compile time so the
# noise walk can index `BandMeasurements` without a runtime lookup.
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

# Pre-built SignalGroup → pass through, unless it is empty and its slot type
# doesn't match the TrackState's estimator; then rebuild the empty dict. The
# group's own `num_ants` wins over the TrackState kwarg.
@inline function _normalize_group_entry(
    g::SignalGroup,
    doppler_estimator::AbstractDopplerEstimator,
    _num_ants::NumAnts,
)
    sats = g.satellites
    if !isempty(sats)
        # Populated: the slot type is fixed by the sats, so it must already match.
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

# Same-band groups must declare identical `num_ants` (one front-end per band).
# O(num_groups²), folded at compile time for a concrete groups type.
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

Internal helper: build a template `TrackedSat` for a group's signal tuple. Its
data is meaningless (PRN 0, zero Doppler); only its type is used, to fix the
dictionary value type at TrackState construction.
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
)
    # `signal` is unused (implied by each sat's driver signal); kept for
    # backward compatibility.
    sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(sats_dict, doppler_estimator)
    groups = (default = _signal_group_from_dict(sats_dict),)
    TrackState(
        groups,
        doppler_estimator,
        _resolve_noise_estimators(noise_estimators, groups),
    )
end

function TrackState(
    tracked_sats::Dictionary{<:Any,<:TrackedSat};
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    noise_estimators::Maybe{NamedTuple} = nothing,
)
    _assert_doppler_estimator_types_match(tracked_sats, doppler_estimator)
    groups = (default = _signal_group_from_dict(tracked_sats),)
    TrackState(
        groups,
        doppler_estimator,
        _resolve_noise_estimators(noise_estimators, groups),
    )
end

# Type-only check that the sats' `doppler_estimator_state` matches what
# `estimator` would produce; throws an ArgumentError naming both types.
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
)
    foreach(
        d -> _assert_doppler_estimator_types_match(d, doppler_estimator),
        Tuple(satellites),
    )
    groups = map(_signal_group_from_dict, satellites)
    TrackState(
        groups,
        doppler_estimator,
        _resolve_noise_estimators(noise_estimators, groups),
    )
end

# Copy-with-overrides constructor; overrides are constrained to the input's
# concrete types so the result's type parameters are preserved.
function TrackState(
    track_state::TrackState{G,DE,NE};
    groups::Maybe{G} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
    noise_estimators::Maybe{NE} = nothing,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    TrackState{G,DE,NE}(
        isnothing(groups) ? track_state.groups : groups,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
        isnothing(noise_estimators) ? track_state.noise_estimators : noise_estimators,
        track_state.noise_descriptor,
    )
end

# Build a SignalGroup from a non-empty satellites dictionary, recovering the
# signal tuple, band and antenna count from its first sat.
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

# Immutable reset — the first copy `track` makes of the caller's state. Detaches
# keys and values (`_detach_groups_slot_vectors`, #123) so later
# `add_satellite!`/`remove_satellite!` on the result cannot corrupt the input.
function reset_start_sample_and_bit_buffer(track_state::TrackState)
    new_groups = _detach_groups_slot_vectors(track_state.groups)
    reset_start_sample_and_bit_buffer!(new_groups)
    TrackState(track_state; groups = new_groups)
end

function reset_start_sample_and_bit_buffer!(track_state::TrackState)
    reset_start_sample_and_bit_buffer!(track_state.groups)
    return track_state
end

"""
$(SIGNATURES)

Get the per-group satellite dictionary for a specific signal-group index.
"""
function get_sat_states(
    satellites::SatelliteDicts{N},
    group_idx::Union{Symbol,Integer,Val},
) where {N}
    _index_group(satellites, group_idx)
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
    group_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sat_states(map(g -> g.satellites, track_state.groups), group_idx)
end

function get_sat_states(track_state::TrackState{<:SignalGroups{1}})
    get_sat_states(map(g -> g.satellites, track_state.groups))
end

"""
$(SIGNATURES)

Return the first (driver) signal of the given group, read off the group's
declared signal tuple, so it also works on groups without satellites.

For a single-group `TrackState` the index can be omitted.
"""
function get_signal(
    track_state::TrackState{<:SignalGroups{N}},
    group_idx::Union{Symbol,Integer,Val},
) where {N}
    first(_index_group(track_state.groups, group_idx).signals)
end

function get_signal(track_state::TrackState{<:SignalGroups{1}})
    get_signal(track_state, 1)
end

function get_sat_state(
    track_state::TrackState{<:SignalGroups{N}},
    group_idx::Union{Symbol,Integer,Val},
    sat_identifier,
) where {N}
    get_sat_state(get_sat_states(track_state, group_idx), sat_identifier)
end

function get_sat_state(track_state::TrackState{<:SignalGroups{1}}, sat_identifier)
    get_sat_state(track_state, 1, sat_identifier)
end

function get_sat_state(track_state::TrackState{<:SignalGroups{1}})
    only(get_sat_states(track_state, 1))
end

"""
$(SIGNATURES)

Merge already-built [`TrackedSat`](@ref)s into the `group_idx` group of
`track_state`, returning a new [`TrackState`](@ref) (the input is left
unchanged). `tracked_sats` may be a single `TrackedSat`, a `Vector`, or a
`Dictionary` keyed by PRN; existing PRNs in the group are overwritten.

Each sat's `doppler_estimator_state` must match the type the
`track_state`'s configured estimator produces (checked up front), and the
sat's signal-tuple shape must match the group's slot type. The estimator's
[`update_estimator_on_handoff`](@ref) hook is invoked once with the incoming
sats so estimators with cross-satellite shared state can grow it.

For a single-group `TrackState` the `group_idx` may be omitted.
"""
function merge_sats(
    track_state::TrackState{G,DE,NE},
    group_idx::Union{Symbol,Integer},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    new_sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(new_sats_dict, track_state.doppler_estimator)
    new_estimator =
        update_estimator_on_handoff(track_state.doppler_estimator, new_sats_dict)
    groups = track_state.groups
    g = groups[group_idx]
    _assert_sats_match_slot_type(g, new_sats_dict, group_idx)
    new_group = SignalGroup(g; satellites = merge(g.satellites, new_sats_dict))
    new_groups = @set groups[group_idx] = new_group
    TrackState{G,DE,NE}(
        new_groups,
        new_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
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

Add (or replace) a satellite in `track_state` in place. Builds a
multi-signal [`TrackedSat`](@ref) for the requested group using the
library default correlator and post-corr filter, with the supplied
acquisition-handoff values (`prn`, `code_phase`, `code_doppler`,
`carrier_phase`, `carrier_doppler`) wired into each `TrackedSignal`.
The per-satellite doppler-estimator state is initialized via
[`init_estimator_state`](@ref) against the `TrackState`'s configured
estimator.

When `group` is omitted, a single-group `TrackState` uses its only group
(whatever it is named); a multi-group `TrackState` requires the key and
otherwise throws an `ArgumentError` naming the available groups. If a
satellite with the same `prn` already exists in that group's dictionary,
it is overwritten.

The satellite dictionary is mutated in place, but keep using the *returned*
`TrackState`: if [`update_estimator_on_handoff`](@ref) rebuilds the estimator,
only the returned (immutable) `TrackState` carries it. Estimators that update in
place get the same `track_state` back.

```julia
track_state = TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
add_satellite!(
    track_state;
    prn = 11,
    group = :modern_gps,
    code_phase = 0.0,
    carrier_doppler = 1234.0Hz,
)
```

To use a non-default correlator or post-corr-filter *type*, construct
the [`TrackedSat`](@ref) yourself and call the
`add_satellite!(track_state, group, sat)` overload — see below.
"""
function add_satellite!(
    track_state::TrackState;
    prn::Int,
    group::Union{Symbol,Nothing} = nothing,
    kwargs...,
)
    resolved = _resolve_group(track_state, group)
    sat = _make_default_tracked_sat_for_group(track_state, resolved; prn, kwargs...)
    add_satellite!(track_state, resolved, sat)
end

"""
$(SIGNATURES)

In-place add (or replace) with a pre-built [`TrackedSat`](@ref), for non-default
correlator or post-corr-filter types. The sat's type must match the group's slot
type fixed at [`TrackState`](@ref) construction, or an `ArgumentError` is thrown.

Like the keyword form, the returned `TrackState` carries the estimator
returned by [`update_estimator_on_handoff`](@ref) — keep using the
return value.
"""
function add_satellite!(
    track_state::TrackState{G,DE,NE},
    group::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    _assert_sat_matches_slot_type(track_state, group, sat)
    insert_or_set!(_dict_for_group(track_state, group), sat.prn, sat)
    new_estimator = update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((sat.prn => sat,)),
    )
    # A rebuilt estimator can only be honored through the return value.
    new_estimator === track_state.doppler_estimator && return track_state
    TrackState{G,DE,NE}(
        track_state.groups,
        new_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
    )
end

# Check `sat` has exactly the group's slot type, with a clear ArgumentError
# instead of Dictionaries.jl's deep `convert` MethodError.
@inline function _assert_sat_matches_slot_type(
    track_state::TrackState,
    group::Symbol,
    sat::TrackedSat,
)
    SlotT = eltype(track_state.groups[group].satellites)
    typeof(sat) === SlotT && return nothing
    throw(
        ArgumentError(
            string(
                "TrackedSat type does not match the `:",
                group,
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

# Dictionary-level variant of `_assert_sat_matches_slot_type`, for `merge_sats`.
@inline function _assert_sats_match_slot_type(
    g::SignalGroup,
    new_sats_dict::Dictionary{<:Any,<:TrackedSat},
    group_idx,
)
    SlotT = eltype(g.satellites)
    T = eltype(new_sats_dict)
    T === SlotT && return nothing
    throw(
        ArgumentError(
            string(
                "TrackedSat type does not match the `",
                group_idx,
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

# Overwrite semantics (matching `merge_sats`); `insert!` would error on a key.
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
    group::Union{Symbol,Nothing} = nothing,
    kwargs...,
)
    resolved = _resolve_group(track_state, group)
    sat = _make_default_tracked_sat_for_group(track_state, resolved; prn, kwargs...)
    add_satellite(track_state, resolved, sat)
end

function add_satellite(
    track_state::TrackState{G,DE,NE},
    group::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    _assert_sat_matches_slot_type(track_state, group, sat)
    new_estimator = update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((sat.prn => sat,)),
    )
    groups = track_state.groups
    g = groups[group]
    new_dict = merge(g.satellites, dictionary((sat.prn => sat,)))
    new_group = SignalGroup(g; satellites = new_dict)
    new_groups = @set groups[group] = new_group
    TrackState{G,DE,NE}(
        new_groups,
        new_estimator,
        track_state.noise_estimators,
        track_state.noise_descriptor,
    )
end

"""
$(SIGNATURES)

Remove a satellite from `track_state` in place. Throws a `KeyError` if no
satellite with the given `prn` exists in the named group (same contract as
the immutable [`remove_satellite`](@ref)).

```julia
remove_satellite!(track_state; prn = 11, group = :modern_gps)
```

When `group` is omitted it is inferred the same way as
[`add_satellite!`](@ref): a single-group TrackState uses its only group;
a multi-group TrackState requires the key.
Returns `track_state` unchanged (the dictionary is mutated in place).
"""
function remove_satellite!(
    track_state::TrackState;
    prn::Int,
    group::Union{Symbol,Nothing} = nothing,
)
    dict = _dict_for_group(track_state, _resolve_group(track_state, group))
    haskey(dict, prn) || throw(KeyError(prn))
    delete!(dict, prn)
    # `delete!` may leave a `#undef` hole (see `_satellites_without`); rebuild in
    # place only when it did. In place, and not by swapping in a fresh dict,
    # because the `SignalGroup.satellites` field is immutable.
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
    group::Union{Symbol,Nothing} = nothing,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    resolved = _resolve_group(track_state, group)
    dict = _dict_for_group(track_state, resolved)
    haskey(dict, prn) || throw(KeyError(prn))
    # `Dictionary` adopts the value vector without copying, so wrapping the
    # surviving entries is a single allocation (no copy-then-delete round trip).
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
    )
end

# Collect the `(keys, values)` of every satellite except `prn` into fresh,
# concretely-typed vectors. Removal must be hole-free (issue #182):
# `Dictionaries.delete!` can leave a `#undef` slot in `values` that the hot paths
# (which iterate `values` directly) hit as an `UndefRefError`.
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

# The satellites dictionary of group `group`; unknown keys raise NamedTuple errors.
@inline _dict_for_group(track_state::TrackState, group::Symbol) =
    track_state.groups[group].satellites

# Resolve the `group=` keyword: `nothing` picks the only group of a single-group
# TrackState (folds to a constant, keeping type stability) and errors otherwise.
@inline _resolve_group(::TrackState, group::Symbol) = group
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
                "; pass `group = :one_of_them` to select one.",
            ),
        ),
    )
end

# Build a default TrackedSat for the group's signals and antenna count. Holds the
# handoff kwarg defaults for `add_satellite!` / `add_satellite`.
function _make_default_tracked_sat_for_group(
    track_state::TrackState,
    group::Symbol;
    prn::Int,
    code_phase = 0.0,
    code_doppler = nothing,
    carrier_phase = 0.0,
    carrier_doppler = 0.0Hz,
)
    g = track_state.groups[group]
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

# Sat-level accessors. These do not vary across signals on a sat (one PRN
# per sat, one shared code/carrier Doppler/phase, one signal_start_sample).
get_prn(s::TrackState, id...) = get_prn(get_sat_state(s, id...))
get_num_ants(s::TrackState, id...) = get_num_ants(get_sat_state(s, id...))
get_code_phase(s::TrackState, id...) = get_code_phase(get_sat_state(s, id...))
get_code_doppler(s::TrackState, id...) = get_code_doppler(get_sat_state(s, id...))
get_carrier_phase(s::TrackState, id...) = get_carrier_phase(get_sat_state(s, id...))
get_carrier_phase_polarity(s::TrackState, id...) =
    get_carrier_phase_polarity(get_sat_state(s, id...))
get_carrier_doppler(s::TrackState, id...) = get_carrier_doppler(get_sat_state(s, id...))
get_signal_start_sample(s::TrackState, id...) =
    get_signal_start_sample(get_sat_state(s, id...))

# Per-signal accessors. Addressing forms (e.g. `get_correlator`):
#   * `get_correlator(track_state)` — 1 group, 1 sat, 1 signal.
#   * `get_correlator(track_state, prn)` — 1 group, 1 signal.
#   * `get_correlator(track_state, group, prn)` — multi-group, 1 signal.
#   * `get_correlator(track_state, group, prn, sig)` — per-signal; `sig` is an
#     `Integer` index or a signal type. Always names the group, even on a
#     single-group TrackState (use `:default` or `1`).
const _SignalSelector = Union{Integer,Type{<:AbstractGNSSSignal}}

# `(s, id...)` forwards to `get_sat_state` (the sat-level forms);
# `(s, group, prn, sig)` is the per-signal form.
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
    :get_group_delay,
)
    @eval begin
        $fn(s::TrackState, id...) = $fn(get_sat_state(s, id...))
        $fn(
            s::TrackState{<:SignalGroups},
            group::Union{Symbol,Integer,Val},
            sat_id,
            sig::_SignalSelector,
        ) = $fn(get_sat_state(s, group, sat_id), sig)
    end
end

"""
$(SIGNATURES)

Append an externally built [`CorrelatorOutput`](@ref) to the buffer of the
addressed signal and return `track_state`. The `output` is followed by the same
satellite/signal addressing forms as the per-signal accessors (e.g.
[`get_correlator_outputs`](@ref)):

  - `append_correlator_output!(track_state, output)` — one group, one sat, one signal.
  - `append_correlator_output!(track_state, output, prn)` — one group, one signal.
  - `append_correlator_output!(track_state, output, group, prn)` — multi-group, one signal.
  - `append_correlator_output!(track_state, output, group, prn, sig)` — per-signal.

The signal-level method (`append_correlator_output!(::TrackedSignal, output)`)
carries the contract; see it and [External correlator producers](@ref).
"""
append_correlator_output!(s::TrackState, output::CorrelatorOutput, id...) =
    (append_correlator_output!(get_sat_state(s, id...), output); s)
append_correlator_output!(
    s::TrackState{<:SignalGroups},
    output::CorrelatorOutput,
    group::Union{Symbol,Integer,Val},
    sat_id,
    sig::_SignalSelector,
) = (append_correlator_output!(get_sat_state(s, group, sat_id), output, sig); s)

"""
$(SIGNATURES)

Append an externally built [`NoiseObservation`](@ref) to the addressed
**signal**'s noise estimator and return `track_state`:

  - `append_noise_observation!(track_state, obs)` — single-signal `TrackState`.
  - `append_noise_observation!(track_state, obs, :GPSL1CA)` — explicit signal id.
  - `append_noise_observation!(track_state, obs, GPSL1CA)` — or the signal type /
    an instance of it, which is what a caller usually has to hand.

This is the hardware/FPGA fill path, used alongside
[`append_correlator_output!`](@ref) before
[`estimate_dopplers_and_filter_prompt!`](@ref) folds; see the estimator-level
method for how the two differ. Keyed per signal (not per band) because the noise
floor is post-correlation; see [`AbstractNoiseEstimator`](@ref).

The signal must have a noise estimator; it has one whenever its C/N₀ estimator
reads a density (see [`requires_noise_density`](@ref)), or whenever you declared
one through `TrackState`'s `noise_estimators` keyword.
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

# Look up a signal's noise estimator, with an error naming what is configured.
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

# Resolve the index of the addressed signal within a sat's signals tuple.
# Config-time only (not the hot path), so plain control flow is fine.
_signal_index(signals::Tuple) = length(signals) == 1 ? 1 : _throw_needs_signal_selector()
_signal_index(::Tuple, i::Integer) = Int(i)
function _signal_index(signals::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    idx = findfirst(s -> s.signal isa T, signals)
    isnothing(idx) && throw(ArgumentError("no signal of type $T on this satellite"))
    # Like `_find_signal_by_type`, a type selector must be unambiguous.
    isnothing(findnext(s -> s.signal isa T, signals, idx + 1)) || throw(
        ArgumentError(
            "more than one signal of type $T on this satellite — " *
            "pass an integer index to address one of them.",
        ),
    )
    idx
end

# Rebuild `sat` with the addressed signal's coherent-integration length set to
# `N`, preserving the satellite's concrete type.
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

Set the preferred coherent-integration length, in primary code blocks, for one
signal on one satellite — the `preferred_num_code_blocks_to_integrate` field of
the addressed [`TrackedSignal`](@ref). The actual length is still capped per
integration by the signal's bit/secondary-code period and held at 1 until
bit/secondary sync (see `calc_num_code_blocks_to_integrate`); with the
conventional estimator the loop bandwidths are capped for stability, so no
re-tuning is needed (see [`ConventionalPLLAndDLL`](@ref)).

For data-bearing signals the length must evenly divide the number of code
blocks that form one bit (e.g. a divisor of 20 for GPS L1 C/A, of 10 for GPS
L5I) so integrations stay aligned to bit boundaries; an `ArgumentError` is
thrown otherwise (issue #128). Pilot signals accept any length of at least
one block.

The satellite is addressed exactly like the per-signal accessors
(e.g. [`estimate_cn0`](@ref)) — in particular, the per-signal form always
names the group explicitly, even on a single-group `TrackState`:

```julia
set_preferred_num_code_blocks_to_integrate!(ts, :gps_l5, 1, GPSL5I, 10)  # (group, prn, signal)
set_preferred_num_code_blocks_to_integrate!(ts, :gps_l5, 1, 10)          # single-signal sat
set_preferred_num_code_blocks_to_integrate!(ts, 1, 10)                   # single-group state
set_preferred_num_code_blocks_to_integrate!(ts, 10)                      # 1 group, 1 sat, 1 signal
```

Mutates `track_state` in place and returns it.
"""
function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    sig::_SignalSelector,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state,
        get_sat_states(track_state, group),
        sat_id,
        num_code_blocks,
        sig,
    )
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state,
        get_sat_states(track_state, group),
        sat_id,
        num_code_blocks,
    )
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    sat_id::Integer,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state,
        get_sat_states(track_state),
        sat_id,
        num_code_blocks,
    )
end

# Shared body for the three addressing overloads above.
@inline function _set_preferred_blocks!(
    track_state::TrackState,
    sats,
    sat_id,
    num_code_blocks::Integer,
    sel...,
)
    sats[sat_id] =
        _set_sat_signal_preferred_blocks(sats[sat_id], Int(num_code_blocks), sel...)
    track_state
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    num_code_blocks::Integer,
)
    sats = get_sat_states(track_state)
    sat_id = only(keys(sats))
    sats[sat_id] = _set_sat_signal_preferred_blocks(sats[sat_id], Int(num_code_blocks))
    track_state
end

# Rebuild `sat` with the addressed signal's group delay set to `delay`, or keep it
# where the signal has that delay already.
function _set_sat_signal_group_delay(sat::TrackedSat, delay, sel...)
    isequal(get_group_delay(sat, sel...), convert(Union{Nothing,typeof(1.0s)}, delay)) &&
        return sat
    idx = _signal_index(sat.signals, sel...)
    idx_tuple = ntuple(identity, length(sat.signals))
    new_signals = map(sat.signals, idx_tuple) do s, i
        i == idx ? TrackedSignal(s; group_delay = Some(delay)) : s
    end
    TrackedSat(sat; signals = new_signals)
end

"""
$(SIGNATURES)

Set the group delay of one signal on one satellite, in units of time, or
`nothing` to mark it unknown again (the default). Read it back with
`get_group_delay`.

A group delay is what [signal combining](@ref "Signal combining") needs to refer
a passenger's code discriminator to the driver's code phase: a signal with the
larger group delay arrives later. Only the difference between a passenger's and
the driver's (`signals[1]`) group delay is used, so any common datum cancels.
Zero is not assumed for an unknown delay: a passenger is combined into the code
loop only where both its own and the driver's group delay are set. Setting the
delay a signal has already leaves its satellite as it is.

Addressed like [`set_preferred_num_code_blocks_to_integrate!`](@ref):

```julia
set_group_delay!(ts, :e1, 11, GalileoE1B, 0.0u"ns")  # (group, prn, signal)
set_group_delay!(ts, :e1, 11, 0.0u"ns")              # single-signal sat
set_group_delay!(ts, 11, 0.0u"ns")                   # single-group state
set_group_delay!(ts, 0.0u"ns")                       # 1 group, 1 sat, 1 signal
```

Mutates `track_state` in place and returns it.
"""
function set_group_delay!(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    sig::_SignalSelector,
    delay,
)
    sats = get_sat_states(track_state, group)
    sats[sat_id] = _set_sat_signal_group_delay(sats[sat_id], delay, sig)
    track_state
end

function set_group_delay!(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    delay,
)
    sats = get_sat_states(track_state, group)
    sats[sat_id] = _set_sat_signal_group_delay(sats[sat_id], delay)
    track_state
end

function set_group_delay!(
    track_state::TrackState{<:SignalGroups{1}},
    sat_id::Integer,
    delay,
)
    sats = get_sat_states(track_state)
    sats[sat_id] = _set_sat_signal_group_delay(sats[sat_id], delay)
    track_state
end

function set_group_delay!(track_state::TrackState{<:SignalGroups{1}}, delay)
    sats = get_sat_states(track_state)
    sat_id = only(keys(sats))
    sats[sat_id] = _set_sat_signal_group_delay(sats[sat_id], delay)
    track_state
end

# Per-satellite body of `reset_loop_filters!`. A zeroed previous prompt makes
# `fll_disc` return 0, so the first post-reset FLL update is skipped.
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
satellite) from its current Doppler, giving each a freshly initialized loop
filter. For the conventional PLL/DLL estimator this zeroes the loop-filter
integrators while keeping `carrier_doppler` / `code_doppler` and any
per-satellite bandwidth override, and restarts the carrier loop's staging on
the FLL-assisted PLL (see [`Tracking.FrequencyLockIndicator`](@ref)). Each signal's
`last_fully_integrated_filtered_prompt` is cleared too, so the first FLL update
doesn't span the old integration interval.

Use it after changing a signal's coherent-integration length mid-track (e.g.
via [`set_preferred_num_code_blocks_to_integrate!`](@ref)): the loop filter's
integrator state is not portable across a change in update interval and would
cause a transient that can break lock.

For [`VectorPLLAndDLL`](@ref) the re-seed additionally zeroes the
externally-supplied NCO corrections (`code_freq_update` / `carrier_freq_update`)
and the discriminator accumulators — the converged Doppler already folds in the
last correction, so keeping it would apply it twice. A navigation filter must
therefore **re-issue** its corrections via [`set_code_freq_updates!`](@ref) /
[`set_carrier_freq_updates!`](@ref) after resetting a vector-loop satellite.
The `vt_on` flag and any per-satellite bandwidth override are preserved.

Addressed like the per-signal accessors — no satellite id resets every
satellite in `track_state`; `(group, prn)` or (single-group) `prn` resets one.
Mutates `track_state` in place and returns it. Works for any
[`AbstractDopplerEstimator`](@ref) through its [`init_estimator_state`](@ref) hook.

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
    group::Union{Symbol,Integer,Val},
    sat_id,
)
    sats = get_sat_states(track_state, group)
    sats[sat_id] = _reset_sat_loop_filters(track_state, sats[sat_id])
    track_state
end

reset_loop_filters!(track_state::TrackState{<:SignalGroups{1}}, sat_id) =
    reset_loop_filters!(track_state, 1, sat_id)

# Tuple walk collecting distinct band ids in first-encounter order; unrolled at
# compile time, allocation-free.
@inline _band_keys_in_groups(::Tuple{}, acc::Tuple{Vararg{Symbol}}) = acc
@inline function _band_keys_in_groups(t::Tuple, acc::Tuple{Vararg{Symbol}})
    k = get_band_id(first(t).band)
    new_acc = k in acc ? acc : (acc..., k)
    _band_keys_in_groups(Base.tail(t), new_acc)
end

"""
$(SIGNATURES)

The set of distinct band ids (`GNSSSignals.get_band_id`) used across all
groups in `track_state`, returned as a tuple of Symbols in first-encounter
order. These are the keys a multi-band measurements NamedTuple must use.
Resolves at compile time when the groups type is known.

```julia
ts = TrackState(; signals = (legacy_gps = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
band_keys(ts) == (:L1, :L5)
```
"""
@inline band_keys(track_state::TrackState) =
    _band_keys_in_groups(Tuple(track_state.groups), ())

# The single band shared by all groups, or an error — gates the bare-buffer
# `track(buf, state, fs)` entry point.
@inline function _single_band(track_state::TrackState)
    keys_tuple = band_keys(track_state)
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

# Allocation-free `Set(a) == Set(b)` for short tuples.
@inline function _tuple_sets_equal(a::Tuple, b::Tuple)
    length(a) == length(b) || return false
    all(x -> x in b, a) && all(y -> y in a, b)
end

# Validate a multi-band measurements NamedTuple against the groups: matching keys,
# per-band antenna shapes, and identical durations. Called once per `track`.
@inline function _validate_measurements(
    track_state::TrackState,
    measurements::BandMeasurements,
)
    expected_keys = band_keys(track_state)
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

# Check each group's band measurement has the declared antenna count (Vector → 1,
# Matrix → columns). Tuple recursion, since a `for` loop over heterogeneous
# groups boxes and allocates.
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

# Exact-equality duration check across all measurements. Cross-multiplies
# promoted rates instead of dividing, so mixed-precision rates with identical
# durations compare equal.
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
