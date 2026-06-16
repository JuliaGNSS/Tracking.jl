"""
$(SIGNATURES)

Construct a fresh `TrackState` from a declaration of which signals each
group tracks. A *group* is a set of satellites that share the same
signal-tuple shape (and therefore the same concrete `TrackedSat` type,
which is what gives the hot loop type stability). Each entry in `signals`
is a tuple of `AbstractGNSSSignal` instances; the first signal is the
estimator-driver signal — the one the Doppler estimator uses to update the
satellite-shared carrier and code Doppler (with the conventional PLL/DLL,
it's the signal the discriminator runs on). The group key
(`:modern_gps`, `:legacy_gps`, …) is what `add_satellite!` later refers
to.

```julia
track_state = TrackState(;
    signals = (
        legacy_gps = (GPSL1CA(),),
        modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
        galileo    = (GalileoE1B(),),
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

`TrackState` is parameterized on the per-group `TrackedSat` value type,
which captures the correlator type, post-corr-filter type, and
doppler-estimator-state type — these are frozen at construction.
`add_satellite!` cannot change these types; it can only fill in
satellites of the already-fixed shape. Power users who need non-default
correlator or PCF *types* should construct `TrackedSat`s themselves and
hand them to the `add_satellite!(track_state, group, sat)` overload.
"""
function TrackState(;
    signal::Maybe{AbstractGNSSSignal} = nothing,
    signals = nothing,
    doppler_estimator::Maybe{AbstractDopplerEstimator} = nothing,
    num_ants::NumAnts = NumAnts(1),
)
    if isnothing(signal) && isnothing(signals)
        throw(ArgumentError(
            "TrackState requires either `signal = <AbstractGNSSSignal>` or " *
            "`signals = (...)`. See the docstring for examples.",
        ))
    end
    if !isnothing(signal) && !isnothing(signals)
        throw(ArgumentError(
            "Pass either `signal` (singular, one AbstractGNSSSignal) or " *
            "`signals` (plural, a NamedTuple of signal tuples) — not both.",
        ))
    end
    sig_groups_nt = isnothing(signal) ?
        _normalize_signal_groups(signals) :
        (default = (signal,),)
    # Default estimator: the auto-bandwidth `ConventionalAssistedPLLAndDLL`.
    # Each group's satellites are then seeded (via `init_estimator_state`)
    # with the loop bandwidth recommended for that group's own driver signal,
    # so every group gets the right bandwidth without a cross-group compromise.
    estimator = isnothing(doppler_estimator) ?
        ConventionalAssistedPLLAndDLL() : doppler_estimator
    # Each entry of `sig_groups_nt` is either:
    #   - a bare `Tuple{Vararg{AbstractGNSSSignal}}` (the common case,
    #     uses the constructor's `num_ants` kwarg); or
    #   - a pre-built `SignalGroup` instance (carries its own band /
    #     num_ants — used when the user wants per-band overrides).
    groups = map(sig_groups_nt) do entry
        _normalize_group_entry(entry, estimator, num_ants)
    end
    _validate_same_band_num_ants(groups)
    TrackState(groups, estimator)
end

# Bare tuple of AbstractGNSSSignal → single :default capability NamedTuple.
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
    sats = Dictionary{Int, typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, sig_tuple, num_ants)
end

# Pre-built SignalGroup → pass through, but if its `satellites` dict
# value type doesn't match the estimator the user passed to TrackState,
# rebuild the template so the slot type lines up. The common case where
# the user built the SignalGroup with `SignalGroup((sigs,); num_ants =
# ..., doppler_estimator = same)` then it just passes through.
@inline function _normalize_group_entry(
    g::SignalGroup,
    doppler_estimator::AbstractDopplerEstimator,
    _num_ants::NumAnts,
)
    # If the existing slot type already matches the estimator, keep it.
    # Otherwise rebuild the empty dict with a fresh template that uses
    # the TrackState's estimator. The user's `num_ants` on the
    # SignalGroup wins — the TrackState's `num_ants` kwarg is only the
    # default for bare-tuple entries.
    sats = g.satellites
    if !isempty(sats)
        # Pre-populated SignalGroup — the slot type is fixed by the
        # existing sats, so it must already match the estimator.
        _assert_doppler_estimator_types_match(sats, doppler_estimator)
        return g
    end
    template = _make_template_tracked_sat(g.signals, doppler_estimator, g.num_ants)
    if typeof(template) === eltype(sats)
        return g
    end
    new_sats = Dictionary{Int, typeof(template)}(Int[], typeof(template)[])
    SignalGroup(g.band, new_sats, g.signals, g.num_ants)
end

# Same-band groups must declare identical `num_ants`. Two groups on the
# same physical band can't be sampled by front-ends with different
# antenna counts. Walked over the concrete-typed groups tuple at
# construction; O(num_groups²) but folded at compile time when the
# groups type is known.
@inline function _validate_same_band_num_ants(groups::NamedTuple)
    _check_same_band_num_ants(Tuple(groups), ())
end

@inline _check_same_band_num_ants(::Tuple{}, ::Tuple) = nothing
@inline function _check_same_band_num_ants(t::Tuple, seen::Tuple)
    g = first(t)
    bk = band_key(g.band)
    for (sk, sna) in seen
        if sk === bk && sna !== g.num_ants
            throw(ArgumentError(string(
                "Two groups on band `:", bk,
                "` declare different antenna counts (",
                sna, " vs ", g.num_ants,
                "). Groups sharing a band must share NumAnts ",
                "— they're sampled by the same front-end.",
            )))
        end
    end
    _check_same_band_num_ants(Base.tail(t), (seen..., (bk, g.num_ants)))
end

"""
$(SIGNATURES)

Internal helper: build a template `TrackedSat` for a capability declared
as `signal_tuple = (GPSL1C_P(), GPSL1C_D(), GPSL1CA())`. The template's
data is meaningless (PRN 0, zero Doppler); only its concrete type is
used to fix the dictionary value type at TrackState construction.
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
)
    # `signal` is implied by each sat's `signals[1].signal` in the new design;
    # the positional argument is kept for backward-compatible construction but
    # is otherwise unused — the default estimator auto-sizes each sat's loop
    # bandwidth from its own driver signal.
    sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(sats_dict, doppler_estimator)
    groups = (default = _signal_group_from_dict(sats_dict),)
    TrackState(groups, doppler_estimator)
end

function TrackState(
    tracked_sats::Dictionary{<:Any,<:TrackedSat};
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    _assert_doppler_estimator_types_match(tracked_sats, doppler_estimator)
    groups = (default = _signal_group_from_dict(tracked_sats),)
    TrackState(groups, doppler_estimator)
end

# Verify every sat in `dict` has a `doppler_estimator_state` matching
# what `estimator` would produce. The check is type-only so it has no
# runtime cost in the typed-correct case. Throws an ArgumentError that
# names the mismatched concrete types if the user built sats with a
# different estimator than the one configured on the TrackState.
@inline function _assert_doppler_estimator_types_match(
    dict::Dictionary{<:Any,<:TrackedSat},
    estimator::AbstractDopplerEstimator,
)
    isempty(dict) && return nothing
    sat = first(dict.values)
    expected_state_type = typeof(init_estimator_state(estimator, sat))
    actual_state_type = typeof(sat.doppler_estimator_state)
    expected_state_type === actual_state_type && return nothing
    throw(ArgumentError(string(
        "TrackedSat has doppler_estimator_state of type ",
        actual_state_type,
        ", but the configured doppler_estimator (",
        typeof(estimator),
        ") would produce ",
        expected_state_type,
        ". Construct each TrackedSat with `doppler_estimator = <same instance>` ",
        "as the one passed here.",
    )))
end

function TrackState(
    satellites::SatelliteDicts;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    foreach(
        d -> _assert_doppler_estimator_types_match(d, doppler_estimator),
        Tuple(satellites),
    )
    groups = map(_signal_group_from_dict, satellites)
    TrackState(groups, doppler_estimator)
end

# Copy-with-overrides constructor: rebuild a TrackState reusing the
# input's groups / estimator unless an override is supplied. The
# overrides are constrained to the input's concrete `G` / `DE` types so
# the result's type parameters are preserved.
function TrackState(
    track_state::TrackState{G,DE};
    groups::Maybe{G} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    TrackState{G,DE}(
        isnothing(groups) ? track_state.groups : groups,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
    )
end

# Build a SignalGroup from a non-empty satellites dictionary by recovering
# the signal-instance tuple, band, and antenna count from any sat. Used by
# the positional `TrackState(signal, sats)` / `TrackState(satellites)`
# constructors. Empty dicts can't be recovered this way (no sats to inspect)
# — requires at least one sat in the dict.
@inline function _signal_group_from_dict(dict::Dictionary{<:Any,<:TrackedSat})
    isempty(dict) && throw(ArgumentError(
        "Cannot recover the signal-instance tuple from an empty " *
        "satellites dictionary. Use the `TrackState(; signals = ...)` " *
        "constructor to declare signal groups before populating sats.",
    ))
    sat = first(dict.values)
    sig_tuple = map(s -> s.signal, sat.signals)
    band = get_band(first(sig_tuple))
    num_ants = NumAnts(get_num_ants(sat))
    SignalGroup(band, dict, sig_tuple, num_ants)
end

# Immutable reset — the first copy `track` makes of the caller's live
# `TrackState`. Detaches the key set (`Indices`) as well as the values
# (`_detach_groups_slot_vectors`, #123) so that a later
# `add_satellite!`/`remove_satellite!` on the returned state (or on `track`'s
# output, which derives from it) cannot corrupt the input's key set. The
# per-iteration loop steps inside `track` reuse this already-detached key set
# via the cheaper, key-sharing `_copy_groups_slot_vectors`.
function reset_start_sample_and_bit_buffer(track_state::TrackState)
    new_groups = _detach_groups_slot_vectors(track_state.groups)
    reset_start_sample_and_bit_buffer!(new_groups)
    TrackState(track_state; groups = new_groups)
end

function reset_start_sample_and_bit_buffer!(track_state::TrackState)
    reset_start_sample_and_bit_buffer!(track_state.groups)
    return track_state
end

# Loop-termination helper for `track`/`track!` lives in `track.jl` as
# `_all_groups_reached_end` — it iterates the per-band measurement
# lengths so each group can terminate against its own band's chunk.

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

Return the first signal of the given capability — useful when the caller
needs the signal instance (for `gen_code`, frequency lookups, …) but
doesn't already have a sat in hand. The dictionary's value type carries
the signal type, so this resolves at compile time when `group_idx` is a
literal `Symbol` / `Integer`.

For a single-capability `TrackState` the index can be omitted.
"""
function get_signal(
    track_state::TrackState{<:SignalGroups{N}},
    group_idx::Union{Symbol,Integer,Val},
) where {N}
    # Read the signal off the group's declared signal tuple (not off a
    # tracked sat) so this also works on declared-but-unpopulated groups.
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

function get_sat_state(
    track_state::TrackState{<:SignalGroups{1}},
    sat_identifier,
)
    get_sat_state(track_state, 1, sat_identifier)
end

function get_sat_state(track_state::TrackState{<:SignalGroups{1}})
    only(get_sat_states(track_state, 1))
end

# `estimate_cn0` follows the same dispatch ladder as the per-signal
# accessors — its `TrackState` overloads are generated alongside them
# in the `@eval` loop further down (search for `:estimate_cn0`).

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
    track_state::TrackState{G,DE},
    group_idx::Union{Symbol,Integer},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    new_sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(new_sats_dict, track_state.doppler_estimator)
    new_estimator =
        update_estimator_on_handoff(track_state.doppler_estimator, new_sats_dict)
    groups = track_state.groups
    g = groups[group_idx]
    _assert_sats_match_slot_type(g, new_sats_dict, group_idx)
    new_group = SignalGroup(g; satellites = merge(g.satellites, new_sats_dict))
    new_groups = @set groups[group_idx] = new_group
    TrackState{G,DE}(new_groups, new_estimator)
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

`group` defaults to `:default` (matching the single-group construction
shortcut). For multi-group `TrackState`s, pass the group key explicitly.
If a satellite with the same `prn` already exists in that group's
dictionary, it is overwritten.

The satellite dictionary is mutated in place, but callers should keep
using the *returned* `TrackState`: when the configured estimator's
[`update_estimator_on_handoff`](@ref) returns a rebuilt estimator (rather
than mutating in place), the rebuilt estimator is carried by the returned
`TrackState` — `track_state.doppler_estimator` cannot be replaced in
place because `TrackState` is immutable. For estimators that update in
place (including the default conventional ones), the very same
`track_state` comes back.

```julia
track_state = TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
add_satellite!(track_state;
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
    group::Symbol = :default,
    kwargs...,
)
    sat = _make_default_tracked_sat_for_group(track_state, group; prn, kwargs...)
    add_satellite!(track_state, group, sat)
end

"""
$(SIGNATURES)

In-place add (or replace) with a pre-built [`TrackedSat`](@ref) — the
escape hatch for power users who need non-default correlator or
post-corr-filter types. The sat's type must match the group's slot type
already fixed at [`TrackState`](@ref) construction; passing a sat of the
wrong type errors at dispatch time.

Like the keyword form, the returned `TrackState` carries the estimator
returned by [`update_estimator_on_handoff`](@ref) — keep using the
return value.
"""
function add_satellite!(
    track_state::TrackState{G,DE},
    group::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    _assert_sat_matches_slot_type(track_state, group, sat)
    insert_or_set!(_dict_for_group(track_state, group), sat.prn, sat)
    new_estimator = update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((sat.prn => sat,)),
    )
    # `TrackState` is immutable, so a rebuilt estimator can only be honored
    # through the return value. Estimators that update in place return the
    # identical object (the contract guarantees the concrete type either
    # way), and then the input `track_state` is handed back unchanged.
    new_estimator === track_state.doppler_estimator && return track_state
    TrackState{G,DE}(track_state.groups, new_estimator)
end

# Verify that `sat` has exactly the concrete type the capability's
# dictionary slot expects. Throws a clear ArgumentError that names the
# mismatching types if not; called from the escape-hatch overloads so
# the user gets a useful message before Dictionaries.jl's `set!` raises
# a deep MethodError about `convert`.
@inline function _assert_sat_matches_slot_type(
    track_state::TrackState, group::Symbol, sat::TrackedSat,
)
    SlotT = eltype(track_state.groups[group].satellites)
    typeof(sat) === SlotT && return nothing
    throw(ArgumentError(string(
        "TrackedSat type does not match the `:", group, "` group's slot. ",
        "Got: ", typeof(sat),
        ". Expected: ", SlotT,
        ". The slot type is fixed at TrackState construction; rebuild the ",
        "sat with the matching correlator, post_corr_filter, and ",
        "doppler_estimator types.",
    )))
end

# Dictionary-level variant used by `merge_sats`: the incoming dict's
# value type must match the group's slot type exactly (not just the
# estimator-state type) so correlator / PCF / signal-shape mismatches
# also get the curated error instead of a confusing MethodError.
@inline function _assert_sats_match_slot_type(
    g::SignalGroup, new_sats_dict::Dictionary{<:Any,<:TrackedSat}, group_idx,
)
    SlotT = eltype(g.satellites)
    T = eltype(new_sats_dict)
    T === SlotT && return nothing
    throw(ArgumentError(string(
        "TrackedSat type does not match the `", group_idx, "` group's slot. ",
        "Got: ", T,
        ". Expected: ", SlotT,
        ". The slot type is fixed at TrackState construction; rebuild the ",
        "sats with the matching signals, correlator, post_corr_filter, and ",
        "doppler_estimator types.",
    )))
end

# Dictionaries.jl: `insert!` errors on existing key; `set!` overwrites.
# We want overwrite semantics, matching `merge_sats`.
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
    group::Symbol = :default,
    kwargs...,
)
    sat = _make_default_tracked_sat_for_group(track_state, group; prn, kwargs...)
    add_satellite(track_state, group, sat)
end

function add_satellite(
    track_state::TrackState{G,DE},
    group::Symbol,
    sat::TrackedSat,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
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
    TrackState{G,DE}(new_groups, new_estimator)
end

"""
$(SIGNATURES)

Remove a satellite from `track_state` in place. Errors if no satellite
with the given `prn` exists in the named group (matches Dictionaries.jl's
`delete!` semantics).

```julia
remove_satellite!(track_state; prn = 11, group = :modern_gps)
```

`group` defaults to `:default` for single-group TrackStates.
Returns `track_state` unchanged (the dictionary is mutated in place).
"""
function remove_satellite!(
    track_state::TrackState;
    prn::Int,
    group::Symbol = :default,
)
    delete!(_dict_for_group(track_state, group), prn)
    track_state
end

"""
$(SIGNATURES)

Immutable variant of [`remove_satellite!`](@ref). Returns a new
[`TrackState`](@ref) with the satellite removed; the input is left
unchanged. Errors if no satellite with the given `prn` exists.
"""
function remove_satellite(
    track_state::TrackState{G,DE};
    prn::Int,
    group::Symbol = :default,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    dict = _dict_for_group(track_state, group)
    haskey(dict, prn) || throw(KeyError(prn))
    # Copy-then-delete preserves the dict's concrete value type, which the
    # `map`/`filter` round-trip would otherwise leave under-specified to
    # `Dictionary{<:Any,<:TrackedSat}` from the type system's perspective.
    new_dict = copy(dict)
    delete!(new_dict, prn)
    groups = track_state.groups
    g = groups[group]
    new_group = SignalGroup(g; satellites = new_dict)
    new_groups = @set groups[group] = new_group
    TrackState{G,DE}(new_groups, track_state.doppler_estimator)
end

# Compile-time dispatch helper: hand back the dictionary slot for the
# given group key. Bounds and existence are checked at TrackState
# construction time (the NamedTuple only contains declared keys), so an
# unknown `group` here triggers the standard NamedTuple KeyErrors.
@inline _dict_for_group(track_state::TrackState, group::Symbol) =
    track_state.groups[group].satellites

# Build a default-correlator, default-PCF TrackedSat whose
# signal-tuple shape matches the group's slot in `track_state`. Reads the
# signal-instance tuple and antenna count straight off the SignalGroup.
# Carries the acquisition-handoff kwarg defaults for both `add_satellite!`
# and `add_satellite`, which forward their kwargs here verbatim.
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
        g.signals, prn, code_phase, carrier_doppler;
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
get_carrier_doppler(s::TrackState, id...) = get_carrier_doppler(get_sat_state(s, id...))
get_signal_start_sample(s::TrackState, id...) =
    get_signal_start_sample(get_sat_state(s, id...))

# Per-signal accessors. Each takes a trailing signal selector
# (`Integer` index or `Type{<:AbstractGNSSSignal}`) to disambiguate
# between the signals tracked on a multi-signal sat. Without a selector
# they fall back through `get_sat_state` → `only(sat.signals)` and so
# are only valid on single-signal sats.
#
# Addressing forms (using `get_correlator` as the example — applies to
# every accessor in this block):
#   * `get_correlator(track_state)` — 1 group, 1 sat, 1 signal.
#   * `get_correlator(track_state, prn)` — 1 group, 1 signal.
#   * `get_correlator(track_state, group, prn)` — multi-group, 1 signal.
#   * `get_correlator(track_state, group, prn, sig)` — per-signal.
#
# The per-signal form always names the group explicitly, even on a
# single-group TrackState, to keep the API unambiguous. Use `:default`
# (or `1`) as the group key in the single-group case.
const _SignalSelector = Union{Integer,Type{<:AbstractGNSSSignal}}

# Each accessor in the loop below gets two `TrackState` overloads:
#   * `(s, id...)` — varargs forward to `get_sat_state`, covering the
#     no-arg, prn-only, and (group, prn) shapes via that function's own
#     dispatch ladder. Sat-level fallback when no signal selector is given.
#   * `(s, group, prn, sig)` — the per-signal form. `sig` is `Integer` or
#     `Type{<:AbstractGNSSSignal}`; the lookup goes via the sat-level
#     accessor (which routes through `_find_signal`).
#
# Generated via `@eval` at module load — methods bake into the precompile
# image, indistinguishable from hand-written ones at runtime.
for fn in (
    :get_integrated_samples, :get_correlator,
    :get_last_fully_integrated_correlator, :get_last_fully_integrated_filtered_prompt,
    :get_filtered_prompts,
    :get_post_corr_filter, :get_cn0_estimator, :get_bit_buffer, :get_bits,
    :get_soft_bits,
    :get_num_bits, :has_bit_or_secondary_code_been_found, :estimate_cn0,
    :get_preferred_num_code_blocks_to_integrate,
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

# Resolve the index of the addressed signal within a sat's signals tuple.
# Config-time only (not the hot path), so plain control flow is fine.
_signal_index(signals::Tuple) =
    length(signals) == 1 ? 1 : _throw_needs_signal_selector()
_signal_index(::Tuple, i::Integer) = Int(i)
function _signal_index(signals::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    idx = findfirst(s -> s.signal isa T, signals)
    isnothing(idx) && throw(ArgumentError("no signal of type $T on this satellite"))
    idx
end

# Rebuild `sat` with the addressed signal's coherent-integration length set
# to `N`; the other signals are left untouched (types unchanged, so the
# satellite's concrete type is preserved).
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
conventional estimator the loop bandwidth auto-scales by `1/N` so the loop
stays stable at any length.

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
        track_state, get_sat_states(track_state, group), sat_id, num_code_blocks, sig,
    )
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state, get_sat_states(track_state, group), sat_id, num_code_blocks,
    )
end

function set_preferred_num_code_blocks_to_integrate!(
    track_state::TrackState{<:SignalGroups{1}},
    sat_id::Integer,
    num_code_blocks::Integer,
)
    _set_preferred_blocks!(
        track_state, get_sat_states(track_state), sat_id, num_code_blocks,
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

# Re-seed one satellite's Doppler-estimator state from its current Doppler via
# the estimator's `_reset_estimator_state` hook (a fresh, zeroed loop filter
# for the conventional estimator, keeping any per-sat bandwidth override),
# preserving `carrier_doppler` / `code_doppler`. Each signal's
# `last_fully_integrated_filtered_prompt` is cleared too: the FLL
# discriminator measures the prompt rotation since the previous integration,
# and after a cadence change that previous prompt belongs to the old
# integration interval — `fll_disc` returns 0 for a zeroed previous prompt,
# so the first post-reset update skips the stale measurement.
@inline function _reset_sat_loop_filters(track_state::TrackState, sat::TrackedSat)
    new_signals = map(
        s -> TrackedSignal(s; last_fully_integrated_filtered_prompt = complex(0.0, 0.0)),
        sat.signals,
    )
    TrackedSat(
        sat;
        signals = new_signals,
        doppler_estimator_state =
            _reset_estimator_state(track_state.doppler_estimator, sat),
    )
end

"""
$(SIGNATURES)

Re-seed the Doppler-estimator state of every satellite (or one addressed
satellite) from its current Doppler, giving each a freshly initialized loop
filter. For the conventional PLL/DLL estimator this zeroes the carrier and code
loop-filter integrators while preserving the converged `carrier_doppler` /
`code_doppler` — and any per-satellite loop-bandwidth override carried on the
`SatConventionalPLLAndDLL` state — so the loop continues from the
converged frequency with a clean filter. Each signal's
`last_fully_integrated_filtered_prompt` is cleared as well, so the first
FLL update after the reset doesn't measure a prompt rotation that spans the
old integration interval.

This is the recommended handoff when a signal's coherent-integration length
changes mid-track — e.g. promoting GPS L5I from 1 ms to 10 ms via
[`set_preferred_num_code_blocks_to_integrate!`](@ref). The bilinear loop
filter's integrator state is not portable across the change in update interval
(`Δt` grows by the integration factor), so resetting it avoids a transient that
can drag the loop out of lock; the converged Doppler is the right seed for the
new, longer integration.

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

# Recursive tuple walk that folds the set of distinct `band_key`s across
# a `groups` tuple. Returns an `NTuple{N,Symbol}` of unique keys, in
# first-encounter order. Concrete-typed input → fully unrolled at compile
# time, no allocation.
@inline _band_keys_in_groups(::Tuple{}, acc::Tuple{Vararg{Symbol}}) = acc
@inline function _band_keys_in_groups(t::Tuple, acc::Tuple{Vararg{Symbol}})
    k = band_key(first(t).band)
    new_acc = k in acc ? acc : (acc..., k)
    _band_keys_in_groups(Base.tail(t), new_acc)
end

"""
$(SIGNATURES)

The set of distinct band keys used across all groups in `track_state`,
returned as a tuple of Symbols in first-encounter order. Resolves at
compile time when the groups type is known.

```julia
ts = TrackState(; signals = (legacy_gps = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
band_keys(ts) == (:l1, :l5)
```
"""
@inline band_keys(track_state::TrackState) =
    _band_keys_in_groups(Tuple(track_state.groups), ())

# Single-band shortcut helper: returns the unique band instance shared by
# all groups when there is exactly one distinct band. Errors otherwise —
# this is what gates whether the bare-buffer `track(buf, state, fs)`
# entry point can route the measurement to a single auto-keyed
# `BandMeasurements` NamedTuple.
@inline function _single_band(track_state::TrackState)
    keys_tuple = band_keys(track_state)
    if length(keys_tuple) != 1
        throw(ArgumentError(string(
            "Bare-buffer `track`/`track!` requires a single-band TrackState, ",
            "but this TrackState spans bands ", keys_tuple,
            ". Pass a NamedTuple of `BandMeasurement`s instead, ",
            "one per band.",
        )))
    end
    # All groups share one band — pull it off the first group.
    first(track_state.groups).band
end

# Two-way membership check on two short symbol tuples. Equivalent to
# `Set(a) == Set(b)` but allocation-free: tuple `in` is unrolled and
# the comparison runs in O(length(a)·length(b)) — fine for band tuples
# which are 1..a few entries.
@inline function _tuple_sets_equal(a::Tuple, b::Tuple)
    length(a) == length(b) || return false
    all(x -> x in b, a) && all(y -> y in a, b)
end

# Validate a multi-band measurements NamedTuple against the TrackState's
# groups: keys must match exactly, antenna shape per band must match,
# and all observation durations must be identical (no tolerance). Called
# once at the top of `track` / `track!` — O(num_bands), irrelevant next
# to the inner loop.
@inline function _validate_measurements(
    track_state::TrackState, measurements::BandMeasurements,
)
    expected_keys = band_keys(track_state)
    got_keys = keys(measurements)
    if !_tuple_sets_equal(expected_keys, got_keys)
        throw(ArgumentError(string(
            "BandMeasurement keys do not match the TrackState's band set. ",
            "Expected: ", expected_keys,
            ". Got: ", got_keys, ".",
        )))
    end
    _validate_antenna_shapes(track_state, measurements)
    _validate_equal_durations(measurements)
    return nothing
end

# For each group, check the measurement at its band has the antenna
# shape the group declares. Vector → 1 antenna; Matrix → cols = num_ants.
# A `for g in track_state.groups` loop boxes when `groups` is a
# heterogeneous tuple (multi-band TrackStates), allocating 2 entries per
# call. Walk via tuple recursion so each step has concrete types and
# inlines cleanly.
@inline function _validate_antenna_shapes(
    track_state::TrackState, measurements::BandMeasurements,
)
    _validate_antenna_shapes_walk(Tuple(track_state.groups), measurements)
end

@inline _validate_antenna_shapes_walk(::Tuple{}, ::BandMeasurements) = nothing
@inline function _validate_antenna_shapes_walk(groups::Tuple, measurements::BandMeasurements)
    g = first(groups)
    m = measurements[band_key(g.band)]
    _assert_antenna_shape(g, m)
    _validate_antenna_shapes_walk(Base.tail(groups), measurements)
end

@inline function _assert_antenna_shape(g::SignalGroup{B,S,Sigs,NumAnts{M}}, m::BandMeasurement) where {B,S,Sigs,M}
    cols = m.samples isa AbstractMatrix ? size(m.samples, 2) : 1
    cols == M && return nothing
    throw(ArgumentError(string(
        "Antenna shape mismatch for band `:", band_key(g.band),
        "`. Group declares NumAnts(", M, ") but the measurement's ",
        "`samples` has ", cols, " column(s).",
    )))
end

# Exact-equality duration check across all measurements — no tolerance.
# Cross-multiply instead of dividing so that mixed-precision sampling
# rates (e.g. Float32 vs Float64) with identical real durations compare
# equal; the rates are promoted first so neither product rounds in a
# narrower type than the comparison.
@inline function _validate_equal_durations(measurements::BandMeasurements)
    ms = Tuple(measurements)
    isempty(ms) && return nothing
    first_m = first(ms)
    ref_num_samples = get_num_samples(first_m)
    for m in Base.tail(ms)
        ref_fs, fs = promote(first_m.sampling_frequency, m.sampling_frequency)
        get_num_samples(m) * ref_fs == ref_num_samples * fs && continue
        throw(ArgumentError(string(
            "BandMeasurement durations must be exactly equal across bands. ",
            "Got ", uconvert(s, get_num_samples(m) / m.sampling_frequency),
            " and ", uconvert(s, ref_num_samples / first_m.sampling_frequency), ". ",
            "Check `num_samples / sampling_frequency` for each band.",
        )))
    end
    return nothing
end
