"""
$(SIGNATURES)

Construct a fresh `TrackState` from a declaration of which signals each
group tracks. A *group* is a set of satellites that share the same
signal-tuple shape (and therefore the same concrete `TrackedSat` type,
which is what gives the hot loop type stability). Each entry in `signals`
is a tuple of `AbstractGNSSSignal` instances; the first signal is the
Doppler source (its correlator drives the PLL/DLL). The group key
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
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
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
    # For each capability, build a template TrackedSat to capture the concrete
    # dict value type, build an empty dict of that type, then wrap with the
    # signal-tuple, band, and antenna count into a SignalGroup.
    groups = map(sig_groups_nt) do sig_tuple
        template = _make_template_tracked_sat(sig_tuple, doppler_estimator, num_ants)
        sats = Dictionary{Int, typeof(template)}(Int[], typeof(template)[])
        band = get_band(first(sig_tuple))
        SignalGroup(band, sats, sig_tuple, num_ants)
    end
    TrackState(groups, doppler_estimator)
end

# Bare tuple of AbstractGNSSSignal → single :default capability NamedTuple.
@inline _normalize_signal_groups(signals::Tuple{Vararg{AbstractGNSSSignal}}) =
    (default = signals,)
# NamedTuple of signal tuples → pass through unchanged.
@inline _normalize_signal_groups(signals::NamedTuple) = signals

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
    tracked_signals = map(signal_tuple) do sig
        TrackedSignal(
            sig;
            num_ants,
            correlator = get_default_correlator(sig, num_ants),
        )
    end
    # signals[1] is the Doppler source. Use it to infer a sensible zero
    # code-doppler value (the math itself doesn't matter for the template).
    bare = TrackedSat(
        0,
        0.0,
        0.0Hz,
        0.0,
        0.0Hz,
        1,
        tracked_signals,
        nothing,
    )
    de_state = init_estimator_state(doppler_estimator, bare)
    TrackedSat(
        bare.prn, bare.code_phase, bare.code_doppler,
        bare.carrier_phase, bare.carrier_doppler,
        bare.signal_start_sample, bare.signals, de_state,
    )
end

function TrackState(
    system::AbstractGNSSSignal,
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}};
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    # `system` is implied by each sat's `signals[1].signal` in the new
    # design; the positional argument is kept for backward-compatible
    # construction but is otherwise unused.
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
    groups = map(_signal_group_from_dict, satellites)
    TrackState(groups, doppler_estimator)
end

# Internal copy-constructor used by the immutable downconvert/estimator path
# to swap in a fresh `satellites` NamedTuple (slot-vector-detached). Rebuilds
# each group around the new per-group dictionary while preserving band,
# signals, and num_ants.
function TrackState(
    track_state::TrackState{G,DE},
    satellites::NamedTuple,
    doppler_estimator::DE,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    new_groups = map(track_state.groups, satellites) do g, sats
        SignalGroup(g; satellites = sats)
    end
    TrackState{G,DE}(new_groups, doppler_estimator)
end

# Be careful when calling this.
# It might lead to types that are inferred at runtime?!
# Tested with 1.11.6
function TrackState(
    track_state::TrackState{G,DE};
    satellites::Maybe{NamedTuple} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    new_groups = isnothing(satellites) ? track_state.groups :
        map(track_state.groups, satellites) do g, sats
            SignalGroup(g; satellites = sats)
        end
    TrackState{G,DE}(
        new_groups,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
    )
end

# Build a SignalGroup from a non-empty satellites dictionary by recovering
# the signal-instance tuple, band, and antenna count from any sat. Used by
# the legacy `TrackState(system, sats)` / `TrackState(satellites)`
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

function reset_start_sample_and_bit_buffer(track_state::TrackState)
    TrackState(
        track_state;
        satellites = reset_start_sample_and_bit_buffer(
            track_state.satellites,
        ),
    )
end

function reset_start_sample_and_bit_buffer!(track_state::TrackState)
    reset_start_sample_and_bit_buffer!(track_state.satellites)
    return track_state
end

function has_integration_reached_signal_end_for_all_satellites(
    track_state::TrackState,
    num_samples::Int,
)
    target = num_samples + 1
    # NamedTuple unwraps to its underlying Tuple via Tuple(...) — concrete and
    # cheap.
    _all_sats_at(Tuple(track_state.satellites), target)
end

# Recursive tuple-walk: each step has fully concrete types, no closure boxes,
# no allocation. Base case is an empty tuple.
@inline _all_sats_at(::Tuple{}, target::Int) = true
@inline function _all_sats_at(t::Tuple, target::Int)
    sats = first(t)
    @inbounds for sat in sats.values
        sat.signal_start_sample == target || return false
    end
    _all_sats_at(Base.tail(t), target)
end

"""
$(SIGNATURES)

Get the per-system satellite dictionary for a specific GNSS system index.
"""
function get_sat_states(
    satellites::SatelliteDicts{N},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    _index_system(satellites, system_idx)
end

@inline _index_system(t::SatelliteDicts, i::Union{Symbol,Integer}) = t[i]
@inline _index_system(t::SatelliteDicts, ::Val{M}) where {M} = t[M]

function get_sat_states(satellites::SatelliteDicts{1})
    get_sat_states(satellites, 1)
end

function get_sat_states(
    track_state::TrackState{<:SignalGroups{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sat_states(track_state.satellites, system_idx)
end

function get_sat_states(track_state::TrackState{<:SignalGroups{1}})
    get_sat_states(track_state.satellites)
end

"""
$(SIGNATURES)

Return the first signal of the given capability — useful when the caller
needs the signal instance (for `gen_code`, frequency lookups, …) but
doesn't already have a sat in hand. The dictionary's value type carries
the signal type, so this resolves at compile time when `system_idx` is a
literal `Symbol` / `Integer`.

For a single-capability `TrackState` the index can be omitted.
"""
function get_system(
    track_state::TrackState{<:SignalGroups{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    sats = get_sat_states(track_state, system_idx)
    get_signal(first(sats.values))
end

function get_system(track_state::TrackState{<:SignalGroups{1}})
    get_system(track_state, 1)
end

function get_sat_state(
    track_state::TrackState{<:SignalGroups{N}},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
) where {N}
    get_sat_state(get_sat_states(track_state, system_idx), sat_identifier)
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

function estimate_cn0(
    track_state::TrackState{<:SignalGroups},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
)
    estimate_cn0(get_sat_state(track_state, system_idx, sat_identifier))
end

function estimate_cn0(track_state::TrackState{<:SignalGroups{1}}, sat_identifier)
    estimate_cn0(track_state, 1, sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:SignalGroups{1}})
    estimate_cn0(get_sat_state(track_state))
end

function merge_sats(
    track_state::TrackState{G,DE},
    system_idx::Union{Symbol,Integer},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
) where {G<:SignalGroups,DE<:AbstractDopplerEstimator}
    new_sats_dict = to_dictionary(tracked_sats)
    _assert_doppler_estimator_types_match(new_sats_dict, track_state.doppler_estimator)
    new_estimator =
        update_estimator_on_handoff(track_state.doppler_estimator, new_sats_dict)
    groups = track_state.groups
    g = groups[system_idx]
    new_group = SignalGroup(g; satellites = merge(g.satellites, new_sats_dict))
    new_groups = @set groups[system_idx] = new_group
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

Returns `track_state` unchanged (the dictionary is mutated in place).

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
    code_phase = 0.0,
    code_doppler::Maybe{typeof(1.0Hz)} = nothing,
    carrier_phase = 0.0,
    carrier_doppler = 0.0Hz,
)
    sat = _make_default_tracked_sat_for_group(
        track_state, group;
        prn, code_phase, code_doppler, carrier_phase, carrier_doppler,
    )
    add_satellite!(track_state, group, sat)
end

"""
$(SIGNATURES)

In-place add (or replace) with a pre-built [`TrackedSat`](@ref) — the
escape hatch for power users who need non-default correlator or
post-corr-filter types. The sat's type must match the group's slot type
already fixed at [`TrackState`](@ref) construction; passing a sat of the
wrong type errors at dispatch time.
"""
function add_satellite!(
    track_state::TrackState,
    group::Symbol,
    sat::TrackedSat,
)
    _assert_sat_matches_slot_type(track_state, group, sat)
    insert_or_set!(_dict_for_group(track_state, group), sat.prn, sat)
    # `update_estimator_on_handoff` is called for its side effect on
    # estimators with cross-sat shared state (the default returns the
    # estimator unchanged). The return value must have the same concrete
    # type as the input — the TrackState is parameterized on the
    # estimator type and `track_state.doppler_estimator` is immutable.
    update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((sat.prn => sat,)),
    )
    track_state
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
    code_phase = 0.0,
    code_doppler::Maybe{typeof(1.0Hz)} = nothing,
    carrier_phase = 0.0,
    carrier_doppler = 0.0Hz,
)
    sat = _make_default_tracked_sat_for_group(
        track_state, group;
        prn, code_phase, code_doppler, carrier_phase, carrier_doppler,
    )
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
    haskey(dict, prn) ||
        throw(KeyError("Dictionary does not contain index: $prn"))
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
function _make_default_tracked_sat_for_group(
    track_state::TrackState,
    group::Symbol;
    prn::Int,
    code_phase,
    code_doppler,
    carrier_phase,
    carrier_doppler,
)
    g = track_state.groups[group]
    sig_tuple = g.signals
    first_signal = first(sig_tuple)
    cd = isnothing(code_doppler) ?
        carrier_doppler * get_code_center_frequency_ratio(first_signal) :
        code_doppler
    num_ants = g.num_ants
    tracked_signals = map(sig_tuple) do sig
        TrackedSignal(
            sig;
            num_ants,
            correlator = get_default_correlator(sig, num_ants),
        )
    end
    bare = TrackedSat(
        prn,
        float(code_phase),
        cd,
        float(carrier_phase) / 2π,
        carrier_doppler,
        1,
        tracked_signals,
        nothing,
    )
    de_state = init_estimator_state(track_state.doppler_estimator, bare)
    TrackedSat(
        bare.prn, bare.code_phase, bare.code_doppler,
        bare.carrier_phase, bare.carrier_doppler,
        bare.signal_start_sample, bare.signals, de_state,
    )
end

# Convenient methods
get_prn(s::TrackState, id...) = get_prn(get_sat_state(s, id...))
get_num_ants(s::TrackState, id...) = get_num_ants(get_sat_state(s, id...))
get_code_phase(s::TrackState, id...) = get_code_phase(get_sat_state(s, id...))
get_code_doppler(s::TrackState, id...) = get_code_doppler(get_sat_state(s, id...))
get_carrier_phase(s::TrackState, id...) = get_carrier_phase(get_sat_state(s, id...))
get_carrier_doppler(s::TrackState, id...) = get_carrier_doppler(get_sat_state(s, id...))
get_integrated_samples(s::TrackState, id...) =
    get_integrated_samples(get_sat_state(s, id...))
get_signal_start_sample(s::TrackState, id...) =
    get_signal_start_sample(get_sat_state(s, id...))
get_correlator(s::TrackState, id...) = get_correlator(get_sat_state(s, id...))
get_last_fully_integrated_correlator(s::TrackState, id...) =
    get_last_fully_integrated_correlator(get_sat_state(s, id...))
get_last_fully_integrated_filtered_prompt(s::TrackState, id...) =
    get_last_fully_integrated_filtered_prompt(get_sat_state(s, id...))
get_post_corr_filter(s::TrackState, id...) = get_post_corr_filter(get_sat_state(s, id...))
get_cn0_estimator(s::TrackState, id...) = get_cn0_estimator(get_sat_state(s, id...))
get_bit_buffer(s::TrackState, id...) = get_bit_buffer(get_sat_state(s, id...))
get_bits(s::TrackState, id...) = get_bits(get_sat_state(s, id...))
get_num_bits(s::TrackState, id...) = get_num_bits(get_sat_state(s, id...))
has_bit_or_secondary_code_been_found(s::TrackState, id...) =
    has_bit_or_secondary_code_been_found(get_sat_state(s, id...))
