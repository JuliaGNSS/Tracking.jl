"""
$(SIGNATURES)

Construct a fresh `TrackState` from a declaration of which signals each
capability (signal-set) tracks. Each entry in `signals` is a tuple of
`AbstractGNSSSignal` instances; the first signal is the Doppler source
(its correlator drives the PLL/DLL). The capability key (`:modern_gps`,
`:legacy_gps`, …) is what `add_satellite!` later refers to.

```julia
track_state = TrackState(;
    signals = (
        legacy_gps = (GPSL1CA(),),
        modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
        galileo    = (GalileoE1B(),),
    ),
)
```

The simplified single-capability form `signals = (GPSL1CA(),)` desugars
to a capability called `:default`. `add_satellite!` may then omit the
`capability=` keyword.

`TrackState` is parameterized on the per-capability `TrackedSat` value
type, which captures the correlator type, post-corr-filter type, and
doppler-estimator-state type — these are frozen at construction.
`add_satellite!` cannot change these types; it can only fill in
satellites of the already-fixed shape. Power users who need non-default
correlator or PCF *types* should construct `TrackedSat`s themselves and
hand them to the `add_satellite!(track_state, capability, sat)`
overload.
"""
function TrackState(;
    signals,
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    num_ants::NumAnts = NumAnts(1),
)
    # Normalize the input: a bare signal tuple gets wrapped under the
    # `:default` capability key.
    signal_groups = _normalize_signal_groups(signals)
    # For each capability, build a template TrackedSat to capture the
    # concrete dict value type, then create an empty dict of that type.
    sat_dicts = map(signal_groups) do sig_tuple
        template = _make_template_tracked_sat(sig_tuple, doppler_estimator, num_ants)
        Dictionary{Int, typeof(template)}(Int[], typeof(template)[])
    end
    TrackState(sat_dicts, signal_groups, doppler_estimator)
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
    reseeded = map(sat -> _reseed_doppler_estimator_state(sat, doppler_estimator), sats_dict)
    satellites = (reseeded,)
    signal_groups = _signal_groups_from_satellites(satellites)
    TrackState(satellites, signal_groups, doppler_estimator)
end

function TrackState(
    tracked_sats::Dictionary{<:Any,<:TrackedSat};
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    reseeded =
        map(sat -> _reseed_doppler_estimator_state(sat, doppler_estimator), tracked_sats)
    satellites = (reseeded,)
    signal_groups = _signal_groups_from_satellites(satellites)
    TrackState(satellites, signal_groups, doppler_estimator)
end

function TrackState(
    satellites::SatelliteDicts;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    signal_groups = _signal_groups_from_satellites(satellites)
    TrackState(satellites, signal_groups, doppler_estimator)
end

function TrackState(
    track_state::TrackState{S,SG,DE},
    satellites::S,
    doppler_estimator::DE,
) where {S<:SatelliteDicts,SG<:TupleLike,DE<:AbstractDopplerEstimator}
    TrackState{S,SG,DE}(satellites, track_state.signal_groups, doppler_estimator)
end

# Be careful when calling this.
# It might lead to types that are inferred at runtime?!
# Tested with 1.11.6
function TrackState(
    track_state::TrackState{S,SG,DE};
    satellites::Maybe{S} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
) where {S<:SatelliteDicts,SG<:TupleLike,DE<:AbstractDopplerEstimator}
    TrackState{S,SG,DE}(
        isnothing(satellites) ? track_state.satellites :
        satellites,
        track_state.signal_groups,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
    )
end

# Recover the per-capability signal-instance tuple from a non-empty
# satellites NamedTuple by pulling `signals[1].signal` (and friends) out
# of any sat in each dictionary. Used by the legacy
# `TrackState(system, sats)` constructor — empty-dict capabilities can't
# be recovered this way (no sats to inspect), so this path requires at
# least one sat per capability.
@inline _signal_groups_from_satellites(satellites::Tuple) =
    map(_signals_in_dict, satellites)
@inline _signal_groups_from_satellites(satellites::NamedTuple) =
    map(_signals_in_dict, satellites)

@inline function _signals_in_dict(dict::Dictionary{<:Any,<:TrackedSat})
    isempty(dict) && throw(ArgumentError(
        "Cannot recover the signal-instance tuple from an empty " *
        "satellites dictionary. Use the `TrackState(; signals = ...)` " *
        "constructor to declare signal groups before populating sats.",
    ))
    sat = first(dict.values)
    map(s -> s.signal, sat.signals)
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
    track_state::TrackState{<:SatelliteDicts{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sat_states(track_state.satellites, system_idx)
end

function get_sat_states(track_state::TrackState{<:SatelliteDicts{1}})
    get_sat_states(track_state.satellites)
end

# Pull the signal type for system `system_idx` out of any sat in its
# dictionary. The dictionary's value type carries the signal type, so this
# resolves at compile time when `system_idx` is a literal Symbol/Int.
function get_system(
    track_state::TrackState{<:SatelliteDicts{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    sats = get_sat_states(track_state, system_idx)
    get_signal(first(sats.values))
end

function get_system(track_state::TrackState{<:SatelliteDicts{1}})
    get_system(track_state, 1)
end

function get_sat_state(
    track_state::TrackState{<:SatelliteDicts{N}},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
) where {N}
    get_sat_state(get_sat_states(track_state, system_idx), sat_identifier)
end

function get_sat_state(
    track_state::TrackState{<:SatelliteDicts{1}},
    sat_identifier,
)
    get_sat_state(track_state, 1, sat_identifier)
end

function get_sat_state(track_state::TrackState{<:SatelliteDicts{1}})
    only(get_sat_states(track_state, 1))
end

function estimate_cn0(
    track_state::TrackState{<:SatelliteDicts},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
)
    estimate_cn0(get_sat_state(track_state, system_idx, sat_identifier))
end

function estimate_cn0(track_state::TrackState{<:SatelliteDicts{1}}, sat_identifier)
    estimate_cn0(track_state, 1, sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:SatelliteDicts{1}})
    estimate_cn0(get_sat_state(track_state))
end

function merge_sats(
    track_state::TrackState{S,SG,DE},
    system_idx::Union{Symbol,Integer},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
) where {S<:SatelliteDicts,SG<:TupleLike,DE<:AbstractDopplerEstimator}
    new_sats_dict = to_dictionary(tracked_sats)
    # Re-seed each incoming sat's per-sat estimator state with the
    # track_state's estimator config — the sat may have been constructed
    # with a different (often default) estimator, but its slot in the
    # tracking dictionary must hold state of the right concrete type.
    reseeded = map(new_sats_dict) do sat
        _reseed_doppler_estimator_state(sat, track_state.doppler_estimator)
    end
    new_estimator =
        update_estimator_on_handoff(track_state.doppler_estimator, reseeded)
    TrackState{S,SG,DE}(
        merge_sats(track_state.satellites, system_idx, reseeded),
        track_state.signal_groups,
        new_estimator,
    )
end

function merge_sats(
    track_state::TrackState{<:SatelliteDicts{1}},
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
)
    merge_sats(track_state, 1, tracked_sats)
end

"""
$(SIGNATURES)

Add (or replace) a satellite in `track_state` in place. Builds a
multi-signal [`TrackedSat`](@ref) for the requested capability using
the library default correlator and post-corr filter, with the supplied
acquisition-handoff values (`prn`, `code_phase`, `code_doppler`,
`carrier_phase`, `carrier_doppler`) wired into each `TrackedSignal`.
The per-satellite doppler-estimator state is initialized via
[`init_estimator_state`](@ref) against the `TrackState`'s configured
estimator.

`capability` defaults to `:default` (matching the single-capability
construction shortcut). For multi-capability `TrackState`s, pass the
capability key explicitly. If a satellite with the same `prn` already
exists in that capability's dictionary, it is overwritten.

Returns `track_state` unchanged (the dictionary is mutated in place).

```julia
track_state = TrackState(; signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),))
add_satellite!(track_state;
    prn = 11,
    capability = :modern_gps,
    code_phase = 0.0,
    carrier_doppler = 1234.0Hz,
)
```

To use a non-default correlator or post-corr-filter *type*, construct
the [`TrackedSat`](@ref) yourself and call the
`add_satellite!(track_state, capability, sat)` overload — see below.
"""
function add_satellite!(
    track_state::TrackState;
    prn::Int,
    capability::Symbol = :default,
    code_phase = 0.0,
    code_doppler::Maybe{typeof(1.0Hz)} = nothing,
    carrier_phase = 0.0,
    carrier_doppler = 0.0Hz,
)
    sat = _make_default_tracked_sat_for_capability(
        track_state, capability;
        prn, code_phase, code_doppler, carrier_phase, carrier_doppler,
    )
    add_satellite!(track_state, capability, sat)
end

"""
$(SIGNATURES)

In-place add (or replace) with a pre-built [`TrackedSat`](@ref) — the
escape hatch for power users who need non-default correlator or
post-corr-filter types. The sat's type must match the capability's
slot type already fixed at [`TrackState`](@ref) construction; passing a
sat of the wrong type errors at dispatch time.
"""
function add_satellite!(
    track_state::TrackState,
    capability::Symbol,
    sat::TrackedSat,
)
    # `_reseed_doppler_estimator_state` upgrades sats built with the
    # library-default estimator to whatever the TrackState was actually
    # configured with, when the per-sat estimator-state type doesn't
    # match what `init_estimator_state(track_state.doppler_estimator, sat)`
    # would produce.
    reseeded = _reseed_doppler_estimator_state(sat, track_state.doppler_estimator)
    insert_or_set!(_dict_for_capability(track_state, capability), reseeded.prn, reseeded)
    # `update_estimator_on_handoff` is called for its side effect on
    # estimators with cross-sat shared state (the default returns the
    # estimator unchanged). The return value must have the same concrete
    # type as the input — the TrackState is parameterized on the
    # estimator type and `track_state.doppler_estimator` is immutable.
    update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((reseeded.prn => reseeded,)),
    )
    track_state
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
    capability::Symbol = :default,
    code_phase = 0.0,
    code_doppler::Maybe{typeof(1.0Hz)} = nothing,
    carrier_phase = 0.0,
    carrier_doppler = 0.0Hz,
)
    sat = _make_default_tracked_sat_for_capability(
        track_state, capability;
        prn, code_phase, code_doppler, carrier_phase, carrier_doppler,
    )
    add_satellite(track_state, capability, sat)
end

function add_satellite(
    track_state::TrackState,
    capability::Symbol,
    sat::TrackedSat,
)
    reseeded = _reseed_doppler_estimator_state(sat, track_state.doppler_estimator)
    new_estimator = update_estimator_on_handoff(
        track_state.doppler_estimator,
        dictionary((reseeded.prn => reseeded,)),
    )
    new_dict = _dict_for_capability(track_state, capability)
    new_dict = merge(new_dict, dictionary((reseeded.prn => reseeded,)))
    new_satellites = Base.setindex(track_state.satellites, new_dict, capability)
    TrackState(new_satellites, track_state.signal_groups, new_estimator)
end

# Compile-time dispatch helper: hand back the dictionary slot for the
# given capability key. Bounds and existence are checked at TrackState
# construction time (the NamedTuple only contains declared keys), so an
# unknown `capability` here triggers the standard NamedTuple
# KeyErrors.
@inline _dict_for_capability(track_state::TrackState, capability::Symbol) =
    track_state.satellites[capability]

# Build a default-correlator, default-PCF TrackedSat whose
# signal-tuple shape matches the capability's slot in `track_state`.
# Reads the signal-instance tuple from `track_state.signal_groups`.
function _make_default_tracked_sat_for_capability(
    track_state::TrackState,
    capability::Symbol;
    prn::Int,
    code_phase,
    code_doppler,
    carrier_phase,
    carrier_doppler,
)
    sig_tuple = track_state.signal_groups[capability]
    first_signal = first(sig_tuple)
    cd = isnothing(code_doppler) ?
        carrier_doppler * get_code_center_frequency_ratio(first_signal) :
        code_doppler
    tracked_signals = map(sig_tuple) do sig
        TrackedSignal(
            sig;
            correlator = get_default_correlator(sig, NumAnts(1)),
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

function filter_out_sats(
    track_state::TrackState{S,SG,DE},
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {S<:SatelliteDicts,SG<:TupleLike,DE<:AbstractDopplerEstimator}
    TrackState{S,SG,DE}(
        filter_out_sats(track_state.satellites, system_idx, identifiers),
        track_state.signal_groups,
        track_state.doppler_estimator,
    )
end

function filter_out_sats(track_state::TrackState{<:SatelliteDicts{1}}, identifiers)
    filter_out_sats(track_state, 1, identifiers)
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
