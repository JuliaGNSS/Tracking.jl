"""
$(SIGNATURES)

Per-signal tracking state. One `TrackedSignal` exists for each signal being
tracked on a satellite — for a satellite tracked on GPS L1 C/A only there is
one; for a satellite tracked on GPS L1 C/A + L1C-D + L1C-P there are three.

Holds the signal-specific state: the correlator, post-correlation filter,
CN0 estimator, bit buffer, and the integration-progress flags. The
per-satellite carrier/code Doppler and phase are shared across signals and
live on the enclosing [`TrackedSat`](@ref).
"""
struct TrackedSignal{
    Sig<:AbstractGNSSSignal,
    C<:AbstractCorrelator,
    PCF<:AbstractPostCorrFilter,
}
    signal::Sig
    integrated_samples::Int
    is_integration_completed::Bool
    correlator::C
    last_fully_integrated_correlator::C
    last_fully_integrated_filtered_prompt::ComplexF64
    cn0_estimator::MomentsCN0Estimator
    bit_buffer::BitBuffer
    post_corr_filter::PCF
    filtered_prompts::Vector{ComplexF64}
end

"""
$(SIGNATURES)

Construct a fresh [`TrackedSignal`](@ref) for `signal`. The correlator and
post-corr filter default to the system's recommended values; pass
`correlator` and `post_corr_filter` explicitly to override.
"""
function TrackedSignal(
    signal::AbstractGNSSSignal;
    num_ants::NumAnts = NumAnts(1),
    correlator::AbstractCorrelator = get_default_correlator(signal, num_ants),
    num_prompts_for_cn0_estimation::Int = 100,
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
)
    TrackedSignal(
        signal,
        0,
        false,
        correlator,
        correlator,
        complex(0.0, 0.0),
        MomentsCN0Estimator(num_prompts_for_cn0_estimation),
        BitBuffer(),
        post_corr_filter,
        ComplexF64[],
    )
end

# Kwarg-update constructor — produces a new TrackedSignal sharing concrete
# correlator and PCF types with `t`. The two constrained fields use
# `Maybe{C}` / `Maybe{PCF}` to keep `nothing` distinguishable from a real
# value, since the user might legitimately want to set them to anything of
# the same type.
function TrackedSignal(
    t::TrackedSignal{Sig,C,PCF};
    signal = nothing,
    integrated_samples = nothing,
    is_integration_completed = nothing,
    correlator::Maybe{C} = nothing,
    last_fully_integrated_correlator::Maybe{C} = nothing,
    last_fully_integrated_filtered_prompt = nothing,
    cn0_estimator = nothing,
    bit_buffer = nothing,
    post_corr_filter::Maybe{PCF} = nothing,
    filtered_prompts::Maybe{Vector{ComplexF64}} = nothing,
) where {Sig<:AbstractGNSSSignal,C<:AbstractCorrelator,PCF<:AbstractPostCorrFilter}
    TrackedSignal{Sig,C,PCF}(
        isnothing(signal) ? t.signal : signal,
        isnothing(integrated_samples) ? t.integrated_samples : integrated_samples,
        isnothing(is_integration_completed) ? t.is_integration_completed :
        is_integration_completed,
        isnothing(correlator) ? t.correlator : correlator,
        isnothing(last_fully_integrated_correlator) ? t.last_fully_integrated_correlator :
        last_fully_integrated_correlator,
        isnothing(last_fully_integrated_filtered_prompt) ?
        t.last_fully_integrated_filtered_prompt : last_fully_integrated_filtered_prompt,
        isnothing(cn0_estimator) ? t.cn0_estimator : cn0_estimator,
        isnothing(bit_buffer) ? t.bit_buffer : bit_buffer,
        isnothing(post_corr_filter) ? t.post_corr_filter : post_corr_filter,
        isnothing(filtered_prompts) ? t.filtered_prompts : filtered_prompts,
    )
end

get_signal(t::TrackedSignal) = t.signal
get_correlator(t::TrackedSignal) = t.correlator
get_last_fully_integrated_correlator(t::TrackedSignal) = t.last_fully_integrated_correlator
get_last_fully_integrated_filtered_prompt(t::TrackedSignal) =
    t.last_fully_integrated_filtered_prompt
get_filtered_prompts(t::TrackedSignal) = t.filtered_prompts
get_post_corr_filter(t::TrackedSignal) = t.post_corr_filter
get_cn0_estimator(t::TrackedSignal) = t.cn0_estimator
get_bit_buffer(t::TrackedSignal) = t.bit_buffer
get_bits(t::TrackedSignal) = get_bits(t.bit_buffer)
get_num_bits(t::TrackedSignal) = length(t.bit_buffer)
has_bit_or_secondary_code_been_found(t::TrackedSignal) =
    has_bit_or_secondary_code_been_found(t.bit_buffer)
get_integrated_samples(t::TrackedSignal) = t.integrated_samples

"""
$(SIGNATURES)

Holds the state of a single satellite being tracked. Carries the satellite-
level carrier/code Doppler and phase (shared across all signals on this
satellite), the per-signal correlator state in `signals::Tuple{Vararg{TrackedSignal}}`,
and the per-satellite Doppler-estimator state in `doppler_estimator_state`.

The first signal in `signals` is the Doppler-source signal — its correlator
is what the PLL/DLL discriminator runs on. Per-satellite carrier/code Doppler
updates therefore happen at the rate of the first signal's integration
boundary. Other signals filter their own prompts and update their own CN0
estimates and bit buffers on their own boundaries.

The shared `code_phase` wraps at the longest code period across all signals
including secondary code (see [`max_code_length`](@ref)).
"""
struct TrackedSat{Signals<:Tuple{Vararg{TrackedSignal}},D}
    prn::Int
    code_phase::Float64
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    carrier_doppler::typeof(1.0Hz)
    signal_start_sample::Int
    signals::Signals
    doppler_estimator_state::D
end

@inline _max_code_length(::Tuple{}) = 0
@inline _max_code_length(t::Tuple) = max(
    get_code_length(first(t).signal) * get_secondary_code_length(first(t).signal),
    _max_code_length(Base.tail(t)),
)

"""
$(SIGNATURES)

The longest code period (in chips, including secondary code) across all
signals in `signals`. Returns the constant the shared `sat.code_phase`
wraps at — for a sat tracking only GPS L1 C/A (1023 chips × 1) this is
1023; for one tracking L1C-P (10230 × 1800 ≈ 18.4 M) it is 18 414 000.

Implemented via tuple recursion (not `@generated`) so the heterogeneous
walk unrolls at type-inference time; the result folds to a literal in
the calling site for any concrete `signals` tuple type.
"""
@inline max_code_length(signals::Tuple{Vararg{TrackedSignal}}) = _max_code_length(signals)

"""
$(SIGNATURES)

Construct a [`TrackedSat`](@ref) from acquisition handoff values plus the
Doppler estimator. The signal is wrapped in a single [`TrackedSignal`](@ref);
multi-signal tracking is wired in a later step. The estimator's per-satellite
state is built via [`init_estimator_state`](@ref).
"""
function TrackedSat(
    signal::AbstractGNSSSignal,
    prn::Int,
    code_phase,
    carrier_doppler;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    num_ants::NumAnts = NumAnts(1),
    correlator::AbstractCorrelator = get_default_correlator(signal, num_ants),
    carrier_phase = 0.0,
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(signal),
    num_prompts_for_cn0_estimation::Int = 100,
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
)
    tracked_signal = TrackedSignal(
        signal;
        num_ants,
        correlator,
        num_prompts_for_cn0_estimation,
        post_corr_filter,
    )
    # Two-stage build so the estimator's `init_estimator_state` sees a real
    # `TrackedSat` (handy for estimators that need carrier/code Doppler at
    # init time). The first stage builds a sat with `D = Nothing`; the
    # second stage rebuilds it with the actual `doppler_estimator_state`.
    bare = TrackedSat(
        prn,
        float(code_phase),
        float(code_doppler),
        float(carrier_phase) / 2π,
        float(carrier_doppler),
        1,
        (tracked_signal,),
        nothing,
    )
    doppler_estimator_state = init_estimator_state(doppler_estimator, bare)
    TrackedSat(
        bare.prn,
        bare.code_phase,
        bare.code_doppler,
        bare.carrier_phase,
        bare.carrier_doppler,
        bare.signal_start_sample,
        bare.signals,
        doppler_estimator_state,
    )
end

# Kwarg-update constructor. `signals` and `doppler_estimator_state` carry
# `Maybe{...}` constraints so that the new value retains the same concrete
# type as the original (preventing accidental type changes that would break
# inference on the enclosing TrackState).
function TrackedSat(
    sat::TrackedSat{Signals,D};
    prn = nothing,
    code_phase = nothing,
    code_doppler = nothing,
    carrier_phase = nothing,
    carrier_doppler = nothing,
    signal_start_sample = nothing,
    signals::Maybe{Signals} = nothing,
    doppler_estimator_state::Maybe{D} = nothing,
) where {Signals<:Tuple{Vararg{TrackedSignal}},D}
    TrackedSat{Signals,D}(
        isnothing(prn) ? sat.prn : prn,
        isnothing(code_phase) ? sat.code_phase : code_phase,
        isnothing(code_doppler) ? sat.code_doppler : code_doppler,
        isnothing(carrier_phase) ? sat.carrier_phase : carrier_phase,
        isnothing(carrier_doppler) ? sat.carrier_doppler : carrier_doppler,
        isnothing(signal_start_sample) ? sat.signal_start_sample : signal_start_sample,
        isnothing(signals) ? sat.signals : signals,
        isnothing(doppler_estimator_state) ? sat.doppler_estimator_state :
        doppler_estimator_state,
    )
end

"""
$(SIGNATURES)

Get the PRN (Pseudo-Random Noise) number of the satellite.
"""
get_prn(s::TrackedSat) = s.prn
get_num_ants(s::TrackedSat{<:Tuple{TrackedSignal{<:Any,<:AbstractCorrelator{M}},Vararg}}) where {M} = M

"""
$(SIGNATURES)

Get the satellite's shared code phase. Wraps at [`max_code_length`](@ref) of
the signals tuple (in chips, including secondary code). To get the
replica-relative phase for a specific signal, mod by that signal's primary
code length.
"""
get_code_phase(s::TrackedSat) = s.code_phase

"""
$(SIGNATURES)

Get the current code Doppler frequency.
"""
get_code_doppler(s::TrackedSat) = s.code_doppler

"""
$(SIGNATURES)

Get the current carrier phase in radians.
"""
get_carrier_phase(s::TrackedSat) = s.carrier_phase * 2π

"""
$(SIGNATURES)

Get the current carrier Doppler frequency.
"""
get_carrier_doppler(s::TrackedSat) = s.carrier_doppler

"""
$(SIGNATURES)

Get the starting sample index in the signal for the next integration.
"""
get_signal_start_sample(s::TrackedSat) = s.signal_start_sample

"""
$(SIGNATURES)

Get the satellite's tuple of [`TrackedSignal`](@ref)s.
"""
get_signals(s::TrackedSat) = s.signals

"""
$(SIGNATURES)

Get the per-satellite Doppler estimator state (e.g. the loop-filter state
for the conventional PLL/DLL).
"""
get_doppler_estimator_state(s::TrackedSat) = s.doppler_estimator_state

# Convenience accessors that forward through to the single-signal case.
# These will need to disambiguate by signal type once multi-signal tracking
# lands in a later step.
get_signal(s::TrackedSat) = get_signal(only(s.signals))
get_correlator(s::TrackedSat) = get_correlator(only(s.signals))
get_last_fully_integrated_correlator(s::TrackedSat) =
    get_last_fully_integrated_correlator(only(s.signals))
get_last_fully_integrated_filtered_prompt(s::TrackedSat) =
    get_last_fully_integrated_filtered_prompt(only(s.signals))
get_filtered_prompts(s::TrackedSat) = get_filtered_prompts(only(s.signals))
get_post_corr_filter(s::TrackedSat) = get_post_corr_filter(only(s.signals))
get_cn0_estimator(s::TrackedSat) = get_cn0_estimator(only(s.signals))
get_bit_buffer(s::TrackedSat) = get_bit_buffer(only(s.signals))
get_bits(s::TrackedSat) = get_bits(only(s.signals))
get_num_bits(s::TrackedSat) = get_num_bits(only(s.signals))
has_bit_or_secondary_code_been_found(s::TrackedSat) =
    has_bit_or_secondary_code_been_found(only(s.signals))
get_integrated_samples(s::TrackedSat) = get_integrated_samples(only(s.signals))

# Reset the satellite's signal-start sample and per-signal bit buffer between
# `track` calls. Per-signal `filtered_prompts` vectors are emptied in place;
# the per-signal `bit_buffer` is reset to its no-sync state; the
# `signal_start_sample` returns to 1 (the first sample of the next buffer).
function reset_start_sample_and_bit_buffer(sat::TrackedSat)
    new_signals = map(s -> _reset_signal(s), sat.signals)
    TrackedSat(sat; signal_start_sample = 1, signals = new_signals)
end

@inline function _reset_signal(t::TrackedSignal)
    empty!(t.filtered_prompts)
    TrackedSignal(t; bit_buffer = reset(t.bit_buffer))
end

"""
$(SIGNATURES)

Build the per-satellite Doppler-estimator state used by `estimator` for the
given satellite. Called once per satellite when a new sat enters the track
set. A custom doppler estimator must define this method for its
[`AbstractDopplerEstimator`](@ref) subtype.
"""
function init_estimator_state end

"""
$(SIGNATURES)

Optionally update an estimator's cross-satellite or cross-system shared state
when new satellites enter the track set. Called once per [`merge_sats`](@ref)
call with the dictionary of incoming satellites, *after* per-sat seeding via
[`init_estimator_state`](@ref).

The default returns `estimator` unchanged, so estimators with no shared state
need not implement it.

The returned estimator must have the same concrete type as the input —
[`TrackState`](@ref) is parameterized on the estimator type, and changing it
would break inference. For growing shared state, hold the storage in a
resizable container (e.g. `Vector`, `Matrix`) on an otherwise immutable
estimator and `push!`/`resize!` it in place; rebuild the estimator with
[`Setfield.@set`](https://jw3126.github.io/Setfield.jl/stable/) or a copying
constructor when fields need replacing.
"""
update_estimator_on_handoff(estimator::AbstractDopplerEstimator, _new_sats) = estimator

"""
Type alias for a tuple or named tuple of per-system satellite dictionaries.
Each entry is `Dictionary{I, <:TrackedSat}` — the keys are satellite
identifiers (PRNs) and the values are the per-sat tracking state. The signal
type for each system lives in the dictionary value type, accessed via
`only(sat.signals).signal` at use sites.
"""
const SatelliteDicts{N} = TupleLike{<:NTuple{N,Dictionary{<:Any,<:TrackedSat}}}

"""
$(SIGNATURES)

A group of satellites that all track the same tuple of GNSS signal types,
on the same RF band, observed by the same antenna array.

Groups are the unit of type stability: every `TrackedSat` inside a
`SignalGroup` shares the same concrete `Tuple{Vararg{TrackedSignal}}`
shape, so the dictionary's value type is concrete and the hot loop sees
no dynamic dispatch.

Two groups may share a band (e.g. `:legacy_gps` tracking `(GPSL1CA(),)` and
`:galileo` tracking `(GalileoE1B(),)` both on `L1()`). The grouping is by
signal-tuple shape, not by band — band is metadata each group carries so
the right measurement is routed to it during `track`.

Fields:
- `band`: an `AbstractGNSSSignal` `Band` instance (`L1()`, `L5()`, …)
- `satellites`: `Dictionary{Int, <:TrackedSat}` keyed by PRN
- `signals`: the signal-instance tuple (e.g. `(GPSL1C_P(), GPSL1C_D(), GPSL1CA())`)
- `num_ants`: the antenna count for this group's band
"""
struct SignalGroup{
    B,                                             # GNSSSignals Band instance
    S<:Dictionary{<:Any,<:TrackedSat},
    Sigs<:Tuple{Vararg{AbstractGNSSSignal}},
    NA<:NumAnts,
}
    band::B
    satellites::S
    signals::Sigs
    num_ants::NA
end

# Kwarg-update constructor — produces a new SignalGroup sharing concrete types
# with `g`. The `satellites` field uses `Maybe{S}` so `nothing` stays
# distinguishable from a real dict.
function SignalGroup(
    g::SignalGroup{B,S,Sigs,NA};
    band::Maybe{B} = nothing,
    satellites::Maybe{S} = nothing,
    signals::Maybe{Sigs} = nothing,
    num_ants::Maybe{NA} = nothing,
) where {B,S<:Dictionary{<:Any,<:TrackedSat},Sigs<:Tuple{Vararg{AbstractGNSSSignal}},NA<:NumAnts}
    SignalGroup{B,S,Sigs,NA}(
        isnothing(band) ? g.band : band,
        isnothing(satellites) ? g.satellites : satellites,
        isnothing(signals) ? g.signals : signals,
        isnothing(num_ants) ? g.num_ants : num_ants,
    )
end

"""
Type alias: NamedTuple of `SignalGroup`s — the storage shape inside
[`TrackState`](@ref). The `N` parameter is the number of groups.
"""
const SignalGroups{N} = NamedTuple{<:Any,<:NTuple{N,SignalGroup}}

"""
$(SIGNATURES)

Merge already-built tracked satellites into the existing tracking state.
Used internally by [`TrackState`](@ref)'s `merge_sats` — external callers
should prefer the `TrackState`-level method.
"""
function merge_sats(
    satellites::SatelliteDicts,
    system_idx,
    new_tracked_sats::Dictionary{<:Any,<:TrackedSat},
)
    sats_dict = satellites[system_idx]
    @set satellites[system_idx] = merge(sats_dict, new_tracked_sats)
end

# Build a `Dictionary` that shares its keys with the original but holds a
# freshly-copied `values::Vector{TrackedSat}`. Used by the immutable variants
# of `downconvert_and_correlate` / `estimate_dopplers_and_filter_prompt` /
# `reset_start_sample_and_bit_buffer` to detach the slot vector before
# delegating to the in-place form — preserves `track`'s "input is not
# mutated" contract while reusing the in-place implementation.
@inline function _copy_slot_vector(sats::Dictionary{<:Any,<:TrackedSat})
    Dictionary(keys(sats), copy(sats.values))
end

@inline function _copy_slot_vectors(
    satellites::SatelliteDicts,
)
    map(_copy_slot_vector, satellites)
end

# Recursive tuple walker: applies `f(sats_dict, args...)` to each per-system
# satellite dictionary in the (named-)tuple. Each step has fully concrete
# types, so no runtime dispatch even when systems differ (heterogeneous
# `Tuple{Dictionary{Int, TrackedSat{Tuple{TrackedSignal{GPSL1CA,...}},...}}, Dictionary{Int, TrackedSat{Tuple{TrackedSignal{GalileoE1B,...}},...}}}`
# would otherwise box each element when iterated with `for s in tuple`).
@inline _foreach_system!(f::F, ::Tuple{}, args::Vararg{Any,N}) where {F,N} = nothing
@inline function _foreach_system!(f::F, t::Tuple, args::Vararg{Any,N}) where {F,N}
    f(first(t), args...)
    _foreach_system!(f, Base.tail(t), args...)
end
# NamedTuple unwraps to its underlying Tuple — concrete and cheap.
@inline _foreach_system!(f::F, nt::NamedTuple, args::Vararg{Any,N}) where {F,N} =
    _foreach_system!(f, Tuple(nt), args...)

# In-place variant: walks each per-system dictionary's `Vector{TrackedSat}`
# slot storage and overwrites each entry with a freshly reset value. The
# vector itself, the Dictionary, and the enclosing tuple are all reused.
# `TrackedSat` itself is immutable, so the slot is reassigned rather than
# mutated; we still call the non-`!` per-`TrackedSat` form.
@inline function _reset_one_system!(sats::Dictionary{<:Any,<:TrackedSat})
    vals = sats.values
    @inbounds for i in eachindex(vals)
        vals[i] = reset_start_sample_and_bit_buffer(vals[i])
    end
    return nothing
end

function reset_start_sample_and_bit_buffer!(
    satellites::SatelliteDicts,
)
    _foreach_system!(_reset_one_system!, satellites)
    return satellites
end

# Immutable variant: copies the slot vectors first, then delegates to the
# in-place form. Same per-sat work, only the storage ownership differs.
function reset_start_sample_and_bit_buffer(
    satellites::SatelliteDicts,
)
    reset_start_sample_and_bit_buffer!(_copy_slot_vectors(satellites))
end

function to_dictionary(tracked_sats::Dictionary{I,<:TrackedSat}) where {I}
    tracked_sats
end

function to_dictionary(tracked_sats::Vector{<:TrackedSat})
    Dictionary(map(get_prn, tracked_sats), tracked_sats)
end

function to_dictionary(t::TrackedSat)
    dictionary((get_prn(t) => t,))
end

"""
$(SIGNATURES)

Get the satellite state for a specific satellite identifier.
"""
get_sat_state(sats::Dictionary{<:Any,<:TrackedSat}, identifier) = sats[identifier]
get_sat_state(sats::Dictionary{<:Any,<:TrackedSat}) = only(sats)

function estimate_cn0(sat::TrackedSat)
    signal = get_signal(sat)
    estimate_cn0(
        get_cn0_estimator(sat),
        get_code_length(signal) / get_code_frequency(signal),
    )
end
