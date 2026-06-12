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
    B<:Unsigned,
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
    bit_buffer::BitBuffer{B}
    post_corr_filter::PCF
    filtered_prompts::Vector{ComplexF64}
    # Preferred coherent-integration length for THIS signal, in primary code
    # blocks. The actual length is capped per integration by the signal's
    # bit/secondary-code period and held at 1 until bit/secondary sync (see
    # `calc_num_code_blocks_to_integrate`). Defaults to 1; change it with
    # [`set_preferred_num_code_blocks_to_integrate!`](@ref).
    preferred_num_code_blocks_to_integrate::Int
end

# Reject a preferred coherent-integration length that cannot work for this
# signal. For data-bearing signals the length must evenly divide the number
# of code blocks that form one bit: a non-divisor length would make the
# post-sync integrations straddle bit boundaries, so the bit buffer's
# blocks-per-bit boundary check would never fire and no bit would ever be
# emitted again (issue #128). Pilot signals (`data_frequency == 0`) have no
# bit boundary to straddle, so any length of at least one block is accepted.
function validate_preferred_num_code_blocks_to_integrate(
    signal::AbstractGNSSSignal,
    preferred_num_code_blocks::Integer,
)
    preferred_num_code_blocks >= 1 || throw(
        ArgumentError(
            "preferred_num_code_blocks_to_integrate must be at least 1, got " *
            "$preferred_num_code_blocks",
        ),
    )
    num_code_blocks_that_form_a_bit = _calc_num_code_blocks_that_form_a_bit(signal)
    num_code_blocks_that_form_a_bit == 0 && return nothing
    if num_code_blocks_that_form_a_bit % preferred_num_code_blocks != 0
        valid = filter(
            d -> num_code_blocks_that_form_a_bit % d == 0,
            1:num_code_blocks_that_form_a_bit,
        )
        throw(
            ArgumentError(
                "preferred_num_code_blocks_to_integrate = $preferred_num_code_blocks " *
                "must evenly divide the $num_code_blocks_that_form_a_bit code blocks " *
                "that form one bit of $(typeof(signal)) — an integration straddling " *
                "a bit boundary would never emit a bit. Valid values: " *
                "$(join(valid, ", ")).",
            ),
        )
    end
    nothing
end

"""
$(SIGNATURES)

Construct a fresh [`TrackedSignal`](@ref) for `signal`. The correlator and
post-corr filter default to the signal's recommended values; pass
`correlator` and `post_corr_filter` explicitly to override.

Throws an `ArgumentError` if `preferred_num_code_blocks_to_integrate` is
invalid for `signal` (see
[`set_preferred_num_code_blocks_to_integrate!`](@ref)).
"""
function TrackedSignal(
    signal::AbstractGNSSSignal;
    num_ants::NumAnts = NumAnts(1),
    correlator::AbstractCorrelator = get_default_correlator(signal, num_ants),
    num_prompts_for_cn0_estimation::Int = 100,
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    preferred_num_code_blocks_to_integrate::Int = 1,
)
    validate_preferred_num_code_blocks_to_integrate(
        signal,
        preferred_num_code_blocks_to_integrate,
    )
    # Per-signal sync-search buffer width — see `get_code_block_buffer_type`.
    # Picking the type here makes `B` concrete in the resulting
    # `TrackedSignal{Sig, B, C, PCF}`.
    B = get_code_block_buffer_type(signal)
    TrackedSignal(
        signal,
        0,
        false,
        correlator,
        correlator,
        complex(0.0, 0.0),
        MomentsCN0Estimator(num_prompts_for_cn0_estimation),
        BitBuffer{B}(),
        post_corr_filter,
        ComplexF64[],
        preferred_num_code_blocks_to_integrate,
    )
end

# Kwarg-update constructor — produces a new TrackedSignal sharing concrete
# correlator and PCF types with `t`. The two constrained fields use
# `Maybe{C}` / `Maybe{PCF}` to keep `nothing` distinguishable from a real
# value, since the user might legitimately want to set them to anything of
# the same type.
function TrackedSignal(
    t::TrackedSignal{Sig,B,C,PCF};
    signal = nothing,
    integrated_samples = nothing,
    is_integration_completed = nothing,
    correlator::Maybe{C} = nothing,
    last_fully_integrated_correlator::Maybe{C} = nothing,
    last_fully_integrated_filtered_prompt = nothing,
    cn0_estimator = nothing,
    bit_buffer::Maybe{BitBuffer{B}} = nothing,
    post_corr_filter::Maybe{PCF} = nothing,
    filtered_prompts::Maybe{Vector{ComplexF64}} = nothing,
    preferred_num_code_blocks_to_integrate = nothing,
) where {Sig<:AbstractGNSSSignal,B<:Unsigned,C<:AbstractCorrelator,PCF<:AbstractPostCorrFilter}
    isnothing(preferred_num_code_blocks_to_integrate) ||
        validate_preferred_num_code_blocks_to_integrate(
            isnothing(signal) ? t.signal : signal,
            preferred_num_code_blocks_to_integrate,
        )
    TrackedSignal{Sig,B,C,PCF}(
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
        isnothing(preferred_num_code_blocks_to_integrate) ?
        t.preferred_num_code_blocks_to_integrate : preferred_num_code_blocks_to_integrate,
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
get_soft_bits(t::TrackedSignal) = get_soft_bits(t.bit_buffer)
get_num_bits(t::TrackedSignal) = length(t.bit_buffer)
has_bit_or_secondary_code_been_found(t::TrackedSignal) =
    has_bit_or_secondary_code_been_found(t.bit_buffer)
get_integrated_samples(t::TrackedSignal) = t.integrated_samples
get_preferred_num_code_blocks_to_integrate(t::TrackedSignal) =
    t.preferred_num_code_blocks_to_integrate

"""
$(SIGNATURES)

Holds the state of a single satellite being tracked. Carries the satellite-
level carrier/code Doppler and phase (shared across all signals on this
satellite), the per-signal correlator state in `signals::Tuple{Vararg{TrackedSignal}}`,
and the per-satellite Doppler-estimator state in `doppler_estimator_state`.

The first signal in `signals` is the **estimator-driver signal** — the one
the Doppler estimator uses to update the satellite-shared carrier and code
Doppler. With the default [`ConventionalPLLAndDLL`](@ref) /
[`ConventionalAssistedPLLAndDLL`](@ref), that means `signals[1]`'s
correlator is what the PLL/DLL discriminator runs on, and per-satellite
Doppler updates happen at the rate of the first signal's integration
boundary; other signals filter their own prompts and update their own CN0
estimates and bit buffers on their own boundaries. A user-supplied
[`AbstractDopplerEstimator`](@ref) is free to use the other signals' state
too — `signals[1]`'s privileged role is a convention of the conventional
estimators, not a structural constraint of the type.

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

# Worst-case wrap *length* contributed by a single signal — used both by
# the compile-time `max_code_length` upper bound and by the runtime
# `current_code_wrap` (which falls back to the primary-only length when
# the signal hasn't synced yet).
#
# For a pilot (`data_frequency == 0`) the long wrap is one full secondary
# code period. For a data-bearing signal it is one full data-bit period,
# which may or may not coincide with the secondary code length:
#   - GPS L5I: 10230 × 10 = 102300 chips (NH10 spans exactly one data
#     bit, so the two notions agree).
#   - GPS L1 C/A: 1023 × 20 = 20460 chips (no secondary code; the long
#     wrap is purely the 20-block bit period).
# We take the `max` of the two so multi-signal callers stay
# upper-bounded by whichever is larger.
@inline function _post_sync_code_length(tsig::TrackedSignal)
    sig = tsig.signal
    primary = get_code_length(sig)
    secondary = get_secondary_code_length(sig)
    df = get_data_frequency(sig)
    if iszero(df)
        primary * secondary
    else
        blocks_per_bit = Int(get_code_frequency(sig) / (primary * df))
        primary * max(secondary, blocks_per_bit)
    end
end

@inline _max_code_length(::Tuple{}) = 0
@inline _max_code_length(t::Tuple) = max(
    _post_sync_code_length(first(t)),
    _max_code_length(Base.tail(t)),
)

"""
$(SIGNATURES)

Upper bound on the shared `sat.code_phase` wrap period, in chips —
the largest value `code_phase` ever takes once *every* signal on the
sat has synced. For a sat tracking only GPS L1 C/A this is 1023 × 20 =
20460 (one full data bit); for one tracking L1C-P this is 10230 × 1800
≈ 18.4 M (one full secondary-code cycle).

This is the *compile-time* bound. The actual runtime wrap shrinks to
`get_code_length(signal) × 1` for any signal whose bit/secondary-code
sync hasn't been found yet — see [`current_code_wrap`](@ref) for the
runtime value used by the inner loop.

Implemented via tuple recursion (not `@generated`) so the heterogeneous
walk unrolls at type-inference time; the result folds to a literal in
the calling site for any concrete `signals` tuple type.
"""
@inline max_code_length(signals::Tuple{Vararg{TrackedSignal}}) = _max_code_length(signals)

# Runtime per-signal wrap contribution, honoring the current sync state:
# a signal that hasn't synced yet wraps at just its primary code length
# (we don't yet know which bit / secondary-chip we're in), so its
# contribution to the shared wrap is `primary × 1`. Once synced its
# contribution widens to `_post_sync_code_length(tsig)`.
@inline function _current_code_length(tsig::TrackedSignal)
    if tsig.bit_buffer.found
        _post_sync_code_length(tsig)
    else
        get_code_length(tsig.signal)
    end
end

@inline _current_code_wrap(::Tuple{}) = 0
@inline _current_code_wrap(t::Tuple) = max(
    _current_code_length(first(t)),
    _current_code_wrap(Base.tail(t)),
)

"""
$(SIGNATURES)

The runtime wrap period for the shared `sat.code_phase`, in chips —
the value the inner loop uses to wrap `code_phase` modulo each
integration step.

Unlike [`max_code_length`](@ref) (which is the worst-case bound), this
honors the current per-signal sync state. For each signal:

- If `bit_buffer.found = true`, the signal contributes
  `primary × max(secondary_code_length, blocks_per_data_bit)` — i.e.
  the full secondary-code period for pilots, or the full data-bit
  period for data-bearing signals.
- If `bit_buffer.found = false`, the signal contributes just its
  primary code length — we don't yet know which bit / secondary chip
  we're in, so wrapping at the primary length is the most we can
  legitimately do.

The shared wrap is the `max` across all signals, so the longest synced
signal pins the wrap to its full sync period and shorter signals just
ride along.
"""
@inline current_code_wrap(signals::Tuple{Vararg{TrackedSignal}}) = _current_code_wrap(signals)

# Detect whether any signal on the sat transitioned `bit_buffer.found`
# from `false` (in `old_signals`) to `true` (in `new_signals`) during
# this estimator iteration. Used to gate the one-time code-phase snap:
# the snap is only valid at the sync transition (when `code_phase` sits
# on a primary-code boundary) and must not run on subsequent iterations
# (see `_update_tracked_sat_doppler` and issue #117). Walks the two
# tuples in lockstep; folds to a compile-time-decided chain of `||`s.
@inline _any_signal_just_synced(::Tuple{}, ::Tuple{}) = false
@inline _any_signal_just_synced(old::Tuple, new::Tuple) =
    (!first(old).bit_buffer.found && first(new).bit_buffer.found) ||
    _any_signal_just_synced(Base.tail(old), Base.tail(new))

# Phase-snap fallback chain: walk `signals` and pick the synced signal
# whose `(primary × secondary)` code length is the largest. That signal's
# `bit_buffer.secondary_phase` carries the secondary-chip offset for the
# next primary-code period, which determines the absolute position in
# `sat.code_phase`'s wrap window. If no signal is synced (or the synced
# ones all have secondary code length 1 and so don't constrain the
# wrap), return the input `code_phase` unchanged.
#
# Walks the tuple recursively; the heterogeneous signal types fold to a
# compile-time-decided sequence of comparisons, so this is type-stable
# and allocation-free.
@inline function _snap_code_phase_from_synced_signal(
    signals::Tuple{Vararg{TrackedSignal}},
    code_phase::Float64,
)
    best_len, best_phase_chips, best_prim = _find_best_secondary_anchor(signals, 0, 0, 1)
    if best_len == 0
        return code_phase
    end
    # `code_phase` already wraps mod `max_code_length(signals)`. The
    # best-anchored signal has wrap = best_len; align the low
    # `best_len` chips of `code_phase` to the synced secondary-chip
    # window (`best_phase_chips`, a multiple of the primary length),
    # leaving the higher-order wrap untouched. Crucially we *preserve*
    # the within-primary-block phase `mod(code_phase, best_prim)`: the
    # snap places `code_phase` into the right secondary window without
    # discarding any partial (chunk-bounded) integration progress within
    # the current primary block — dropping it would inject a one-block
    # phase error at sync and force the loops to re-converge (issue #117).
    base = floor(Int, code_phase / best_len) * best_len
    within_primary = code_phase - floor(code_phase / best_prim) * best_prim
    Float64(base + best_phase_chips) + within_primary
end

# Tuple walker — finds the synced signal with the largest
# `(primary × secondary)` length and returns
# `(length, secondary_phase_in_chips, primary_length)` for it.
# `(0, 0, 1)` if no signal is synced.
@inline _find_best_secondary_anchor(::Tuple{}, best_len::Int, best_chips::Int, best_prim::Int) =
    (best_len, best_chips, best_prim)
@inline function _find_best_secondary_anchor(
    t::Tuple,
    best_len::Int,
    best_chips::Int,
    best_prim::Int,
)
    head = first(t)
    sig = head.signal
    bb = head.bit_buffer
    if bb.found
        prim = get_code_length(sig)
        sec  = get_secondary_code_length(sig)
        total = prim * sec
        # Only signals with a non-trivial secondary code contribute
        # information about the wrap-window offset. Signals with
        # secondary_code_length == 1 (e.g. GPS L1 C/A) lock only the
        # bit-edge inside the primary period and don't pin the secondary
        # phase; skip them.
        if sec > 1 && total > best_len
            best_len = total
            best_chips = bb.secondary_phase * prim
            best_prim = prim
        end
    end
    _find_best_secondary_anchor(Base.tail(t), best_len, best_chips, best_prim)
end

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
    doppler_estimator::AbstractDopplerEstimator =
        ConventionalAssistedPLLAndDLL(;
            carrier_loop_filter_bandwidth = default_carrier_loop_filter_bandwidth(signal),
            code_loop_filter_bandwidth = default_code_loop_filter_bandwidth(signal),
        ),
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
get_num_ants(s::TrackedSat{<:Tuple{TrackedSignal{<:Any,<:Any,<:AbstractCorrelator{M}},Vararg}}) where {M} = M

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

# Per-signal accessors on a TrackedSat.
#
# Three selector forms (in increasing specificity), all routed through
# `_find_signal` so every per-signal accessor is one line:
#   * no selector — only valid for single-signal sats; falls back to
#     `only(s.signals)`. Errors on multi-signal sats.
#   * `Integer` index — canonical: picks `s.signals[i]`. Unambiguous even
#     when the same signal type appears twice in the tuple.
#   * `Type{<:AbstractGNSSSignal}` — sugar: picks the unique signal of
#     that type. Errors if zero or >1 matches.
#
# Type-based selection walks the signals tuple recursively and folds at
# compile time when the sat's `Signals` type is concrete.
@inline _find_signal(s::Tuple) = only(s)
@inline _find_signal(s::Tuple, i::Integer) = s[i]
@inline _find_signal(s::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal} =
    _find_signal_by_type(s, T)

@inline function _find_signal_by_type(::Tuple{}, ::Type{T}) where {T<:AbstractGNSSSignal}
    throw(ArgumentError("no signal of type $T on this satellite"))
end
@inline function _find_signal_by_type(t::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    head = first(t)
    if head.signal isa T
        _assert_no_more_of_type(Base.tail(t), T)
        return head
    end
    _find_signal_by_type(Base.tail(t), T)
end
@inline _assert_no_more_of_type(::Tuple{}, ::Type{T}) where {T<:AbstractGNSSSignal} = nothing
@inline function _assert_no_more_of_type(t::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    first(t).signal isa T && throw(ArgumentError(
        "signal type $T matches more than one signal on this satellite — " *
        "use an integer index to disambiguate (e.g. `get_signals(sat)[i]`)."))
    _assert_no_more_of_type(Base.tail(t), T)
end

get_signal(s::TrackedSat, sel...) = get_signal(_find_signal(s.signals, sel...))
get_correlator(s::TrackedSat, sel...) = get_correlator(_find_signal(s.signals, sel...))
get_last_fully_integrated_correlator(s::TrackedSat, sel...) =
    get_last_fully_integrated_correlator(_find_signal(s.signals, sel...))
get_last_fully_integrated_filtered_prompt(s::TrackedSat, sel...) =
    get_last_fully_integrated_filtered_prompt(_find_signal(s.signals, sel...))
get_filtered_prompts(s::TrackedSat, sel...) = get_filtered_prompts(_find_signal(s.signals, sel...))
get_post_corr_filter(s::TrackedSat, sel...) = get_post_corr_filter(_find_signal(s.signals, sel...))
get_cn0_estimator(s::TrackedSat, sel...) = get_cn0_estimator(_find_signal(s.signals, sel...))
get_bit_buffer(s::TrackedSat, sel...) = get_bit_buffer(_find_signal(s.signals, sel...))
get_bits(s::TrackedSat, sel...) = get_bits(_find_signal(s.signals, sel...))
get_soft_bits(s::TrackedSat, sel...) = get_soft_bits(_find_signal(s.signals, sel...))
get_num_bits(s::TrackedSat, sel...) = get_num_bits(_find_signal(s.signals, sel...))
has_bit_or_secondary_code_been_found(s::TrackedSat, sel...) =
    has_bit_or_secondary_code_been_found(_find_signal(s.signals, sel...))
get_integrated_samples(s::TrackedSat, sel...) =
    get_integrated_samples(_find_signal(s.signals, sel...))
get_preferred_num_code_blocks_to_integrate(s::TrackedSat, sel...) =
    get_preferred_num_code_blocks_to_integrate(_find_signal(s.signals, sel...))

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
Type alias for a tuple or named tuple of per-group satellite dictionaries.
Each entry is `Dictionary{I, <:TrackedSat}` — the keys are satellite
identifiers (PRNs) and the values are the per-sat tracking state. The signal
type for each group lives in the dictionary value type, accessed via
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
$(SIGNATURES)

User-facing outer constructor: build a fresh `SignalGroup` from a signal
tuple with band and antenna count as kwargs. `band` defaults to
`get_band(first(signals))` (so users only override for the rare case of
naming a band differently); `num_ants` defaults to `NumAnts(1)`.

The `satellites` dictionary is left empty — populate it via
[`add_satellite!`](@ref) after the enclosing [`TrackState`](@ref) is
built. The signal-tuple shape determines the dict's concrete value type,
so the slot is type-stable.

```julia
SignalGroup((GPSL1CA(),))                              # band L1(), 1 antenna
SignalGroup((GPSL5I(),); num_ants = NumAnts(2))        # 2-antenna L5
```
"""
function SignalGroup(
    signals::Tuple{Vararg{AbstractGNSSSignal}};
    band = get_band(first(signals)),
    num_ants::NumAnts = NumAnts(1),
    doppler_estimator::AbstractDopplerEstimator =
        ConventionalAssistedPLLAndDLL(;
            carrier_loop_filter_bandwidth =
                default_carrier_loop_filter_bandwidth(first(signals)),
            code_loop_filter_bandwidth =
                default_code_loop_filter_bandwidth(first(signals)),
        ),
)
    # Build a template TrackedSat so the dict's value type is concrete.
    # Reuses the existing helper from tracking_state.jl, which is fine
    # because doppler_estimator only affects per-sat state type, not the
    # storage's outer shape.
    template = _make_template_tracked_sat(signals, doppler_estimator, num_ants)
    sats = Dictionary{Int, typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, signals, num_ants)
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
    group_idx,
    new_tracked_sats::Dictionary{<:Any,<:TrackedSat},
)
    sats_dict = satellites[group_idx]
    @set satellites[group_idx] = merge(sats_dict, new_tracked_sats)
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

# Groups-shape variant: produce a fresh `SignalGroups` where each `SignalGroup`
# reuses the original's band / signals / num_ants but holds a Dictionary whose
# `values` vector is freshly copied. Used by the immutable
# downconvert/estimator/reset paths to detach storage from the input
# `TrackState` before delegating to the in-place form.
@inline _copy_group_slot_vectors(g::SignalGroup) =
    SignalGroup(g; satellites = _copy_slot_vector(g.satellites))

@inline _copy_groups_slot_vectors(groups::SignalGroups) =
    map(_copy_group_slot_vectors, groups)

# Recursive tuple walker: applies `f(sats_dict, args...)` to each per-group
# satellite dictionary in the (named-)tuple. Each step has fully concrete
# types, so no runtime dispatch even when groups differ (heterogeneous
# `Tuple{Dictionary{Int, TrackedSat{Tuple{TrackedSignal{GPSL1CA,...}},...}}, Dictionary{Int, TrackedSat{Tuple{TrackedSignal{GalileoE1B,...}},...}}}`
# would otherwise box each element when iterated with `for s in tuple`).
@inline _foreach_group!(f::F, ::Tuple{}, args::Vararg{Any,N}) where {F,N} = nothing
@inline function _foreach_group!(f::F, t::Tuple, args::Vararg{Any,N}) where {F,N}
    f(first(t), args...)
    _foreach_group!(f, Base.tail(t), args...)
end
# NamedTuple unwraps to its underlying Tuple — concrete and cheap.
@inline _foreach_group!(f::F, nt::NamedTuple, args::Vararg{Any,N}) where {F,N} =
    _foreach_group!(f, Tuple(nt), args...)

# In-place variant: walks each group's `Vector{TrackedSat}` slot storage
# and overwrites each entry with a freshly reset value. The vector
# itself, the Dictionary, the SignalGroup, and the enclosing NamedTuple
# are all reused. `TrackedSat` itself is immutable, so the slot is
# reassigned rather than mutated; we still call the non-`!` per-
# `TrackedSat` form.
@inline function _reset_one_group!(g::SignalGroup)
    vals = g.satellites.values
    @inbounds for i in eachindex(vals)
        vals[i] = reset_start_sample_and_bit_buffer(vals[i])
    end
    return nothing
end

function reset_start_sample_and_bit_buffer!(
    groups::SignalGroups,
)
    _foreach_group!(_reset_one_group!, groups)
    return groups
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

function estimate_cn0(tsig::TrackedSignal)
    signal = get_signal(tsig)
    estimate_cn0(
        get_cn0_estimator(tsig),
        get_code_length(signal) / get_code_frequency(signal),
    )
end

# Per-signal selectors handled by `_find_signal` (no selector → only-fold;
# Integer → indexed; Type{T} → unique-type match).
estimate_cn0(sat::TrackedSat, sel...) = estimate_cn0(_find_signal(sat.signals, sel...))
