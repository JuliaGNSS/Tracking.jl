"""
$(SIGNATURES)

Per-signal tracking state, one per signal tracked on a satellite: correlator,
post-correlation filter, CN0 estimator, bit buffer and integration progress. The
shared carrier/code Doppler and phase live on the enclosing [`TrackedSat`](@ref).
"""
struct TrackedSignal{
    Sig<:AbstractGNSSSignal,
    B<:Unsigned,
    C<:AbstractCorrelator,
    PCF<:AbstractPostCorrFilter,
    CN0<:AbstractCN0Estimator,
}
    signal::Sig
    integrated_samples::Int
    correlator::C
    last_fully_integrated_correlator::C
    last_fully_integrated_filtered_prompt::ComplexF64
    cn0_estimator::CN0
    bit_buffer::BitBuffer{B}
    post_corr_filter::PCF
    filtered_prompts::Vector{ComplexF64}
    # Records completed within the current chunk; preallocated, filled by the
    # correlate phase and emptied by the Doppler estimator (see
    # `get_correlator_outputs`).
    correlator_outputs::Vector{CorrelatorOutput{C}}
    # Preferred coherent-integration length in primary code blocks; see
    # `set_preferred_num_code_blocks_to_integrate!`.
    preferred_num_code_blocks_to_integrate::Int
    # Primary-code blocks the most recently folded record actually spanned, set
    # by `_apply_correlator_output` from the record's sample count (so it also
    # holds for external producers). `estimate_cn0` needs it: the buffered prompts
    # are sample-normalized, so ignoring N over-reports C/N₀ by 10·log₁₀(N).
    last_fully_integrated_num_code_blocks::Int
end

# Reject a preferred coherent-integration length that cannot work for this
# signal. For data-bearing signals it must evenly divide the blocks per bit, or
# integrations straddle bit boundaries and no bit is ever emitted (issue #128).
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
                "that form one bit of $(get_signal_name(signal)) — an integration " *
                "straddling a bit boundary would never emit a bit. Valid values: " *
                "$(join(valid, ", ")).",
            ),
        )
    end
    nothing
end

"""
$(SIGNATURES)

Construct a fresh [`TrackedSignal`](@ref) for `signal`. The correlator,
post-corr filter and CN0 estimator default to the signal's recommended values;
pass `correlator`, `post_corr_filter` and `cn0_estimator` explicitly to
override.

`cn0_estimator` accepts any [`AbstractCN0Estimator`](@ref) (it is a type
parameter of `TrackedSignal`). The default is [`default_cn0_estimator`](@ref);
the alternatives are on the [CN0 Estimator](@ref) page. On a correlator-ingest path feed it with
[`append_noise_observation!`](@ref), or C/N₀ stays at `-Inf dB-Hz`. Each signal
needs its **own** estimator instance: estimators buffer into a vector, so sharing
one corrupts both.

`preferred_num_code_blocks_to_integrate` defaults to
[`default_num_code_blocks_to_integrate`](@ref); an invalid value throws an
`ArgumentError` (see [`set_preferred_num_code_blocks_to_integrate!`](@ref)).
"""
function TrackedSignal(
    signal::AbstractGNSSSignal;
    num_ants::NumAnts = NumAnts(1),
    correlator::AbstractCorrelator = get_default_correlator(signal, num_ants),
    num_prompts_for_cn0_estimation::Int = 100,
    cn0_estimator::AbstractCN0Estimator = default_cn0_estimator(
        signal,
        num_prompts_for_cn0_estimation,
    ),
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    preferred_num_code_blocks_to_integrate::Int = default_num_code_blocks_to_integrate(
        signal,
    ),
)
    validate_preferred_num_code_blocks_to_integrate(
        signal,
        preferred_num_code_blocks_to_integrate,
    )
    # Per-signal sync-search buffer width (see `get_code_block_buffer_type`).
    B = get_code_block_buffer_type(signal)
    # A default chunk yields at most one record; a small sizehint keeps `push!`
    # from growing it even for a moderately enlarged `doppler_update_interval`.
    correlator_outputs = CorrelatorOutput{typeof(correlator)}[]
    sizehint!(correlator_outputs, 4)
    TrackedSignal(
        signal,
        0,
        correlator,
        correlator,
        complex(0.0, 0.0),
        cn0_estimator,
        BitBuffer{B}(),
        post_corr_filter,
        ComplexF64[],
        correlator_outputs,
        preferred_num_code_blocks_to_integrate,
        1,
    )
end

# Kwarg-update constructor; keeps `t`'s concrete correlator and PCF types.
# `cn0_estimator` is deliberately not pinned, so a custom estimator can be swapped
# in; call sites pass `nothing` or a concrete estimator, so the type still infers
# and `_apply_correlator_output` stays allocation-free.
function TrackedSignal(
    t::TrackedSignal{Sig,B,C,PCF};
    signal = nothing,
    integrated_samples = nothing,
    correlator::Maybe{C} = nothing,
    last_fully_integrated_correlator::Maybe{C} = nothing,
    last_fully_integrated_filtered_prompt = nothing,
    cn0_estimator::Maybe{AbstractCN0Estimator} = nothing,
    bit_buffer::Maybe{BitBuffer{B}} = nothing,
    post_corr_filter::Maybe{PCF} = nothing,
    filtered_prompts::Maybe{Vector{ComplexF64}} = nothing,
    correlator_outputs::Maybe{Vector{CorrelatorOutput{C}}} = nothing,
    preferred_num_code_blocks_to_integrate = nothing,
    last_fully_integrated_num_code_blocks = nothing,
) where {
    Sig<:AbstractGNSSSignal,
    B<:Unsigned,
    C<:AbstractCorrelator,
    PCF<:AbstractPostCorrFilter,
}
    isnothing(preferred_num_code_blocks_to_integrate) ||
        validate_preferred_num_code_blocks_to_integrate(
            isnothing(signal) ? t.signal : signal,
            preferred_num_code_blocks_to_integrate,
        )
    new_cn0_estimator = isnothing(cn0_estimator) ? t.cn0_estimator : cn0_estimator
    TrackedSignal{Sig,B,C,PCF,typeof(new_cn0_estimator)}(
        isnothing(signal) ? t.signal : signal,
        isnothing(integrated_samples) ? t.integrated_samples : integrated_samples,
        isnothing(correlator) ? t.correlator : correlator,
        isnothing(last_fully_integrated_correlator) ? t.last_fully_integrated_correlator :
        last_fully_integrated_correlator,
        isnothing(last_fully_integrated_filtered_prompt) ?
        t.last_fully_integrated_filtered_prompt : last_fully_integrated_filtered_prompt,
        new_cn0_estimator,
        isnothing(bit_buffer) ? t.bit_buffer : bit_buffer,
        isnothing(post_corr_filter) ? t.post_corr_filter : post_corr_filter,
        isnothing(filtered_prompts) ? t.filtered_prompts : filtered_prompts,
        isnothing(correlator_outputs) ? t.correlator_outputs : correlator_outputs,
        isnothing(preferred_num_code_blocks_to_integrate) ?
        t.preferred_num_code_blocks_to_integrate : preferred_num_code_blocks_to_integrate,
        isnothing(last_fully_integrated_num_code_blocks) ?
        t.last_fully_integrated_num_code_blocks : last_fully_integrated_num_code_blocks,
    )
end

get_signal(t::TrackedSignal) = t.signal
get_correlator(t::TrackedSignal) = t.correlator
get_last_fully_integrated_correlator(t::TrackedSignal) = t.last_fully_integrated_correlator
get_last_fully_integrated_filtered_prompt(t::TrackedSignal) =
    t.last_fully_integrated_filtered_prompt
get_last_fully_integrated_num_code_blocks(t::TrackedSignal) =
    t.last_fully_integrated_num_code_blocks

"""
$(SIGNATURES)

The integration time of the most recently completed record: its primary-code
block count times one code period.

This is the `T` [`estimate_cn0`](@ref) divides by, and what turns C/N₀ into the
record's post-integration SNR (`C/N₀ · T`). Unlike `get_integrated_samples`, it
describes the last *completed* record, not the one being accumulated.
"""
get_last_fully_integrated_integration_time(t::TrackedSignal) =
    get_last_fully_integrated_num_code_blocks(t) * get_code_length(get_signal(t)) /
    get_code_frequency(get_signal(t))

get_filtered_prompts(t::TrackedSignal) = t.filtered_prompts

"""
$(SIGNATURES)

The [`CorrelatorOutput`](@ref)s this signal completed during the most recent
processing chunk, in order. Populated by the correlate phase and consumed +
cleared by the Doppler estimator after each chunk, so it is empty between
`track!` calls; read it inside a custom estimator, or right after a bare
`downconvert_and_correlate!`. See [Chunked Doppler updates](@ref).
"""
get_correlator_outputs(t::TrackedSignal) = t.correlator_outputs

"""
$(SIGNATURES)

Append an externally built [`CorrelatorOutput`](@ref) to `signal`'s
per-chunk `correlator_outputs` buffer and return `signal`.

This is the ingest path for an **external correlator producer** (e.g. an FPGA):
append records per signal in `sample_index` order, then run
[`estimate_dopplers_and_filter_prompt!`](@ref), which folds and clears them.
Prefer it over mutating [`get_correlator_outputs`](@ref) directly; it is
type-checked. See [External correlator producers](@ref) for the full contract.
"""
append_correlator_output!(t::TrackedSignal, output::CorrelatorOutput) =
    (push!(t.correlator_outputs, output); t)

get_post_corr_filter(t::TrackedSignal) = t.post_corr_filter
get_cn0_estimator(t::TrackedSignal) = t.cn0_estimator
get_bit_buffer(t::TrackedSignal) = t.bit_buffer
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
The first signal in `signals` is the **estimator-driver signal**, whose
correlator feeds the Doppler estimator (see [Estimator-driver signal](@ref)).

The shared `code_phase` wraps at the least common multiple of the signals'
code periods including secondary code (see [`max_code_length`](@ref)).
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

# Post-sync wrap length of one signal: a full secondary-code period for a pilot;
# for a data signal the larger of the secondary-code and data-bit periods (e.g.
# GPS L1 C/A: 1023 × 20 chips, no secondary code).
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

# Folded with `lcm`, not `max` — see `current_code_wrap` (issue #129).
@inline _max_code_length(::Tuple{}) = 1
@inline _max_code_length(t::Tuple) =
    lcm(_post_sync_code_length(first(t)), _max_code_length(Base.tail(t)))

"""
$(SIGNATURES)

Upper bound on the shared `sat.code_phase` wrap period, in chips —
the least common multiple of the per-signal wrap periods once *every*
signal on the sat has synced (e.g. 1023 × 20 = 20460 for GPS L1 C/A alone).

This is the compile-time bound (it folds to a literal for a concrete `signals`
type); see [`current_code_wrap`](@ref) for the runtime value, which honors sync
state.
"""
@inline max_code_length(signals::Tuple{Vararg{TrackedSignal}}) = _max_code_length(signals)

# Runtime per-signal wrap contribution; see `current_code_wrap`.
@inline function _current_code_length(tsig::TrackedSignal)
    if tsig.bit_buffer.found
        _post_sync_code_length(tsig)
    else
        get_code_length(tsig.signal)
    end
end

@inline _current_code_wrap(::Tuple{}) = 1
@inline _current_code_wrap(t::Tuple) =
    lcm(_current_code_length(first(t)), _current_code_wrap(Base.tail(t)))

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
    primary code length, since the bit / secondary chip is still unknown.

The shared wrap is the least common multiple of the per-signal
contributions, not their maximum: per-signal phases are re-derived as
`mod(code_phase, replica_wrap)`, which a non-common-multiple wrap would
silently corrupt (issue #129). For shipped pairings the two coincide.
"""
@inline current_code_wrap(signals::Tuple{Vararg{TrackedSignal}}) =
    _current_code_wrap(signals)

# Whether any signal's `bit_buffer.found` went `false` → `true` this iteration.
# Gates the one-time code-phase snap, valid only at the sync transition (see
# `_update_tracked_sat_doppler` and issue #117).
@inline _any_signal_just_synced(::Tuple{}, ::Tuple{}) = false
@inline _any_signal_just_synced(old::Tuple, new::Tuple) =
    (!first(old).bit_buffer.found && first(new).bit_buffer.found) ||
    _any_signal_just_synced(Base.tail(old), Base.tail(new))

# Snap `code_phase` into the secondary-code window given by the synced signal
# with the longest `primary × secondary` length (its `secondary_phase` fixes the
# absolute position). Unchanged if no synced signal has a secondary code.
@inline function _snap_code_phase_from_synced_signal(
    signals::Tuple{Vararg{TrackedSignal}},
    code_phase::Float64,
)
    best_len, best_phase_chips, best_prim = _find_best_secondary_anchor(signals, 0, 0, 1)
    if best_len == 0
        return code_phase
    end
    # Align the low `best_len` chips to the synced secondary window, keeping the
    # higher-order wrap and the within-primary-block phase — dropping the latter
    # would inject a phase error at sync (issue #117).
    base = floor(Int, code_phase / best_len) * best_len
    within_primary = code_phase - floor(code_phase / best_prim) * best_prim
    Float64(base + best_phase_chips) + within_primary
end

# Returns `(length, secondary_phase_in_chips, primary_length)` of the best anchor,
# or `(0, 0, 1)` if none.
@inline _find_best_secondary_anchor(
    ::Tuple{},
    best_len::Int,
    best_chips::Int,
    best_prim::Int,
) = (best_len, best_chips, best_prim)
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
        sec = get_secondary_code_length(sig)
        total = prim * sec
        # Signals without a secondary code (e.g. GPS L1 C/A) don't pin the window.
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

Core constructor: build a [`TrackedSat`](@ref) from pre-built
[`TrackedSignal`](@ref)s plus acquisition handoff values. The first signal in
the tuple is the estimator-driver signal — it scales the default
`code_doppler` and sizes the default loop bandwidths. The estimator's
per-satellite state is built via [`init_estimator_state`](@ref).

`code_doppler = nothing` (the default) derives the code Doppler from
`carrier_doppler` and the driver signal's code/center frequency ratio.
`carrier_phase` is in radians.

Use this form when a signal needs a non-default correlator or post-corr
filter; otherwise prefer the signal-instance forms
`TrackedSat(signal, prn, ...)` / `TrackedSat((sig_a, sig_b, ...), prn, ...)`,
which build the `TrackedSignal`s with the library defaults.
"""
function TrackedSat(
    tracked_signals::Tuple{TrackedSignal,Vararg{TrackedSignal}},
    prn::Int,
    code_phase,
    carrier_doppler;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
    carrier_phase = 0.0,
    code_doppler = nothing,
)
    # Float-ize so an integer `200Hz` fits the `typeof(1.0Hz)` fields.
    cdop = float(carrier_doppler)
    cd =
        isnothing(code_doppler) ?
        cdop * get_code_center_frequency_ratio(first(tracked_signals).signal) :
        float(code_doppler)
    # Two-stage build so `init_estimator_state` sees a real `TrackedSat`.
    bare = TrackedSat(
        prn,
        float(code_phase),
        cd,
        float(carrier_phase) / 2π,
        cdop,
        1,
        tracked_signals,
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

"""
$(SIGNATURES)

Public multi-signal constructor: build a [`TrackedSat`](@ref) tracking the
given tuple of signal instances, e.g.
`TrackedSat((GPSL1C_P(), GPSL1C_D(), GPSL1CA()), prn, code_phase, carrier_doppler)`.
Each signal is wrapped in a [`TrackedSignal`](@ref) with its recommended
default correlator. The first signal is the estimator-driver signal. Further
kwargs (`doppler_estimator`, `carrier_phase`, `code_doppler`) forward to the
`TrackedSignal`-tuple core constructor above.

`cn0_estimator` takes a **tuple** of estimators, one per signal in the same
order (each signal needs its own instance); a single estimator is accepted only
for a single-signal tuple.
"""
function TrackedSat(
    signals::Tuple{AbstractGNSSSignal,Vararg{AbstractGNSSSignal}},
    prn::Int,
    code_phase,
    carrier_doppler;
    num_ants::NumAnts = NumAnts(1),
    num_prompts_for_cn0_estimation::Int = 100,
    cn0_estimator = nothing,
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    kwargs...,
)
    cn0_estimators =
        _per_signal_cn0_estimators(signals, cn0_estimator, num_prompts_for_cn0_estimation)
    tracked_signals = map(signals, cn0_estimators) do sig, estimator
        TrackedSignal(
            sig;
            num_ants,
            correlator = get_default_correlator(sig, num_ants),
            cn0_estimator = estimator,
            post_corr_filter,
        )
    end
    TrackedSat(tracked_signals, prn, code_phase, carrier_doppler; kwargs...)
end

# Resolve the multi-signal `cn0_estimator` kwarg into one estimator per signal.
@inline _per_signal_cn0_estimators(
    signals::Tuple{Vararg{AbstractGNSSSignal}},
    ::Nothing,
    num_prompts_for_cn0_estimation::Int,
) = map(sig -> default_cn0_estimator(sig, num_prompts_for_cn0_estimation), signals)

@inline _per_signal_cn0_estimators(
    ::Tuple{AbstractGNSSSignal},
    cn0_estimator::AbstractCN0Estimator,
    ::Int,
) = (cn0_estimator,)

@noinline _per_signal_cn0_estimators(
    signals::Tuple{Vararg{AbstractGNSSSignal}},
    ::AbstractCN0Estimator,
    ::Int,
) = throw(
    ArgumentError(
        "a single cn0_estimator cannot be shared by the $(length(signals)) signals " *
        "of this satellite — the estimators buffer into a vector they would then " *
        "both write to. Pass one estimator per signal as a tuple, e.g. " *
        "`cn0_estimator = ($(join(fill("MomentsCN0Estimator(100)", length(signals)), ", ")))`.",
    ),
)

@inline function _per_signal_cn0_estimators(
    signals::Tuple{Vararg{AbstractGNSSSignal}},
    cn0_estimators::Tuple{Vararg{AbstractCN0Estimator}},
    ::Int,
)
    length(cn0_estimators) == length(signals) || throw(
        ArgumentError(
            "got $(length(cn0_estimators)) cn0_estimators for $(length(signals)) " *
            "signals — pass exactly one estimator per signal, in the same order.",
        ),
    )
    cn0_estimators
end

"""
$(SIGNATURES)

Construct a single-signal [`TrackedSat`](@ref) from acquisition handoff values
plus the Doppler estimator. The signal is wrapped in a single
[`TrackedSignal`](@ref); pass a tuple of signals for multi-signal tracking.
"""
function TrackedSat(
    signal::AbstractGNSSSignal,
    prn::Int,
    code_phase,
    carrier_doppler;
    num_ants::NumAnts = NumAnts(1),
    correlator::AbstractCorrelator = get_default_correlator(signal, num_ants),
    num_prompts_for_cn0_estimation::Int = 100,
    cn0_estimator::AbstractCN0Estimator = default_cn0_estimator(
        signal,
        num_prompts_for_cn0_estimation,
    ),
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    kwargs...,
)
    tracked_signal =
        TrackedSignal(signal; num_ants, correlator, cn0_estimator, post_corr_filter)
    TrackedSat((tracked_signal,), prn, code_phase, carrier_doppler; kwargs...)
end

# Kwarg-update constructor; `signals` and `doppler_estimator_state` keep their
# concrete types so the enclosing TrackState stays inferable.
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
get_num_ants(
    s::TrackedSat{<:Tuple{TrackedSignal{<:Any,<:Any,<:AbstractCorrelator{M}},Vararg}},
) where {M} = M

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

# Per-signal accessors on a TrackedSat, routed through `_find_signal`:
#   * no selector — single-signal sats only.
#   * `Integer` index — always unambiguous.
#   * signal type — the unique signal of that type; errors on zero or >1
#     matches. Folds at compile time for a concrete `Signals` type.
@noinline _throw_needs_signal_selector() = throw(
    ArgumentError(
        "satellite tracks multiple signals — pass a signal selector " *
        "(integer index or signal type) to address one of them.",
    ),
)
@inline _find_signal(s::Tuple{TrackedSignal}) = s[1]
@inline _find_signal(::Tuple) = _throw_needs_signal_selector()
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
@inline _assert_no_more_of_type(::Tuple{}, ::Type{T}) where {T<:AbstractGNSSSignal} =
    nothing
@inline function _assert_no_more_of_type(t::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    first(t).signal isa T && throw(
        ArgumentError(
            "signal type $T matches more than one signal on this satellite — " *
            "use an integer index to disambiguate (e.g. `get_signals(sat)[i]`).",
        ),
    )
    _assert_no_more_of_type(Base.tail(t), T)
end

get_signal(s::TrackedSat, sel...) = get_signal(_find_signal(s.signals, sel...))
get_correlator(s::TrackedSat, sel...) = get_correlator(_find_signal(s.signals, sel...))
get_last_fully_integrated_correlator(s::TrackedSat, sel...) =
    get_last_fully_integrated_correlator(_find_signal(s.signals, sel...))
get_last_fully_integrated_filtered_prompt(s::TrackedSat, sel...) =
    get_last_fully_integrated_filtered_prompt(_find_signal(s.signals, sel...))
get_last_fully_integrated_num_code_blocks(s::TrackedSat, sel...) =
    get_last_fully_integrated_num_code_blocks(_find_signal(s.signals, sel...))
get_last_fully_integrated_integration_time(s::TrackedSat, sel...) =
    get_last_fully_integrated_integration_time(_find_signal(s.signals, sel...))
get_filtered_prompts(s::TrackedSat, sel...) =
    get_filtered_prompts(_find_signal(s.signals, sel...))
get_post_corr_filter(s::TrackedSat, sel...) =
    get_post_corr_filter(_find_signal(s.signals, sel...))
get_cn0_estimator(s::TrackedSat, sel...) =
    get_cn0_estimator(_find_signal(s.signals, sel...))
get_bit_buffer(s::TrackedSat, sel...) = get_bit_buffer(_find_signal(s.signals, sel...))
get_soft_bits(s::TrackedSat, sel...) = get_soft_bits(_find_signal(s.signals, sel...))
get_num_bits(s::TrackedSat, sel...) = get_num_bits(_find_signal(s.signals, sel...))
has_bit_or_secondary_code_been_found(s::TrackedSat, sel...) =
    has_bit_or_secondary_code_been_found(_find_signal(s.signals, sel...))
get_integrated_samples(s::TrackedSat, sel...) =
    get_integrated_samples(_find_signal(s.signals, sel...))
get_correlator_outputs(s::TrackedSat, sel...) =
    get_correlator_outputs(_find_signal(s.signals, sel...))
get_preferred_num_code_blocks_to_integrate(s::TrackedSat, sel...) =
    get_preferred_num_code_blocks_to_integrate(_find_signal(s.signals, sel...))

# `output` comes first so a trailing signal selector works as for the accessors.
append_correlator_output!(s::TrackedSat, output::CorrelatorOutput, sel...) =
    (append_correlator_output!(_find_signal(s.signals, sel...), output); s)

# Reset between `track` calls: `signal_start_sample` back to 1, per-signal
# buffers emptied in place and bit buffers reset.
function reset_start_sample_and_bit_buffer(sat::TrackedSat)
    new_signals = map(s -> _reset_signal(s), sat.signals)
    TrackedSat(sat; signal_start_sample = 1, signals = new_signals)
end

@inline function _reset_signal(t::TrackedSignal)
    empty!(t.filtered_prompts)
    empty!(t.correlator_outputs)
    TrackedSignal(t; bit_buffer = reset(t.bit_buffer))
end

"""
$(SIGNATURES)

Build the per-satellite Doppler-estimator state used by `estimator` for the
given satellite. A custom doppler estimator must define this method for its
[`AbstractDopplerEstimator`](@ref) subtype.

This function must be **pure**: it is also called on throwaway template sats
(to fix slot types at [`TrackState`](@ref) construction), as a type probe, and
by [`reset_loop_filters!`](@ref). Register satellites in cross-satellite shared
state in [`update_estimator_on_handoff`](@ref) instead, which runs exactly once
per handoff with the real incoming satellites.
"""
function init_estimator_state end

"""
$(SIGNATURES)

Optionally update an estimator's cross-satellite or cross-system shared state
when new satellites enter the track set. Called once per handoff entry point
([`merge_sats`](@ref), [`add_satellite`](@ref), [`add_satellite!`](@ref))
with the dictionary of incoming satellites, *after* per-sat seeding via
[`init_estimator_state`](@ref).

The default returns `estimator` unchanged, so estimators with no shared state
need not implement it.

The returned estimator must have the same concrete type as the input
([`TrackState`](@ref) is parameterized on it). For growing shared state, hold
resizable storage (e.g. `Vector`) on the estimator and `push!`/`resize!` it in
place, or rebuild it with
[`Setfield.@set`](https://jw3126.github.io/Setfield.jl/stable/) or a copying
constructor.

Every entry point's returned `TrackState` carries the returned estimator, so
callers of the in-place [`add_satellite!`](@ref) must keep using its return value.
"""
update_estimator_on_handoff(estimator::AbstractDopplerEstimator, _new_sats) = estimator

"""
Type alias for a tuple or named tuple of per-group satellite dictionaries,
each a `Dictionary{I, <:TrackedSat}` keyed by satellite identifier (PRN).
"""
const SatelliteDicts{N} = TupleLike{<:NTuple{N,Dictionary{<:Any,<:TrackedSat}}}

"""
$(SIGNATURES)

A group of satellites that all track the same tuple of GNSS signal types,
on the same RF band, observed by the same antenna array.

Groups are the unit of type stability (one concrete `TrackedSat` type per
group). Several groups may share a band; the band only routes the right
measurement to the group during `track`.

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

    # Unchecked, for rebuilding a group from one already validated; the positional
    # `SignalGroup(band, satellites, signals, num_ants)` validates.
    SignalGroup{B,S,Sigs,NA}(band, satellites, signals, num_ants) where {B,S,Sigs,NA} =
        new{B,S,Sigs,NA}(band, satellites, signals, num_ants)
end

# Every group formed from its parts is validated here, also one built positionally
# rather than through the keyword constructor.
function SignalGroup(
    band::B,
    satellites::S,
    signals::Sigs,
    num_ants::NA,
) where {
    B,
    S<:Dictionary{<:Any,<:TrackedSat},
    Sigs<:Tuple{Vararg{AbstractGNSSSignal}},
    NA<:NumAnts,
}
    _validate_signal_group(signals, band)
    SignalGroup{B,S,Sigs,NA}(band, satellites, signals, num_ants)
end

# Kwarg-update constructor; keeps `g`'s concrete types.
function SignalGroup(
    g::SignalGroup{B,S,Sigs,NA};
    band::Maybe{B} = nothing,
    satellites::Maybe{S} = nothing,
    signals::Maybe{Sigs} = nothing,
    num_ants::Maybe{NA} = nothing,
) where {
    B,
    S<:Dictionary{<:Any,<:TrackedSat},
    Sigs<:Tuple{Vararg{AbstractGNSSSignal}},
    NA<:NumAnts,
}
    SignalGroup{B,S,Sigs,NA}(
        isnothing(band) ? g.band : band,
        isnothing(satellites) ? g.satellites : satellites,
        isnothing(signals) ? g.signals : signals,
        isnothing(num_ants) ? g.num_ants : num_ants,
    )
end

# A group's signals share one band (#129), one chip rate (#129) and one
# constellation (#224), and are all different; the error messages say why. The
# constellation stands in for the PRN namespace: move to a PRN-namespace trait
# once QZSS/SBAS share GPS's code space. Bands compare by id, so a user-defined
# band aliasing an existing measurement key still validates.
# The instance method is generic and forwards to the type, so look for the latter.
_has_constellation_id(s::AbstractGNSSSignal) =
    hasmethod(get_constellation_id, Tuple{Type{typeof(s)}})

function _require_constellation_id(s::AbstractGNSSSignal)
    _has_constellation_id(s) || throw(
        ArgumentError(
            string(
                "Every signal in a SignalGroup of more than one signal must ",
                "declare its constellation, so that the group cannot mix ",
                "constellations, but `",
                get_signal_id(s),
                "` has no `get_constellation_id` method. Define ",
                "`GNSSSignals.get_constellation_id(::Type{",
                nameof(typeof(s)),
                "})`, or subtype `AbstractGPSSignal`, ",
                "`AbstractGalileoSignal` or `AbstractBeiDouSignal`.",
            ),
        ),
    )
    nothing
end

function _validate_signal_group(signals::Tuple{Vararg{AbstractGNSSSignal}}, band)
    driver = first(signals)
    for (i, s) in enumerate(signals)
        any(t -> typeof(t) === typeof(s), signals[1:(i-1)]) && throw(
            ArgumentError(
                string(
                    "A SignalGroup lists `",
                    get_signal_id(s),
                    "` twice. Each signal is tracked once per satellite, and a ",
                    "signal-type selector must name exactly one of them: list it once.",
                ),
            ),
        )
    end
    foreach(signals) do s
        if get_band_id(get_band(s)) !== get_band_id(band)
            throw(
                ArgumentError(
                    string(
                        "All signals in a SignalGroup must be on the same RF band: `",
                        get_signal_id(s),
                        "` is on band `:",
                        get_band_id(get_band(s)),
                        "` but the group is on band `:",
                        get_band_id(band),
                        "`. Put signals on different bands into separate groups, e.g. ",
                        "`signals = (gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),))`.",
                    ),
                ),
            )
        end
        if get_code_frequency(s) != get_code_frequency(driver)
            throw(
                ArgumentError(
                    string(
                        "All signals in a SignalGroup must share one chip rate: the ",
                        "shared code phase advances at the first signal's code ",
                        "frequency (`",
                        get_signal_id(driver),
                        "`: ",
                        get_code_frequency(driver),
                        "), but `",
                        get_signal_id(s),
                        "` has ",
                        get_code_frequency(s),
                        ". Track it in its own group instead.",
                    ),
                ),
            )
        end
    end
    length(signals) > 1 && _require_constellation_id(driver)
    foreach(Base.tail(signals)) do s
        _require_constellation_id(s)
        if get_constellation_id(s) !== get_constellation_id(driver)
            throw(
                ArgumentError(
                    string(
                        "All signals in a SignalGroup must come from one ",
                        "constellation: `",
                        get_signal_id(s),
                        "` is a `:",
                        get_constellation_id(s),
                        "` signal but `",
                        get_signal_id(driver),
                        "` a `:",
                        get_constellation_id(driver),
                        "` one. A group's satellites carry one PRN and one ",
                        "Doppler, so a second constellation's signal would track a ",
                        "different satellite. Put them into separate groups, e.g. ",
                        "`signals = (gps_l1 = (GPSL1CA(),), galileo_e1 = (GalileoE1B(),))`.",
                    ),
                ),
            )
        end
    end
    nothing
end

"""
$(SIGNATURES)

User-facing outer constructor: build a fresh `SignalGroup` from a signal
tuple with band and antenna count as kwargs. `band` defaults to
`get_band(first(signals))` (so users only override for the rare case of
naming a band differently); `num_ants` defaults to `NumAnts(1)`.

All signals must live on `band`, share one chip rate and come from one
constellation (a group's satellites carry one PRN and one Doppler); the
constructor throws an `ArgumentError` otherwise. In a group of more than one
signal, a user-defined signal type therefore needs a
`GNSSSignals.get_constellation_id(::Type{MySignal})` method, or to subtype
`AbstractGPSSignal`, `AbstractGalileoSignal` or `AbstractBeiDouSignal`.

The `satellites` dictionary is left empty — populate it via
[`add_satellite!`](@ref) after the enclosing [`TrackState`](@ref) is
built. The signal-tuple shape determines the dict's concrete value type,
so the slot is type-stable.

The `doppler_estimator` kwarg only shapes the dictionary's slot type
(via the per-sat estimator-state type); the estimator that actually runs
is the one configured on the enclosing `TrackState`. When the two
disagree, `TrackState` rebuilds the empty slot to match its own
estimator, so pass the same estimator instance to both — or simply leave
this kwarg at its default and configure the estimator on `TrackState`
alone.

```julia
SignalGroup((GPSL1CA(),))                              # band L1(), 1 antenna
SignalGroup((GPSL5I(),); num_ants = NumAnts(2))        # 2-antenna L5
```
"""
function SignalGroup(
    signals::Tuple{Vararg{AbstractGNSSSignal}};
    band = get_band(first(signals)),
    num_ants::NumAnts = NumAnts(1),
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    # Before the template satellite, which an invalid signal would fail with a
    # less telling error; the positional constructor checks again.
    _validate_signal_group(signals, band)
    # Template TrackedSat fixes the dict's concrete value type.
    template = _make_template_tracked_sat(signals, doppler_estimator, num_ants)
    sats = Dictionary{Int,typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, signals, num_ants)
end

"""
Type alias: NamedTuple of `SignalGroup`s — the storage shape inside
[`TrackState`](@ref). The `N` parameter is the number of groups.
"""
const SignalGroups{N} = NamedTuple{<:Any,<:NTuple{N,SignalGroup}}

# Copy that shares the keys (`Indices`) but copies `values`. Used by the
# immutable per-iteration steps (`downconvert_and_correlate`,
# `estimate_dopplers_and_filter_prompt`), which never change the key set, to
# avoid copying the hash table every iteration. Keys are detached once, at the
# `track` boundary (`_detach_slot_vector`, #123); calling these steps directly
# and then `add_satellite!`/`remove_satellite!` on the result would corrupt the
# input's keys.
@inline function _copy_slot_vector(sats::Dictionary{<:Any,<:TrackedSat})
    Dictionary(keys(sats), copy(sats.values))
end

# Fully-detached copy (keys and values, #123), used once per `track` call by
# `reset_start_sample_and_bit_buffer`.
@inline function _detach_slot_vector(sats::Dictionary{<:Any,<:TrackedSat})
    copy(sats)
end

# Groups-shape variant of `_copy_slot_vector`.
@inline _copy_group_slot_vectors(g::SignalGroup) =
    SignalGroup(g; satellites = _copy_slot_vector(g.satellites))

@inline _copy_groups_slot_vectors(groups::SignalGroups) =
    map(_copy_group_slot_vectors, groups)

# Groups-shape variant of `_detach_slot_vector`.
@inline _detach_group_slot_vectors(g::SignalGroup) =
    SignalGroup(g; satellites = _detach_slot_vector(g.satellites))

@inline _detach_groups_slot_vectors(groups::SignalGroups) =
    map(_detach_group_slot_vectors, groups)

# Apply `f(group, args...)` to each element of a (named) tuple via recursion, so
# heterogeneous groups get concrete types instead of being boxed by a `for` loop.
@inline _foreach_group!(f::F, ::Tuple{}, args::Vararg{Any,N}) where {F,N} = nothing
@inline function _foreach_group!(f::F, t::Tuple, args::Vararg{Any,N}) where {F,N}
    f(first(t), args...)
    _foreach_group!(f, Base.tail(t), args...)
end
# NamedTuple unwraps to its underlying Tuple — concrete and cheap.
@inline _foreach_group!(f::F, nt::NamedTuple, args::Vararg{Any,N}) where {F,N} =
    _foreach_group!(f, Tuple(nt), args...)

# In-place reset: overwrites each (immutable) `TrackedSat` slot in the group's
# values vector, reusing all containers.
@inline function _reset_one_group!(g::SignalGroup)
    vals = g.satellites.values
    @inbounds for i in eachindex(vals)
        vals[i] = reset_start_sample_and_bit_buffer(vals[i])
    end
    return nothing
end

function reset_start_sample_and_bit_buffer!(groups::SignalGroups)
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
    # Divide by the record's integration time, not one code period (see the
    # `last_fully_integrated_num_code_blocks` field).
    estimate_cn0(get_cn0_estimator(tsig), get_last_fully_integrated_integration_time(tsig))
end

# Per-signal selectors handled by `_find_signal`.
estimate_cn0(sat::TrackedSat, sel...) = estimate_cn0(_find_signal(sat.signals, sel...))
