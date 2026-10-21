"""
$(SIGNATURES)

Per-signal tracking state, one per signal tracked on a satellite: the
correlator, the integration progress and the per-record loop state, a
TrackingLoops [`SignalLoopState`](@extref TrackingLoops.SignalLoopState) (post-corr
filter, C/N₀ estimator, bit buffer, and what the last record left for the next).
The carrier/code Doppler and phase are shared and live on the [`TrackedSat`](@ref).
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
    # Advanced per record by TrackingLoops' `apply_record`. Its
    # `last_num_code_blocks` is what the record actually spanned, which the C/N₀
    # estimate divides by.
    loop_state::SignalLoopState{B,PCF,CN0}
    filtered_prompts::Vector{ComplexF64}
    # Records completed in the current chunk: pushed by the correlate phase,
    # folded and emptied by the Doppler estimator.
    correlator_outputs::Vector{CorrelatorOutput{C}}
    # In primary code blocks; see `set_preferred_num_code_blocks_to_integrate!`.
    preferred_num_code_blocks_to_integrate::Int
end

# For data-bearing signals the length must divide the code blocks per bit, or
# integrations straddle bit boundaries and no bit is ever emitted again (#128).
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
post-corr filter and C/N₀ estimator default to the signal's recommended values.

`cn0_estimator` accepts any
[`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator). The default,
[`default_cn0_estimator`](@extref TrackingLoops.default_cn0_estimator), is a
[`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator) over
`num_prompts_for_cn0_estimation` records against the signal's **measured** noise
density: `track!` fills that on the sample-driven path, but a
**correlator-ingest** path must feed [`append_noise_observation!`](@ref) per
signal, or C/N₀ stays at `-Inf dB-Hz`. `NWPRCN0Estimator()` needs no noise
observation; `MomentsCN0Estimator(num_prompts_for_cn0_estimation)` is the plain
moment-ratio estimator. Each signal needs its **own** estimator instance: they
buffer into a vector, so sharing one corrupts both.

`preferred_num_code_blocks_to_integrate` defaults to
[`default_num_code_blocks_to_integrate`](@extref TrackingLoops.default_num_code_blocks_to_integrate)
and must be valid for `signal` (see
[`set_preferred_num_code_blocks_to_integrate!`](@ref)).
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
    # A default chunk yields at most one record; the hint covers a moderately
    # enlarged `doppler_update_interval` without regrowing.
    correlator_outputs = CorrelatorOutput{typeof(correlator)}[]
    sizehint!(correlator_outputs, 4)
    TrackedSignal(
        signal,
        0,
        correlator,
        correlator,
        SignalLoopState(signal; cn0_estimator, post_corr_filter),
        ComplexF64[],
        correlator_outputs,
        preferred_num_code_blocks_to_integrate,
    )
end

# Kwarg-update constructor keeping `t`'s correlator and PCF types. The loop
# state's parts can be replaced one by one; `last_fully_integrated_*` map to its
# `last_filtered_prompt` / `last_num_code_blocks`. `cn0_estimator` may change type,
# so a custom estimator can be swapped in; call sites pass `nothing` or a concrete
# estimator, so the result infers and the per-record rebuild stays allocation-free.
function TrackedSignal(
    t::TrackedSignal{Sig,B,C,PCF};
    signal = nothing,
    integrated_samples = nothing,
    correlator::Maybe{C} = nothing,
    last_fully_integrated_correlator::Maybe{C} = nothing,
    loop_state::Maybe{SignalLoopState{B,PCF}} = nothing,
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
    new_loop_state = SignalLoopState(
        isnothing(loop_state) ? t.loop_state : loop_state;
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        last_filtered_prompt = last_fully_integrated_filtered_prompt,
        last_num_code_blocks = last_fully_integrated_num_code_blocks,
    )
    TrackedSignal{Sig,B,C,PCF,typeof(new_loop_state.cn0_estimator)}(
        isnothing(signal) ? t.signal : signal,
        isnothing(integrated_samples) ? t.integrated_samples : integrated_samples,
        isnothing(correlator) ? t.correlator : correlator,
        isnothing(last_fully_integrated_correlator) ? t.last_fully_integrated_correlator :
        last_fully_integrated_correlator,
        new_loop_state,
        isnothing(filtered_prompts) ? t.filtered_prompts : filtered_prompts,
        isnothing(correlator_outputs) ? t.correlator_outputs : correlator_outputs,
        isnothing(preferred_num_code_blocks_to_integrate) ?
        t.preferred_num_code_blocks_to_integrate : preferred_num_code_blocks_to_integrate,
    )
end

get_signal(t::TrackedSignal) = t.signal
get_correlator(t::TrackedSignal) = t.correlator
get_last_fully_integrated_correlator(t::TrackedSignal) = t.last_fully_integrated_correlator
get_last_fully_integrated_filtered_prompt(t::TrackedSignal) =
    t.loop_state.last_filtered_prompt
get_last_fully_integrated_num_code_blocks(t::TrackedSignal) =
    t.loop_state.last_num_code_blocks

"""
$(SIGNATURES)

The integration time `T` of the most recently completed record (its primary-code
blocks times one code period), which
[`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{NoiseRefCN0Estimator, Any})
divides by. The record's post-integration SNR is `C/N₀ · T`.

Not to be confused with `get_integrated_samples`, the record in progress.
"""
get_last_fully_integrated_integration_time(t::TrackedSignal) =
    get_last_fully_integrated_num_code_blocks(t) * get_code_length(get_signal(t)) /
    get_code_frequency(get_signal(t))

get_filtered_prompts(t::TrackedSignal) = t.filtered_prompts

"""
$(SIGNATURES)

The [`CorrelatorOutput`](@extref TrackingLoops.CorrelatorOutput)s this signal
completed in the current chunk, in order. The Doppler estimator consumes and
clears them, so read them inside a custom estimator or right after a bare
`downconvert_and_correlate!`. See [Chunked Doppler updates](@ref).
"""
get_correlator_outputs(t::TrackedSignal) = t.correlator_outputs

"""
$(SIGNATURES)

Append an externally built
[`CorrelatorOutput`](@extref TrackingLoops.CorrelatorOutput) to `signal`'s per-chunk
buffer and return `signal`: the ingest path for an **external correlator
producer** (e.g. FPGA dumps). Append per signal in `sample_index` order, then run
[`estimate_dopplers_and_filter_prompt!`](@ref), which folds and clears the buffer.
Type-checked against the signal's correlator, unlike pushing to
[`get_correlator_outputs`](@ref). See [External correlator producers](@ref).
"""
append_correlator_output!(t::TrackedSignal, output::CorrelatorOutput) =
    (push!(t.correlator_outputs, output); t)

get_post_corr_filter(t::TrackedSignal) = t.loop_state.post_corr_filter
get_cn0_estimator(t::TrackedSignal) = t.loop_state.cn0_estimator
get_bit_buffer(t::TrackedSignal) = t.loop_state.bit_buffer
get_soft_bits(t::TrackedSignal) = get_soft_bits(get_bit_buffer(t))
get_num_bits(t::TrackedSignal) = length(get_bit_buffer(t))
has_bit_or_secondary_code_been_found(t::TrackedSignal) =
    has_bit_or_secondary_code_been_found(get_bit_buffer(t))
get_integrated_samples(t::TrackedSignal) = t.integrated_samples
get_preferred_num_code_blocks_to_integrate(t::TrackedSignal) =
    t.preferred_num_code_blocks_to_integrate

"""
$(SIGNATURES)

The state of one tracked satellite: the carrier/code Doppler and phase shared by
its signals, the per-signal state in `signals`, and the Doppler-estimator state
in `doppler_estimator_state`.

`signals[1]` is the **estimator-driver signal**. The conventional estimators
([`ConventionalPLLAndDLL`](@extref TrackingLoops.ConventionalPLLAndDLL),
[`ConventionalAssistedPLLAndDLL`](@extref TrackingLoops.ConventionalAssistedPLLAndDLL))
run their discriminators on it, at its integration boundaries; the other signals
update their own filters, C/N₀ and bit buffers on theirs. A custom
[`AbstractDopplerEstimator`](@extref TrackingLoops.AbstractDopplerEstimator) may
use any signal.

The shared `code_phase` wraps at [`max_code_length`](@ref).
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

# A signal's post-sync wrap in chips: the secondary-code period for a pilot, and
# the longer of secondary-code and data-bit period otherwise (GPS L1 C/A:
# 1023 × 20, no secondary code).
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

# `lcm`, not `max`: per-signal phases are `mod(code_phase, _replica_code_wrap(...))`,
# correct only if the shared wrap is a multiple of every signal's (#129).
@inline _max_code_length(::Tuple{}) = 1
@inline _max_code_length(t::Tuple) =
    lcm(_post_sync_code_length(first(t)), _max_code_length(Base.tail(t)))

"""
$(SIGNATURES)

Upper bound on the shared `sat.code_phase` wrap period, in chips: the least
common multiple of the per-signal wraps once *every* signal has synced (GPS
L1 C/A: 1023 × 20 = 20460, one data bit; L1C-P: 10230 × 1800, one
secondary-code cycle). Folds to a literal for a concrete `signals` type.

The runtime wrap, [`current_code_wrap`](@ref), is shorter for unsynced signals.
"""
@inline max_code_length(signals::Tuple{Vararg{TrackedSignal}}) = _max_code_length(signals)

# Primary code length until synced, then `_post_sync_code_length`.
@inline function _current_code_length(tsig::TrackedSignal)
    if get_bit_buffer(tsig).found
        _post_sync_code_length(tsig)
    else
        get_code_length(tsig.signal)
    end
end

# `lcm` for the same reason as `_max_code_length`.
@inline _current_code_wrap(::Tuple{}) = 1
@inline _current_code_wrap(t::Tuple) =
    lcm(_current_code_length(first(t)), _current_code_wrap(Base.tail(t)))

"""
$(SIGNATURES)

The runtime wrap period of the shared `sat.code_phase`, in chips: like
[`max_code_length`](@ref), but a signal not yet synced (`bit_buffer.found`
false) contributes only its primary code length, since its bit or secondary
chip is unknown. The least common multiple keeps the wrap a multiple of every
signal's replica wrap (#129).
"""
@inline current_code_wrap(signals::Tuple{Vararg{TrackedSignal}}) =
    _current_code_wrap(signals)

# Did any signal's `bit_buffer.found` go from false to true? Gates the one-time
# code-phase snap, valid only at the sync transition (#117).
@inline _any_signal_just_synced(::Tuple{}, ::Tuple{}) = false
@inline _any_signal_just_synced(old::Tuple, new::Tuple) =
    (!get_bit_buffer(first(old)).found && get_bit_buffer(first(new)).found) ||
    _any_signal_just_synced(Base.tail(old), Base.tail(new))

# Snap `code_phase` into the secondary-code window of the synced signal with the
# longest `primary × secondary` length, whose `bit_buffer.secondary_phase` pins
# the position. Unchanged if no synced signal has a secondary code.
@inline function _snap_code_phase_from_synced_signal(
    signals::Tuple{Vararg{TrackedSignal}},
    code_phase::Float64,
)
    best_len, best_phase_chips, best_prim = _find_best_secondary_anchor(signals, 0, 0, 1)
    if best_len == 0
        return code_phase
    end
    # Align the low `best_len` chips to the synced window, keeping the higher
    # wrap and the phase within the primary block: dropping that would inject a
    # phase error at sync (#117).
    base = floor(Int, code_phase / best_len) * best_len
    within_primary = code_phase - floor(code_phase / best_prim) * best_prim
    Float64(base + best_phase_chips) + within_primary
end

# `(length, secondary_phase_in_chips, primary_length)` of that signal, or
# `(0, 0, 1)`.
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
    bb = get_bit_buffer(head)
    if bb.found
        prim = get_code_length(sig)
        sec = get_secondary_code_length(sig)
        total = prim * sec
        # Without a secondary code (e.g. GPS L1 C/A) nothing is pinned.
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

Build a [`TrackedSat`](@ref) from pre-built [`TrackedSignal`](@ref)s and
acquisition handoff values, e.g. for a non-default correlator or post-corr
filter; otherwise use the signal-instance forms. The first signal is the
estimator driver. The estimator state comes from
[`init_estimator_state`](@ref).

`code_doppler = nothing` derives it from `carrier_doppler` and the driver's
code/center frequency ratio. `carrier_phase` is in radians.
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
    # Accept integer Dopplers such as `200Hz`.
    cdop = float(carrier_doppler)
    cd =
        isnothing(code_doppler) ?
        cdop * get_code_center_frequency_ratio(first(tracked_signals).signal) :
        float(code_doppler)
    # Built twice so `init_estimator_state` sees a real `TrackedSat`.
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

Build a [`TrackedSat`](@ref) tracking a tuple of signals, e.g.
`TrackedSat((GPSL1C_P(), GPSL1C_D(), GPSL1CA()), prn, code_phase, carrier_doppler)`,
each with its default correlator; the first is the estimator driver.
`doppler_estimator`, `carrier_phase` and `code_doppler` forward to the
`TrackedSignal`-tuple constructor.

`cn0_estimator` is a **tuple**, one per signal in order (each signal needs its
own instance); a single estimator only for a single signal.
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

# One estimator per signal: `nothing` → defaults; a tuple as is; a single
# estimator only for a single signal.
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

Construct a single-signal [`TrackedSat`](@ref) from acquisition handoff values.
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

# Kwarg-update constructor; `signals` and `doppler_estimator_state` keep their types.
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

Get the satellite's shared code phase in chips, wrapping at
[`max_code_length`](@ref). Mod by a signal's primary code length for its replica
phase.
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

# Per-signal accessors on a TrackedSat, through `_find_signal`: no selector
# (single-signal sats only), an `Integer` index, or a signal type (must match
# exactly one signal).
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

# `output` first, so a trailing signal selector works as for the accessors.
append_correlator_output!(s::TrackedSat, output::CorrelatorOutput, sel...) =
    (append_correlator_output!(_find_signal(s.signals, sel...), output); s)

# Between `track` calls: start sample back to 1, per-signal prompts and outputs
# emptied in place, bit buffers `reset`.
function reset_start_sample_and_bit_buffer(sat::TrackedSat)
    new_signals = map(s -> _reset_signal(s), sat.signals)
    TrackedSat(sat; signal_start_sample = 1, signals = new_signals)
end

@inline function _reset_signal(t::TrackedSignal)
    empty!(t.filtered_prompts)
    empty!(t.correlator_outputs)
    TrackedSignal(t; bit_buffer = reset(get_bit_buffer(t)))
end

"""
    init_estimator_state(estimator, sat::TrackedSat)

Build the per-satellite Doppler-estimator state for `sat`: TrackingLoops'
[`init_estimator_state`](@ref) on the satellite's estimator-driver signal
(`signals[1]`) and its handover Dopplers. A custom estimator defines either
this form or TrackingLoops' four-argument one.

It must be **pure**: it also builds the PRN-0 template sat that fixes a group's
slot type, probes types when validating pre-built sats, and re-seeds in
[`reset_loop_filters!`](@ref). Register satellites in shared state in
[`update_estimator_on_handoff`](@ref) instead.
"""
init_estimator_state(estimator::AbstractDopplerEstimator, sat::TrackedSat) =
    init_estimator_state(
        estimator,
        first(sat.signals).signal,
        sat.carrier_doppler,
        sat.code_doppler,
    )

"""
$(SIGNATURES)

Update an estimator's cross-satellite shared state when satellites enter the
track set. Called once per handoff ([`merge_sats`](@ref),
[`add_satellite`](@ref), [`add_satellite!`](@ref)) with the incoming satellites,
after [`init_estimator_state`](@ref). The default returns `estimator` unchanged.

The result must have the input's concrete type ([`TrackState`](@ref) is
parameterized on it): grow resizable containers in place, or rebuild with
[`Setfield.@set`](https://jw3126.github.io/Setfield.jl/stable/) or a copying
constructor. The returned `TrackState` carries the result, so callers of the
in-place [`add_satellite!`](@ref) too must keep using the returned state.
"""
update_estimator_on_handoff(estimator::AbstractDopplerEstimator, _new_sats) = estimator

"""
A tuple or named tuple of per-group satellite dictionaries, each a
`Dictionary{I, <:TrackedSat}` keyed by PRN.
"""
const SatelliteDicts{N} = TupleLike{<:NTuple{N,Dictionary{<:Any,<:TrackedSat}}}

"""
$(SIGNATURES)

A group of satellites tracking the same tuple of signals, on one RF band and
antenna array: the unit of type stability (one concrete `TrackedSat` type).
Groups may share a band (e.g. `(GPSL1CA(),)` and `(GalileoE1B(),)` on `L1()`);
the band routes the right measurement to the group.

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

# Kwarg-update constructor keeping `g`'s types.
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

# A group's signals must (#129) be on its band, whose measurement the whole group
# correlates against, and share one chip rate, since the shared `code_phase`
# advances at the driver's (`update` in downconvert_and_correlate.jl). Bands
# compare by id, so a user-defined alias of a band validates.
function _validate_signal_group(signals::Tuple{Vararg{AbstractGNSSSignal}}, band)
    driver = first(signals)
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
    nothing
end

"""
$(SIGNATURES)

Build an empty `SignalGroup` from a signal tuple; populate it with
[`add_satellite!`](@ref) once the [`TrackState`](@ref) is built. `band` defaults to
`get_band(first(signals))`, `num_ants` to `NumAnts(1)`. All signals must be on
`band` and share one chip rate (`ArgumentError` otherwise).

`doppler_estimator` only shapes the slot type; the estimator that runs is the
`TrackState`'s, which rebuilds a mismatched empty slot. Leave it at its default
and configure the estimator on `TrackState`.

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
    _validate_signal_group(signals, band)
    template = _make_template_tracked_sat(signals, doppler_estimator, num_ants)
    sats = Dictionary{Int,typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, signals, num_ants)
end

"""
A NamedTuple of `N` `SignalGroup`s, as stored in [`TrackState`](@ref).
"""
const SignalGroups{N} = NamedTuple{<:Any,<:NTuple{N,SignalGroup}}

# Copy the values but share the keys: for the immutable per-iteration steps,
# which never change the key set. `track` detaches the keys once up front
# (`_detach_slot_vector`, #123); called directly, a later `add_satellite!` /
# `remove_satellite!` on the result would corrupt the input's keys.
@inline function _copy_slot_vector(sats::Dictionary{<:Any,<:TrackedSat})
    Dictionary(keys(sats), copy(sats.values))
end

# Keys and values copied, at the `track` boundary (#123).
@inline function _detach_slot_vector(sats::Dictionary{<:Any,<:TrackedSat})
    copy(sats)
end

@inline _copy_group_slot_vectors(g::SignalGroup) =
    SignalGroup(g; satellites = _copy_slot_vector(g.satellites))

@inline _copy_groups_slot_vectors(groups::SignalGroups) =
    map(_copy_group_slot_vectors, groups)

@inline _detach_group_slot_vectors(g::SignalGroup) =
    SignalGroup(g; satellites = _detach_slot_vector(g.satellites))

@inline _detach_groups_slot_vectors(groups::SignalGroups) =
    map(_detach_group_slot_vectors, groups)

# `f(group, args...)` for each group, by tuple recursion: a `for` loop over
# heterogeneous groups would box.
@inline _foreach_group!(f::F, ::Tuple{}, args::Vararg{Any,N}) where {F,N} = nothing
@inline function _foreach_group!(f::F, t::Tuple, args::Vararg{Any,N}) where {F,N}
    f(first(t), args...)
    _foreach_group!(f, Base.tail(t), args...)
end
@inline _foreach_group!(f::F, nt::NamedTuple, args::Vararg{Any,N}) where {F,N} =
    _foreach_group!(f, Tuple(nt), args...)

# Overwrite each (immutable) `TrackedSat` slot with its reset value.
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
    # The buffered prompts' SNR is a whole record's, so divide by its time.
    estimate_cn0(get_cn0_estimator(tsig), get_last_fully_integrated_integration_time(tsig))
end

estimate_cn0(sat::TrackedSat, sel...) = estimate_cn0(_find_signal(sat.signals, sel...))
