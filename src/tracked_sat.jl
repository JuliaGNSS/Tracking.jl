"""
$(SIGNATURES)

Per-signal tracking state, one per signal tracked on a satellite: the
correlator, the integration progress and the post-correlation state, a
TrackingLoops [`PostCorrelationState`](@extref TrackingLoops.PostCorrelationState) (the signal and PRN, post-correlation filter, C/N₀ estimator, bit synchronizer, and
the last predetection integration).
The carrier/code Doppler and phase are shared and live on the [`TrackedSat`](@ref).
"""
struct TrackedSignal{
    S<:AbstractGNSSSignal,
    C<:AbstractCorrelatorBank,
    P<:PostCorrelationState{S},
}
    num_integrated_samples::Int
    correlator_bank::C
    last_correlator_bank::C
    # Advanced per integration by TrackingLoops' `process_correlator_bank_output!`.
    post_correlation_state::P
    filtered_prompts::Vector{ComplexF64}
    # Integrations completed in the current NCO update span: pushed by the correlate phase,
    # stepped through the tracking loop and emptied by `step_tracking_loops!`.
    correlator_bank_outputs::Vector{CorrelatorBankOutput{C}}
    # In primary code periods; see `set_preferred_num_code_periods_to_integrate!`.
    preferred_num_code_periods_to_integrate::Int
end

# For data-bearing signals the length must divide the code periods per bit, or
# integrations straddle bit boundaries and no bit is ever emitted again (#128).
function validate_preferred_num_code_periods_to_integrate(
    signal::AbstractGNSSSignal,
    preferred_num_code_periods::Integer,
)
    preferred_num_code_periods >= 1 || throw(
        ArgumentError(
            "preferred_num_code_periods_to_integrate must be at least 1, got " *
            "$preferred_num_code_periods",
        ),
    )
    num_code_periods_per_bit = get_code_periods_per_bit(signal)
    num_code_periods_per_bit == 0 && return nothing
    if num_code_periods_per_bit % preferred_num_code_periods != 0
        valid = filter(d -> num_code_periods_per_bit % d == 0, 1:num_code_periods_per_bit)
        throw(
            ArgumentError(
                "preferred_num_code_periods_to_integrate = $preferred_num_code_periods " *
                "must evenly divide the $num_code_periods_per_bit code periods " *
                "of one bit of $(get_signal_name(signal)) — an integration " *
                "straddling a bit boundary would never emit a bit. Valid values: " *
                "$(join(valid, ", ")).",
            ),
        )
    end
    nothing
end

"""
$(SIGNATURES)

Construct a fresh [`TrackedSignal`](@ref) for `signal`. The correlator bank,
post-correlation filter and C/N₀ estimator default to the signal's recommended values.
The [`TrackedSat`](@ref) it joins sets its PRN and its driver's carrier-phase offset.

`cn0_estimator` accepts any
[`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator). The default,
[`default_cn0_estimator`](@extref TrackingLoops.default_cn0_estimator-Tuple{AbstractGNSSSignal, Int64}), is a
[`NoiseReferencedCN0Estimator`](@extref TrackingLoops.NoiseReferencedCN0Estimator) over
`num_averaged_integrations` integrations against the signal's **measured**
noise density: `track!` fills that on the sample-driven path, but a
**correlator-ingest** path must feed [`append_noise_observation!`](@ref) per signal, or
C/N₀ stays at `-Inf dB-Hz`. `NWPRCN0Estimator()` needs no noise observation;
`MomentsCN0Estimator(; num_averaged_integrations)` is the
plain moment-ratio estimator. Each signal needs its **own** estimator instance: they
buffer into a vector, so sharing one corrupts both.

`preferred_num_code_periods_to_integrate` defaults to
[`default_num_code_periods_to_integrate`](@extref TrackingLoops.default_num_code_periods_to_integrate-Tuple{AbstractGNSSSignal})
and must be valid for `signal` (see
[`set_preferred_num_code_periods_to_integrate!`](@ref)).
"""
function TrackedSignal(
    signal::AbstractGNSSSignal;
    num_ants::NumAnts = NumAnts(1),
    correlator_bank::AbstractCorrelatorBank = get_default_correlator_bank(signal, num_ants),
    num_averaged_integrations::Int = 100,
    cn0_estimator::AbstractCN0Estimator = default_cn0_estimator(
        signal,
        num_averaged_integrations,
    ),
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    preferred_num_code_periods_to_integrate::Int = default_num_code_periods_to_integrate(
        signal,
    ),
)
    validate_preferred_num_code_periods_to_integrate(
        signal,
        preferred_num_code_periods_to_integrate,
    )
    # A default NCO update span yields at most one integration; the hint covers a moderately
    # enlarged `nco_update_interval` without regrowing.
    correlator_bank_outputs = CorrelatorBankOutput{typeof(correlator_bank)}[]
    sizehint!(correlator_bank_outputs, 4)
    TrackedSignal(
        0,
        correlator_bank,
        correlator_bank,
        PostCorrelationState(signal, 0; correlator_bank, cn0_estimator, post_corr_filter),
        ComplexF64[],
        correlator_bank_outputs,
        preferred_num_code_periods_to_integrate,
    )
end

# Kwarg-update constructor keeping `t`'s correlator bank type. `cn0_estimator` may
# change type, so a custom estimator can be swapped in; call sites pass `nothing` or a
# concretely typed one, so the result infers and the per-integration rebuild stays
# allocation-free.
function TrackedSignal(
    t::TrackedSignal{S,C};
    num_integrated_samples = nothing,
    correlator_bank::Maybe{C} = nothing,
    last_correlator_bank::Maybe{C} = nothing,
    post_correlation_state::Maybe{PostCorrelationState{S}} = nothing,
    cn0_estimator::Maybe{AbstractCN0Estimator} = nothing,
    preferred_num_code_periods_to_integrate = nothing,
) where {S<:AbstractGNSSSignal,C<:AbstractCorrelatorBank}
    isnothing(preferred_num_code_periods_to_integrate) ||
        validate_preferred_num_code_periods_to_integrate(
            get_signal(t),
            preferred_num_code_periods_to_integrate,
        )
    new_post_correlation_state = PostCorrelationState(
        something(post_correlation_state, t.post_correlation_state);
        cn0_estimator,
    )
    TrackedSignal{S,C,typeof(new_post_correlation_state)}(
        something(num_integrated_samples, t.num_integrated_samples),
        something(correlator_bank, t.correlator_bank),
        something(last_correlator_bank, t.last_correlator_bank),
        new_post_correlation_state,
        t.filtered_prompts,
        t.correlator_bank_outputs,
        something(
            preferred_num_code_periods_to_integrate,
            t.preferred_num_code_periods_to_integrate,
        ),
    )
end

get_signal(t::TrackedSignal) = t.post_correlation_state.signal

"""
$(SIGNATURES)

The working correlator bank: the correlator outputs of the predetection integration in progress.
"""
get_correlator_bank(t::TrackedSignal) = t.correlator_bank

"""
$(SIGNATURES)

The correlator bank of the most recently completed predetection integration.
"""
get_last_correlator_bank(t::TrackedSignal) = t.last_correlator_bank

"""
$(SIGNATURES)

The filtered prompt of the most recently completed predetection
integration.
"""
get_last_filtered_prompt(t::TrackedSignal) =
    get_prompt(t.post_correlation_state.last_filtered_correlator_bank)

"""
$(SIGNATURES)

The number of primary code periods the most recently completed predetection
integration spanned.
"""
get_last_num_code_periods(t::TrackedSignal) =
    get_last_num_code_periods(t.post_correlation_state)

"""
$(SIGNATURES)

The coherent integration time `T` of the latest integration, in seconds; the
post-integration SNR is `C/N₀ · T`. `NaN` s before the first integration.

Not to be confused with `get_num_integrated_samples`, the integration in progress.
"""
get_last_coherent_integration_time(t::TrackedSignal) =
    get_last_coherent_integration_time(t.post_correlation_state)

"""
$(SIGNATURES)

The filtered prompts of the predetection integrations the most recent `track` call
completed, in order.
"""
get_filtered_prompts(t::TrackedSignal) = t.filtered_prompts

"""
$(SIGNATURES)

The [`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput)s this signal
completed in the current NCO update span, in order. [`step_tracking_loops!`](@ref) consumes
and clears them, so read them inside a custom tracking loop or right after a bare
`downconvert_and_correlate!`. See [NCO updates](@ref).
"""
get_correlator_bank_outputs(t::TrackedSignal) = t.correlator_bank_outputs

"""
$(SIGNATURES)

Append an externally built
[`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput) to `signal`'s per-span
buffer and return `signal`: the ingest path for an **external correlator
producer** (e.g. FPGA dumps). Append per signal in `end_sample_index` order, then run
[`step_tracking_loops!`](@ref), which steps the tracking loop with them and clears the buffer.
Type-checked against the signal's correlator bank, unlike pushing to
[`get_correlator_bank_outputs`](@ref). See [External correlator producers](@ref).
"""
append_correlator_bank_output!(t::TrackedSignal, output::CorrelatorBankOutput) =
    (push!(t.correlator_bank_outputs, output); t)

get_post_corr_filter(t::TrackedSignal) = t.post_correlation_state.post_corr_filter
get_cn0_estimator(t::TrackedSignal) = t.post_correlation_state.cn0_estimator

"""
$(SIGNATURES)

The signal's [`BitSynchronizer`](@extref TrackingLoops.BitSynchronizer): bit or secondary-code
sync and the soft bits for the decoder.
"""
get_bit_synchronizer(t::TrackedSignal) = t.post_correlation_state.bit_synchronizer
get_soft_bits(t::TrackedSignal) = get_soft_bits(get_bit_synchronizer(t))

"""
$(SIGNATURES)

The number of soft bits in the bit synchronizer (see `get_soft_bits`).
"""
get_num_bits(t::TrackedSignal) = length(get_bit_synchronizer(t))
is_bit_or_secondary_code_synced(t::TrackedSignal) =
    is_bit_or_secondary_code_synced(get_bit_synchronizer(t))

"""
$(SIGNATURES)

The number of samples accumulated into the predetection integration in progress.
"""
get_num_integrated_samples(t::TrackedSignal) = t.num_integrated_samples

"""
$(SIGNATURES)

The number of primary code periods the signal integrates coherently per predetection
integration once synced (see
[`set_preferred_num_code_periods_to_integrate!`](@ref)).
"""
get_preferred_num_code_periods_to_integrate(t::TrackedSignal) =
    t.preferred_num_code_periods_to_integrate

"""
$(SIGNATURES)

The state of one tracked satellite: the carrier/code Doppler and phase shared by
its signals, the per-signal state in `signals`, and the tracking loop's
per-satellite state in `loop_state`.

`signals[1]` is the **driver signal**. The scalar tracking loops
([`ScalarTrackingLoop`](@extref TrackingLoops.ScalarTrackingLoop),
[`DelayCompensatedTrackingLoop`](@extref TrackingLoops.DelayCompensatedTrackingLoop))
run their discriminators on it, at its integration boundaries; the other signals
update their own filters, C/N₀ and bit synchronizers on theirs. A custom
[`AbstractTrackingLoop`](@extref TrackingLoops.AbstractTrackingLoop) may
use any signal.

The shared `code_phase` wraps at [`current_code_phase_period`](@ref Tracking.current_code_phase_period).
"""
struct TrackedSat{Signals<:Tuple{Vararg{TrackedSignal}},D}
    prn::Int
    code_phase::Float64
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    carrier_doppler::typeof(1.0Hz)
    next_sample_index::Int
    signals::Signals
    loop_state::D
end

# A signal's post-sync period in chips: the secondary-code period for a pilot, and
# the longer of secondary-code and data-bit period otherwise (GPS L1 C/A:
# 1023 × 20, no secondary code).
@inline function _post_sync_code_phase_period(tsig::TrackedSignal)
    sig = get_signal(tsig)
    get_code_length(sig) *
    max(get_secondary_code_length(sig), get_code_periods_per_bit(sig))
end

# Primary code length until synced, then `_post_sync_code_phase_period`.
@inline function _signal_code_phase_period(tsig::TrackedSignal)
    if is_bit_or_secondary_code_synced(tsig)
        _post_sync_code_phase_period(tsig)
    else
        get_code_length(get_signal(tsig))
    end
end

"""
$(SIGNATURES)

The period of the shared `sat.code_phase`, in chips: the least common multiple of the
signals' periods. A signal not yet synced (`is_bit_or_secondary_code_synced` false)
contributes its primary code length, since its bit or secondary chip is unknown; a
synced one its secondary-code period for a pilot, else the longer of secondary-code
and data-bit period (GPS L1 C/A: 1023 × 20 = 20460, one data bit). The least common
multiple, not the maximum, keeps the period a multiple of every signal's replica period
(#129).
"""
@inline current_code_phase_period(::Tuple{}) = 1
@inline current_code_phase_period(signals::Tuple{Vararg{TrackedSignal}}) = lcm(
    _signal_code_phase_period(first(signals)),
    current_code_phase_period(Base.tail(signals)),
)

# Did any signal go from unsynced to synced? Gates the one-time
# code-phase snap, valid only at the sync transition (#117).
@inline _any_signal_just_synced(::Tuple{}, ::Tuple{}) = false
@inline _any_signal_just_synced(old::Tuple, new::Tuple) =
    (
        !is_bit_or_secondary_code_synced(first(old)) &&
        is_bit_or_secondary_code_synced(first(new))
    ) || _any_signal_just_synced(Base.tail(old), Base.tail(new))

# Snap `code_phase` into the secondary-code window of the synced signal with the
# longest `primary × secondary` length, whose `get_secondary_code_phase` pins the
# position. Unchanged if no synced signal has a secondary code.
@inline function _snap_code_phase_from_synced_signal(
    signals::Tuple{Vararg{TrackedSignal}},
    code_phase::Float64,
)
    best_len, best_phase_chips, best_prim = _find_best_secondary_anchor(signals, 0, 0, 1)
    if best_len == 0
        return code_phase
    end
    # Align the low `best_len` chips to the synced window, keeping the higher
    # wrap and the phase within the primary code period: dropping that would inject a
    # phase error at sync (#117).
    base = floor(Int, code_phase / best_len) * best_len
    within_primary = code_phase - floor(code_phase / best_prim) * best_prim
    Float64(base + best_phase_chips) + within_primary
end

# `(length, secondary_code_phase_in_chips, primary_length)` of that signal, or
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
    sig = get_signal(head)
    post_correlation_state = head.post_correlation_state
    if is_bit_or_secondary_code_synced(post_correlation_state)
        prim = get_code_length(sig)
        sec = get_secondary_code_length(sig)
        total = prim * sec
        # Without a secondary code (e.g. GPS L1 C/A) nothing is pinned.
        if sec > 1 && total > best_len
            best_len = total
            best_chips = get_secondary_code_phase(post_correlation_state) * prim
            best_prim = prim
        end
    end
    _find_best_secondary_anchor(Base.tail(t), best_len, best_chips, best_prim)
end

"""
$(SIGNATURES)

Build a [`TrackedSat`](@ref) from pre-built [`TrackedSignal`](@ref)s and
acquisition handoff values, e.g. for a non-default correlator or post-correlation filter; otherwise use the signal-instance forms. The first signal is the
driver signal. The loop state comes from [`init_loop_state`](@ref).

`code_doppler = nothing` derives it from `carrier_doppler` and the driver's
code/center frequency ratio. `carrier_phase` is in cycles.
"""
function TrackedSat(
    tracked_signals::Tuple{TrackedSignal,Vararg{TrackedSignal}},
    prn::Int,
    code_phase,
    carrier_doppler;
    tracking_loop::AbstractTrackingLoop = _default_tracking_loop(),
    carrier_phase = 0.0,
    code_doppler = nothing,
)
    _assert_unique_signal_ids(map(get_signal, tracked_signals))
    # Accept integer Dopplers such as `200Hz`.
    carrier_doppler = float(carrier_doppler)
    code_doppler =
        isnothing(code_doppler) ?
        carrier_doppler *
        get_chip_rate_to_center_frequency_ratio(get_signal(first(tracked_signals))) :
        float(code_doppler)
    driver_carrier_phase_offset =
        get_carrier_phase_offset(get_signal(first(tracked_signals)))
    # Built twice so `init_loop_state` sees a real `TrackedSat`.
    bare = TrackedSat(
        prn,
        float(code_phase),
        code_doppler,
        float(carrier_phase),
        carrier_doppler,
        1,
        map(s -> _join_satellite(s, prn, driver_carrier_phase_offset), tracked_signals),
        nothing,
    )
    loop_state = init_loop_state(tracking_loop, bare)
    TrackedSat(
        bare.prn,
        bare.code_phase,
        bare.code_doppler,
        bare.carrier_phase,
        bare.carrier_doppler,
        bare.next_sample_index,
        bare.signals,
        loop_state,
    )
end

# A signal's post-correlation state on its satellite: its PRN, and the carrier-phase
# offset of the driver, against which its bit synchronizer de-rotates the prompt so a
# quadrature component (GPS L5, Galileo E5a) decodes off the right axis.
_join_satellite(tsig::TrackedSignal, prn::Int, driver_carrier_phase_offset::Real) =
    TrackedSignal(
        tsig;
        post_correlation_state = PostCorrelationState(
            tsig.post_correlation_state;
            prn,
            driver_carrier_phase_offset,
        ),
    )

"""
$(SIGNATURES)

Build a [`TrackedSat`](@ref) tracking a tuple of signals, e.g.
`TrackedSat((GPSL1C_P(), GPSL1C_D(), GPSL1CA()), prn, code_phase, carrier_doppler)`,
each with its default correlator bank; the first is the driver signal.
`tracking_loop`, `carrier_phase` and `code_doppler` forward to the
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
    num_averaged_integrations::Int = 100,
    cn0_estimator = nothing,
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    kwargs...,
)
    cn0_estimators =
        _per_signal_cn0_estimators(signals, cn0_estimator, num_averaged_integrations)
    tracked_signals = map(signals, cn0_estimators) do sig, estimator
        TrackedSignal(
            sig;
            num_ants,
            correlator_bank = get_default_correlator_bank(sig, num_ants),
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
    num_averaged_integrations::Int,
) = map(sig -> default_cn0_estimator(sig, num_averaged_integrations), signals)

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
        "`cn0_estimator = ($(join(fill("MomentsCN0Estimator(; num_averaged_integrations = 100)", length(signals)), ", ")))`.",
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
    correlator_bank::AbstractCorrelatorBank = get_default_correlator_bank(signal, num_ants),
    num_averaged_integrations::Int = 100,
    cn0_estimator::AbstractCN0Estimator = default_cn0_estimator(
        signal,
        num_averaged_integrations,
    ),
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
    kwargs...,
)
    tracked_signal =
        TrackedSignal(signal; num_ants, correlator_bank, cn0_estimator, post_corr_filter)
    TrackedSat((tracked_signal,), prn, code_phase, carrier_doppler; kwargs...)
end

# Kwarg-update constructor; `signals` and `loop_state` keep their types.
function TrackedSat(
    sat::TrackedSat{Signals,D};
    code_phase = nothing,
    code_doppler = nothing,
    carrier_phase = nothing,
    carrier_doppler = nothing,
    next_sample_index = nothing,
    signals::Maybe{Signals} = nothing,
    loop_state::Maybe{D} = nothing,
) where {Signals<:Tuple{Vararg{TrackedSignal}},D}
    TrackedSat{Signals,D}(
        sat.prn,
        isnothing(code_phase) ? sat.code_phase : code_phase,
        isnothing(code_doppler) ? sat.code_doppler : code_doppler,
        isnothing(carrier_phase) ? sat.carrier_phase : carrier_phase,
        isnothing(carrier_doppler) ? sat.carrier_doppler : carrier_doppler,
        isnothing(next_sample_index) ? sat.next_sample_index : next_sample_index,
        isnothing(signals) ? sat.signals : signals,
        isnothing(loop_state) ? sat.loop_state : loop_state,
    )
end

"""
$(SIGNATURES)

Get the PRN (Pseudo-Random Noise) number of the satellite.
"""
get_prn(sat::TrackedSat) = sat.prn
get_num_ants(
    sat::TrackedSat{<:Tuple{TrackedSignal{<:Any,<:AbstractCorrelatorBank{M}},Vararg}},
) where {M} = M

"""
$(SIGNATURES)

Get the satellite's shared code phase in chips, wrapping at
[`current_code_phase_period`](@ref Tracking.current_code_phase_period). Mod by a signal's primary code length
for its replica phase.
"""
get_code_phase(sat::TrackedSat) = sat.code_phase

"""
$(SIGNATURES)

Get the current code Doppler frequency.
"""
get_code_doppler(sat::TrackedSat) = sat.code_doppler

"""
$(SIGNATURES)

Get the carrier phase of the replica in cycles, folded to [-0.5, 0.5): not the
accumulated carrier-phase measurement, for which a host adds
[`calc_carrier_phase_correction`](@ref TrackingLoops.calc_carrier_phase_correction(::TrackState, ::Union{Integer, Symbol, Val}, ::Any))
and counts the cycles itself.
"""
get_carrier_phase(sat::TrackedSat) = sat.carrier_phase

"""
$(SIGNATURES)

Get the current carrier Doppler frequency.
"""
get_carrier_doppler(sat::TrackedSat) = sat.carrier_doppler

"""
$(SIGNATURES)

Get the index of the next sample of the current samples this satellite
processes, mid-integration if one is in flight; one past the end once the
buffer is done. `track` resets it to 1 for each new sample buffer.
"""
get_next_sample_index(sat::TrackedSat) = sat.next_sample_index

"""
$(SIGNATURES)

Get the satellite's tuple of [`TrackedSignal`](@ref)s.
"""
get_signals(sat::TrackedSat) = sat.signals

"""
$(SIGNATURES)

Get the tracking loop's per-satellite state (e.g. the loop-filter state of a
scalar tracking loop, a `TrackingLoops.ScalarLoopState`).
"""
get_loop_state(sat::TrackedSat) = sat.loop_state

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
    if get_signal(head) isa T
        _assert_no_more_of_type(Base.tail(t), T)
        return head
    end
    _find_signal_by_type(Base.tail(t), T)
end
@inline _assert_no_more_of_type(::Tuple{}, ::Type{T}) where {T<:AbstractGNSSSignal} =
    nothing
@inline function _assert_no_more_of_type(t::Tuple, ::Type{T}) where {T<:AbstractGNSSSignal}
    get_signal(first(t)) isa T && throw(
        ArgumentError(
            "signal type $T matches more than one signal on this satellite — " *
            "use an integer index to disambiguate (e.g. `get_signals(sat)[i]`).",
        ),
    )
    _assert_no_more_of_type(Base.tail(t), T)
end

get_signal(sat::TrackedSat, sel...) = get_signal(_find_signal(sat.signals, sel...))
get_correlator_bank(sat::TrackedSat, sel...) =
    get_correlator_bank(_find_signal(sat.signals, sel...))
get_last_correlator_bank(sat::TrackedSat, sel...) =
    get_last_correlator_bank(_find_signal(sat.signals, sel...))
get_last_filtered_prompt(sat::TrackedSat, sel...) =
    get_last_filtered_prompt(_find_signal(sat.signals, sel...))
get_last_num_code_periods(sat::TrackedSat, sel...) =
    get_last_num_code_periods(_find_signal(sat.signals, sel...))
get_last_coherent_integration_time(sat::TrackedSat, sel...) =
    get_last_coherent_integration_time(_find_signal(sat.signals, sel...))
get_filtered_prompts(sat::TrackedSat, sel...) =
    get_filtered_prompts(_find_signal(sat.signals, sel...))
get_post_corr_filter(sat::TrackedSat, sel...) =
    get_post_corr_filter(_find_signal(sat.signals, sel...))
get_cn0_estimator(sat::TrackedSat, sel...) =
    get_cn0_estimator(_find_signal(sat.signals, sel...))
get_bit_synchronizer(sat::TrackedSat, sel...) =
    get_bit_synchronizer(_find_signal(sat.signals, sel...))
get_soft_bits(sat::TrackedSat, sel...) = get_soft_bits(_find_signal(sat.signals, sel...))
get_num_bits(sat::TrackedSat, sel...) = get_num_bits(_find_signal(sat.signals, sel...))
is_bit_or_secondary_code_synced(sat::TrackedSat, sel...) =
    is_bit_or_secondary_code_synced(_find_signal(sat.signals, sel...))
get_num_integrated_samples(sat::TrackedSat, sel...) =
    get_num_integrated_samples(_find_signal(sat.signals, sel...))
get_correlator_bank_outputs(sat::TrackedSat, sel...) =
    get_correlator_bank_outputs(_find_signal(sat.signals, sel...))
get_preferred_num_code_periods_to_integrate(sat::TrackedSat, sel...) =
    get_preferred_num_code_periods_to_integrate(_find_signal(sat.signals, sel...))

# `output` first, so a trailing signal selector works as for the accessors.
append_correlator_bank_output!(sat::TrackedSat, output::CorrelatorBankOutput, sel...) =
    (append_correlator_bank_output!(_find_signal(sat.signals, sel...), output); sat)

# Between `track` calls: next sample index back to 1, per-signal prompts and
# outputs emptied in place, soft bits drained.
function _prepare_for_new_buffer(sat::TrackedSat)
    new_signals = map(s -> _reset_signal(s), sat.signals)
    TrackedSat(sat; next_sample_index = 1, signals = new_signals)
end

@inline function _reset_signal(t::TrackedSignal)
    empty!(t.filtered_prompts)
    empty!(t.correlator_bank_outputs)
    drain_soft_bits!(get_bit_synchronizer(t))
    t
end

"""
    init_loop_state(tracking_loop, sat::TrackedSat)

Build the tracking loop's per-satellite state for `sat`: TrackingLoops'
[`init_loop_state`](@ref) on the satellite's driver signal (`signals[1]`) and its
handover Dopplers. A custom tracking loop defines either this form or
TrackingLoops' four-argument one.

It must be **pure**: it also builds the PRN-0 template sat that fixes a group's
slot type, probes types when validating pre-built sats, and re-seeds in
[`reset_loop_filters!`](@ref). Register satellites in shared state in
[`update_tracking_loop_on_new_sats`](@ref Tracking.update_tracking_loop_on_new_sats) instead.
"""
init_loop_state(tracking_loop::AbstractTrackingLoop, sat::TrackedSat) = init_loop_state(
    tracking_loop,
    get_signal(first(sat.signals)),
    sat.carrier_doppler,
    sat.code_doppler,
)

"""
$(SIGNATURES)

Update a tracking loop's cross-satellite shared state when satellites enter the
track set. Called once per batch ([`merge_sats`](@ref),
[`add_satellite`](@ref), [`add_satellite!`](@ref)) with the incoming satellites,
after [`init_loop_state`](@ref). The default returns `tracking_loop` unchanged.

The result must have the input's concrete type ([`TrackState`](@ref) is
parameterized on it): grow resizable containers in place, or rebuild with
[`Setfield.@set`](https://jw3126.github.io/Setfield.jl/stable/) or a copying
constructor. The returned `TrackState` carries the result, so callers of the
in-place [`add_satellite!`](@ref) too must keep using the returned state.
"""
update_tracking_loop_on_new_sats(tracking_loop::AbstractTrackingLoop, _new_sats) =
    tracking_loop

"""
A tuple or named tuple of per-group [`TrackedSat`](@ref) dictionaries, each a
`Dictionary{I, <:TrackedSat}` keyed by PRN: what `TrackState(tracked_sats)` builds
groups from and [`get_tracked_sats`](@ref) reads one group of.
"""
const TrackedSatDicts{N} = TupleLike{<:NTuple{N,Dictionary{<:Any,<:TrackedSat}}}

"""
$(SIGNATURES)

A group of satellites tracking the same tuple of signals, on one RF band and
antenna array: the unit of type stability (one concrete `TrackedSat` type).
Groups may share a band (e.g. `(GPSL1CA(),)` and `(GalileoE1B(),)` on `L1()`);
the band routes the right samples to the group.

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

# TrackingLoops tells a satellite's driver integrations from its other signals' by
# signal id, so a satellite tracks each signal once: a second GPS L1 C/A would close
# the loop a second time per integration.
function _assert_unique_signal_ids(signals::Tuple{Vararg{AbstractGNSSSignal}})
    ids = map(get_signal_id, signals)
    allunique(ids) && return nothing
    duplicate = first(id for id in ids if count(==(id), ids) > 1)
    throw(
        ArgumentError(
            string(
                "Each signal can be tracked only once per satellite, but `",
                duplicate,
                "` is given more than once: TrackingLoops tells the driver from ",
                "the other signals by signal id.",
            ),
        ),
    )
end

# A group's signals must (#129) be on its band, whose samples the whole group
# correlates against, and share one chip rate, since the shared `code_phase`
# advances at the driver's (`_advance_sat` in downconvert_and_correlate.jl). Bands
# compare by id, so a user-defined alias of a band validates.
function _validate_signal_group(signals::Tuple{Vararg{AbstractGNSSSignal}}, band)
    _assert_unique_signal_ids(signals)
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
        if get_chip_rate(s) != get_chip_rate(driver)
            throw(
                ArgumentError(
                    string(
                        "All signals in a SignalGroup must share one chip rate: the ",
                        "shared code phase advances at the first signal's code ",
                        "frequency (`",
                        get_signal_id(driver),
                        "`: ",
                        get_chip_rate(driver),
                        "), but `",
                        get_signal_id(s),
                        "` has ",
                        get_chip_rate(s),
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
`band`, share one chip rate and differ in signal id (`ArgumentError` otherwise).

`tracking_loop` only shapes the slot type; the tracking loop that runs is the
`TrackState`'s, which rebuilds a mismatched empty slot. Leave it at its default
and configure the tracking loop on `TrackState`.

```julia
SignalGroup((GPSL1CA(),))                              # band L1(), 1 antenna
SignalGroup((GPSL5I(),); num_ants = NumAnts(2))  # 2-antenna L5
```
"""
function SignalGroup(
    signals::Tuple{Vararg{AbstractGNSSSignal}};
    band = get_band(first(signals)),
    num_ants::NumAnts = NumAnts(1),
    tracking_loop::AbstractTrackingLoop = _default_tracking_loop(),
)
    _validate_signal_group(signals, band)
    template = _make_template_tracked_sat(signals, tracking_loop, num_ants)
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
        vals[i] = _prepare_for_new_buffer(vals[i])
    end
    return nothing
end

function _prepare_for_new_buffer!(groups::SignalGroups)
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

Get the [`TrackedSat`](@ref) of `prn` from a group's dictionary, or its only one.
"""
get_tracked_sat(sats::Dictionary{<:Any,<:TrackedSat}, prn) = sats[prn]
get_tracked_sat(sats::Dictionary{<:Any,<:TrackedSat}) = only(sats)

estimate_cn0(tsig::TrackedSignal) = estimate_cn0(tsig.post_correlation_state)

estimate_cn0(sat::TrackedSat, sel...) = estimate_cn0(_find_signal(sat.signals, sel...))
