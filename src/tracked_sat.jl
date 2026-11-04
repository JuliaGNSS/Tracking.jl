"""
$(SIGNATURES)

Per-signal tracking state, one per signal tracked on a satellite: the
correlator, the integration progress and the post-correlation state, a
TrackingLoops [`PostCorrelationState`](@extref TrackingLoops.PostCorrelationState) (beamformer, C/N₀ estimator, bit synchronizer, and what the last predetection integration
left for the next).
The carrier/code Doppler and phase are shared and live on the [`TrackedSat`](@ref).
"""
struct TrackedSignal{
    Sig<:AbstractGNSSSignal,
    B<:Unsigned,
    C<:AbstractCorrelatorBank,
    BF<:AbstractBeamformer,
    CN0<:AbstractCN0Estimator,
}
    signal::Sig
    num_integrated_samples::Int
    correlator_bank::C
    last_integration_correlator_bank::C
    # Advanced per integration by TrackingLoops' `process_correlator_bank_output!`. Its
    # `last_num_code_periods` is what the integration actually spanned, which the
    # C/N₀ estimate divides by.
    post_correlation_state::PostCorrelationState{B,BF,CN0}
    beamformed_prompts::Vector{ComplexF64}
    # Integrations completed in the current chunk: pushed by the correlate phase,
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

Construct a fresh [`TrackedSignal`](@ref) for `signal`. The correlator,
beamformer and C/N₀ estimator default to the signal's recommended values.

`cn0_estimator` accepts any
[`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator). The default,
[`default_cn0_estimator`](@extref TrackingLoops.default_cn0_estimator-Tuple{AbstractGNSSSignal, Int64}), is a
[`NoiseReferencedCN0Estimator`](@extref TrackingLoops.NoiseReferencedCN0Estimator) over
`num_integrations_for_cn0_estimation` integrations against the signal's **measured** noise
density: `track!` fills that on the sample-driven path, but a
**correlator-ingest** path must feed [`append_noise_observation!`](@ref) per
signal, or C/N₀ stays at `-Inf dB-Hz`. `NWPRCN0Estimator()` needs no noise
observation; `MomentsCN0Estimator(; num_integrations = num_integrations_for_cn0_estimation)`
is the plain moment-ratio estimator. Each signal needs its **own** estimator instance: they
buffer into a vector, so sharing one corrupts both.

`preferred_num_code_periods_to_integrate` defaults to
[`default_num_code_periods_to_integrate`](@extref TrackingLoops.default_num_code_periods_to_integrate-Tuple{AbstractGNSSSignal})
and must be valid for `signal` (see
[`set_preferred_num_code_periods_to_integrate!`](@ref)).
"""
function TrackedSignal(
    signal::AbstractGNSSSignal;
    num_antennas::NumAntennas = NumAntennas(1),
    correlator_bank::AbstractCorrelatorBank = get_default_correlator_bank(
        signal,
        num_antennas,
    ),
    num_integrations_for_cn0_estimation::Int = 100,
    cn0_estimator::AbstractCN0Estimator = default_cn0_estimator(
        signal,
        num_integrations_for_cn0_estimation,
    ),
    beamformer::AbstractBeamformer = DefaultBeamformer(),
    preferred_num_code_periods_to_integrate::Int = default_num_code_periods_to_integrate(
        signal,
    ),
)
    validate_preferred_num_code_periods_to_integrate(
        signal,
        preferred_num_code_periods_to_integrate,
    )
    # A default chunk yields at most one integration; the hint covers a moderately
    # enlarged `nco_update_interval` without regrowing.
    correlator_bank_outputs = CorrelatorBankOutput{typeof(correlator_bank)}[]
    sizehint!(correlator_bank_outputs, 4)
    TrackedSignal(
        signal,
        0,
        correlator_bank,
        correlator_bank,
        PostCorrelationState(signal; cn0_estimator, beamformer),
        ComplexF64[],
        correlator_bank_outputs,
        preferred_num_code_periods_to_integrate,
    )
end

# Kwarg-update constructor keeping `t`'s correlator bank and beamformer types. The
# post-correlation state's parts can be replaced one by one; `last_integration_*` map
# to its `last_beamformed_prompt` / `last_num_code_periods`. `cn0_estimator` may change type,
# so a custom estimator can be swapped in; call sites pass `nothing` or a concrete
# estimator, so the result infers and the per-integration rebuild stays allocation-free.
function TrackedSignal(
    t::TrackedSignal{Sig,B,C,BF};
    signal = nothing,
    num_integrated_samples = nothing,
    correlator_bank::Maybe{C} = nothing,
    last_integration_correlator_bank::Maybe{C} = nothing,
    post_correlation_state::Maybe{PostCorrelationState{B,BF}} = nothing,
    last_integration_beamformed_prompt = nothing,
    cn0_estimator::Maybe{AbstractCN0Estimator} = nothing,
    bit_synchronizer::Maybe{BitSynchronizer{B}} = nothing,
    beamformer::Maybe{BF} = nothing,
    beamformed_prompts::Maybe{Vector{ComplexF64}} = nothing,
    correlator_bank_outputs::Maybe{Vector{CorrelatorBankOutput{C}}} = nothing,
    preferred_num_code_periods_to_integrate = nothing,
    last_integration_num_code_periods = nothing,
) where {
    Sig<:AbstractGNSSSignal,
    B<:Unsigned,
    C<:AbstractCorrelatorBank,
    BF<:AbstractBeamformer,
}
    isnothing(preferred_num_code_periods_to_integrate) ||
        validate_preferred_num_code_periods_to_integrate(
            isnothing(signal) ? t.signal : signal,
            preferred_num_code_periods_to_integrate,
        )
    new_post_correlation_state = PostCorrelationState(
        isnothing(post_correlation_state) ? t.post_correlation_state :
        post_correlation_state;
        bit_synchronizer,
        cn0_estimator,
        beamformer,
        last_beamformed_prompt = last_integration_beamformed_prompt,
        last_num_code_periods = last_integration_num_code_periods,
    )
    TrackedSignal{Sig,B,C,BF,typeof(new_post_correlation_state.cn0_estimator)}(
        isnothing(signal) ? t.signal : signal,
        isnothing(num_integrated_samples) ? t.num_integrated_samples :
        num_integrated_samples,
        isnothing(correlator_bank) ? t.correlator_bank : correlator_bank,
        isnothing(last_integration_correlator_bank) ? t.last_integration_correlator_bank :
        last_integration_correlator_bank,
        new_post_correlation_state,
        isnothing(beamformed_prompts) ? t.beamformed_prompts : beamformed_prompts,
        isnothing(correlator_bank_outputs) ? t.correlator_bank_outputs :
        correlator_bank_outputs,
        isnothing(preferred_num_code_periods_to_integrate) ?
        t.preferred_num_code_periods_to_integrate : preferred_num_code_periods_to_integrate,
    )
end

get_signal(t::TrackedSignal) = t.signal

"""
$(SIGNATURES)

The working correlator bank: the correlator outputs of the predetection integration in progress.
"""
get_correlator_bank(t::TrackedSignal) = t.correlator_bank

"""
$(SIGNATURES)

The correlator bank of the most recently completed predetection integration.
"""
get_last_integration_correlator_bank(t::TrackedSignal) = t.last_integration_correlator_bank

"""
$(SIGNATURES)

The beamformed prompt of the most recently completed predetection
integration.
"""
get_last_integration_beamformed_prompt(t::TrackedSignal) =
    t.post_correlation_state.last_beamformed_prompt

"""
$(SIGNATURES)

The number of primary code periods the most recently completed predetection
integration spanned.
"""
get_last_integration_num_code_periods(t::TrackedSignal) =
    t.post_correlation_state.last_num_code_periods

"""
$(SIGNATURES)

The predetection integration time `T` of the most recently completed integration
(its primary code periods times one code period), which
[`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{NoiseReferencedCN0Estimator, Any})
divides by. The integration's post-integration SNR is `C/N₀ · T`.

Not to be confused with `get_num_integrated_samples`, the integration in progress.
"""
get_last_integration_time(t::TrackedSignal) =
    get_last_integration_num_code_periods(t) * get_code_length(get_signal(t)) /
    get_chip_rate(get_signal(t))

"""
$(SIGNATURES)

The beamformed prompts of the predetection integrations the most recent `track` call
completed, in order.
"""
get_beamformed_prompts(t::TrackedSignal) = t.beamformed_prompts

"""
$(SIGNATURES)

The [`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput)s this signal
completed in the current chunk, in order. [`step_tracking_loops!`](@ref) consumes
and clears them, so read them inside a custom tracking loop or right after a bare
`downconvert_and_correlate!`. See [Chunked NCO updates](@ref).
"""
get_correlator_bank_outputs(t::TrackedSignal) = t.correlator_bank_outputs

"""
$(SIGNATURES)

Append an externally built
[`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput) to `signal`'s per-chunk
buffer and return `signal`: the ingest path for an **external correlator
producer** (e.g. FPGA dumps). Append per signal in `sample_index` order, then run
[`step_tracking_loops!`](@ref), which steps the tracking loop with them and clears the buffer.
Type-checked against the signal's correlator bank, unlike pushing to
[`get_correlator_bank_outputs`](@ref). See [External correlator producers](@ref).
"""
append_correlator_bank_output!(t::TrackedSignal, output::CorrelatorBankOutput) =
    (push!(t.correlator_bank_outputs, output); t)

get_beamformer(t::TrackedSignal) = t.post_correlation_state.beamformer
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

The shared `code_phase` wraps at [`max_code_phase_period`](@ref Tracking.max_code_phase_period).
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
    sig = tsig.signal
    get_code_length(sig) *
    max(get_secondary_code_length(sig), get_code_periods_per_bit(sig))
end

# `lcm`, not `max`: per-signal phases are `mod(code_phase, _replica_code_phase_period(...))`,
# correct only if the shared period is a multiple of every signal's (#129).
@inline _max_code_phase_period(::Tuple{}) = 1
@inline _max_code_phase_period(t::Tuple) =
    lcm(_post_sync_code_phase_period(first(t)), _max_code_phase_period(Base.tail(t)))

"""
$(SIGNATURES)

Upper bound on the shared `sat.code_phase` period, in chips: the least
common multiple of the per-signal periods once *every* signal has synced (GPS
L1 C/A: 1023 × 20 = 20460, one data bit; L1C-P: 10230 × 1800, one
secondary-code cycle). Folds to a literal for a concrete `signals` type.

The runtime period, [`current_code_phase_period`](@ref Tracking.current_code_phase_period), is shorter for unsynced signals.
"""
@inline max_code_phase_period(signals::Tuple{Vararg{TrackedSignal}}) =
    _max_code_phase_period(signals)

# Primary code length until synced, then `_post_sync_code_phase_period`.
@inline function _signal_code_phase_period(tsig::TrackedSignal)
    if is_bit_or_secondary_code_synced(tsig)
        _post_sync_code_phase_period(tsig)
    else
        get_code_length(tsig.signal)
    end
end

# `lcm` for the same reason as `_max_code_phase_period`.
@inline _current_code_phase_period(::Tuple{}) = 1
@inline _current_code_phase_period(t::Tuple) =
    lcm(_signal_code_phase_period(first(t)), _current_code_phase_period(Base.tail(t)))

"""
$(SIGNATURES)

The runtime period of the shared `sat.code_phase`, in chips: like
[`max_code_phase_period`](@ref Tracking.max_code_phase_period), but a signal not yet synced (`is_bit_or_secondary_code_synced`
false) contributes only its primary code length, since its bit or secondary
chip is unknown. The least common multiple keeps the wrap a multiple of every
signal's replica wrap (#129).
"""
@inline current_code_phase_period(signals::Tuple{Vararg{TrackedSignal}}) =
    _current_code_phase_period(signals)

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
    sig = head.signal
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
acquisition handoff values, e.g. for a non-default correlator or beamformer; otherwise use the signal-instance forms. The first signal is the
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
    # Accept integer Dopplers such as `200Hz`.
    carrier_doppler = float(carrier_doppler)
    code_doppler =
        isnothing(code_doppler) ?
        carrier_doppler *
        get_chip_rate_to_center_frequency_ratio(first(tracked_signals).signal) :
        float(code_doppler)
    # Built twice so `init_loop_state` sees a real `TrackedSat`.
    bare = TrackedSat(
        prn,
        float(code_phase),
        code_doppler,
        float(carrier_phase),
        carrier_doppler,
        1,
        tracked_signals,
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
    num_antennas::NumAntennas = NumAntennas(1),
    num_integrations_for_cn0_estimation::Int = 100,
    cn0_estimator = nothing,
    beamformer::AbstractBeamformer = DefaultBeamformer(),
    kwargs...,
)
    cn0_estimators = _per_signal_cn0_estimators(
        signals,
        cn0_estimator,
        num_integrations_for_cn0_estimation,
    )
    tracked_signals = map(signals, cn0_estimators) do sig, estimator
        TrackedSignal(
            sig;
            num_antennas,
            correlator_bank = get_default_correlator_bank(sig, num_antennas),
            cn0_estimator = estimator,
            beamformer,
        )
    end
    TrackedSat(tracked_signals, prn, code_phase, carrier_doppler; kwargs...)
end

# One estimator per signal: `nothing` → defaults; a tuple as is; a single
# estimator only for a single signal.
@inline _per_signal_cn0_estimators(
    signals::Tuple{Vararg{AbstractGNSSSignal}},
    ::Nothing,
    num_integrations_for_cn0_estimation::Int,
) = map(sig -> default_cn0_estimator(sig, num_integrations_for_cn0_estimation), signals)

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
        "`cn0_estimator = ($(join(fill("MomentsCN0Estimator(; num_integrations = 100)", length(signals)), ", ")))`.",
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
    num_antennas::NumAntennas = NumAntennas(1),
    correlator_bank::AbstractCorrelatorBank = get_default_correlator_bank(
        signal,
        num_antennas,
    ),
    num_integrations_for_cn0_estimation::Int = 100,
    cn0_estimator::AbstractCN0Estimator = default_cn0_estimator(
        signal,
        num_integrations_for_cn0_estimation,
    ),
    beamformer::AbstractBeamformer = DefaultBeamformer(),
    kwargs...,
)
    tracked_signal =
        TrackedSignal(signal; num_antennas, correlator_bank, cn0_estimator, beamformer)
    TrackedSat((tracked_signal,), prn, code_phase, carrier_doppler; kwargs...)
end

# Kwarg-update constructor; `signals` and `loop_state` keep their types.
function TrackedSat(
    sat::TrackedSat{Signals,D};
    prn = nothing,
    code_phase = nothing,
    code_doppler = nothing,
    carrier_phase = nothing,
    carrier_doppler = nothing,
    next_sample_index = nothing,
    signals::Maybe{Signals} = nothing,
    loop_state::Maybe{D} = nothing,
) where {Signals<:Tuple{Vararg{TrackedSignal}},D}
    TrackedSat{Signals,D}(
        isnothing(prn) ? sat.prn : prn,
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
get_num_antennas(
    sat::TrackedSat{<:Tuple{TrackedSignal{<:Any,<:Any,<:AbstractCorrelatorBank{M}},Vararg}},
) where {M} = M

"""
$(SIGNATURES)

Get the satellite's shared code phase in chips, wrapping at
[`max_code_phase_period`](@ref Tracking.max_code_phase_period). Mod by a signal's primary code length for its replica
phase.
"""
get_code_phase(sat::TrackedSat) = sat.code_phase

"""
$(SIGNATURES)

Get the current code Doppler frequency.
"""
get_code_doppler(sat::TrackedSat) = sat.code_doppler

"""
$(SIGNATURES)

Get the current carrier phase in cycles; tracking wraps it to [-0.5, 0.5).
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

get_signal(sat::TrackedSat, sel...) = get_signal(_find_signal(sat.signals, sel...))
get_correlator_bank(sat::TrackedSat, sel...) =
    get_correlator_bank(_find_signal(sat.signals, sel...))
get_last_integration_correlator_bank(sat::TrackedSat, sel...) =
    get_last_integration_correlator_bank(_find_signal(sat.signals, sel...))
get_last_integration_beamformed_prompt(sat::TrackedSat, sel...) =
    get_last_integration_beamformed_prompt(_find_signal(sat.signals, sel...))
get_last_integration_num_code_periods(sat::TrackedSat, sel...) =
    get_last_integration_num_code_periods(_find_signal(sat.signals, sel...))
get_last_integration_time(sat::TrackedSat, sel...) =
    get_last_integration_time(_find_signal(sat.signals, sel...))
get_beamformed_prompts(sat::TrackedSat, sel...) =
    get_beamformed_prompts(_find_signal(sat.signals, sel...))
get_beamformer(sat::TrackedSat, sel...) = get_beamformer(_find_signal(sat.signals, sel...))
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
    empty!(t.beamformed_prompts)
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
    first(sat.signals).signal,
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
the band routes the right measurement to the group.

Fields:

  - `band`: an `AbstractGNSSSignal` `Band` instance (`L1()`, `L5()`, …)
  - `satellites`: `Dictionary{Int, <:TrackedSat}` keyed by PRN
  - `signals`: the signal-instance tuple (e.g. `(GPSL1C_P(), GPSL1C_D(), GPSL1CA())`)
  - `num_antennas`: the antenna count for this group's band
"""
struct SignalGroup{
    B,                                             # GNSSSignals Band instance
    S<:Dictionary{<:Any,<:TrackedSat},
    Sigs<:Tuple{Vararg{AbstractGNSSSignal}},
    NA<:NumAntennas,
}
    band::B
    satellites::S
    signals::Sigs
    num_antennas::NA
end

# Kwarg-update constructor keeping `g`'s types.
function SignalGroup(
    g::SignalGroup{B,S,Sigs,NA};
    band::Maybe{B} = nothing,
    satellites::Maybe{S} = nothing,
    signals::Maybe{Sigs} = nothing,
    num_antennas::Maybe{NA} = nothing,
) where {
    B,
    S<:Dictionary{<:Any,<:TrackedSat},
    Sigs<:Tuple{Vararg{AbstractGNSSSignal}},
    NA<:NumAntennas,
}
    SignalGroup{B,S,Sigs,NA}(
        isnothing(band) ? g.band : band,
        isnothing(satellites) ? g.satellites : satellites,
        isnothing(signals) ? g.signals : signals,
        isnothing(num_antennas) ? g.num_antennas : num_antennas,
    )
end

# A group's signals must (#129) be on its band, whose measurement the whole group
# correlates against, and share one chip rate, since the shared `code_phase`
# advances at the driver's (`_advance_sat` in downconvert_and_correlate.jl). Bands
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
`get_band(first(signals))`, `num_antennas` to `NumAntennas(1)`. All signals must be on
`band` and share one chip rate (`ArgumentError` otherwise).

`tracking_loop` only shapes the slot type; the tracking loop that runs is the
`TrackState`'s, which rebuilds a mismatched empty slot. Leave it at its default
and configure the tracking loop on `TrackState`.

```julia
SignalGroup((GPSL1CA(),))                              # band L1(), 1 antenna
SignalGroup((GPSL5I(),); num_antennas = NumAntennas(2))  # 2-antenna L5
```
"""
function SignalGroup(
    signals::Tuple{Vararg{AbstractGNSSSignal}};
    band = get_band(first(signals)),
    num_antennas::NumAntennas = NumAntennas(1),
    tracking_loop::AbstractTrackingLoop = _default_tracking_loop(),
)
    _validate_signal_group(signals, band)
    template = _make_template_tracked_sat(signals, tracking_loop, num_antennas)
    sats = Dictionary{Int,typeof(template)}(Int[], typeof(template)[])
    SignalGroup(band, sats, signals, num_antennas)
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

function estimate_cn0(tsig::TrackedSignal)
    # The buffered prompts' SNR is a whole integration's, so divide by its time.
    estimate_cn0(get_cn0_estimator(tsig), get_last_integration_time(tsig))
end

estimate_cn0(sat::TrackedSat, sel...) = estimate_cn0(_find_signal(sat.signals, sel...))
