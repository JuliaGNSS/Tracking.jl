# Largest `BL · Δt` per loop. The default FLL-assisted carrier filter diverges
# at ≈ 0.4 and runs 25 % wider than configured at 0.09; the 0.018 code cap is
# conservative for the second-order filter (stable to ≈ 0.4, Stephens & Thomas
# 1995). See the "Stability cap" section of docs/src/loop_filter.md.
const MAX_CARRIER_LOOP_BANDWIDTH_TIME_PRODUCT = 0.09
const MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT = 0.018

"""
$(SIGNATURES)

Recommended carrier-loop-filter bandwidth for `signal`: a flat 18 Hz, the
third-order PLL bandwidth of the literature. It is the loop's one-sided noise
bandwidth `BL`, seeded for every satellite whose estimator leaves
`carrier_loop_filter_bandwidth` as `nothing`, and capped at filter time by
[`effective_carrier_loop_filter_bandwidth`](@ref).

Override by defining a method for your signal type, or pass
`carrier_loop_filter_bandwidth =` to the estimator.
"""
function default_carrier_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    18.0Hz
end

"""
$(SIGNATURES)

Recommended code-loop-filter (DLL) bandwidth for `signal`: a flat 1 Hz, inside
the 0.25–2 Hz of the reference software receivers. Carrier-aided (see
`aid_dopplers`), the DLL has almost no dynamics to track, so the bandwidth is a
thermal-noise-versus-pull-in trade independent of the signal. Capped at filter
time by [`effective_code_loop_filter_bandwidth`](@ref).

Override by defining a method for your signal type.
"""
function default_code_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    1.0Hz
end

"""
$(SIGNATURES)

The configured carrier bandwidth, capped at `0.09 / integration_time` to keep
the loop stable on long integrations (see the loop-filter docs). The default
18 Hz runs unchanged up to 5 ms of integration, 9 Hz at 10 ms, 4.5 Hz at 20 ms.
"""
@inline function effective_carrier_loop_filter_bandwidth(bandwidth, integration_time)
    min(bandwidth, uconvert(Hz, MAX_CARRIER_LOOP_BANDWIDTH_TIME_PRODUCT / integration_time))
end

"""
$(SIGNATURES)

The configured code bandwidth, capped at `0.018 / integration_time` for
stability. Longer integration must not otherwise narrow the DLL: neither its
dynamics nor its noise floor depend on it. The cap binds only past 18 ms at the
default 1 Hz.
"""
@inline function effective_code_loop_filter_bandwidth(bandwidth, integration_time)
    min(bandwidth, uconvert(Hz, MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT / integration_time))
end

"""
$(SIGNATURES)

Window of the frequency lock indicator for `signal` at the record length
`integration_time`: the FLL is dropped once its mean reading over this window
stays below [`frequency_lock_threshold`](@ref). By default 0.5 s, and at least
four records. Override by defining a method for your signal type.
"""
frequency_lock_window(signal::AbstractGNSSSignal, integration_time) =
    max(0.5s, uconvert(s, 4 * integration_time))

"""
$(SIGNATURES)

Threshold of the frequency lock indicator for `signal` at the record length
`integration_time`, see [`frequency_lock_window`](@ref). By default 3 Hz, and at
most 1/(16T), a quarter of the two-quadrant FLL's range and an eighth of the
four-quadrant one's (0.04 Hz at GPS L2 CL's 1.5 s). Override by defining a method
for your signal type.
"""
frequency_lock_threshold(signal::AbstractGNSSSignal, integration_time) =
    min(3.0Hz, uconvert(Hz, 1 / (16 * integration_time)))

"""
$(TYPEDEF)

The frequency lock indicator of a satellite's carrier loop: FLL-assisted until
`locked` latches, a pure PLL after, until [`reset_loop_filters!`](@ref).
See [Carrier loop staging](@ref).

$(TYPEDFIELDS)
"""
struct FrequencyLockIndicator
    """
    FLL readings integrated over the current window
    """
    integrated_frequency_error::typeof(1.0Hz * 1.0s)
    """
    Length of the current window
    """
    window_time::typeof(1.0s)
    """
    Whether frequency lock has been declared (latched)
    """
    locked::Bool
end

FrequencyLockIndicator() = FrequencyLockIndicator(0.0Hz * 0.0s, 0.0s, false)

# Advance the frequency lock indicator by one record's FLL reading, in windows
# of `frequency_lock_window`. A record without a previous prompt has no FLL
# reading and is left out; a latched lock is kept as is.
@inline function _update_frequency_lock(
    indicator::FrequencyLockIndicator,
    signal::AbstractGNSSSignal,
    fll_discriminator,
    previous_prompt::Complex,
    integration_time,
)
    (indicator.locked || iszero(previous_prompt)) && return indicator
    dt = uconvert(s, integration_time)
    integrated_error = indicator.integrated_frequency_error + fll_discriminator * dt
    window_time = indicator.window_time + dt
    window_time < frequency_lock_window(signal, integration_time) &&
        return FrequencyLockIndicator(integrated_error, window_time, false)
    locked =
        abs(integrated_error / window_time) <
        frequency_lock_threshold(signal, integration_time)
    FrequencyLockIndicator(0.0Hz * 0.0s, 0.0s, locked)
end

"""
Per-satellite state for the conventional PLL and DLL Doppler estimator.
Holds initial Doppler values, loop filter states and the carrier loop's
[`FrequencyLockIndicator`](@ref).
"""
@kwdef struct SatConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
    frequency_lock::FrequencyLockIndicator = FrequencyLockIndicator()
end

function SatConventionalPLLAndDLL(
    sat::TrackedSat,
    carrier_loop_filter::CA,
    code_loop_filter::CO;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        FrequencyLockIndicator(),
    )
end

function SatConventionalPLLAndDLL(
    sat_conventional_pll_and_dll::SatConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    frequency_lock::Maybe{FrequencyLockIndicator} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL{CA,CO}(
        sat_conventional_pll_and_dll.init_carrier_doppler,
        sat_conventional_pll_and_dll.init_code_doppler,
        isnothing(carrier_loop_filter) ? sat_conventional_pll_and_dll.carrier_loop_filter :
        carrier_loop_filter,
        isnothing(code_loop_filter) ? sat_conventional_pll_and_dll.code_loop_filter :
        code_loop_filter,
        isnothing(carrier_loop_filter_bandwidth) ?
        sat_conventional_pll_and_dll.carrier_loop_filter_bandwidth :
        carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ?
        sat_conventional_pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
        isnothing(frequency_lock) ? sat_conventional_pll_and_dll.frequency_lock :
        frequency_lock,
    )
end

"""
$(SIGNATURES)

Conventional Phase-Locked Loop (PLL) and Delay-Locked Loop (DLL) Doppler
estimator. Configuration-only — per-satellite state lives in each
[`TrackedSat`](@ref) wrapper, produced via [`init_estimator_state`](@ref).

Type parameters `CA` and `CO` select the carrier and code loop filter types.
A `nothing` bandwidth (the default) means **auto**: [`init_estimator_state`](@ref)
sizes it per satellite from the sat's estimator-driver signal (`signals[1]`) via
[`default_carrier_loop_filter_bandwidth`](@ref) /
[`default_code_loop_filter_bandwidth`](@ref); an explicit bandwidth applies to
every satellite. At filter time both are capped against the record's integration
time ([`effective_carrier_loop_filter_bandwidth`](@ref),
[`effective_code_loop_filter_bandwidth`](@ref)), so lengthening the coherent
integration with [`set_preferred_num_code_blocks_to_integrate!`](@ref) needs no
re-tuning.
"""
struct ConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
       AbstractDopplerEstimator
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
end

function ConventionalPLLAndDLL(
    ::Type{CA} = ThirdOrderBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{CA,CO}(carrier_loop_filter_bandwidth, code_loop_filter_bandwidth)
end

"""
$(SIGNATURES)

Create a ConventionalPLLAndDLL with FLL-assisted carrier tracking. This is the
default Doppler estimator used by TrackState. Its `ThirdOrderAssistedBilinearLF`
carrier filter combines PLL and FLL discriminators for high dynamics. Bandwidths
default to auto, see [`ConventionalPLLAndDLL`](@ref).
"""
function ConventionalAssistedPLLAndDLL(
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL(
        ThirdOrderAssistedBilinearLF,
        CO;
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Kwarg-update constructor for tweaking bandwidths in place.
function ConventionalPLLAndDLL(
    pll_and_dll::ConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ? pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Build the per-satellite estimator state stored in a [`TrackedSat`](@ref) for a
satellite tracked under [`ConventionalPLLAndDLL`](@ref).

Auto bandwidths (`nothing` on the estimator) are resolved here from the sat's
estimator-driver signal (`signals[1]`), so a multi-group [`TrackState`](@ref)
gets the right bandwidth per group from one shared estimator.
"""
function init_estimator_state(
    estimator::ConventionalPLLAndDLL{CA,CO},
    sat::TrackedSat,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    carrier_loop_filter = constructorof(CA)()
    code_loop_filter = constructorof(CO)()
    driver_signal = first(sat.signals).signal
    carrier_loop_filter_bandwidth =
        isnothing(estimator.carrier_loop_filter_bandwidth) ?
        default_carrier_loop_filter_bandwidth(driver_signal) :
        estimator.carrier_loop_filter_bandwidth
    code_loop_filter_bandwidth =
        isnothing(estimator.code_loop_filter_bandwidth) ?
        default_code_loop_filter_bandwidth(driver_signal) :
        estimator.code_loop_filter_bandwidth
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        FrequencyLockIndicator(),
    )
end

# Re-seed hook used by `reset_loop_filters!`. The fallback rebuilds the per-sat
# state via `init_estimator_state`; estimators may specialize to preserve
# per-sat configuration.
_reset_estimator_state(estimator::AbstractDopplerEstimator, sat::TrackedSat) =
    init_estimator_state(estimator, sat)

# Zero the loop-filter integrators, restart the carrier loop's staging, and
# re-seed the init Dopplers from the sat's current Dopplers, but keep the
# existing per-sat bandwidths so a per-sat override survives the reset.
function _reset_estimator_state(
    ::ConventionalPLLAndDLL,
    sat::TrackedSat{<:Tuple{Vararg{TrackedSignal}},<:SatConventionalPLLAndDLL},
)
    state = sat.doppler_estimator_state
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
        constructorof(typeof(state.carrier_loop_filter))(),
        constructorof(typeof(state.code_loop_filter))(),
        state.carrier_loop_filter_bandwidth,
        state.code_loop_filter_bandwidth,
        FrequencyLockIndicator(),
    )
end

"""
$(SIGNATURES)

Aid dopplers. That is velocity aiding for the carrier doppler and carrier aiding
for the code doppler.
"""
function aid_dopplers(
    signal::AbstractGNSSSignal,
    init_carrier_doppler,
    init_code_doppler,
    carrier_freq_update,
    code_freq_update,
)
    carrier_doppler = carrier_freq_update
    code_doppler =
        code_freq_update + carrier_doppler * get_code_center_frequency_ratio(signal)
    init_carrier_doppler + carrier_doppler, init_code_doppler + code_doppler
end

# Pure per-sat update, shared by every estimator whose per-sat state closes its
# loops through `_close_loops`. Every signal's completed records are folded into
# its prompt filter, CN0 estimator and bit buffer; `signals[1]` (the estimator
# driver) additionally runs the loops and sets the sat-shared Dopplers.
#
# `noise` is one `(density, ready)` pair per signal, in `sat.signals` order (see
# `_signal_noise_densities`).
function _update_tracked_sat_doppler(sat::TrackedSat, sampling_frequency, noise::Tuple)
    estimator_state = sat.doppler_estimator_state
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    # See `_carrier_phase_derotation`.
    driver_carrier_phase_offset = get_carrier_phase_offset(head.signal)

    driver_noise_density, driver_noise_density_ready = first(noise)
    new_head, new_doppler_estimator_state, new_carrier_doppler, new_code_doppler =
        _process_estimator_driver_signal(
            head,
            sat,
            estimator_state,
            sampling_frequency,
            driver_noise_density,
            driver_noise_density_ready,
            driver_carrier_phase_offset,
        )

    new_tail = _process_passenger_signals(
        tail_signals,
        sat.prn,
        sampling_frequency,
        Base.tail(noise),
        driver_carrier_phase_offset,
    )

    # One-time phase snap on the iteration a signal first syncs: anchor
    # `sat.code_phase` to the secondary-chip window of the synced signal with the
    # longest `(primary × secondary)` code, keeping the within-primary-block
    # phase. Re-running it later would wedge the satellite (issue #117); after
    # sync, `update`'s `mod(…, current_code_wrap)` keeps the alignment.
    #
    # Sync is detected after the chunk was correlated, so every signal's
    # in-flight partial was accumulated at the pre-snap phase and would be
    # inconsistent with the new alignment (a flipped NH chip can even cancel it
    # to zero, a 0/0 in the discriminators). It is therefore dropped and
    # re-integrated from the snapped phase.
    new_signals = (new_head, new_tail...)
    just_synced = _any_signal_just_synced(sat.signals, new_signals)
    snapped_code_phase =
        just_synced ? _snap_code_phase_from_synced_signal(new_signals, sat.code_phase) :
        sat.code_phase
    final_signals =
        just_synced ? map(_reset_on_snap, sat.signals, new_signals) : new_signals

    TrackedSat(
        sat;
        code_phase = snapped_code_phase,
        carrier_doppler = new_carrier_doppler,
        code_doppler = new_code_doppler,
        signals = final_signals,
        doppler_estimator_state = new_doppler_estimator_state,
    )
end

# Drop an in-flight (partial) integration at the sync phase snap (see
# `_update_tracked_sat_doppler`), leaving all other per-signal state intact.
@inline _reset_inflight_integration(s::TrackedSignal) =
    TrackedSignal(s; correlator = zero(s.correlator), integrated_samples = 0)

# The reset at the sync phase snap, given a signal before and after the fold. A
# signal whose prompt this sync wipes off (a pilot synced to its secondary code)
# also drops its previous prompt: it was correlated without the wipe-off and
# may differ in sign from the next record, which the four-quadrant FLL would
# read as half a cycle of rotation. A zero prompt makes the FLL skip that
# record, as after `reset_loop_filters!`.
@inline function _reset_on_snap(old::TrackedSignal, new::TrackedSignal)
    s = _reset_inflight_integration(new)
    newly_wiped_off =
        _is_wiped_off(new.signal, new.bit_buffer.found) &&
        !_is_wiped_off(old.signal, old.bit_buffer.found)
    newly_wiped_off ?
    TrackedSignal(
        s;
        last_fully_integrated_filtered_prompt = zero(
            s.last_fully_integrated_filtered_prompt,
        ),
    ) : s
end

# The loops lock the driver onto the real axis; this rotation brings a
# component's bit-buffer prompt back onto it, so a quadrature component (GPS L5 /
# Galileo E5a I vs Q) does not decode off the collapsed real part. It is
# `cis(0) === 1 + 0im`, a bit-identical no-op, for the driver and co-phased
# components; the residual ±90° sign is resolved by the navigation decoder's
# preamble.
@inline _carrier_phase_derotation(driver_carrier_phase_offset::Real, signal) =
    cis(driver_carrier_phase_offset - get_carrier_phase_offset(signal))

# Apply one completed `CorrelatorOutput` record to a signal, shared by the
# driver and passenger folds (issue #133): normalize the record's correlator,
# update/apply the post-corr filter, advance the CN0 estimator and bit buffer,
# and move the record to `last_fully_integrated_*`. Returns the rebuilt signal
# and the filtered correlator for the driver's loops. The live accumulator is
# not touched; the correlate phase already reset it.
#
# The bit accumulator is credited with the blocks *actually* integrated,
# recovered from the record's sample count: the first post-sync integration is
# truncated to the data-bit boundary (issue #125).
#
# `correlated_pre_sync = true` marks a record correlated with a pre-sync replica
# (sync was detected earlier in the same fold). Its prompt lacks the
# secondary-code wipe-off, so it is dropped from the coherent bit sum for
# secondary-coded signals only. Its blocks are always credited: dropping them
# slides the bit window off the navigation-bit grid for the rest of the run
# (issue #219).
@inline function _apply_correlator_output(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase_offset::Real = 0.0;
    correlated_pre_sync::Bool = false,
)
    signal = tracked_signal.signal
    normalized_correlator =
        normalize(output.correlator, output.integrated_samples, get_code_amplitude(signal))
    post_corr_filter =
        update(tracked_signal.post_corr_filter, get_prompt(normalized_correlator))
    # Used both to combine the antennas and to reduce the noise covariance to this
    # satellite's floor below: both sides of the C/N₀ ratio must share one `w`.
    weights = get_weights(post_corr_filter, _num_ants_val(normalized_correlator))
    filtered_correlator = _combine_correlator(normalized_correlator, weights)
    prompt = get_prompt(filtered_correlator)
    push!(tracked_signal.filtered_prompts, prompt)
    bit_block_count = calc_num_code_blocks_for_bit_buffer(
        signal,
        output.integrated_samples,
        sampling_frequency,
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    )
    # Floored at 1 for the fractional-block record after a sync phase-snap reset.
    integrated_code_blocks = max(1, bit_block_count)
    # De-rotated before both the sync search and the coherent bit sum in `buffer`
    # (see `_carrier_phase_derotation`).
    bit_prompt = prompt * _carrier_phase_derotation(driver_carrier_phase_offset, signal)
    drop_prompt = correlated_pre_sync && get_secondary_code_length(signal) > 1
    # The CN0 context carries the bit state from *before* this record, because
    # `NWPRCN0Estimator` sums prompts coherently over exactly one data bit. Where
    # `drop_prompt` holds, the context reports "no bit grid", which keeps the
    # sign-corrupted prompt out of any coherent window, as before sync.
    #
    # It also carries this signal's noise density and this record's integration
    # time for `NoiseRefCN0Estimator`. A not-ready density skips the update only
    # for estimators that read it (`requires_noise_density`, compile-time):
    # skipping the whole record would make a co-resident NWPR drop its open window.
    # The density (scalar, or covariance for an array) is reduced through the
    # prompt's weights to this combiner's scalar floor; `nothing` passes through
    # (see `_reduce_noise_density`).
    scalar_noise_density = _reduce_noise_density(noise_density, weights)
    cn0_estimator = _update_cn0_estimator(
        get_cn0_estimator(tracked_signal),
        prompt,
        signal,
        tracked_signal.bit_buffer,
        bit_block_count,
        !drop_prompt,
        scalar_noise_density,
        noise_density_ready,
        output.integrated_samples / sampling_frequency,
    )
    bit_buffer = buffer(
        signal,
        prn,
        tracked_signal.bit_buffer,
        bit_block_count,
        drop_prompt ? zero(bit_prompt) : bit_prompt,
    )
    # A pre-sync-correlated record also moves the secondary-code anchor: the
    # phase snap after this fold aligns the upcoming integration to
    # `bit_buffer.secondary_phase`, reported for the block after the syncing one.
    if correlated_pre_sync
        bit_buffer = _advance_secondary_phase(signal, bit_buffer, bit_block_count)
    end
    new_signal = TrackedSignal(
        tracked_signal;
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        last_fully_integrated_correlator = output.correlator,
        # What the CN0 estimator's newest prompt was integrated over.
        last_fully_integrated_num_code_blocks = integrated_code_blocks,
    )
    return new_signal, filtered_correlator
end

# Build the context and fold the record into one signal's CN0 estimator, or skip
# it where the estimator needs a density that signal has not measured yet. Both
# branches return the same concrete estimator type, so this stays type-stable
# and allocation-free; the `requires_noise_density` half of the condition folds
# away at compile time.
@inline function _update_cn0_estimator(
    estimator::AbstractCN0Estimator,
    prompt,
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer,
    num_code_blocks::Integer,
    bit_sync_usable::Bool,
    noise_density,
    noise_density_ready::Bool,
    integration_time,
)
    if !noise_density_ready && requires_noise_density(estimator)
        return estimator
    end
    # Positional, not keyword: this is the per-record site, and a keyword call
    # allocates per call on Julia 1.10 (see `CN0UpdateContext`).
    update(
        estimator,
        prompt,
        CN0UpdateContext(
            signal,
            bit_buffer,
            num_code_blocks,
            bit_sync_usable,
            noise_density,
            integration_time,
        ),
    )
end

# Fold the estimator-driver signal's (signals[1]) chunk of `CorrelatorOutput`s
# in order, closing the loops per record and threading the per-sat state and the
# FLL `previous_prompt` across records. Returns the new per-sat state and the
# *last* record's Dopplers (the NCO is written once per chunk); with no outputs
# the Doppler holds. How a record closes the loops is `_close_loops`, dispatched
# on the per-sat state. A custom AbstractDopplerEstimator may use any signal.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    estimator_state,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase_offset::Real = 0.0,
)
    outputs = tracked_signal.correlator_outputs
    if isempty(outputs)
        return tracked_signal, estimator_state, sat.carrier_doppler, sat.code_doppler
    end
    signal = tracked_signal.signal
    ts = tracked_signal
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    found_before_fold = has_bit_or_secondary_code_been_found(ts.bit_buffer)
    # Every record of the chunk was correlated before this fold, so the sync
    # state from before the fold decides which carrier discriminators apply.
    wiped_off = _is_wiped_off(signal, found_before_fold)
    polarity = _sync_polarity(signal, ts.bit_buffer, sat.prn)
    @inbounds for k in eachindex(outputs)
        output = outputs[k]
        # The FLL's previous prompt (the previous chunk's last for the first
        # record); read it BEFORE the advance overwrites it.
        previous_prompt = get_last_fully_integrated_filtered_prompt(ts)
        # Per-record integration time — the block time, NOT the chunk time.
        integration_time = output.integrated_samples / sampling_frequency
        # See `correlated_pre_sync` in `_apply_correlator_output`.
        synced_earlier_in_fold =
            !found_before_fold && has_bit_or_secondary_code_been_found(ts.bit_buffer)
        ts, filtered_correlator = _apply_correlator_output(
            ts,
            output,
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase_offset;
            correlated_pre_sync = synced_earlier_in_fold,
        )

        # Capped against the record's actual integration time, not the intended
        # one: records folded after a mid-fold sync, or the truncated first
        # post-sync integration, are still short.
        carrier_bandwidth = effective_carrier_loop_filter_bandwidth(
            estimator_state.carrier_loop_filter_bandwidth,
            integration_time,
        )
        code_bandwidth = effective_code_loop_filter_bandwidth(
            estimator_state.code_loop_filter_bandwidth,
            integration_time,
        )

        # `dll_disc` gets the chunk-fixed `sat.code_doppler` that generated this
        # chunk's replicas; only the per-sat state threads across records.
        carrier_freq_update, code_freq_update, estimator_state = _close_loops(
            estimator_state,
            signal,
            filtered_correlator,
            previous_prompt,
            sat.code_doppler,
            sampling_frequency,
            integration_time,
            carrier_bandwidth,
            code_bandwidth,
            wiped_off,
            polarity,
        )
        carrier_doppler, code_doppler = aid_dopplers(
            signal,
            estimator_state.init_carrier_doppler,
            estimator_state.init_code_doppler,
            carrier_freq_update,
            code_freq_update,
        )
    end
    empty!(outputs)
    return ts, estimator_state, carrier_doppler, code_doppler
end

# How one record closes the loops, dispatched on the per-sat state: the scalar
# closure unless the state defines its own method.
@inline _close_loops(state, args::Vararg{Any,N}) where {N} =
    _close_scalar_loops(state, args...)

# Close the scalar loops for one record: form the discriminators the filters
# read, run them through the loop filters held in the per-sat state, and return
# the two NCO updates and the state carrying the advanced filters and frequency
# lock indicator. Any per-sat state with the `carrier_loop_filter`,
# `code_loop_filter` and `frequency_lock` fields and a `_with_loop_state` method
# will do.
#
# The carrier loop is staged (see `FrequencyLockIndicator`); `wiped_off` and
# `polarity` pick its discriminators (see `pll_disc`, `fll_disc`).
@inline function _close_scalar_loops(
    state,
    signal::AbstractGNSSSignal,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    code_doppler,
    sampling_frequency,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
    wiped_off::Bool,
    polarity::Int8,
)
    pll_discriminator = pll_disc(signal, correlator; polarity)
    # Only formed for a carrier filter that reads it, until frequency lock.
    frequency_lock = state.frequency_lock
    fll_discriminator = 0.0Hz
    if _uses_fll(state.carrier_loop_filter) && !frequency_lock.locked
        fll_discriminator = fll_disc(
            signal,
            correlator,
            previous_prompt,
            integration_time;
            four_quadrant = wiped_off,
        )
        frequency_lock = _update_frequency_lock(
            frequency_lock,
            signal,
            fll_discriminator,
            previous_prompt,
            integration_time,
        )
    end
    dll_discriminator = dll_disc(signal, correlator, code_doppler, sampling_frequency)
    carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
        state.carrier_loop_filter,
        pll_discriminator,
        fll_discriminator,
        integration_time,
        carrier_bandwidth,
    )
    code_freq_update, code_loop_filter = calculate_code_frequency_update(
        state.code_loop_filter,
        dll_discriminator,
        integration_time,
        code_bandwidth,
    )
    carrier_freq_update,
    code_freq_update,
    _with_loop_state(state; carrier_loop_filter, code_loop_filter, frequency_lock)
end

# A per-sat state with some of its fields replaced, through the state type's
# keyword-update constructor.
@inline _with_loop_state(state::SatConventionalPLLAndDLL; kwargs...) =
    SatConventionalPLLAndDLL(state; kwargs...)

# Process the non-driver signals (signals[2:end]): the per-signal advance only,
# no loop filtering. Recursive over the tuple for type stability, stepping the
# `(density, ready)` tuple in lockstep.
@inline _process_passenger_signals(::Tuple{}, ::Integer, _, ::Tuple{}, ::Real) = ()
@inline function _process_passenger_signals(
    signals::Tuple,
    prn::Integer,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
)
    noise_density, noise_density_ready = first(noise)
    new_head = _process_one_passenger_signal(
        first(signals),
        prn,
        sampling_frequency,
        noise_density,
        noise_density_ready,
        driver_carrier_phase_offset,
    )
    (
        new_head,
        _process_passenger_signals(
            Base.tail(signals),
            prn,
            sampling_frequency,
            Base.tail(noise),
            driver_carrier_phase_offset,
        )...,
    )
end

@inline function _process_one_passenger_signal(
    tracked_signal::TrackedSignal,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase_offset::Real = 0.0,
)
    outputs = tracked_signal.correlator_outputs
    isempty(outputs) && return tracked_signal
    ts = tracked_signal
    found_before_fold = has_bit_or_secondary_code_been_found(ts.bit_buffer)
    @inbounds for k in eachindex(outputs)
        # See `correlated_pre_sync` in `_apply_correlator_output`.
        synced_earlier_in_fold =
            !found_before_fold && has_bit_or_secondary_code_been_found(ts.bit_buffer)
        ts = first(
            _apply_correlator_output(
                ts,
                outputs[k],
                prn,
                sampling_frequency,
                noise_density,
                noise_density_ready,
                driver_carrier_phase_offset;
                correlated_pre_sync = synced_earlier_in_fold,
            ),
        )
    end
    empty!(outputs)
    ts
end

"""
$(SIGNATURES)

Estimate Dopplers and filter prompts for all satellites where the correlation has reached
the end of the code or multiples of that, using the conventional PLL and DLL. The
Dopplers drive the next replicas; the prompts go through the configured post
correlation filter. Satellites without a completed integration are passed through
unchanged.

The second argument is either a [`BandMeasurements`](@ref) NamedTuple (the `track`
path) or a bare per-band sampling-frequency source, a `NamedTuple`/`Dict` keyed by
`get_band_id`. The latter is the entry point for an **external correlator
producer** (e.g. an FPGA), which needs no sample buffer. See
[External correlator producers](@ref).
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    # Copy the slot *values* but share the key set, which this step never
    # changes (it is detached once at the `track` boundary, #123), then delegate
    # to the in-place form.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, sampling_frequencies)
end

# Per-band sampling frequency for a group, looked up by its band id in either
# a `BandMeasurements` NamedTuple or a bare per-band rate source.
@inline _band_sampling_frequency(m::BandMeasurements, key) = m[key].sampling_frequency
@inline _band_sampling_frequency(fs::NamedTuple, key) = fs[key]
@inline _band_sampling_frequency(fs::AbstractDict, key) = fs[key]

# Per-signal noise density as a union-free `(density, ready)` pair. Keyed by
# signal id: the post-correlation floor depends on the despreading modulation,
# not the RF band (see `AbstractNoiseEstimator`).
#
#   * no entry (no noise estimator configured, a static wiring mistake):
#     `(nothing, true)`, so `NoiseRefCN0Estimator` throws at the first record
#     instead of being skipped forever.
#   * an empty window (a runtime condition): a dummy density of the right type
#     with `ready = false`, so the fold skips the estimators that read it.
#   * otherwise the measured density.
@inline _signal_noise_density(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? _noise_density_and_ready(noise_estimators[K]) :
    (nothing, true)

# One `(density, ready)` pair per signal, in the order of every satellite's
# `signals` tuple, so the fold pairs them positionally without a per-sat lookup.
# The shape is compile-time (signal ids come from the type parameters).
@inline _signal_noise_densities(
    noise_estimators::NamedTuple,
    ::Type{<:TrackedSat{Signals}},
) where {Signals} = _signal_noise_densities(noise_estimators, Signals)
@inline _signal_noise_densities(::NamedTuple, ::Type{Tuple{}}) = ()
@inline _signal_noise_densities(noise_estimators::NamedTuple, ::Type{T}) where {T<:Tuple} =
    (
        _signal_noise_density(
            noise_estimators,
            Val(get_signal_id(_signal_type(Base.tuple_type_head(T)))),
        ),
        _signal_noise_densities(noise_estimators, Base.tuple_type_tail(T))...,
    )

# Per-group body for the doppler estimator, a named function so `_foreach_group!`
# can call it without boxing over heterogeneous groups (e.g. GPS L1 + Galileo
# E1B). The sampling rate is per band, the noise density per signal.
@inline function _est_one_group!(
    g::SignalGroup,
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
    noise_estimators::NamedTuple,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    sampling_frequency = _band_sampling_frequency(sampling_frequencies, get_band_id(g.band))
    noise = _signal_noise_densities(noise_estimators, eltype(g.satellites))
    _warn_noise_density_missing(eltype(g.satellites), noise, noise_estimators)
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_doppler(vals[i], sampling_frequency, noise)
    end
    return nothing
end

# Warn when a signal whose CN0 estimator needs a noise density has none usable
# (it would report `-Inf dB-Hz`). The fold cannot tell the two causes apart, so
# the message names both: an empty window (e.g. `append_noise_observation!`
# never called) or a zero floor (a dead input or underrun, see
# `_noise_density_and_ready`). This is the only point that knows a fold ran
# with an unusable density.
#
# Warn rather than throw: unlike the static no-estimator case, this has
# legitimate transient causes (a producer that folds before it appends, a short
# buffer). A multi-antenna window still filling to its dimension count is a
# normal startup transient and is not warned about (`_noise_window_filling`).
# `maxlog` is per callsite, so `_id` is made signal-specific to keep one
# signal's warning from silencing another's.
@inline _warn_noise_density_missing(
    ::Type{<:TrackedSat{Signals}},
    noise::Tuple,
    noise_estimators::NamedTuple,
) where {Signals} = _warn_noise_density_missing(Signals, noise, noise_estimators)
@inline _warn_noise_density_missing(::Type{Tuple{}}, ::Tuple{}, ::NamedTuple) = nothing
@inline function _warn_noise_density_missing(
    ::Type{T},
    noise::Tuple,
    noise_estimators::NamedTuple,
) where {T<:Tuple}
    head = Base.tuple_type_head(T)
    if !last(first(noise)) && requires_noise_density(_cn0_estimator_type(head))
        signal_id = get_signal_id(_signal_type(head))
        _noise_window_still_filling(noise_estimators, Val(signal_id)) ||
            _emit_noise_density_warning(signal_id)
    end
    _warn_noise_density_missing(Base.tuple_type_tail(T), Base.tail(noise), noise_estimators)
end

# `haskey` on a NamedTuple is compile-time, so the no-estimator case (which must
# warn) folds away.
@inline _noise_window_still_filling(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? _noise_window_filling(noise_estimators[K]) : false

@noinline function _emit_noise_density_warning(signal_id::Symbol)
    @warn """
          Signal `:$signal_id` has a noise estimator but no usable noise density, so \
          `NoiseRefCN0Estimator` will not update and C/N₀ stays at `-Inf dB-Hz`. \
          Either its window is still empty — call `append_noise_observation!` before \
          `estimate_dopplers_and_filter_prompt!`, or use a noise estimator that \
          measures from the band's samples — or the measured floor is zero, which \
          means the samples carry no power at all (a dead input or an underrun).""" _id =
        Symbol(:no_noise_density_, signal_id) maxlog = 1
    nothing
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref). Walks each
group's `Vector{TrackedSat}` backing storage and overwrites slots with the
new immutable `TrackedSat` value. Returns the same `track_state` object —
allocation-free in steady state when [`track!`](@ref)'s preconditions are met.

The second argument is as for the immutable form; only the per-band sampling
frequency is read from it. Each signal's `correlator_outputs` are consumed and
**cleared**, whether produced by the software correlate phase or appended via
[`append_correlator_output!`](@ref).
"""
function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    _foreach_group!(
        _est_one_group!,
        track_state.groups,
        sampling_frequencies,
        track_state.noise_estimators,
    )
    return track_state
end

# The filter half of a loop update. The filter's coefficients assume consistent
# units, so it is fed the phase error in cycles (`pll_disc`) and the FLL error in
# Hz to output a Doppler in Hz. Fed radians, the loop gain was 2π too high
# (#244). `fll_discriminator` is ignored by a filter that does not read it (see
# `_uses_fll`).
calculate_carrier_frequency_update(
    carrier_loop_filter::AbstractLoopFilter,
    pll_discriminator,
    fll_discriminator,
    integration_time,
    loop_bandwidth,
) = filter_loop(
    carrier_loop_filter,
    _uses_fll(carrier_loop_filter) ? (pll_discriminator, fll_discriminator) :
    pll_discriminator,
    integration_time,
    loop_bandwidth,
)

calculate_code_frequency_update(
    code_loop_filter::AbstractLoopFilter,
    dll_discriminator,
    integration_time,
    loop_bandwidth,
) = filter_loop(code_loop_filter, dll_discriminator, integration_time, loop_bandwidth)

# Whether the carrier loop filter reads an FLL discriminator: the FLL-assisted
# one does, any other carrier filter reads the PLL one alone.
_uses_fll(::AbstractLoopFilter) = false
_uses_fll(::ThirdOrderAssistedBilinearLF) = true
