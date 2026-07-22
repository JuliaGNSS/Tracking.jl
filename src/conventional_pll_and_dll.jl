"""
$(SIGNATURES)

Recommended carrier-loop-filter bandwidth for `signal`'s primary integration
period. Sized so that the PLL time-bandwidth product `BL * T` lands at
about 0.018 (≈10× margin from the 0.18 stability edge of the bilinear
third-order filter). Used by [`TrackState(; signal=…)`](@ref) when the
user doesn't pass an explicit `doppler_estimator`.

Override by defining a method for your signal type, or by constructing
[`ConventionalAssistedPLLAndDLL`](@ref) yourself with explicit
`carrier_loop_filter_bandwidth =` / `code_loop_filter_bandwidth =` kwargs.

```julia
T = get_code_length(signal) / get_code_frequency(signal)   # primary period
BL = 0.018 / T                                              # this default
```

`T` here is the **primary**-code period, not the chosen coherent
integration length. For GPS L1 C/A (T = 1 ms) and GPS L5I (T = 1 ms, a
10230-chip code at 10.23 MHz) this returns 18 Hz — matching the historical
hand-picked default. For L1C-D / L1C-P (T = 10 ms) it returns 1.8 Hz, and
for Galileo E1B (T = 4 ms) 4.5 Hz — the well-inside-stability values the
multi-signal flagship use case needs.

This value is the **reference** bandwidth for a one-primary-code-period
integration; it is not the bandwidth that ends up in the loop when you
integrate longer. Coherently integrating `N` primary blocks grows the loop
update interval to `N·T`, which would push `BL·N·T` toward the ~0.18
stability edge of the bilinear filter. To avoid that, the conventional
estimator **automatically scales the effective loop bandwidth by
`1/N`** at filter time (see [`ConventionalPLLAndDLL`](@ref)), holding the
`BL·Δt` stability product fixed at its single-period value. So you set this
reference bandwidth once and the loop stays stable at any integration length
— no manual `1/N` adjustment is needed.
"""
function default_carrier_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    # T = primary_code_period_seconds. The estimator's bandwidth fields are
    # typed `typeof(1.0Hz)`, so explicitly land on Hz (otherwise `1/s`
    # propagates and trips the typed field assignment).
    primary_period = get_code_length(signal) / get_code_frequency(signal)
    uconvert(Hz, 0.018 / primary_period)
end

"""
$(SIGNATURES)

Recommended code-loop-filter (DLL) bandwidth for `signal`'s primary
integration period. Picks 1/18 of [`default_carrier_loop_filter_bandwidth`](@ref)
— the historical 18:1 carrier:code-bandwidth ratio that gives the DLL good
noise rejection without lagging the PLL.

For L1 C/A and L5I (T = 1 ms) this returns 1 Hz; for L1C-D / L1C-P
(T = 10 ms) it returns 0.1 Hz.
Override by defining a method for your signal type.
"""
function default_code_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    default_carrier_loop_filter_bandwidth(signal) / 18
end

"""
Per-satellite state for the conventional PLL and DLL Doppler estimator.
Holds initial Doppler values and loop filter states.
"""
@kwdef struct SatConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
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
    )
end

function SatConventionalPLLAndDLL(
    sat_conventional_pll_and_dll::SatConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
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
    )
end

"""
$(SIGNATURES)

Conventional Phase-Locked Loop (PLL) and Delay-Locked Loop (DLL) Doppler
estimator. Configuration-only — per-satellite state lives in each
[`TrackedSat`](@ref) wrapper, produced via [`init_estimator_state`](@ref).

Type parameters `CA` and `CO` select the carrier and code loop filter types;
the bandwidth fields configure the loop bandwidths used when seeding new
satellites. Each bandwidth field is `Maybe{typeof(1.0Hz)}`: a `nothing`
field (the default) means **auto** — [`init_estimator_state`](@ref) sizes the
bandwidth per satellite from that sat's estimator-driver signal (`signals[1]`)
via [`default_carrier_loop_filter_bandwidth`](@ref) /
[`default_code_loop_filter_bandwidth`](@ref), so each signal gets a loop sized
for its own integration period (18 Hz for GPS L1 C/A, 4.5 Hz for Galileo E1B,
1.8 Hz for L1C-D / L1C-P, …). Pass an explicit bandwidth to override the
auto-sizing for every satellite this estimator seeds.

The bandwidth is referenced to a **one-primary-code-period** integration.
When a signal coherently integrates `N` primary blocks (its
per-[`TrackedSignal`](@ref) `preferred_num_code_blocks_to_integrate`, set via
[`set_preferred_num_code_blocks_to_integrate!`](@ref)), the effective loop
bandwidth is automatically scaled to `BL/N` at filter time so the loop's
`BL·Δt` stability product stays at its single-period value. This keeps the
loop stable across integration lengths without the caller re-tuning the
bandwidth — e.g. a 1 ms→10 ms switch needs no bandwidth change.
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
default Doppler estimator used by TrackState. Uses a ThirdOrderAssistedBilinearLF
for the carrier loop filter which combines PLL and FLL discriminators for
improved tracking under high dynamics.

Bandwidths default to `nothing` (auto): each satellite is seeded with the
loop bandwidth recommended for its own estimator-driver signal — see
[`ConventionalPLLAndDLL`](@ref). Pass explicit bandwidths to override.
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

Auto bandwidths (`nothing` on the estimator) are resolved here, per satellite,
from the sat's estimator-driver signal (`signals[1]`): each sat gets the loop
bandwidth recommended for the signal that actually drives its loop, so a
multi-group / multi-constellation [`TrackState`](@ref) ends up with the right
bandwidth per group even though it carries one shared estimator. An explicit
bandwidth on the estimator is used verbatim for every satellite.
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
    )
end

# Re-seed hook used by `reset_loop_filters!`. The generic fallback simply
# rebuilds the per-sat state from scratch via `init_estimator_state`; custom
# estimators may specialize to preserve per-sat configuration across the
# reset.
_reset_estimator_state(estimator::AbstractDopplerEstimator, sat::TrackedSat) =
    init_estimator_state(estimator, sat)

# Conventional PLL/DLL: zero the loop-filter integrators and re-seed the
# init Dopplers from the sat's current (converged) Dopplers, but keep the
# bandwidths from the EXISTING per-sat state — a per-sat
# `SatConventionalPLLAndDLL` bandwidth override must survive the reset
# (going through `init_estimator_state` would silently revert it to the
# estimator-level defaults).
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

# Per-sat update for the conventional PLL/DLL estimator. Pure: takes a
# TrackedSat and returns the updated TrackedSat. Shared by the immutable
# `estimate_dopplers_and_filter_prompt` and the in-place
# `estimate_dopplers_and_filter_prompt!` so the two cannot drift.
function _update_tracked_sat_doppler(sat::TrackedSat, sampling_frequency)
    # Walk all signals. For each one whose integration completed this
    # iteration, normalize/filter its prompt, advance CN0 and bit buffer,
    # and move its correlator to `last_fully_integrated_*`. Additionally,
    # for `signals[1]` (the estimator-driver signal), run PLL/DLL and
    # update the sat-shared carrier/code Doppler. Each signal's coherent-
    # integration length comes from its own `preferred_num_code_blocks_to_integrate`.
    pll_and_dll_state = sat.doppler_estimator_state
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    # The loops lock the driver (`signals[1]`) onto the real axis; every signal's
    # bit-buffer prompt is de-rotated by its carrier-phase offset from the driver
    # so a quadrature component (QPSK data/pilot, e.g. GPS L5 / Galileo E5a) does
    # not decode off the collapsed real part. The per-signal carrier phase comes
    # from `get_carrier_phase_offset`.
    driver_carrier_phase = get_carrier_phase_offset(head.signal)

    new_head, new_doppler_estimator_state, new_carrier_doppler, new_code_doppler =
        _process_estimator_driver_signal(
            head,
            sat,
            pll_and_dll_state,
            sampling_frequency,
            driver_carrier_phase,
        )

    new_tail = _process_passenger_signals(
        tail_signals,
        sat.prn,
        sampling_frequency,
        driver_carrier_phase,
    )

    # Phase-snap fallback chain. Picks the synced signal with the
    # longest `(primary × secondary)` code length, and uses its
    # secondary-code phase to anchor `sat.code_phase` to the right
    # secondary-chip window.
    #
    # This is a *one-time* anchoring applied only on the iteration a
    # signal transitions `found == false → true`. It preserves the
    # within-primary-block phase (`mod(code_phase, primary)`) so the loop
    # keeps the current chunk-bounded position; re-running it on later
    # iterations would wedge the satellite (see issue #117). After sync,
    # `update`'s `mod(…, current_code_wrap)` maintains the alignment.
    #
    # Because sync is detected in this estimate pass — *after* the whole
    # chunk was correlated — any in-flight partial integration for this
    # chunk was accumulated at the pre-snap code phase (and, pre-sync, with
    # no secondary-code overlay). Once the snap jumps `code_phase` into the
    # secondary window that partial's data is phase-inconsistent with the
    # new alignment (for a flipped NH chip it can even cancel to zero and
    # feed a 0/0 into the discriminators). So on the snap we also reset every
    # signal's in-flight accumulator: the phase bookkeeping is kept, but the
    # next chunk re-integrates cleanly from the snapped phase to the next
    # boundary. Block-aligned starts have no residue, so this is a no-op there.
    new_signals = (new_head, new_tail...)
    just_synced = _any_signal_just_synced(sat.signals, new_signals)
    snapped_code_phase =
        just_synced ? _snap_code_phase_from_synced_signal(new_signals, sat.code_phase) :
        sat.code_phase
    final_signals =
        just_synced ? map(_reset_inflight_integration, new_signals) : new_signals

    TrackedSat(
        sat;
        code_phase = snapped_code_phase,
        carrier_doppler = new_carrier_doppler,
        code_doppler = new_code_doppler,
        signals = final_signals,
        doppler_estimator_state = new_doppler_estimator_state,
    )
end

# Drop an in-flight (partial) integration: zero the accumulator and its sample
# counter, leaving all other per-signal state intact. Used at the sync-transition
# phase snap, where the shared `code_phase` moves and any partial accumulated at
# the old phase must not be carried into the re-anchored window.
@inline _reset_inflight_integration(s::TrackedSignal) =
    TrackedSignal(s; correlator = zero(s.correlator), integrated_samples = 0)

# De-rotation applied to a component's bit-buffer prompt so its own energy is
# real again, given the loops lock the driver onto the real axis. The rotation
# is `cis(driver_carrier_phase − get_carrier_phase(signal))`, where the
# per-signal carrier phase (radians, relative to the band's in-phase reference)
# comes from `get_carrier_phase_offset`.
# For an in-phase component (co-phased with the driver, or the driver itself)
# the difference is 0 and `cis(0) === 1 + 0im`, a bit-identical no-op; a
# quadrature component (GPS L5 / Galileo E5a I-vs-Q) rotates by `±90°` onto the
# real axis, where the navigation decoder resolves the residual sign via its
# preamble.
@inline _carrier_phase_derotation(driver_carrier_phase::Real, signal) =
    cis(driver_carrier_phase - get_carrier_phase_offset(signal))

# Apply one completed `CorrelatorOutput` record to a signal — shared by the
# estimator-driver and passenger folds so they cannot drift (issue #133):
# normalize the record's (raw) correlator by its sample count, update/apply the
# post-corr filter, record the filtered prompt, advance the CN0 estimator and
# bit buffer, and rebuild the `TrackedSignal` with the record moved to
# `last_fully_integrated_*`. Returns the rebuilt signal plus the intermediate
# values the driver's loop-filter section needs (`filtered_correlator`,
# `integrated_code_blocks`).
#
# Unlike the old per-integration advance, this does NOT reset the live
# accumulator or `integrated_samples`: the correlate phase already reset them
# when it snapshotted this record and began (or is carrying) the next
# integration. It consumes only the record's stored correlator.
#
# The bit accumulator is credited with the blocks *actually* integrated
# (`calc_num_code_blocks_for_bit_buffer`), not the intended
# `integrated_code_blocks`: post-sync the first integration is truncated to
# land on the data-bit boundary, so crediting the intended length would
# misalign the decoded bits (issue #125). The intended count is still
# returned for the driver's `1/N` loop-bandwidth scaling.
# `skip_bit_buffer = true` applies everything EXCEPT the bit-buffer update.
# Used for records that follow a bit/secondary sync detected earlier in the
# same fold: those records were correlated with pre-sync replicas (no
# secondary-code wipe-off), so accumulating them as if wiped would feed
# sign-corrupted prompts into the first post-sync bits. They re-enter cleanly
# next chunk, correlated at the snapped phase with the overlay in the replica.
@inline function _apply_correlator_output(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    prn::Integer,
    sampling_frequency,
    driver_carrier_phase::Real = 0.0;
    skip_bit_buffer::Bool = false,
)
    signal = tracked_signal.signal
    integrated_code_blocks = calc_num_code_blocks_to_integrate(
        signal,
        tracked_signal.preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    )
    normalized_correlator =
        normalize(output.correlator, output.integrated_samples, get_code_amplitude(signal))
    post_corr_filter =
        update(tracked_signal.post_corr_filter, get_prompt(normalized_correlator))
    filtered_correlator = apply(post_corr_filter, normalized_correlator)
    prompt = get_prompt(filtered_correlator)
    push!(tracked_signal.filtered_prompts, prompt)
    cn0_estimator = update(get_cn0_estimator(tracked_signal), prompt)
    bit_block_count = calc_num_code_blocks_for_bit_buffer(
        signal,
        output.integrated_samples,
        sampling_frequency,
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    )
    # De-rotate the prompt onto the driver's (real) phase frame before both the
    # secondary/bit sync search and the coherent bit accumulation inside
    # `buffer`, so a quadrature component's data lands on the real axis it is
    # decided on. No-op for the driver and for co-phased pairs.
    bit_prompt = prompt * _carrier_phase_derotation(driver_carrier_phase, signal)
    bit_buffer =
        skip_bit_buffer ? tracked_signal.bit_buffer :
        buffer(signal, prn, tracked_signal.bit_buffer, bit_block_count, bit_prompt)
    new_signal = TrackedSignal(
        tracked_signal;
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        last_fully_integrated_correlator = output.correlator,
    )
    return new_signal, filtered_correlator, integrated_code_blocks
end

# Process the estimator-driver signal (signals[1]): fold over every
# `CorrelatorOutput` collected during this chunk, in order — running the PLL/DLL
# plus prompt filter / CN0 / bit-buffer update per record and threading the loop
# filters and FLL `previous_prompt` across them — then return the new
# doppler_estimator_state and the *last* record's carrier/code Doppler (the NCO
# is written once per chunk). With no outputs the Doppler holds. This is where
# ConventionalPLLAndDLL hard-codes the "signals[1] drives the loop filter" rule —
# a custom AbstractDopplerEstimator may use any/all signals' state.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    pll_and_dll_state::SatConventionalPLLAndDLL,
    sampling_frequency,
    driver_carrier_phase::Real = 0.0,
)
    outputs = tracked_signal.correlator_outputs
    if isempty(outputs)
        return tracked_signal, pll_and_dll_state, sat.carrier_doppler, sat.code_doppler
    end
    signal = tracked_signal.signal
    ts = tracked_signal
    carrier_loop_filter = pll_and_dll_state.carrier_loop_filter
    code_loop_filter = pll_and_dll_state.code_loop_filter
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    found_before_fold = has_bit_or_secondary_code_been_found(ts.bit_buffer)
    @inbounds for k in eachindex(outputs)
        output = outputs[k]
        # FLL needs the previous record's filtered prompt; the first record of
        # the chunk chains from the sat's carried-over
        # `last_fully_integrated_filtered_prompt` (the previous chunk's last).
        # Read it off `ts` BEFORE the advance overwrites it.
        previous_prompt = get_last_fully_integrated_filtered_prompt(ts)
        # Per-record integration time — the block time, NOT the chunk time.
        integration_time = output.integrated_samples / sampling_frequency
        # A record that follows a sync detected earlier in THIS fold was
        # correlated with pre-sync replicas — keep it out of the bit buffer
        # (see `_apply_correlator_output`).
        synced_earlier_in_fold =
            !found_before_fold && has_bit_or_secondary_code_been_found(ts.bit_buffer)
        # The driver de-rotates against itself (offset 0), so the derotation is a
        # no-op for it; passed for symmetry with the passenger path.
        ts, filtered_correlator, integrated_code_blocks = _apply_correlator_output(
            ts,
            output,
            sat.prn,
            sampling_frequency,
            driver_carrier_phase;
            skip_bit_buffer = synced_earlier_in_fold,
        )

        # The configured bandwidths are referenced to a one-primary-code-period
        # integration. Coherently integrating `integrated_code_blocks` periods
        # grows the loop update interval by that factor, so scale the effective
        # bandwidth by 1/integrated_code_blocks to hold the loop's BL·Δt
        # stability product at its single-period value. For the N=1 path this
        # divides by 1 and is bit-identical to before.
        carrier_bandwidth =
            pll_and_dll_state.carrier_loop_filter_bandwidth / integrated_code_blocks
        code_bandwidth =
            pll_and_dll_state.code_loop_filter_bandwidth / integrated_code_blocks

        carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
            signal,
            carrier_loop_filter,
            filtered_correlator,
            previous_prompt,
            integration_time,
            carrier_bandwidth,
        )
        # `dll_disc` is fed the chunk-fixed `sat.code_doppler` — the code Doppler
        # that actually generated this chunk's replicas — for every record;
        # only the loop-filter *state* threads across records.
        code_freq_update, code_loop_filter = calculate_code_frequency_update(
            signal,
            code_loop_filter,
            filtered_correlator,
            sat.code_doppler,
            sampling_frequency,
            integration_time,
            code_bandwidth,
        )
        carrier_doppler, code_doppler = aid_dopplers(
            signal,
            pll_and_dll_state.init_carrier_doppler,
            pll_and_dll_state.init_code_doppler,
            carrier_freq_update,
            code_freq_update,
        )
    end
    empty!(outputs)
    new_doppler_estimator_state =
        SatConventionalPLLAndDLL(pll_and_dll_state; carrier_loop_filter, code_loop_filter)
    return ts, new_doppler_estimator_state, carrier_doppler, code_doppler
end

# Process the non-driver signals (signals[2:end]): the shared per-signal
# advance only — no loop-filter work. Walks the tuple recursively to keep
# type-stability and avoid boxing.
@inline _process_passenger_signals(::Tuple{}, ::Integer, _, ::Real) = ()
@inline function _process_passenger_signals(
    signals::Tuple,
    prn::Integer,
    sampling_frequency,
    driver_carrier_phase::Real,
)
    head = first(signals)
    new_head =
        _process_one_passenger_signal(head, prn, sampling_frequency, driver_carrier_phase)
    (
        new_head,
        _process_passenger_signals(
            Base.tail(signals),
            prn,
            sampling_frequency,
            driver_carrier_phase,
        )...,
    )
end

@inline function _process_one_passenger_signal(
    tracked_signal::TrackedSignal,
    prn::Integer,
    sampling_frequency,
    driver_carrier_phase::Real = 0.0,
)
    outputs = tracked_signal.correlator_outputs
    isempty(outputs) && return tracked_signal
    ts = tracked_signal
    found_before_fold = has_bit_or_secondary_code_been_found(ts.bit_buffer)
    @inbounds for k in eachindex(outputs)
        # Same rule as the driver fold: records after a sync detected earlier
        # in this fold stay out of the bit buffer.
        synced_earlier_in_fold =
            !found_before_fold && has_bit_or_secondary_code_been_found(ts.bit_buffer)
        ts = first(
            _apply_correlator_output(
                ts,
                outputs[k],
                prn,
                sampling_frequency,
                driver_carrier_phase;
                skip_bit_buffer = synced_earlier_in_fold,
            ),
        )
    end
    empty!(outputs)
    ts
end

"""
$(SIGNATURES)

Estimate Dopplers and filter prompts for all satellites where the correlation has reached
the end of the code or multiples of that. This function uses the
conventional PLL and DLL implementation to estimate Dopplers for
carrier and code. Those Doppler estimations will be used to create the next
replicas to downconvert and decode the incoming signal. In addition to the
Doppler estimation it will also filter the prompt with the configured
post correlation filter.
In the case that the that the correlation hasn't reached the end, e.g. in the case
the incoming signal did not provide enough samples, it will return struct with
zeroed values.
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    measurements::BandMeasurements,
)
    # Detach the slot *values* from the input (sharing the key set), then
    # delegate to the in-place form. This step never changes the key set, so
    # sharing the `Indices` is safe and avoids copying the hash table every
    # `track` loop iteration; the key set is detached once at the `track`
    # boundary (`reset_start_sample_and_bit_buffer`, #123). The per-sat
    # doppler update is identical between the two forms — only the storage
    # ownership differs.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, measurements)
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref). Walks each
group's `Vector{TrackedSat}` backing storage and overwrites slots with the
new immutable `TrackedSat` value. Returns the same `track_state` object —
allocation-free in steady state when [`track!`](@ref)'s preconditions are met.
"""
# Per-group body for the doppler estimator. Pulled out so
# `_foreach_group!` can call it without boxing when the groups tuple
# is heterogeneous (e.g. GPS L1 + Galileo E1B). The per-signal type
# is recovered from each signal inside `_update_tracked_sat_doppler`.
# Routes to this group's band's `BandMeasurement` for sampling frequency.
@inline function _est_one_group!(g::SignalGroup, measurements::BandMeasurements)
    vals = g.satellites.values
    isempty(vals) && return nothing
    sampling_frequency = measurements[get_band_id(g.band)].sampling_frequency
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_doppler(vals[i], sampling_frequency)
    end
    return nothing
end

function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    measurements::BandMeasurements,
)
    _foreach_group!(_est_one_group!, track_state.groups, measurements)
    return track_state
end

function calculate_carrier_frequency_update(
    signal::AbstractGNSSSignal,
    carrier_loop_filter::ThirdOrderAssistedBilinearLF,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(signal, correlator)
    fll_discriminator = fll_disc(signal, correlator, previous_prompt, integration_time)
    filter_loop(
        carrier_loop_filter,
        (pll_discriminator, fll_discriminator),
        integration_time,
        loop_bandwidth,
    )
end

function calculate_carrier_frequency_update(
    signal::AbstractGNSSSignal,
    carrier_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(signal, correlator)
    filter_loop(carrier_loop_filter, pll_discriminator, integration_time, loop_bandwidth)
end

function calculate_code_frequency_update(
    signal::AbstractGNSSSignal,
    code_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    code_doppler,
    sampling_frequency,
    integration_time,
    loop_bandwidth,
)
    dll_discriminator = dll_disc(signal, correlator, code_doppler, sampling_frequency)
    filter_loop(code_loop_filter, dll_discriminator, integration_time, loop_bandwidth)
end
