# The TrackState-level plumbing of the per-record Doppler estimators. The
# estimators themselves — configuration, per-satellite state and the per-record
# `step` — are TrackingLoops'; what this file adds is the walk over a
# `TrackedSat`'s signals and over a `TrackState`'s groups, and the mapping of a
# `TrackedSignal` onto the bare per-record state the shared fold operates on.

# Re-seed hook used by `reset_loop_filters!`. The generic fallback simply
# rebuilds the per-sat state from scratch via `init_estimator_state`; custom
# estimators may specialize to preserve per-sat configuration across the reset.
_reset_estimator_state(estimator::AbstractDopplerEstimator, sat::TrackedSat) =
    init_estimator_state(estimator, sat)

# The conventional and the NCO-referenced loops: zero the integrators and
# re-seed the init Dopplers from the sat's current (converged) Dopplers, but
# keep the bandwidths from the EXISTING per-sat state — a per-sat bandwidth
# override must survive the reset.
_reset_estimator_state(
    estimator::Union{ConventionalPLLAndDLL,NCOReferencedPLLAndDLL},
    sat::TrackedSat,
) = reset_estimator_state(
    estimator,
    sat.doppler_estimator_state,
    sat.carrier_doppler,
    sat.code_doppler,
)

# The per-satellite state built from a live `TrackedSat` and explicit filters —
# the form Tracking's own tests and a hand-built satellite use.
function TrackingLoops.SatConventionalPLLAndDLL(
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

# The words a software correlator's replicas ran on within a chunk: the
# satellite's own Dopplers, regenerated every chunk.
@inline _software_words(sat::TrackedSat) = FixedNCOWord(
    ustrip(Hz, uconvert(Hz, sat.carrier_doppler)),
    ustrip(Hz, uconvert(Hz, sat.code_doppler)),
)

# Per-sat update for the per-record estimators. Pure: takes a TrackedSat and
# returns the updated TrackedSat. Shared by the immutable
# `estimate_dopplers_and_filter_prompt` and the in-place
# `estimate_dopplers_and_filter_prompt!` so the two cannot drift.
#
# `noise` is one `(density, ready)` pair per signal, in `sat.signals` order (see
# `_signal_noise_densities`) — so the driver takes `first(noise)` and the
# passengers `Base.tail(noise)`, and no signal is ever handed another's floor.
#
# `words` and `landing_sample` are what a hardware correlator's link adds: the
# NCO words each record really ran on and the device sample the command
# computed from this fold lands at. Through the software receiver they are the
# satellite's own Doppler and `NO_LANDING_SAMPLE`.
function _update_tracked_sat_doppler(
    sat::TrackedSat,
    estimator::AbstractDopplerEstimator,
    sampling_frequency,
    noise::Tuple,
    words = _software_words(sat),
    landing_sample::Int64 = NO_LANDING_SAMPLE,
)
    # Walk all signals. For each one whose integration completed this
    # iteration, normalize/filter its prompt, advance CN0 and bit buffer, and
    # move its correlator to `last_fully_integrated_*`. Additionally, for
    # `signals[1]` (the estimator-driver signal), run PLL/DLL and update the
    # sat-shared carrier/code Doppler.
    pll_and_dll_state = sat.doppler_estimator_state
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    # The loops lock the driver (`signals[1]`) onto the real axis; every signal's
    # bit-buffer prompt is de-rotated by its carrier-phase offset from the driver
    # so a quadrature component (QPSK data/pilot, e.g. GPS L5 / Galileo E5a) does
    # not decode off the collapsed real part.
    driver_carrier_phase = get_carrier_phase_offset(head.signal)

    driver_noise_density, driver_noise_density_ready = first(noise)
    new_head, new_doppler_estimator_state, new_carrier_doppler, new_code_doppler =
        _process_estimator_driver_signal(
            head,
            sat,
            estimator,
            pll_and_dll_state,
            sampling_frequency,
            driver_noise_density,
            driver_noise_density_ready,
            driver_carrier_phase,
            words,
            landing_sample,
        )

    new_tail = _process_passenger_signals(
        tail_signals,
        sat.prn,
        sampling_frequency,
        Base.tail(noise),
        driver_carrier_phase,
    )

    # Phase-snap fallback chain. Picks the synced signal with the longest
    # `(primary × secondary)` code length, and uses its secondary-code phase to
    # anchor `sat.code_phase` to the right secondary-chip window.
    #
    # This is a *one-time* anchoring applied only on the iteration a signal
    # transitions `found == false → true`. It preserves the within-primary-block
    # phase (`mod(code_phase, primary)`) so the loop keeps the current
    # chunk-bounded position; re-running it on later iterations would wedge the
    # satellite (see issue #117). After sync, `update`'s
    # `mod(…, current_code_wrap)` maintains the alignment.
    #
    # Because sync is detected in this estimate pass — *after* the whole chunk
    # was correlated — any in-flight partial integration for this chunk was
    # accumulated at the pre-snap code phase (and, pre-sync, with no
    # secondary-code overlay). So on the snap every signal's in-flight
    # accumulator is reset too: the next chunk re-integrates cleanly from the
    # snapped phase to the next boundary.
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

# The estimator is implied by the state for the shipped estimators — the step
# reads only the state's filters and bandwidths — so a caller holding a bare
# satellite (a test, an external producer) may omit it. The default-constructed
# estimator stands in for whichever one built the state: its filter types and
# bandwidths are not read, and `predict_landing` only matters with a landing
# sample, which this path never has.
_update_tracked_sat_doppler(
    sat::TrackedSat{<:Tuple{Vararg{TrackedSignal}},<:SatConventionalPLLAndDLL},
    sampling_frequency,
    noise::Tuple,
) = _update_tracked_sat_doppler(
    sat,
    ConventionalAssistedPLLAndDLL(),
    sampling_frequency,
    noise,
)
_update_tracked_sat_doppler(
    sat::TrackedSat{<:Tuple{Vararg{TrackedSignal}},<:SatNCOReferencedPLLAndDLL},
    sampling_frequency,
    noise::Tuple,
) = _update_tracked_sat_doppler(sat, NCOReferencedPLLAndDLL(), sampling_frequency, noise)

# Drop an in-flight (partial) integration: zero the accumulator and its sample
# counter, leaving all other per-signal state intact. Used at the sync-transition
# phase snap, where the shared `code_phase` moves and any partial accumulated at
# the old phase must not be carried into the re-anchored window.
@inline _reset_inflight_integration(s::TrackedSignal) =
    TrackedSignal(s; correlator = zero(s.correlator), integrated_samples = 0)

# Apply one completed `CorrelatorOutput` record to a signal — shared by the
# estimator-driver and passenger folds so they cannot drift (issue #133). The
# arithmetic is TrackingLoops' `fold_record`, run on the `TrackedSignal`'s bare
# per-record state; this wraps it: records the filtered prompt in the signal's
# per-chunk `filtered_prompts`, and rebuilds the `TrackedSignal` with the record
# moved to `last_fully_integrated_*`. Returns the rebuilt signal plus the
# intermediate values the driver's loop-filter section needs.
#
# Unlike the old per-integration advance, this does NOT reset the live
# accumulator or `integrated_samples`: the correlate phase already reset them
# when it snapshotted this record and began (or is carrying) the next
# integration. It consumes only the record's stored correlator.
@inline function _apply_correlator_output(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase::Real = 0.0;
    correlated_pre_sync::Bool = false,
)
    bit_buffer,
    cn0_estimator,
    post_corr_filter,
    prompt,
    filtered_correlator,
    bit_block_count,
    integrated_code_blocks,
    overshoot = fold_record(
        tracked_signal.signal,
        prn,
        tracked_signal.bit_buffer,
        get_cn0_estimator(tracked_signal),
        tracked_signal.post_corr_filter,
        output,
        sampling_frequency,
        noise_density,
        noise_density_ready,
        driver_carrier_phase,
        correlated_pre_sync,
    )
    overshoot && _warn_bit_boundary_overshoot(
        get_signal_id(tracked_signal.signal),
        prn,
        tracked_signal.bit_buffer.prompt_accumulator_integrated_code_blocks +
        bit_block_count,
        _calc_num_code_blocks_that_form_a_bit(tracked_signal.signal),
    )
    push!(tracked_signal.filtered_prompts, prompt)
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
    return new_signal, filtered_correlator, integrated_code_blocks
end

# Process the estimator-driver signal (signals[1]): fold over every
# `CorrelatorOutput` collected during this chunk, in order — running the
# per-record advance and the estimator's `step` per record, threading the
# loop-filter state and the FLL `previous_prompt` across them — then return the
# new doppler_estimator_state and the *last* record's carrier/code Doppler (the
# NCO is written once per chunk). With no outputs the Doppler holds.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    estimator::Union{ConventionalPLLAndDLL,NCOReferencedPLLAndDLL},
    pll_and_dll_state::Union{SatConventionalPLLAndDLL,SatNCOReferencedPLLAndDLL},
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase::Real,
    words,
    landing_sample::Int64,
)
    outputs = tracked_signal.correlator_outputs
    if isempty(outputs)
        return tracked_signal, pll_and_dll_state, sat.carrier_doppler, sat.code_doppler
    end
    signal = tracked_signal.signal
    ts = tracked_signal
    state = pll_and_dll_state
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    found_before_fold = has_bit_or_secondary_code_been_found(ts.bit_buffer)
    # The command this fold produces is computed after its last record; a
    # delay-aware estimator maps every record of the fold onto the delay-free
    # loop's record the same distance ahead.
    fold_end = last(outputs).sample_index
    @inbounds for k in eachindex(outputs)
        output = outputs[k]
        # FLL needs the previous record's filtered prompt; the first record of
        # the chunk chains from the sat's carried-over
        # `last_fully_integrated_filtered_prompt` (the previous chunk's last).
        # Read it off `ts` BEFORE the advance overwrites it.
        previous_prompt = get_last_fully_integrated_filtered_prompt(ts)
        # A record that follows a sync detected earlier in THIS fold was
        # correlated with pre-sync replicas — its blocks still count towards the
        # bit, only its prompt may have to be dropped (see `fold_record`).
        synced_earlier_in_fold =
            !found_before_fold && has_bit_or_secondary_code_been_found(ts.bit_buffer)
        ts, filtered_correlator, integrated_code_blocks = _apply_correlator_output(
            ts,
            output,
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase;
            correlated_pre_sync = synced_earlier_in_fold,
        )
        record = LoopRecord(
            signal,
            filtered_correlator,
            previous_prompt,
            output,
            integrated_code_blocks,
            sampling_frequency;
            fold_end,
        )
        state, carrier_doppler, code_doppler =
            step_loop(estimator, state, record, words, landing_sample)
    end
    empty!(outputs)
    return ts, state, carrier_doppler, code_doppler
end

# Process the non-driver signals (signals[2:end]): the shared per-signal
# advance only — no loop-filter work. Walks the tuple recursively to keep
# type-stability and avoid boxing, stepping the per-signal `(density, ready)`
# tuple in lockstep so each passenger divides by its own noise floor.
@inline _process_passenger_signals(::Tuple{}, ::Integer, _, ::Tuple{}, ::Real) = ()
@inline function _process_passenger_signals(
    signals::Tuple,
    prn::Integer,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase::Real,
)
    noise_density, noise_density_ready = first(noise)
    new_head = _process_one_passenger_signal(
        first(signals),
        prn,
        sampling_frequency,
        noise_density,
        noise_density_ready,
        driver_carrier_phase,
    )
    (
        new_head,
        _process_passenger_signals(
            Base.tail(signals),
            prn,
            sampling_frequency,
            Base.tail(noise),
            driver_carrier_phase,
        )...,
    )
end

@inline function _process_one_passenger_signal(
    tracked_signal::TrackedSignal,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
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
                noise_density,
                noise_density_ready,
                driver_carrier_phase;
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
the end of the code or multiples of that, with the conventional PLL and DLL or the
delay-aware NCO-referenced loop (through the software receiver the latter runs as the
conventional loop). Those Doppler estimations will be used to create the next
replicas to downconvert and decode the incoming signal. In addition to the
Doppler estimation it will also filter the prompt with the configured
post correlation filter.

The sampling-frequency argument may be either a [`BandMeasurements`](@ref)
NamedTuple (the `track` path, from which the per-band rate is read) or a bare
per-band sampling-frequency source — a `NamedTuple`/`Dict` keyed by
`get_band_id` mapping each band to its sampling frequency. The latter is
the entry point for an **external correlator producer** (e.g. an FPGA): it needs
no sample buffer, only the per-signal `correlator_outputs` and the rate that
maps `integrated_samples` to an integration time. See
[External correlator producers](@ref).
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{
        <:SignalGroups,
        <:Union{ConventionalPLLAndDLL,NCOReferencedPLLAndDLL},
    },
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    # Detach the slot *values* from the input (sharing the key set), then
    # delegate to the in-place form.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, sampling_frequencies)
end

# Per-band sampling frequency for a group, from either a `BandMeasurements`
# NamedTuple (read the rate off the band's `BandMeasurement`) or a bare
# per-band rate source keyed by `get_band_id` (a `NamedTuple`/`Dict`).
@inline _band_sampling_frequency(m::BandMeasurements, key) = m[key].sampling_frequency
@inline _band_sampling_frequency(fs::NamedTuple, key) = fs[key]
@inline _band_sampling_frequency(fs::AbstractDict, key) = fs[key]

# Per-signal noise density, as `(density, ready)` — the union-free pair the fold
# threads down to `_apply_correlator_output`. Looked up by signal id, because the
# floor a record divides by is the *post-correlation* one and that is a property
# of the despreading modulation, not of the RF band.
#
# Three cases: no entry for the signal (no estimator configured — a static
# property; density `nothing`, `ready` true so the wiring mistake surfaces at
# the first record); an entry whose window is still empty (a dummy scalar and
# `ready` false); an entry with a density.
@inline _signal_noise_density(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? _noise_density_and_ready(noise_estimators[K]) :
    (nothing, true)

# One `(density, ready)` pair per signal of a group's slot type, in the slot
# type's own order — which is the order of every satellite's `signals` tuple.
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

# Per-group body for the doppler estimator. Pulled out so `_foreach_group!`
# can call it without boxing when the groups tuple is heterogeneous.
@inline function _est_one_group!(
    g::SignalGroup,
    estimator::AbstractDopplerEstimator,
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
    noise_estimators::NamedTuple,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    sampling_frequency = _band_sampling_frequency(sampling_frequencies, get_band_id(g.band))
    noise = _signal_noise_densities(noise_estimators, eltype(g.satellites))
    _warn_noise_density_missing(eltype(g.satellites), noise, noise_estimators)
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_doppler(vals[i], estimator, sampling_frequency, noise)
    end
    return nothing
end

# A signal with a configured noise estimator that has no *usable* density is a
# loud *symptom* — every satellite reports `-Inf dB-Hz` on it — but not a loud
# diagnosis. There are two causes, and the message names both: the window is
# still empty (an `append_noise_observation!` that is never called), or the
# window holds a floor of zero (a dead input). Warn rather than throw: "no usable
# density at this instant" has legitimate transient readings.
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

# `haskey` on a NamedTuple is compile-time, so the no-estimator case (a static
# wiring mistake, which must warn) folds away rather than being tested per fold.
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

As with the immutable form, the second argument is either a
[`BandMeasurements`](@ref) NamedTuple or a bare per-band sampling-frequency
source keyed by `get_band_id`. The estimator only reads the per-band
sampling frequency (to turn each output's `integrated_samples` into an
integration time and to normalize the DLL discriminator); it consumes each
signal's `correlator_outputs` and **clears** them, whether they were produced
by the software correlate phase or appended by an external producer via
[`append_correlator_output!`](@ref).
"""
function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{
        <:SignalGroups,
        <:Union{ConventionalPLLAndDLL,NCOReferencedPLLAndDLL},
    },
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    _foreach_group!(
        _est_one_group!,
        track_state.groups,
        track_state.doppler_estimator,
        sampling_frequencies,
        track_state.noise_estimators,
    )
    return track_state
end

# `maxlog` is keyed per callsite, so the `_id` is made signal- and PRN-specific:
# a record-sizing bug in an external producer hits every satellite it feeds, and
# silencing all but the first would hide how wide the problem is. Warn rather
# than throw — one producer's off-by-one should cost a bit-sync, not the
# receiver — but warn loudly, because the resync also discards the partial bit
# and the straddling record. The detection itself is TrackingLoops' (`fold_record`
# reports the overshoot); the logging is this package's, so the loop process
# stays free of it.
@noinline function _warn_bit_boundary_overshoot(
    signal_id::Symbol,
    prn::Integer,
    accumulated_blocks::Integer,
    blocks_per_bit::Integer,
)
    @warn """
          Signal `:$signal_id` PRN $prn: a correlator record carried the bit \
          accumulator past the navigation-bit boundary ($accumulated_blocks of \
          $blocks_per_bit code blocks), which only a record not aligned to that \
          boundary can do. Dropping bit sync and re-running the detector; the \
          partial bit is discarded, the bits decoded so far are kept. Check the \
          record sizing of the external `CorrelatorOutput` producer — it must \
          not straddle a bit boundary.""" _id =
        Symbol(:bit_boundary_overshoot_, signal_id, :_, prn) maxlog = 1
    nothing
end
