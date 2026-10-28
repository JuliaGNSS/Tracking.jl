# The TrackState-level walk of the Doppler estimators over groups, satellites and
# signals. The estimators themselves are TrackingLoops'; any
# `AbstractDopplerEstimator` implementing `init_estimator_state`,
# `reset_estimator_state` and `step_loop` is driven by `track!` as it is.

# Re-seed hook for `reset_loop_filters!`: the estimator's `reset_estimator_state`
# (which keeps what it chooses, e.g. bandwidths and carrier-loop stage), else a
# fresh `init_estimator_state`.
function _reset_estimator_state(estimator::AbstractDopplerEstimator, sat::TrackedSat)
    state = sat.doppler_estimator_state
    applicable(
        reset_estimator_state,
        estimator,
        state,
        sat.carrier_doppler,
        sat.code_doppler,
    ) || return init_estimator_state(estimator, sat)
    reset_estimator_state(estimator, state, sat.carrier_doppler, sat.code_doppler)
end

# The NCO words a software correlator ran on: the satellite's own Dopplers.
@inline _software_words(sat::TrackedSat) = FixedNCOWord(
    ustrip(Hz, uconvert(Hz, sat.carrier_doppler)),
    ustrip(Hz, uconvert(Hz, sat.code_doppler)),
)

# Pure per-satellite update: folds every completed record of every signal through
# the estimator and returns the updated `TrackedSat`.
#
# `noise` is one `(density, ready)` pair per signal, in `sat.signals` order (see
# `_signal_noise_densities`). `words` and `landing_sample` come from a hardware
# correlator's link (the NCO words each record ran on and the device sample the
# command lands at); the software receiver passes the satellite's own Doppler and
# `NO_LANDING_SAMPLE`. `sample_offset` moves records onto the band's running
# sample grid (`TrackState`'s `sample_offsets`).
function _update_tracked_sat_doppler(
    sat::TrackedSat,
    estimator::AbstractDopplerEstimator,
    sampling_frequency,
    noise::Tuple,
    words = _software_words(sat),
    landing_sample::Int64 = NO_LANDING_SAMPLE,
    sample_offset::Int = 0,
)
    # Driver (`signals[1]`) first, then the passengers, all through the one
    # estimator state; the satellite's Doppler is the driver's command.
    pll_and_dll_state = sat.doppler_estimator_state
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    # The loops lock the driver onto the real axis; each signal's bit-buffer
    # prompt is de-rotated by its phase offset from it, so a quadrature component
    # (e.g. GPS L5, Galileo E5a) decodes off the right axis.
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
            sample_offset,
        )

    new_tail, new_doppler_estimator_state = _process_passenger_signals(
        tail_signals,
        sat.prn,
        sampling_frequency,
        Base.tail(noise),
        driver_carrier_phase,
        estimator,
        new_doppler_estimator_state,
        words,
        landing_sample,
        sample_offset,
    )

    # On the fold a signal first syncs, anchor `sat.code_phase` to the secondary
    # chip of the synced signal with the longest full code, keeping the
    # within-primary phase. Only then (#117); afterwards `update`'s wrap keeps it.
    # Sync is found after the chunk was correlated, so the in-flight partial
    # integrations, accumulated at the old phase, are dropped too.
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

# Drop an in-flight partial integration (at the sync phase snap).
@inline _reset_inflight_integration(s::TrackedSignal) =
    TrackedSignal(s; correlator = zero(s.correlator), integrated_samples = 0)

# Apply one completed record to a signal, shared by the driver and passenger
# folds (#133): TrackingLoops' `apply_record` on the signal's `loop_state`, plus
# pushing the filtered prompt and keeping the record as
# `last_fully_integrated_correlator`. The live accumulator is left alone; the
# correlate phase already reset it.
#
# Without `overlay_wiped` the record came from the pre-sync replica, which
# GNSSSignals overlays with secondary chip 0 on every block. TrackingLoops reads a
# dataless, secondary-coded signal's pre-sync prompts as primary-only, so such a
# record is multiplied by chip 0 first; otherwise a PRN with a negative chip 0
# gives the four-quadrant PLL a half-cycle error.
@inline function _apply_correlator_output(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase::Real;
    overlay_wiped::Bool,
)
    if !overlay_wiped && _pre_sync_replica_chip(tracked_signal.signal, prn) < 0
        output = _negated(output)
    end
    loop_state, prompt, filtered_correlator, integrated_code_blocks, overshoot =
        apply_record(
            tracked_signal.loop_state,
            tracked_signal.signal,
            prn,
            output,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase;
            overlay_wiped,
        )
    overshoot && _warn_bit_boundary_overshoot(
        get_signal_id(tracked_signal.signal),
        prn,
        get_bit_buffer(tracked_signal).prompt_accumulator_integrated_code_blocks +
        integrated_code_blocks,
        get_code_blocks_per_bit(tracked_signal.signal),
    )
    push!(tracked_signal.filtered_prompts, prompt)
    new_signal = TrackedSignal(
        tracked_signal;
        loop_state,
        last_fully_integrated_correlator = output.correlator,
    )
    return new_signal, filtered_correlator, integrated_code_blocks
end

# Secondary chip 0 for a dataless, secondary-coded signal (TrackingLoops'
# predicate), else `1`; see `_apply_correlator_output`.
@inline _pre_sync_replica_chip(signal::AbstractGNSSSignal, prn::Integer) =
    is_dataless_with_secondary_code(signal) ?
    Int(GNSSSignals.secondary_value(get_secondary_code(signal), prn, 0)) : 1

@inline _negated(output::CorrelatorOutput) = CorrelatorOutput(
    TrackingLoops.apply(-, output.correlator),
    output.integrated_samples,
    output.sample_index,
    output.code_phase,
)

# Fold the driver's records of this chunk through the estimator; return the new
# state and the last record's Dopplers (the NCO is written once per chunk). With
# no records the Doppler holds.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    estimator::AbstractDopplerEstimator,
    pll_and_dll_state,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase::Real,
    words,
    landing_sample::Int64,
    sample_offset::Int = 0,
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
    found_before_fold = has_bit_or_secondary_code_been_found(ts)
    # The fold's command is computed after its last record; a delay-aware
    # estimator maps each record onto the delay-free loop's one as far ahead.
    fold_end = last(outputs).sample_index
    @inbounds for k in eachindex(outputs)
        output = outputs[k]
        # A fold that started unsynced was correlated with the pre-sync replica
        # throughout, even past a sync found in it (see `apply_record`).
        ts, filtered_correlator, _ = _apply_correlator_output(
            ts,
            output,
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase;
            overlay_wiped = found_before_fold,
        )
        record = LoopRecord(
            signal,
            filtered_correlator,
            output,
            ts.loop_state,
            sampling_frequency;
            fold_end,
            prn = sat.prn,
            sample_offset,
        )
        state, carrier_doppler, code_doppler =
            step_loop(estimator, state, record, words, landing_sample)
    end
    empty!(outputs)
    return ts, state, carrier_doppler, code_doppler
end

# The passengers (`signals[2:end]`): records are stepped through the estimator
# for its state (the vector loop decodes bits from them) and the returned command
# is dropped. Recursive over the tuple for type stability, with `noise` in lockstep.
@inline _process_passenger_signals(
    ::Tuple{},
    ::Integer,
    _,
    ::Tuple{},
    ::Real,
    ::AbstractDopplerEstimator,
    state,
    words,
    ::Int64,
    ::Int,
) = (), state
@inline function _process_passenger_signals(
    signals::Tuple,
    prn::Integer,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase::Real,
    estimator::AbstractDopplerEstimator,
    state,
    words,
    landing_sample::Int64,
    sample_offset::Int,
)
    noise_density, noise_density_ready = first(noise)
    new_head, state = _process_one_passenger_signal(
        first(signals),
        prn,
        sampling_frequency,
        noise_density,
        noise_density_ready,
        driver_carrier_phase,
        estimator,
        state,
        words,
        landing_sample,
        sample_offset,
    )
    new_tail, state = _process_passenger_signals(
        Base.tail(signals),
        prn,
        sampling_frequency,
        Base.tail(noise),
        driver_carrier_phase,
        estimator,
        state,
        words,
        landing_sample,
        sample_offset,
    )
    (new_head, new_tail...), state
end

@inline function _process_one_passenger_signal(
    tracked_signal::TrackedSignal,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase::Real,
    estimator::AbstractDopplerEstimator,
    state,
    words,
    landing_sample::Int64,
    sample_offset::Int,
)
    outputs = tracked_signal.correlator_outputs
    isempty(outputs) && return tracked_signal, state
    ts = tracked_signal
    found_before_fold = has_bit_or_secondary_code_been_found(ts)
    fold_end = last(outputs).sample_index
    @inbounds for k in eachindex(outputs)
        output = outputs[k]
        # Same rule as the driver fold.
        ts, filtered_correlator, _ = _apply_correlator_output(
            ts,
            output,
            prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase;
            overlay_wiped = found_before_fold,
        )
        record = LoopRecord(
            ts.signal,
            filtered_correlator,
            output,
            ts.loop_state,
            sampling_frequency;
            fold_end,
            prn,
            sample_offset,
        )
        state, _, _ = step_loop(estimator, state, record, words, landing_sample)
    end
    empty!(outputs)
    ts, state
end

"""
$(SIGNATURES)

Fold every completed correlator record of every satellite through the estimator's
TrackingLoops `step_loop` (any estimator implementing it works), giving the Dopplers for
the next replicas and the post-correlation-filtered prompts.

The second argument is a [`BandMeasurements`](@ref) NamedTuple (the `track` path) or a
`NamedTuple`/`Dict` mapping each `get_band_id` to its sampling frequency: the entry
point for an **external correlator producer** (e.g. an FPGA), which needs no sample
buffer. See [External correlator producers](@ref).
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:Any,<:AbstractDopplerEstimator},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    # Detach the slot values (sharing the key set), then run in place.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, sampling_frequencies)
end

@inline _band_sampling_frequency(m::BandMeasurements, key) = m[key].sampling_frequency
@inline _band_sampling_frequency(fs::NamedTuple, key) = fs[key]
@inline _band_sampling_frequency(fs::AbstractDict, key) = fs[key]

# Per-signal noise density as a union-free `(density, ready)` pair, keyed by signal
# id (see the noise estimator manual page for why per signal). No estimator gives
# `(nothing, true)`, so the wiring mistake surfaces at the first record.
@inline _signal_noise_density(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? noise_density_and_ready(noise_estimators[K]) :
    (nothing, true)

# One `(density, ready)` pair per signal, in `signals` tuple order.
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

# Per-group body, a function so `_foreach_group!` does not box.
@inline function _est_one_group!(
    g::SignalGroup,
    estimator::AbstractDopplerEstimator,
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
    noise_estimators::NamedTuple,
    sample_offsets::Dict{Symbol,Int},
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    band_id = get_band_id(g.band)
    sampling_frequency = _band_sampling_frequency(sampling_frequencies, band_id)
    sample_offset = get(sample_offsets, band_id, 0)
    noise = _signal_noise_densities(noise_estimators, eltype(g.satellites))
    _warn_noise_density_missing(eltype(g.satellites), noise, noise_estimators)
    @inbounds for i in eachindex(vals)
        sat = vals[i]
        vals[i] = _update_tracked_sat_doppler(
            sat,
            estimator,
            sampling_frequency,
            noise,
            _software_words(sat),
            NO_LANDING_SAMPLE,
            sample_offset,
        )
    end
    return nothing
end

# A configured noise estimator with no usable density shows only as `-Inf dB-Hz`,
# so name the causes: an empty window or a zero floor. Warn, not throw: it can be
# transient.
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

# `haskey` on a NamedTuple folds at compile time.
@inline _noise_window_still_filling(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? noise_window_filling(noise_estimators[K]) : false

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

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref), taking the same
second argument: overwrites each group's `TrackedSat` slots and returns `track_state`,
allocation-free in steady state under [`track!`](@ref)'s preconditions.

It reads only the per-band sampling frequency (turning `integrated_samples` into an
integration time), and consumes and **clears** each
signal's `correlator_outputs`, whether from the software correlate phase or from
[`append_correlator_output!`](@ref).
"""
function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:Any,<:AbstractDopplerEstimator},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    _foreach_group!(
        _est_one_group!,
        track_state.groups,
        track_state.doppler_estimator,
        sampling_frequencies,
        track_state.noise_estimators,
        track_state.sample_offsets,
    )
    return track_state
end

# TrackingLoops' `apply_record` detects the overshoot; the logging lives here. The
# `_id` is per signal and PRN so `maxlog` shows how many satellites a producer's
# record-sizing bug hits. Warn, not throw: it costs a bit sync, not the receiver.
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
