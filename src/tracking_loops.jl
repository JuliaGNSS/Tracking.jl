# The TrackState-level walk of the tracking loops over groups, satellites and
# signals. The tracking loops themselves are TrackingLoops'; any
# `AbstractTrackingLoop` implementing `init_loop_state`, `reseed_loop_state` and
# `step_loop!` is driven by `track!` as it is.

# Re-seed hook for `reset_loop_filters!`: the tracking loop's `reseed_loop_state`
# (which keeps what it chooses, e.g. bandwidths and carrier-loop stage), else a
# fresh `init_loop_state`.
function _reseed_loop_state(tracking_loop::AbstractTrackingLoop, sat::TrackedSat)
    loop_state = sat.loop_state
    applicable(
        reseed_loop_state,
        tracking_loop,
        loop_state,
        sat.carrier_doppler,
        sat.code_doppler,
    ) || return init_loop_state(tracking_loop, sat)
    reseed_loop_state(tracking_loop, loop_state, sat.carrier_doppler, sat.code_doppler)
end

# The NCO command a software correlator ran on: the satellite's own Dopplers.
@inline _software_nco_command(sat::TrackedSat) = ConstantNCOCommand(
    ustrip(Hz, uconvert(Hz, sat.carrier_doppler)),
    ustrip(Hz, uconvert(Hz, sat.code_doppler)),
)

# Pure per-satellite update: steps every completed predetection integration of
# every signal through the tracking loop and returns the updated `TrackedSat`.
#
# `noise` is one `(density, ready)` pair per signal, in `sat.signals` order (see
# `_signal_noise_densities`). Each integration ran on the satellite's own Doppler
# and its command takes effect at the integration's end (application sample
# `nothing`): a hardware correlator's loop, which needs another NCO history and an
# application sample, calls TrackingLoops' `step_loop!` itself.
# `sample_index_offset` moves the integrations onto the band's running sample grid
# (`TrackState`'s `sample_index_offsets`).
function _step_sat_tracking_loop(
    sat::TrackedSat,
    tracking_loop::AbstractTrackingLoop,
    sampling_frequency,
    noise::Tuple,
    sample_index_offset::Int = 0,
)
    # Driver (`signals[1]`) first, then the passengers, all through the one
    # loop state; the satellite's Doppler is the driver's command.
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    # The loops lock the driver onto the real axis; each signal's bit-synchronizer
    # prompt is de-rotated by its phase offset from it, so a quadrature component
    # (e.g. GPS L5, Galileo E5a) decodes off the right axis.
    driver_carrier_phase_offset = get_carrier_phase_offset(head.signal)
    nco_command = _software_nco_command(sat)

    driver_noise_density, driver_noise_density_ready = first(noise)
    new_head, new_loop_state, new_carrier_doppler, new_code_doppler =
        _process_driver_signal(
            head,
            sat,
            tracking_loop,
            sat.loop_state,
            sampling_frequency,
            driver_noise_density,
            driver_noise_density_ready,
            driver_carrier_phase_offset,
            nco_command,
            sample_index_offset,
        )

    new_tail, new_loop_state = _process_passenger_signals(
        tail_signals,
        sat.prn,
        sampling_frequency,
        Base.tail(noise),
        driver_carrier_phase_offset,
        tracking_loop,
        new_loop_state,
        nco_command,
        sample_index_offset,
    )

    # On the batch a signal first syncs, anchor `sat.code_phase` to the secondary
    # chip of the synced signal with the longest full code, keeping the
    # within-primary phase. Only then (#117); afterwards the replica's wrap keeps it.
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
        loop_state = new_loop_state,
    )
end

# Drop an in-flight partial integration (at the sync phase snap).
@inline _reset_inflight_integration(tsig::TrackedSignal) = TrackedSignal(
    tsig;
    correlator_bank = zero(tsig.correlator_bank),
    num_integrated_samples = 0,
)

# Process one completed correlator output of a signal into its predetection
# integration, shared by the driver and the passengers (#133): TrackingLoops'
# `process_correlator_bank_output!` on the signal's `post_correlation_state`, plus
# pushing the beamformed prompt and keeping the output as
# `last_integration_correlator_bank`. The live correlator bank is left alone; the
# correlate phase already reset it.
#
# Without `secondary_code_wiped` the output came from the pre-sync replica, which
# GNSSSignals overlays with secondary chip 0 on every code period. TrackingLoops
# reads a dataless, secondary-coded signal's pre-sync prompts as primary-only, so
# such an output is multiplied by chip 0 first; otherwise a PRN with a negative
# chip 0 gives the four-quadrant PLL a half-cycle error.
@inline function _process_correlator_bank_output(
    tracked_signal::TrackedSignal,
    output::CorrelatorBankOutput,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase_offset::Real;
    secondary_code_wiped::Bool,
    batch_end_sample_index::Integer,
    sample_index_offset::Integer,
)
    if !secondary_code_wiped && _pre_sync_replica_chip(tracked_signal.signal, prn) < 0
        output = _negated(output)
    end
    num_code_periods_into_bit_before =
        get_num_code_periods_into_bit(tracked_signal.post_correlation_state)
    post_correlation_state, integration = process_correlator_bank_output!(
        tracked_signal.post_correlation_state,
        tracked_signal.signal,
        prn,
        output,
        sampling_frequency,
        noise_density,
        noise_density_ready,
        driver_carrier_phase_offset;
        secondary_code_wiped,
        batch_end_sample_index,
        sample_index_offset,
    )
    integration.sync_change == TrackingLoops.SYNC_LOST && _warn_bit_boundary_overshoot(
        get_signal_id(tracked_signal.signal),
        prn,
        num_code_periods_into_bit_before + integration.num_code_periods,
        get_code_periods_per_bit(tracked_signal.signal),
    )
    push!(
        tracked_signal.beamformed_prompts,
        get_prompt(integration.beamformed_correlator_bank),
    )
    new_signal = TrackedSignal(
        tracked_signal;
        post_correlation_state,
        last_integration_correlator_bank = output.correlator_bank,
    )
    return new_signal, integration
end

# Secondary chip 0 for a dataless, secondary-coded signal (TrackingLoops'
# predicate), else `1`; see `_process_correlator_bank_output`.
@inline _pre_sync_replica_chip(signal::AbstractGNSSSignal, prn::Integer) =
    is_dataless_with_secondary_code(signal) ?
    Int(GNSSSignals.secondary_value(get_secondary_code(signal), prn, 0)) : 1

@inline _negated(output::CorrelatorBankOutput) = CorrelatorBankOutput(
    map_correlator_outputs(-, output.correlator_bank),
    output.num_integrated_samples,
    output.sample_index,
    output.code_phase,
)

# Step the driver's predetection integrations of this chunk through the tracking
# loop; return the new loop state and the last integration's Dopplers (the NCO is
# written once per chunk). With no integrations the Doppler holds.
@inline function _process_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    tracking_loop::AbstractTrackingLoop,
    loop_state,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase_offset::Real,
    nco_command,
    sample_index_offset::Int = 0,
)
    outputs = tracked_signal.correlator_bank_outputs
    if isempty(outputs)
        return tracked_signal, loop_state, sat.carrier_doppler, sat.code_doppler
    end
    tsig = tracked_signal
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    synced_before_batch = is_bit_or_secondary_code_synced(tsig)
    # The batch's command is computed after its last integration; a
    # delay-compensated loop maps each integration onto the delay-free loop's one
    # as far ahead.
    batch_end_sample_index = last(outputs).sample_index
    @inbounds for k in eachindex(outputs)
        # A batch that started unsynced was correlated with the pre-sync replica
        # throughout, even past a sync found in it (see `process_correlator_bank_output!`).
        tsig, integration = _process_correlator_bank_output(
            tsig,
            outputs[k],
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase_offset;
            secondary_code_wiped = synced_before_batch,
            batch_end_sample_index,
            sample_index_offset,
        )
        loop_state, carrier_doppler, code_doppler =
            step_loop!(tracking_loop, loop_state, integration, nco_command, nothing)
    end
    empty!(outputs)
    return tsig, loop_state, carrier_doppler, code_doppler
end

# The passengers (`signals[2:end]`): their integrations are stepped through the
# tracking loop for its state (the vector loop decodes bits from them) and the
# returned command is dropped. Recursive over the tuple for type stability, with
# `noise` in lockstep.
@inline _process_passenger_signals(
    ::Tuple{},
    ::Integer,
    _,
    ::Tuple{},
    ::Real,
    ::AbstractTrackingLoop,
    loop_state,
    nco_command,
    ::Int,
) = (), loop_state
@inline function _process_passenger_signals(
    signals::Tuple,
    prn::Integer,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
    tracking_loop::AbstractTrackingLoop,
    loop_state,
    nco_command,
    sample_index_offset::Int,
)
    noise_density, noise_density_ready = first(noise)
    new_head, loop_state = _process_one_passenger_signal(
        first(signals),
        prn,
        sampling_frequency,
        noise_density,
        noise_density_ready,
        driver_carrier_phase_offset,
        tracking_loop,
        loop_state,
        nco_command,
        sample_index_offset,
    )
    new_tail, loop_state = _process_passenger_signals(
        Base.tail(signals),
        prn,
        sampling_frequency,
        Base.tail(noise),
        driver_carrier_phase_offset,
        tracking_loop,
        loop_state,
        nco_command,
        sample_index_offset,
    )
    (new_head, new_tail...), loop_state
end

@inline function _process_one_passenger_signal(
    tracked_signal::TrackedSignal,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    driver_carrier_phase_offset::Real,
    tracking_loop::AbstractTrackingLoop,
    loop_state,
    nco_command,
    sample_index_offset::Int,
)
    outputs = tracked_signal.correlator_bank_outputs
    isempty(outputs) && return tracked_signal, loop_state
    tsig = tracked_signal
    synced_before_batch = is_bit_or_secondary_code_synced(tsig)
    batch_end_sample_index = last(outputs).sample_index
    @inbounds for k in eachindex(outputs)
        # Same rule as the driver's.
        tsig, integration = _process_correlator_bank_output(
            tsig,
            outputs[k],
            prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase_offset;
            secondary_code_wiped = synced_before_batch,
            batch_end_sample_index,
            sample_index_offset,
        )
        loop_state, _, _ =
            step_loop!(tracking_loop, loop_state, integration, nco_command, nothing)
    end
    empty!(outputs)
    tsig, loop_state
end

"""
$(SIGNATURES)

Step the tracking loop with every completed predetection integration of every
satellite, through TrackingLoops' `process_correlator_bank_output!` and the tracking
loop's `step_loop!` (any `AbstractTrackingLoop` implementing it works), giving the
Dopplers for the next replicas and the beamformed prompts.

The second argument is a [`BandMeasurements`](@ref Tracking.BandMeasurements) NamedTuple (the `track` path) or a
`NamedTuple`/`Dict` mapping each `get_band_id` to its sampling frequency: the entry
point for an **external correlator producer** (e.g. an FPGA), which needs no sample
buffer. See [External correlator producers](@ref).
"""
function step_tracking_loops(
    track_state::TrackState{<:Any,<:AbstractTrackingLoop},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    # Detach the slot values (sharing the key set), then run in place.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    step_tracking_loops!(new_track_state, sampling_frequencies)
end

@inline _band_sampling_frequency(m::BandMeasurements, key) = m[key].sampling_frequency
@inline _band_sampling_frequency(fs::NamedTuple, key) = fs[key]
@inline _band_sampling_frequency(fs::AbstractDict, key) = fs[key]

# Per-signal noise density as a union-free `(density, ready)` pair, keyed by signal
# id (see the noise estimator manual page for why per signal). No estimator gives
# `(nothing, true)`, so the wiring mistake surfaces at the first integration.
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
@inline function _step_group_tracking_loops!(
    g::SignalGroup,
    tracking_loop::AbstractTrackingLoop,
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
    noise_estimators::NamedTuple,
    sample_index_offsets::Dict{Symbol,Int},
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    band_id = get_band_id(g.band)
    sampling_frequency = _band_sampling_frequency(sampling_frequencies, band_id)
    sample_index_offset = get(sample_index_offsets, band_id, 0)
    noise = _signal_noise_densities(noise_estimators, eltype(g.satellites))
    _warn_noise_density_missing(eltype(g.satellites), noise, noise_estimators)
    @inbounds for i in eachindex(vals)
        sat = vals[i]
        vals[i] = _step_sat_tracking_loop(
            sat,
            tracking_loop,
            sampling_frequency,
            noise,
            sample_index_offset,
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
          `NoiseReferencedCN0Estimator` will not update and C/N₀ stays at `-Inf dB-Hz`. \
          Either its window is still empty — call `append_noise_observation!` before \
          `step_tracking_loops!`, or use a noise estimator that \
          measures from the band's samples — or the measured floor is zero, which \
          means the samples carry no power at all (a dead input or an underrun).""" _id =
        Symbol(:no_noise_density_, signal_id) maxlog = 1
    nothing
end

"""
$(SIGNATURES)

In-place version of [`step_tracking_loops`](@ref), taking the same
second argument: overwrites each group's `TrackedSat` slots and returns `track_state`,
allocation-free in steady state under [`track!`](@ref)'s preconditions.

It reads only the per-band sampling frequency (turning `num_integrated_samples` into an
integration time), and consumes and **clears** each
signal's `correlator_bank_outputs`, whether from the software correlate phase or from
[`append_correlator_bank_output!`](@ref).
"""
function step_tracking_loops!(
    track_state::TrackState{<:Any,<:AbstractTrackingLoop},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    _foreach_group!(
        _step_group_tracking_loops!,
        track_state.groups,
        track_state.tracking_loop,
        sampling_frequencies,
        track_state.noise_estimators,
        track_state.sample_index_offsets,
    )
    return track_state
end

# TrackingLoops' `process_correlator_bank_output!` detects the overshoot; the logging lives here. The
# `_id` is per signal and PRN so `maxlog` shows how many satellites a producer's
# integration-sizing bug hits. Warn, not throw: it costs a bit sync, not the receiver.
@noinline function _warn_bit_boundary_overshoot(
    signal_id::Symbol,
    prn::Integer,
    accumulated_code_periods::Integer,
    code_periods_per_bit::Integer,
)
    @warn """
          Signal `:$signal_id` PRN $prn: a correlator output carried the bit \
          correlator bank past the navigation-bit boundary ($accumulated_code_periods of \
          $code_periods_per_bit code periods), which only an integration not aligned to that \
          boundary can do. Dropping bit sync and re-running the detector; the \
          partial bit is discarded, the bits decoded so far are kept. Check the \
          integration sizing of the external `CorrelatorBankOutput` producer — it must \
          not straddle a bit boundary.""" _id =
        Symbol(:bit_boundary_overshoot_, signal_id, :_, prn) maxlog = 1
    nothing
end
