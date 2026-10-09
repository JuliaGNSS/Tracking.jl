# The passengers, the host side: the driver fold walks the passengers'
# (`signals[2:end]`) records alongside its own and hands each to the estimator's
# `step_loop` before the driver record it ends within, so an estimator that
# combines signals combines each into the driver step it belongs to, and the
# vector loop decodes its decoding signal's. Nothing here asks which estimator it
# is. The combining itself — weights, gating, which loops — is TrackingLoops'. See the
# "Signal combining" section of docs/src/loop_filter.md.

# Per-passenger constants of one fold. `differential_group_delay_chips`, the
# passenger's group delay minus the driver's in chips at the driver's code
# frequency, refers the passenger's DLL discriminator to the driver's code phase;
# `NaN` where either group delay is unknown, which leaves the passenger out of
# the code loop: unknown is not zero.
@inline function _passenger_context(
    passenger::TrackedSignal,
    (noise_density, noise_density_ready),
    driver::TrackedSignal,
    code_frequency,
)
    known = !isnothing(passenger.group_delay) && !isnothing(driver.group_delay)
    differential_group_delay_chips =
        known ?
        Float64(
            ustrip(NoUnits, (passenger.group_delay - driver.group_delay) * code_frequency),
        ) : NaN
    (;
        found_before_fold = has_bit_or_secondary_code_been_found(passenger.bit_buffer),
        noise_density,
        noise_density_ready,
        differential_group_delay_chips,
        # The fold the records belong to ends at the last of them.
        fold_end = isempty(passenger.correlator_outputs) ? 0 :
                   last(passenger.correlator_outputs).sample_index,
    )
end

# Apply each passenger's records that end by sample `until`, from its cursor on,
# and combine them into the estimator state. Recursive over the tuple for type
# stability.
@inline _advance_passengers(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Integer,
    state,
    ::Vararg{Any,N},
) where {N} = ((), (), state)
@inline function _advance_passengers(
    passengers::Tuple,
    cursors::Tuple,
    contexts::Tuple,
    until::Integer,
    state,
    args::Vararg{Any,N},
) where {N}
    passenger, cursor, state = _advance_passenger(
        first(passengers),
        first(cursors),
        first(contexts),
        until,
        state,
        args...,
    )
    rest, rest_cursors, state = _advance_passengers(
        Base.tail(passengers),
        Base.tail(cursors),
        Base.tail(contexts),
        until,
        state,
        args...,
    )
    (passenger, rest...), (cursor, rest_cursors...), state
end

# One passenger's records up to `until`: the shared per-record advance, then the
# record — with the passenger's own previous prompt, on its own record sequence —
# through the estimator's `step_loop`. Its Dopplers are the command in force and
# are not taken: the chunk's command is the driver's.
@inline function _advance_passenger(
    tracked_signal::TrackedSignal,
    cursor::Int,
    context,
    until::Integer,
    state,
    estimator::AbstractDopplerEstimator,
    driver_signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency,
    driver_carrier_phase::Real,
    words,
    landing_sample::Int64,
    sample_offset::Int,
)
    ts = tracked_signal
    outputs = ts.correlator_outputs
    @inbounds while cursor <= length(outputs) && outputs[cursor].sample_index <= until
        output = outputs[cursor]
        previous_prompt = _fll_previous_prompt(ts, output, sampling_frequency)
        # See `correlated_pre_sync` in `_apply_correlator_output`.
        correlated_pre_sync =
            !context.found_before_fold &&
            has_bit_or_secondary_code_been_found(ts.bit_buffer)
        ts, filtered_correlator, integrated_code_blocks, signal_state = _apply_correlator_output(
            ts,
            output,
            prn,
            sampling_frequency,
            context.noise_density,
            context.noise_density_ready,
            driver_carrier_phase;
            correlated_pre_sync,
        )
        record = LoopRecord(
            ts.signal,
            filtered_correlator,
            previous_prompt,
            output,
            integrated_code_blocks,
            sampling_frequency;
            context.fold_end,
            prn,
            sample_offset,
            signal_state,
            context.differential_group_delay_chips,
        )
        state, _, _ = step_loop(estimator, state, record, words, landing_sample)
        cursor += 1
    end
    ts, cursor, state
end
