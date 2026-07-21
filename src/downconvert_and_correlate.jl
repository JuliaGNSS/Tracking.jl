# Per-sat update applied after a single (downconvert + correlate) sub-step over
# `integrated_samples` samples. Advances the shared carrier and code phase and
# each per-signal correlator accumulator. `signal_start_sample` advances on the
# satellite. Returns the new `TrackedSat`.
#
# When a signal's coherent integration completes on this sub-step, its raw
# accumulator is snapshotted into a `CorrelatorOutput` (appended to the reused
# per-signal `correlator_outputs` buffer) and the accumulator is reset — so the
# next integration within the same processing chunk starts fresh. The Doppler
# estimator later folds over `correlator_outputs`; see
# `estimate_dopplers_and_filter_prompt!`.
#
# `new_signals_data` is a tuple of `(new_correlator, completed)`
# pairs, one per element of `sat.signals`. Built via tuple recursion in
# `_update_tracked_sat_correlator` so the heterogeneous walk stays
# allocation-free and inferable.
function update(
    sat::TrackedSat,
    integrated_samples::Int,
    intermediate_frequency,
    sampling_frequency,
    new_signals_data::Tuple,
)
    # Carrier and code phase wrap shared across all signals on this sat.
    # The chip rate is shared across signals within a single band — see the
    # design doc; we drive the code-phase advance from signals[1].
    driver = first(sat.signals).signal
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    code_frequency = sat.code_doppler + get_code_frequency(driver)
    carrier_phase = update_carrier_phase(
        integrated_samples,
        carrier_frequency,
        sampling_frequency,
        sat.carrier_phase,
    )
    # Runtime wrap honors per-signal sync state: an L1-C/A-only sat
    # wraps at 1023 before bit-edge sync and at 20460 (= 20 × 1023)
    # after, so downstream consumers (e.g. PositionVelocityTime.jl) can
    # distinguish which primary-code-period of the 20-block bit they're
    # in. Pilots widen to `primary × secondary_code_length`. See
    # [`current_code_wrap`](@ref).
    code_length = current_code_wrap(sat.signals)
    code_phase = mod(
        code_frequency * integrated_samples / sampling_frequency + sat.code_phase,
        code_length,
    )
    new_signal_start_sample = sat.signal_start_sample + integrated_samples

    new_signals = _build_new_signals(
        sat.signals,
        new_signals_data,
        integrated_samples,
        code_phase,
        new_signal_start_sample,
    )

    TrackedSat(
        sat;
        code_phase,
        carrier_phase,
        signal_start_sample = new_signal_start_sample,
        signals = new_signals,
    )
end

# Recursive tuple walk that rebuilds each TrackedSignal from its
# `(new_correlator, completed)` in `new_signals_data`.
# Allocation-free and inference-friendly.
#
# On completion, snapshot the (raw) accumulator into the shared
# `correlator_outputs` vector — tagged with `sample_index` (the end sample of
# this integration) and the boundary `code_phase` — and reset the accumulator.
# The `push!` mutates the same vector the copy-update constructor threads
# through unchanged, so it stays allocation-free after the buffer's capacity is
# seated. On a partial (chunk/buffer-bounded) sub-step, carry the accumulator.
@inline _build_new_signals(::Tuple{}, ::Tuple{}, ::Int, _, ::Int) = ()
@inline function _build_new_signals(
    signals::Tuple,
    new_data::Tuple,
    integrated_samples::Int,
    code_phase,
    signal_start_sample::Int,
)
    s = first(signals)
    (corr, completed) = first(new_data)
    total_integrated = s.integrated_samples + integrated_samples
    if completed
        push!(
            s.correlator_outputs,
            CorrelatorOutput(corr, total_integrated, signal_start_sample - 1, code_phase),
        )
        new_s = TrackedSignal(s; integrated_samples = 0, correlator = zero(corr))
    else
        new_s = TrackedSignal(s; integrated_samples = total_integrated, correlator = corr)
    end
    (
        new_s,
        _build_new_signals(
            Base.tail(signals),
            Base.tail(new_data),
            integrated_samples,
            code_phase,
            signal_start_sample,
        )...,
    )
end
