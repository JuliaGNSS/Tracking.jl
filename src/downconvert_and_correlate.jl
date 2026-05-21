# Per-sat update applied after a single (downconvert + correlate) step over
# `integrated_samples` samples. Advances the shared carrier and code phase,
# each per-signal correlator accumulator, integration counter, and the
# `is_integration_completed` flag (per signal). `signal_start_sample`
# advances on the satellite. Returns the new `TrackedSat`.
#
# `new_signals_data` is a tuple of `(new_correlator, is_integration_completed)`
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

    new_signals = _build_new_signals(sat.signals, new_signals_data, integrated_samples)

    TrackedSat(
        sat;
        code_phase,
        carrier_phase,
        signal_start_sample = sat.signal_start_sample + integrated_samples,
        signals = new_signals,
    )
end

# Recursive tuple walk that rebuilds each TrackedSignal with its
# corresponding (new_correlator, is_integration_completed) from
# `new_signals_data`. Allocation-free and inference-friendly.
@inline _build_new_signals(::Tuple{}, ::Tuple{}, ::Int) = ()
@inline function _build_new_signals(
    signals::Tuple,
    new_data::Tuple,
    integrated_samples::Int,
)
    s = first(signals)
    (corr, completed) = first(new_data)
    new_s = TrackedSignal(
        s;
        integrated_samples = s.integrated_samples + integrated_samples,
        is_integration_completed = completed,
        correlator = corr,
    )
    (new_s, _build_new_signals(Base.tail(signals), Base.tail(new_data), integrated_samples)...)
end
