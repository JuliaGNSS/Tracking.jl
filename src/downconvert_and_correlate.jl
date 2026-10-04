# Per-sat update after one downconvert+correlate sub-step over
# `integrated_samples` samples: advances the shared carrier and code phase and
# `signal_start_sample`, and rebuilds each signal from `new_signals_data`, a tuple
# of `(new_correlator, completed)` per element of `sat.signals` (see
# `_build_new_signals`). Returns the new `TrackedSat`.
function update(
    sat::TrackedSat,
    integrated_samples::Int,
    intermediate_frequency,
    sampling_frequency,
    new_signals_data::Tuple,
)
    # The chip rate is shared by all signals of a sat, so signals[1] drives the
    # code-phase advance.
    driver = first(sat.signals).signal
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    code_frequency = sat.code_doppler + get_code_frequency(driver)
    carrier_phase = update_carrier_phase(
        integrated_samples,
        carrier_frequency,
        sampling_frequency,
        sat.carrier_phase,
    )
    # Wrap depends on sync state (e.g. L1 C/A: 1023 chips before bit sync, 20460
    # after); see `current_code_wrap`.
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

# Rebuild each `TrackedSignal` by tuple recursion (allocation-free, inferable).
# On completion, push the raw accumulator as a `CorrelatorOutput` (ending at
# `signal_start_sample - 1`) onto the shared `correlator_outputs` buffer
# (allocation-free once its capacity is seated) and reset it, so the next
# integration in the chunk starts fresh; on a partial sub-step, carry it.
@inline _build_new_signals(::Tuple{}, ::Tuple{}, ::Int, ::Int) = ()
@inline function _build_new_signals(
    signals::Tuple,
    new_data::Tuple,
    integrated_samples::Int,
    signal_start_sample::Int,
)
    s = first(signals)
    (corr, completed) = first(new_data)
    total_integrated = s.integrated_samples + integrated_samples
    if completed
        push!(
            s.correlator_outputs,
            CorrelatorOutput(corr, total_integrated, signal_start_sample - 1),
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
            signal_start_sample,
        )...,
    )
end
