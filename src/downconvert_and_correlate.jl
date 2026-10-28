# Advances a sat by one downconvert + correlate sub-step of `integrated_samples`
# samples: the shared carrier and code phase, `signal_start_sample`, and each
# signal's accumulator. `new_signals_data` holds one `(new_correlator, completed)`
# pair per element of `sat.signals`; see `_build_new_signals` for completion.
function update(
    sat::TrackedSat,
    integrated_samples::Int,
    intermediate_frequency,
    sampling_frequency,
    new_signals_data::Tuple,
)
    # All signals of a sat share one band, so one chip rate: drive it from signals[1].
    driver = first(sat.signals).signal
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    code_frequency = sat.code_doppler + get_chip_rate(driver)
    carrier_phase = update_carrier_phase(
        integrated_samples,
        carrier_frequency,
        sampling_frequency,
        sat.carrier_phase,
    )
    # Wraps at the primary code before sync and at the bit / secondary-code
    # period after (GPS L1 C/A: 1023, then 20460), so consumers such as
    # PositionVelocityTime.jl know the code block; see `current_code_wrap`.
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
        Float64(code_phase),
    )

    TrackedSat(
        sat;
        code_phase,
        carrier_phase,
        signal_start_sample = new_signal_start_sample,
        signals = new_signals,
    )
end

# Rebuilds each TrackedSignal by tuple recursion (allocation-free, inferable). A
# completed integration pushes its raw accumulator, tagged with its end sample
# and code phase, onto the reused `correlator_outputs` (which the Doppler
# estimator folds over) and resets it; a partial one carries it.
@inline _build_new_signals(::Tuple{}, ::Tuple{}, ::Int, ::Int, ::Float64) = ()
@inline function _build_new_signals(
    signals::Tuple,
    new_data::Tuple,
    integrated_samples::Int,
    signal_start_sample::Int,
    code_phase::Float64,
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
            signal_start_sample,
            code_phase,
        )...,
    )
end
