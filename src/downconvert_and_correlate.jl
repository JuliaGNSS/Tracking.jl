# Advances a sat by one downconvert + correlate sub-step of `num_integrated_samples`
# samples: the shared carrier and code phase, `next_sample_index`, and each
# signal's correlator bank. `new_signals_data` holds one `(new_correlator_bank, completed)`
# pair per element of `sat.signals`; see `_build_new_signals` for completion.
function _advance_sat(
    sat::TrackedSat,
    num_integrated_samples::Int,
    intermediate_frequency,
    sampling_frequency,
    new_signals_data::Tuple,
)
    # All signals of a sat share one band, so one chip rate: drive it from signals[1].
    driver = first(sat.signals).signal
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    code_frequency = sat.code_doppler + get_chip_rate(driver)
    carrier_phase = update_carrier_phase(
        num_integrated_samples,
        carrier_frequency,
        sampling_frequency,
        sat.carrier_phase,
    )
    # Wraps at the primary code before sync and at the bit / secondary-code
    # period after (GPS L1 C/A: 1023, then 20460), so consumers such as
    # PositionVelocityTime.jl know the code period; see `current_code_phase_period`.
    code_length = current_code_phase_period(sat.signals)
    code_phase = mod(
        code_frequency * num_integrated_samples / sampling_frequency + sat.code_phase,
        code_length,
    )
    new_next_sample_index = sat.next_sample_index + num_integrated_samples

    new_signals = _build_new_signals(
        sat.signals,
        new_signals_data,
        num_integrated_samples,
        new_next_sample_index,
        Float64(code_phase),
    )

    TrackedSat(
        sat;
        code_phase,
        carrier_phase,
        next_sample_index = new_next_sample_index,
        signals = new_signals,
    )
end

# Rebuilds each TrackedSignal by tuple recursion (allocation-free, inferable). A
# completed integration pushes its raw correlator bank, tagged with its end sample
# and code phase, onto the reused `correlator_bank_outputs` (which
# `step_tracking_loops!` steps the tracking loop with) and resets it; a partial
# one carries it.
@inline _build_new_signals(::Tuple{}, ::Tuple{}, ::Int, ::Int, ::Float64) = ()
@inline function _build_new_signals(
    signals::Tuple,
    new_data::Tuple,
    num_integrated_samples::Int,
    next_sample_index::Int,
    code_phase::Float64,
)
    tsig = first(signals)
    (corr, completed) = first(new_data)
    total_integrated = tsig.num_integrated_samples + num_integrated_samples
    if completed
        push!(
            tsig.correlator_bank_outputs,
            CorrelatorBankOutput(corr, total_integrated, next_sample_index - 1, code_phase),
        )
        new_tsig =
            TrackedSignal(tsig; num_integrated_samples = 0, correlator_bank = zero(corr))
    else
        new_tsig = TrackedSignal(
            tsig;
            num_integrated_samples = total_integrated,
            correlator_bank = corr,
        )
    end
    (
        new_tsig,
        _build_new_signals(
            Base.tail(signals),
            Base.tail(new_data),
            num_integrated_samples,
            next_sample_index,
            code_phase,
        )...,
    )
end
