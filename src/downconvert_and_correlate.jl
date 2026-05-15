# Per-sat update applied after a single (downconvert + correlate) step over
# `integrated_samples` samples. Advances the shared carrier and code phase,
# the per-signal correlator accumulator, the per-signal integration counter,
# and the `is_integration_completed` flag. `signal_start_sample` advances on
# the satellite. Returns the new `TrackedSat`.
#
# In step 2 each satellite carries exactly one TrackedSignal, so the single
# new correlator and the single is_integration_completed argument suffice.
# Step 3 generalises this to a tuple of (correlator, is_integration_completed)
# pairs.
function update(
    system::AbstractGNSSSignal,
    sat::TrackedSat,
    integrated_samples::Int,
    intermediate_frequency,
    sampling_frequency,
    correlator,
    is_integration_completed::Bool,
)
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    code_frequency = sat.code_doppler + get_code_frequency(system)
    carrier_phase = update_carrier_phase(
        integrated_samples,
        carrier_frequency,
        sampling_frequency,
        sat.carrier_phase,
    )
    code_phase = update_code_phase(
        system,
        integrated_samples,
        code_frequency,
        sampling_frequency,
        sat.code_phase,
        has_bit_or_secondary_code_been_found(sat),
    )
    old_signal = only(sat.signals)
    new_signal = TrackedSignal(
        old_signal;
        integrated_samples = old_signal.integrated_samples + integrated_samples,
        is_integration_completed,
        correlator,
    )

    TrackedSat(
        sat;
        code_phase,
        carrier_phase,
        signal_start_sample = sat.signal_start_sample + integrated_samples,
        signals = (new_signal,),
    )
end
