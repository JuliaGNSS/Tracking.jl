# The software fill path of the noise reference: TrackingLoops' `despread_noise!`
# hook, despreading an untracked PRN through this package's kernels. The
# estimator, window and observation builders are TrackingLoops'.

"""
$(SIGNATURES)

Measure this signal's noise over samples `first_sample:last_sample` of
`measurement` and append the resulting observations, returning `estimator`.

The slice is split into `num_sub` **equal** sub-integrations,

```julia
num_sub = max(1, round(Int, slice_duration / code_period))
```

with `code_period` the primary code period of the signal being measured (see
[`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator) for why it
is derived, not configured). Equal slices waste no remainder; the `max(1, …)` keeps a
code longer than the chunk (10 ms code, 1 ms chunks) from never yielding an
observation. Each observation is pushed individually, so the window's time span does
not follow the Doppler update rate.

Each sub-integration draws a fresh random code phase and a fresh random carrier offset
within ±`carrier_dither` of the band's IF, so every window entry is an independent trial
(see [`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator)). It
follows no code-period boundary: with a wrong PRN at a wrong phase there is nothing to
accumulate coherently, so `N̂₀ = |B|²/(N·A_c²·f_s)` is unbiased for any `N` and start
phase, and no partial accumulator has to be carried across calls.

All taps, `tap_code_shift` chips apart, are pooled into one power sum, despread on the
caller's backend kernel (the one the prompt goes through), which keeps the measurement
model-free. With `num_ants > 1` every antenna column is despread and the taps pool into
the spatial covariance `R̂ = Σ_k b_k·b_kᴴ`; each satellite reduces it to its own floor
`wᴴR̂w` at C/N₀ time.

The reference is **open-loop** and reads no satellite state. The PRN rotates over the
whole family, tracked or not; the random phase makes avoiding tracked ones unnecessary.
"""
function TrackingLoops.despread_noise!(
    backend::AbstractDownconvertAndCorrelator,
    estimator::CorrelatorNoiseEstimator,
    measurement::BandMeasurement,
    first_sample::Integer,
    last_sample::Integer,
    context::NoiseUpdateContext,
)
    num_samples = Int(last_sample) - Int(first_sample) + 1
    num_samples > 0 || return estimator
    signal_type = context.signal
    sampling_frequency = measurement.sampling_frequency
    # Nominal chip rate: the reference runs at zero Doppler.
    code_frequency = get_code_frequency(signal_type)
    code_period = get_code_length(signal_type) / code_frequency
    slice_duration = num_samples / sampling_frequency
    num_sub = max(1, round(Int, uconvert(NoUnits, slice_duration / code_period)))
    num_sub = min(num_sub, num_samples)          # never fewer than one sample per slice

    # Every antenna column, so the window holds the spatial covariance.
    samples = measurement.samples
    num_ants = _num_ants(estimator)
    correlator = EarlyPromptLateCorrelator(;
        num_ants,
        preferred_early_late_to_prompt_code_shift = estimator.tap_code_shift,
    )
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    num_taps = length(sample_shifts)
    code_amplitude = get_code_amplitude(signal_type)

    prn = _next_noise_prn(estimator, signal_type)
    # `gen_code_replica!` writes *at* `start_sample` while the kernel reads it
    # offset by `start_sample - 1`, so the buffer spans the whole measurement.
    # Backends that pack the code plane in-kernel ignore it.
    replica_size =
        get_num_samples(measurement) + maximum(sample_shifts) - minimum(sample_shifts)

    rng = estimator.rng
    code_length = get_code_length(signal_type)
    carrier_dither = estimator.carrier_dither
    intermediate_frequency = measurement.intermediate_frequency

    for sub = 1:num_sub
        slice_start = Int(first_sample) + div((sub - 1) * num_samples, num_sub)
        slice_stop = Int(first_sample) + div(sub * num_samples, num_sub) - 1
        slice_samples = slice_stop - slice_start + 1
        slice_samples > 0 || continue
        # Fresh draws per sub-integration (see the docstring).
        code_phase = rand(rng) * code_length
        carrier_frequency = intermediate_frequency + (2 * rand(rng) - 1) * carrier_dither
        # The per-satellite despread: `carrier_phase = 0.0` (open-loop, no NCO),
        # `use_band_cache = false` (the item rides some group's loop, whose packed
        # sign planes may belong to another band or buffer).
        accumulators = get_accumulators(
            _despread_one_signal!(
                backend,
                zero(correlator),
                samples,
                signal_type,
                prn,
                sample_shifts,
                code_phase,
                0.0,
                code_frequency,
                carrier_frequency,
                sampling_frequency,
                slice_start,
                slice_samples,
                replica_size,
                false,
            ),
        )
        power = _pool_taps(accumulators, num_ants)
        # The builders' core rather than `noise_observation_from_correlator`, whose
        # keywords allocate per call on Julia 1.10. The duration is `slice_samples`,
        # not `num_taps` times it: the taps integrate the same samples.
        append_noise_observation!(
            estimator,
            _noise_observation(
                power,
                num_taps,
                num_taps * slice_samples,
                sampling_frequency,
                code_amplitude,
                prn,
                slice_samples / sampling_frequency,
            ),
        )
    end
    estimator
end
