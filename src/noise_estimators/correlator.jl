# The software fill path of the noise reference. The estimator, its sliding
# window and the observation builders are TrackingLoops'; what stays here is the
# one thing only this package can do — despread an untracked PRN through its own
# downconvert-and-correlate kernels — as the method of TrackingLoops'
# `despread_noise!` hook that `update_noise!` forwards to.

"""
$(SIGNATURES)

Measure this signal's noise over samples `first_sample:last_sample` of
`measurement` and append the resulting observations, returning `estimator`.

The slice is split into `num_sub` **equal** sub-integrations,

```julia
num_sub = max(1, round(Int, slice_duration / code_period))
```

with `code_period` the primary code period of the signal being measured — see
[`CorrelatorNoiseEstimator`](@ref) for why the sub-integration length is derived
rather than configured. Equal slices rather than fixed-length ones so no
remainder is wasted and every window entry is statistically identical; the
`max(1, …)` matters, because a signal whose code period exceeds the
chunk (a 10 ms code with 1 ms chunks) would otherwise yield no observations at
all, forever. Each observation is pushed **individually**: pre-averaging would
put one entry per chunk in the window, re-welding its time span to the Doppler
update rate.

Each sub-integration draws a **fresh random code phase** and a **fresh random
carrier offset** within ±`carrier_dither` of the band's nominal IF, then runs for
its slice. **Nothing here follows a code-period boundary**, and that is
deliberate: the reference despreads with a wrong PRN at a wrong phase, so there
is no signal to accumulate coherently and `N̂₀ = |B|²/(N·A_c²·f_s)` is unbiased
for any `N` and any starting phase. Following boundaries would be actively
harmful — a chunk rarely holds a whole number of code periods, so a
boundary-aligned reference would have to carry a partial accumulator across
calls, which is exactly the scalar state that has no per-signal anchor.
Successive observations are independent because their *sample* ranges are
disjoint; the replica repeating is irrelevant.

The draws are per sub-integration rather than per chunk so that every window
entry is an independent trial — see [`CorrelatorNoiseEstimator`](@ref) for why a
*stationary* phase and Doppler turn a chance alignment with a present-but-
untracked signal into a permanent bias, and the randomisation turns it back into
an occasional single-observation outlier.

All of the correlator's taps are used and pooled into one power sum, at
`tap_code_shift` chips apart so they are genuinely independent looks (see the
type's docstring). The despread runs on the caller's backend kernel — the same
one the prompt goes through — which is what makes the measurement model-free.

With `num_ants > 1` the reference despreads **every** antenna column and pools
the taps into the array's spatial covariance `R̂ = Σ_k b_k·b_kᴴ` rather than a
scalar power. One window per signal still, shared by every satellite on it and
just as satellite-agnostic as the scalar case: the reduction to a per-satellite
floor happens at C/N₀ time, where `wᴴR̂w` under that satellite's own weights is
the exact noise scale of the very prompt being divided. Measuring one column
instead — as this did before — and assuming unity noise gain silently
invalidated the ratio for any real beamformer.

The reference is **open-loop**: no discriminator, no loop filter, no NCO update,
and nothing read from satellite state at all. The PRN rotates over the whole
family, tracked or not, because a random phase makes avoiding the tracked ones
unnecessary.
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

    # Every antenna column, despread through a correlator of the estimator's own
    # width. The window then holds the array's spatial covariance rather than one
    # column's scalar power, which is what lets each satellite ask for the floor
    # *its own* beamformer sees — `wᴴR̂w` — instead of forcing every satellite to
    # share one antenna's answer and assume unity noise gain on top.
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
    # `gen_code_replica!` writes *at* `start_sample` while the kernel reads the
    # replica offset by `start_sample - 1`, so the buffer spans the whole
    # measurement rather than one slice. Only the software backends read it at
    # all; `_despread_one_signal!` sizes it from this on those that do and ignores
    # it on the ones that pack the code plane in-kernel.
    replica_size =
        get_num_samples(measurement) + maximum(sample_shifts) - minimum(sample_shifts)

    rng = estimator.rng
    code_length = get_code_length(signal_type)
    carrier_dither = estimator.carrier_dither
    intermediate_frequency = measurement.intermediate_frequency

    for sub = 1:num_sub
        # Equal slices: the `sub`-th covers samples `⌊(sub−1)N/num_sub⌋+1 … ⌊subN/num_sub⌋`.
        slice_start = Int(first_sample) + div((sub - 1) * num_samples, num_sub)
        slice_stop = Int(first_sample) + div(sub * num_samples, num_sub) - 1
        slice_samples = slice_stop - slice_start + 1
        slice_samples > 0 || continue
        # One draw each, per sub-integration. Both are free — the kernels already
        # take a code phase and a carrier frequency — and both are what keep a
        # chance alignment with a present-but-untracked signal from freezing into
        # a standing bias. See the type's docstring.
        code_phase = rand(rng) * code_length
        carrier_frequency = intermediate_frequency + (2 * rand(rng) - 1) * carrier_dither
        # The same despread primitive the per-satellite path runs, which is what
        # keeps the measurement model-free: `carrier_phase = 0.0` because the
        # reference is open-loop and has no NCO to continue from, and
        # `use_band_cache = false` because this runs before any group has packed
        # the bit-wise backends' shared sign planes.
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
        # Straight to the builders' shared core rather than through
        # `noise_observation_from_correlator`'s keyword form: the keywords cost a
        # per-call allocation on Julia 1.10 that 1.11+ elides, and this is the one
        # place it runs per chunk. The arguments are the same ones that builder
        # would forward — note the duration, which is `slice_samples` and not
        # `num_taps` times it, because the taps are simultaneous rather than
        # consecutive: they all integrate the same samples.
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
