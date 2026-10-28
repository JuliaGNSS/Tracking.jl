# The software fill path of the noise reference: TrackingLoops' `despread_noise!`
# hook, running the despreads TrackingLoops plans through this package's kernels.
# The sampling policy, the window and the observations are TrackingLoops'.

"""
$(SIGNATURES)

Measure this signal's noise over samples `first_sample:last_sample` of
`measurement` and append the resulting observations, returning `estimator`.

Runs each of the slice's
[`noise_despreads`](@extref TrackingLoops.noise_despreads-Tuple{CorrelatorNoiseEstimator, AbstractGNSSSignal, Integer, Integer, Any, Any}) (their split into
sub-integrations, the PRN rotation and each one's random code phase and carrier are
TrackingLoops') with the taps of
[`noise_correlator`](@extref TrackingLoops.noise_correlator-Tuple{CorrelatorNoiseEstimator}) on `backend`'s kernel, the
one the prompt goes through, which keeps the measurement model-free, and appends it with
[`append_despread_noise!`](@extref TrackingLoops.append_despread_noise!-Tuple{CorrelatorNoiseEstimator, TrackingLoops.NoiseDespread, Any}). With
`num_ants > 1` every antenna column is despread, so the window holds the spatial
covariance. The reference is **open-loop** and reads no satellite state.
"""
function TrackingLoops.despread_noise!(
    backend::AbstractDownconvertAndCorrelator,
    estimator::CorrelatorNoiseEstimator,
    measurement::BandMeasurement,
    first_sample::Integer,
    last_sample::Integer,
    context::NoiseUpdateContext,
)
    signal = context.signal
    sampling_frequency = measurement.sampling_frequency
    correlator = noise_correlator(estimator)
    # The reference runs at the nominal chip rate, which every despread shares.
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, get_chip_rate(signal))
    # `gen_code_replica!` writes *at* `start_sample` while the kernel reads it
    # offset by `start_sample - 1`, so the buffer spans the whole measurement.
    # Backends that pack the code plane in-kernel ignore it.
    replica_size =
        get_num_samples(measurement) + maximum(sample_shifts) - minimum(sample_shifts)
    for despread in noise_despreads(
        estimator,
        signal,
        first_sample,
        last_sample,
        sampling_frequency,
        measurement.intermediate_frequency,
    )
        # The per-satellite despread: `carrier_phase = 0.0` (open-loop, no NCO),
        # `use_band_cache = false` (the item rides some group's loop, whose packed
        # sign planes may belong to another band or buffer).
        accumulators = get_accumulators(
            _despread_one_signal!(
                backend,
                zero(correlator),
                measurement.samples,
                despread.signal,
                despread.prn,
                sample_shifts,
                despread.code_phase,
                0.0,
                despread.code_frequency,
                despread.carrier_frequency,
                despread.sampling_frequency,
                despread.first_sample,
                despread.num_samples,
                replica_size,
                false,
            ),
        )
        append_despread_noise!(estimator, despread, accumulators)
    end
    estimator
end
