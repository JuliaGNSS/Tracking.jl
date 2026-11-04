# Software fill path of the noise reference: TrackingLoops' `despread_noise!` hook,
# running the despreads TrackingLoops plans through this package's kernels. The
# sampling policy, the window and the observations are TrackingLoops'.

"""
$(SIGNATURES)

Measure this signal's noise over samples `sample_range` of
`band_samples` and append the resulting observations, returning `estimator`.

Runs each of the slice's
[`noise_despreads`](@extref TrackingLoops.noise_despreads-Tuple{CorrelatorNoiseEstimator, AbstractGNSSSignal, AbstractUnitRange{<:Integer}, Any, Any}) (their split into
sub-integrations, the PRN rotation and each one's random code phase and carrier are
TrackingLoops') with the correlators of
[`noise_correlator_bank`](@extref TrackingLoops.noise_correlator_bank-Tuple{CorrelatorNoiseEstimator}) on `backend`'s kernel, the
one the prompt goes through, which keeps the measurement model-free, and appends it with
[`append_despread_noise!`](@extref TrackingLoops.append_despread_noise!-Tuple{CorrelatorNoiseEstimator, TrackingLoops.NoiseDespread, Any}). With
`num_antennas > 1` every antenna column is despread, so the window holds the spatial
covariance. The reference is **open-loop** and reads no satellite state.
"""
function TrackingLoops.despread_noise!(
    backend::AbstractDownconvertAndCorrelator,
    estimator::CorrelatorNoiseEstimator,
    band_samples::BandSamples,
    sample_range::AbstractUnitRange{<:Integer},
    context::NoiseUpdateContext,
)
    signal = context.signal
    sampling_frequency = band_samples.sampling_frequency
    correlator_bank = noise_correlator_bank(estimator)
    # The reference runs at the nominal chip rate, which every despread shares.
    sample_shifts = calc_correlator_sample_offsets(
        correlator_bank,
        sampling_frequency,
        get_chip_rate(signal),
    )
    # `gen_code_replica!` writes *at* `first(sample_range)` while the kernel reads it
    # offset by `first(sample_range) - 1`, so the buffer spans all of the band's samples.
    # Backends that pack the code plane in-kernel ignore it.
    replica_size =
        get_num_samples(band_samples) + maximum(sample_shifts) - minimum(sample_shifts)
    for despread in noise_despreads(
        estimator,
        signal,
        sample_range,
        sampling_frequency,
        band_samples.intermediate_frequency,
    )
        probe = despread.probe
        # The per-satellite despread: `carrier_phase = 0.0` (open-loop, no NCO),
        # `use_band_cache = false` (the item rides some group's loop, whose packed
        # sign planes may belong to another band or buffer).
        correlator_outputs = get_correlator_outputs(
            _despread_one_signal!(
                backend,
                zero(correlator_bank),
                band_samples.samples,
                probe.signal,
                probe.prn,
                sample_shifts,
                probe.code_phase,
                0.0,
                probe.code_frequency,
                probe.carrier_frequency,
                despread.sampling_frequency,
                despread.sample_range,
                replica_size,
                false,
            ),
        )
        append_despread_noise!(estimator, despread, correlator_outputs)
    end
    estimator
end
