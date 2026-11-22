# Noise Estimator

[`NoiseReferencedCN0Estimator`](@extref TrackingLoops.NoiseReferencedCN0Estimator), the
default C/N₀ estimator, divides each integration's prompt power by a **measured noise
floor**. This page is about where that floor comes from.

The floor is estimated **once per signal** and shared by every satellite on it, which
makes the reference's own variance negligible. The estimators live in
[`TrackState`](@ref)'s `noise_estimators` field, a NamedTuple keyed by signal id
(`GNSSSignals.get_signal_id`: `:GPSL1CA`, `:GalileoE1B`, …). On the sample-driven path
nothing needs configuring: `TrackState` provisions a
[`CorrelatorNoiseDensityEstimator`](@extref TrackingLoops.CorrelatorNoiseDensityEstimator)
for every signal asking for a density, and `track!` fills it. The rest of this page is
for hosts feeding correlator outputs from hardware.

## Why per signal and not per RF band

An integration divides by the **post-correlation** floor, weighted by the despreading
modulation's spectrum, so two signals on one front end see different floors as soon as
the interference is not white (see
[`AbstractNoiseDensityEstimator`](@extref TrackingLoops.AbstractNoiseDensityEstimator)).
Through `track!` on one set of samples with a 13 dB J/N tone, densities relative to the
thermal floor (4 MHz, pure BOC(1,1) E1B):

| tone      | GPS L1 C/A | Galileo E1B | E1B / C/A |
|-----------|-----------:|------------:|----------:|
| none      | 1.0        | 1.0         | 1.0       |
| 1.023 MHz | 1.4        | 39          | **28×**   |
| 0.4 MHz   | 35         | 26          | **0.75×** |

C/A's figure would overstate E1B's C/N₀ by ≈14.5 dB in the middle row, and the
reversal in the last row makes this a spectral effect, not a scale factor. The real
CBOC E1B shows the same reversal about four times larger. The *samples* remain a band
property; only the despread, and so the floor, is per signal.

The floor is a **density**: white noise gives `E|P|² = N₀/T`, so `N₀ = σ²/f_s` does not
depend on the integration time, and one figure serves 1- and 20-code-period
integrations alike. It carries the unit `1/Hz`.

## The two fill paths

|                                 | software                                            | hardware                              |
|:------------------------------- |:--------------------------------------------------- |:------------------------------------- |
| fills the window via            | `track!` → `downconvert_and_correlate!` → [`update_noise!`](@extref TrackingLoops.update_noise!-Tuple{TrackingLoops.AbstractNoiseDensityEstimator, Any, Any, Any}) | [`append_noise_observation!`](@ref)   |
| needs sample buffers?           | yes                                                 | no                                    |

Both fill the same
[`CorrelatorNoiseDensityEstimator`](@extref TrackingLoops.CorrelatorNoiseDensityEstimator).

```@docs
TrackingLoops.append_noise_observation!(::TrackState, ::TrackingLoops.NoiseObservation, ::Symbol)
```

## The software source

```@docs
TrackingLoops.despread_noise!
```

The reference is open-loop, with a randomised code phase, carrier and PRN per
sub-integration (see
[`CorrelatorNoiseDensityEstimator`](@extref TrackingLoops.CorrelatorNoiseDensityEstimator)).
It runs the prompt's own quantise → downconvert → despread → accumulate path on the
same kernel, so `N̂₀` already contains that chain's losses: the same 45 dB-Hz sky reads
≈45 dB-Hz on the float and Int16 backends, ≈44 on two-bit and ≈43 on one-bit, each
backend's own quantisation loss.

## Hardware producers (FPGA / ASIC)

Allocate a correlator channel as the noise channel rather than building a separate
power meter: it is one more slot on a time-multiplexed engine and cheaper than a
tracking channel on a parallel array, being open-loop. Reuse also keeps the
measurement model-free, with the same quantisation, input scaling and code generator
as the prompt; a path tapping at another point (before a decimation filter, at
another bit width) breaks that silently. The producer must:

1. **Report raw integer accumulations** with their weight and sample count, through
   one of the builders below; all division and scaling happen on the host.
2. **Randomise the code phase per dump and rotate the PRN** (a free-running code
   generator does this); dithering the carrier a few kHz strengthens it.

The weight `M` is the number of independent squared magnitudes, not the sample count:
`M` of them give a relative variance of `1/M`, which makes observations from different
producers combinable. A producer dumping on the code epoch matches the software
default; `≈1000` squared magnitudes per second per correlator is as good as finer below
40 dB-Hz.

- [`NoiseObservation`](@extref TrackingLoops.NoiseObservation)
- [`calc_noise_observation`](@extref TrackingLoops.calc_noise_observation-Tuple{Complex, Any, Any})
- [`calc_noise_observation_from_correlator_powers`](@extref TrackingLoops.calc_noise_observation_from_correlator_powers-NTuple{4, Any})
- [`calc_noise_observation_from_samples`](@extref TrackingLoops.calc_noise_observation_from_samples-Tuple{Any, Any, Any})

A worked ingest, per loop update:

```julia
# correlator bank outputs, per signal, exactly as before
append_correlator_bank_output!(track_state, CorrelatorBankOutput(correlator_bank, n, sample_index), prn)

# the noise channel, per signal: `accumulated_power` is the raw Σ|B|² pooled over
# `num_correlators` correlator outputs, so the sample count is per correlator times that
append_noise_observation!(
    track_state,
    calc_noise_observation_from_correlator_powers(
        accumulated_power,
        num_correlators,
        num_correlators * samples_per_dump,
        sampling_frequency;
        prn = noise_prn,
        # the correlators integrate the same samples, so the span counts once
        duration = samples_per_dump / sampling_frequency,
    ),
    :GPSL1CA,
)

step_tracking_loops!(track_state, (L1 = sampling_frequency,))
```

If a signal's noise estimator is never filled, the loop update **warns once per signal**
and every C/N₀ on it reads `-Inf dB-Hz`.

Unlike [`append_correlator_bank_output!`](@ref), which carries one satellite's whole
correlator bank with a load-bearing `end_sample_index` and is drained every NCO update
span, [`append_noise_observation!`](@ref) serves every satellite on the signal, pools
the correlator outputs into one power sum (one spatial covariance `R̂ = Σ_k b_k·b_kᴴ`
for an array, as a beamformer's floor is `wᴴR̂w`), needs no timing beyond "recent" and
lands in a sliding window bounded in time.

## Writing your own source

A source that is neither shipped path (a front-end power monitor, say) subtypes
[`AbstractNoiseDensityEstimator`](@extref TrackingLoops.AbstractNoiseDensityEstimator):

```julia
struct MyPowerMonitor <: TrackingLoops.AbstractNoiseDensityEstimator
    densities::Vector{typeof(1.0 / 1.0Hz)}
end

# fed from outside; `update_noise!` defaults to a no-op
TrackingLoops.append_noise_observation!(e::MyPowerMonitor, obs) =
    (push!(e.densities, obs.noise_density); e)
TrackingLoops.get_noise_density(e::MyPowerMonitor) =
    isempty(e.densities) ? nothing : sum(e.densities) / length(e.densities)
```

Pass it in through `TrackState`'s `noise_estimators` keyword:

```julia
TrackState(; signal = GPSL1CA(), noise_estimators = (GPSL1CA = MyPowerMonitor(...),))
```

Mutate the window **in place** and return the struct (`TrackState` is never rebuilt for
a noise update), and have
[`TrackingLoops.noise_density_type`](@extref TrackingLoops.noise_density_type-Tuple{TrackingLoops.AbstractNoiseDensityEstimator})
name the concrete type `get_noise_density` returns (default `typeof(1.0/1.0Hz)`), so
the loop update stays monomorphic.
