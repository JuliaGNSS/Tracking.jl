# Noise Estimator

[`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator), the default C/N₀
estimator, divides each record's prompt power by a **measured noise floor**. This page
is about where that floor comes from.

The floor is estimated **once per signal**, shared by every satellite on it, which
makes the reference's own variance negligible. The estimators live in
[`TrackState`](@ref)'s `noise_estimators` field, a NamedTuple keyed by signal id
(`GNSSSignals.get_signal_id`: `:GPSL1CA`, `:GalileoE1B`, …).

On the sample-driven path nothing needs configuring: `TrackState` provisions a
[`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator) for every
signal asking for a density, and `track!` fills it. Read on if you feed correlator
outputs from hardware.

In the TrackingLoops manual:

- [`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)

## Why per signal and not per RF band

What a record divides by is the **post-correlation** floor, weighted by the
despreading modulation's own spectrum, so two signals on one front end see different
floors as soon as the interference is not white (see
[`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)). Measured
through `track!` on one shared set of samples, 13 dB J/N tone, densities relative to
the thermal floor (4 MHz, pure BOC(1,1) E1B):

| tone      | GPS L1 C/A | Galileo E1B | E1B / C/A |
|-----------|-----------:|------------:|----------:|
| none      | 1.0        | 1.0         | 1.0       |
| 1.023 MHz | 1.4        | 39          | **28×**   |
| 0.4 MHz   | 35         | 26          | **0.75×** |

Reporting C/A's figure to E1B in the middle row overstates E1B's C/N₀ by ≈14.5 dB;
the **reversal** in the last row makes this a spectral effect, not a scale factor.
The real CBOC E1B shows the same reversal about four times larger. Despreading with
the consumer's own code measures this weighting rather than modelling it.

The *samples* remain a band property (`BandMeasurements` is keyed by band id); only
the despread, and so the floor, is per signal.

## Why a density and not a power

With Tracking's normalisation, white noise gives `E|P|² = N₀/T`, so `N₀ = σ²/f_s` is
**independent of the integration time** (see
[`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)). One
per-signal figure therefore serves a satellite on 1-block pre-sync records and one on
20-block post-sync records alike, and the reference needs no alignment to correlator
records: it runs on any grid and averages over a much longer window than any prompt
buffer. `N₀` carries Unitful dimensions of `1/Hz`, so the arithmetic is
dimension-checked.

## The two fill paths

|                                 | software                                            | hardware                              |
|:------------------------------- |:--------------------------------------------------- |:------------------------------------- |
| fills the window via            | `track!` → `downconvert_and_correlate!` → [`update_noise!`](@extref TrackingLoops.update_noise!-Tuple{TrackingLoops.AbstractNoiseEstimator, Vararg{Any, 4}}) | [`append_noise_observation!`](@ref)   |
| enters `downconvert_and_correlate!`? | yes                                            | **no**                                |
| needs sample buffers?           | yes                                                  | no                                    |

The paths never meet (a hardware producer never calls `downconvert_and_correlate!`),
so you configure
[`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator) either way
and only fill it differently.

```@docs
append_noise_observation!
```

In the TrackingLoops manual:

- [`update_noise!`](@extref TrackingLoops.update_noise!-Tuple{TrackingLoops.AbstractNoiseEstimator, Vararg{Any, 4}})
- [`get_noise_density`](@extref TrackingLoops.get_noise_density-Tuple{TrackingLoops.AbstractNoiseEstimator})

## The software source

```@docs
TrackingLoops.despread_noise!
```

In the TrackingLoops manual:

- [`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator)
- [`CorrelatorNoiseEstimator()`](@extref TrackingLoops.CorrelatorNoiseEstimator-Tuple{})
- [`update_noise!(::CorrelatorNoiseEstimator, …)`](@extref TrackingLoops.update_noise!-Tuple{CorrelatorNoiseEstimator, Any, Integer, Integer, TrackingLoops.NoiseUpdateContext})

Its design is documented on
[`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator): it is
**open-loop** and reads no satellite state, so it works with zero satellites tracked;
its code phase and carrier are **randomised** per sub-integration, so a chance
alignment with an untracked signal is a one-observation outlier rather than a standing
bias, which also lets the PRN rotate over the whole family; its sub-integration length
is the primary code period, **derived, not configured**, with the variance bought back
by `window_duration`.

It is **model-free**: the reference runs the prompt's own quantise → downconvert →
despread → accumulate path on the same kernel, so `N̂₀` already contains that chain's
losses with no per-backend correction. The same 45 dB-Hz sky reports ≈45 dB-Hz on the
float and Int16 backends, ≈44 on two-bit and ≈43 on one-bit: each backend's own
quantisation loss.

## Hardware producers (FPGA / ASIC)

**Reuse the downconvert-and-correlate datapath: allocate a channel.** On the
architectures FPGA correlators use, that is also the cheaper option:

- On a **time-multiplexed engine** (one datapath across N channels, per-channel
  state in BRAM) a noise channel is one more slot in the schedule; a dedicated path
  is a second datapath.
- On a **fully parallel array** it costs less than a tracking channel, being
  open-loop: fixed NCO increments, no discriminator, no loop-filter interface, no
  dump-on-command synchronisation. A dedicated path would re-implement the carrier
  NCO, code generator and complex multiplier.
- **No divider and no float on chip.** Report the raw integer accumulation with its
  weight and sample count; all division, `A_c²` and `f_s` scaling and averaging
  happen in Tracking.

Reuse is also what makes the measurement model-free: identical quantisation (the
1-bit/2-bit loss cancels in the ratio), identical input scaling (an AGC step moves
numerator and denominator together), identical code generator (`A_c` matches). A
path tapping at a different point, before a decimation filter or at a different bit
width, breaks all three silently.

### What the producer must match

1. **Report raw accumulations** through one of the builders below, on the
   correlator outputs' return path: `append_noise_observation!` sits beside
   [`append_correlator_output!`](@ref) on the host.
2. **Randomise the code phase per dump, and rotate the PRN** (a free-running code
   generator does this; phase 0 would need a reset). Dithering the carrier a few kHz
   around the nominal IF strengthens it further. The noise channel then needs no
   input from the tracking channels.

A producer that dumps on the code epoch matches the software default exactly. One
that differs still lands in the same `M`-weighted, time-bounded window, on the same
scale; only the precision differs, and `K_n ≈ 1000` looks per second per tap is
indistinguishable from anything finer below 40 dB-Hz.

### The observation and its builders

In the TrackingLoops manual:

- [`NoiseObservation`](@extref TrackingLoops.NoiseObservation)
- [`noise_observation`](@extref TrackingLoops.noise_observation-Tuple{Complex, Any, Any})
- [`noise_observation_from_correlator`](@extref TrackingLoops.noise_observation_from_correlator-NTuple{4, Any})
- [`noise_observation_from_samples`](@extref TrackingLoops.noise_observation_from_samples-Tuple{Any, Any, Any})

All three reduce to the same `N₀` on the same white input. The weight is `M`, **the
number of independent looks, not the sample count**: `M` sub-integrations give a
relative variance of `1/M`, so one 1 ms despread and 64 16-chip despreads cover the
same samples at very different precision. Weighting by `M` makes observations from
different producers combinable.

A worked ingest, per fold:

```julia
# taps, per signal, exactly as before
append_correlator_output!(track_state, CorrelatorOutput(correlator, n, sample_index), prn)

# the noise channel, per signal: `accumulated_power` is the raw Σ|B|² pooled over
# `num_taps` accumulators, so the sample count is per tap times that
append_noise_observation!(
    track_state,
    noise_observation_from_correlator(
        accumulated_power,
        num_taps,
        num_taps * samples_per_dump,
        sampling_frequency;
        prn = noise_prn,
        # the taps integrate the same samples, so the span counts once
        duration = samples_per_dump / sampling_frequency,
    ),
    :GPSL1CA,
)

estimate_dopplers_and_filter_prompt!(track_state, (L1 = sampling_frequency,))
```

If a signal's noise estimator is never filled, the fold **warns once per signal**
and every C/N₀ on it reads `-Inf dB-Hz`.

### How `append_noise_observation!` differs from `append_correlator_output!`

|            | [`append_correlator_output!`](@ref)                                    | [`append_noise_observation!`](@ref)                    |
|:---------- |:---------------------------------------------------------------------- |:------------------------------------------------------ |
| scope      | one **signal** of one satellite                                        | one **signal**, across every satellite tracking it     |
| payload    | the whole correlator: every tap kept **separate**, complex, per antenna | the taps **pooled** — into one power sum, or into one spatial covariance for an array — plus `M`, the span and the PRN |
| drives     | discriminators, both loop filters, bit sync, the NCO update, the C/N₀ prompt | the noise floor only                              |
| timing     | `sample_index` is load-bearing — vector tracking needs a common grid    | none: `N₀` is slowly varying, so only "recent" matters |
| lifetime   | **drained** by the fold each chunk, buffer reused                      | **sliding window** bounded in time, never drained      |
| count      | one per completed coherent integration per signal                      | any granularity; `M` carries the weight                |

The prompt path keeps the taps apart because their *differences* are the
discriminants; the reference pools them, as at ≥1 chip spacing they are independent
looks at the same noise.

It keeps the **antenna** dimension: with `num_ants > 1` the payload is the `M×M`
spatial covariance `R̂ = Σ_k b_k·b_kᴴ`, since a beamformer's floor is `wᴴR̂w` under its
own weights. Every builder takes the per-antenna form of its input.

## Writing your own source

A source that is neither shipped path (a front-end power monitor, say) subtypes
[`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator):

```julia
struct MyPowerMonitor <: TrackingLoops.AbstractNoiseEstimator
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

Two contracts: mutate the window **in place** and return the struct (`TrackState` is
never rebuilt for a noise update), and have
[`TrackingLoops.noise_density_type`](@extref TrackingLoops.noise_density_type-Tuple{TrackingLoops.AbstractNoiseEstimator}) name the
concrete type `get_noise_density` returns (default `typeof(1.0/1.0Hz)`), so the fold
stays monomorphic.

In the TrackingLoops manual:

- [`noise_density_type`](@extref TrackingLoops.noise_density_type-Tuple{TrackingLoops.AbstractNoiseEstimator})
- [`NoiseUpdateContext`](@extref TrackingLoops.NoiseUpdateContext)
- [`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:TrackingLoops.AbstractCN0Estimator}})
