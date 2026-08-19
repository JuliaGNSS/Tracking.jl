# Noise Estimator

[`NoiseRefCN0Estimator`](@ref) — the default C/N₀ estimator — divides each
record's prompt power by a **measured noise floor** instead of inferring one from
the prompt's own statistics. This page is about where that floor comes from.

The floor is estimated **once per signal**, never per satellite: every satellite
tracking a signal shares its floor, and averaging once per signal is what makes
the reference's own variance negligible against the per-record prompt statistics.
The estimators live in [`TrackState`](@ref)'s `noise_estimators` field, a
NamedTuple keyed by signal id (`GNSSSignals.get_signal_id` — `:GPSL1CA`,
`:GalileoE1B`, …).

On the sample-driven path none of this needs configuring — `TrackState`
provisions a [`CorrelatorNoiseEstimator`](@ref) for every signal asking for a
density, and `track!` fills it. Read on if you feed correlator outputs from
hardware, or if you want to know what the measurement actually does.

```@docs
AbstractNoiseEstimator
```

## Why per signal and not per RF band

Because what a record divides by is the **post-correlation** floor. For an
interferer of PSD `S_I(f)` on top of the thermal floor that is

```
N₀,eff = N₀,thermal + ∫ S_I(f) · |G(f)|² df
```

— the spectral separation coefficient — weighted by the *despreading
modulation's* own spectrum `|G(f)|²`. So two signals sharing one band, one
antenna and one front end do not share a noise floor the moment the interference
is not white, and the gap is not subtle. BPSK(1) peaks at band centre and nulls
at ±1.023 MHz; BOC(1,1) is the reverse. A CW tone at the chip rate is rejected by
GPS L1 C/A and lands squarely on Galileo E1B; move it into C/A's main lobe and
the order reverses. Measured through `track!` on one shared set of samples,
13 dB J/N tone, densities relative to the thermal floor (4 MHz, pure BOC(1,1)
E1B):

| tone      | GPS L1 C/A | Galileo E1B | E1B / C/A |
|-----------|-----------:|------------:|----------:|
| none      | 1.0        | 1.0         | 1.0       |
| 1.023 MHz | 1.4        | 39          | **28×**   |
| 0.4 MHz   | 35         | 26          | **0.75×** |

Reporting C/A's figure to E1B in the middle row overstates E1B's C/N₀ by
≈14.5 dB; the **reversal** in the last row is what makes this a spectral result
and not a scale factor. The real CBOC E1B shows the same reversal about four
times larger. Front-end tilt, filter roll-off at the band edge and adjacent-band
leakage colour the floor the same way, more quietly.

The payoff of despreading rather than metering power is that this integral is
**measured, never modelled**: the reference runs the consumer's own code, so its
spectral weighting is the consumer's by construction, with no jammer model and no
spectrum analysis anywhere. It is the same argument that makes the measurement
backend-free, extended from the quantiser to the interference environment.

It also matches the hardware: a noise channel is a tracking channel with a wrong
PRN, and a tracking channel is configured with a code.

The *samples* remain a band property — one front end feeds every signal on it,
and `BandMeasurements` stays keyed by band id. Only the despread, and therefore
the measured floor, is per signal.

## Why a density and not a power

For a correlator normalised the way Tracking normalises — by
`integrated_samples * code_amplitude` — white input noise of per-sample variance
`σ²` gives

```
E|P|² = σ²/N = N₀/T
```

so `N₀ = σ²/f_s` is **independent of the integration time**. That is the whole
reason a density is stored rather than a power, and it buys two things:

- **One per-signal figure serves records of any length.** A satellite still on
  1-block records pre-sync and one promoted to 20-block records post-sync divide
  by the same `N̂₀`, each subtracting its own `1/T`. Mixed record lengths are
  handled by construction rather than by guarding.
- **The reference needs no alignment to correlator records.** It can run on
  whatever grid is convenient and be averaged over a much longer window than any
  satellite's prompt buffer, because `N₀` is a slowly-varying front-end property.

`N₀` carries Unitful dimensions of `1/Hz`, so `⟨|P|²⟩/N₀` and `1/T` are both in
`Hz` and the arithmetic is dimension-checked rather than trusted.

## The two fill paths

|                                 | software                                            | hardware                              |
|:------------------------------- |:--------------------------------------------------- |:------------------------------------- |
| fills the window via            | `track!` → `downconvert_and_correlate!` → [`update_noise!`](@ref) | [`append_noise_observation!`](@ref)   |
| enters `downconvert_and_correlate!`? | yes                                            | **no**                                |
| needs sample buffers?           | yes                                                  | no                                    |

They live on **disjoint call graphs** — a hardware producer never calls
`downconvert_and_correlate!`, it injects correlator outputs and folds — so one
concrete type serves both and nothing needs to distinguish them. You configure
[`CorrelatorNoiseEstimator`](@ref) either way; you simply fill it differently.

```@docs
update_noise!
append_noise_observation!
get_noise_density
```

## The software source

```@docs
CorrelatorNoiseEstimator
CorrelatorNoiseEstimator()
Tracking.update_noise!(::CorrelatorNoiseEstimator, ::Tracking.BandMeasurement, ::Integer, ::Integer, ::Tracking.NoiseUpdateContext)
```

Four properties are worth knowing about it.

**It is open-loop.** The reference runs at its band's nominal IF, dithered, with
a random code phase and a rotating PRN. There is no discriminator, no loop filter
and no NCO update ever written to it, and it reads nothing at all from satellite
state — so it is correct with **zero satellites tracked**, which is also what
would make it usable as an acquisition CFAR floor.

**Its code phase and Doppler are randomised**, drawn afresh for every
sub-integration. The reason is not an adversary; it is that a *stationary*
reference makes a chance alignment permanent. The replica is re-anchored to phase
0 on a grid running at the nominal chip rate, so against any incoming signal at
zero Doppler the relative phase is frozen wherever it lands — and a
present-but-untracked signal (a spoofed PRN, or a visible satellite the receiver
has not acquired) has a `3 × 1.5 / 1023 ≈ 0.44 %` chance per PRN of landing
inside a tap, worth `T·(C/N₀)/3 ≈ 10.5·N₀` at 45 dB-Hz on every observation with
that code, forever. The geometry is perverse too: relative phase drifts only at
the signal's own code Doppler, so *low* Doppler both enables a hit and makes it
stick. Drawing the phase and the carrier per sub-integration turns that into an
independent trial per observation — a ≈0.07 dB residual instead of a ≈1.5 dB
standing bias — and since a hit needs both draws, the two multiply.

This is also why the rotation covers the **whole PRN family**, tracked codes
included. Skipping tracked PRNs was what made a fixed phase 0 safe; a random
phase makes it unnecessary, and dropping it removes the reference's only
dependency on satellite state while keeping the pool at a constant 32 — under the
old rule the pool shrank, and diluted a bad draw less well, exactly as the
receiver acquired more satellites. Landing within ±1 chip of a tracked peak now
costs ≈0.045 dB for the single observation it touches.

**It is model-free**, and that is the reason it despreads rather than measuring
`Σ|x|²`. The reference traverses the *identical* quantise → downconvert →
despread → accumulate path as the prompt, on the same kernel, so the measured
`N̂₀` already contains every imperfection of that chain — the quantisation loss,
the quantiser's operating point under load, the input scaling an AGC step moves,
a CBOC replica's code amplitude — with no per-backend analysis and no closed-form
correction. It shows up directly: the same 45 dB-Hz sky reports ≈45 dB-Hz on the
float and Int16 backends, ≈44 on two-bit and ≈43 on one-bit, i.e. each backend's
own quantisation loss, without a line of per-backend correction anywhere.

**Its sub-integration length is derived, not configured.** Coherent integration
buys *nothing* for noise-power estimation: for a complex-Gaussian accumulation
`Var(|B|²) = (E|B|²)²`, so one dump carries 100 % relative error however long it
is, and all the information about `σ²` lives in the number of independent looks.
Shortening the dump is therefore the only way to buy more of them — and three
things say not to. SIMD (64 kernel calls per 1 ms chunk instead of one, and a
shorter dump would need a new kernel variant, forfeiting the bit-identical
arithmetic the approach rests on). DC balance (a full-period Gold code is
balanced, `E[(Σc)²] ≈ 1`; a 16-chip window behaves like iid signs at ≈16, so an
ADC offset biases a short despread and not a full-period one). And cadence parity
with a hardware correlator, which naturally dumps on the code epoch. The variance
is bought back with `window_duration` instead, which is nearly free — a 1 s
window is ~1000 `Float64` densities per signal.

## Hardware producers (FPGA / ASIC)

The question a hardware implementer asks first is whether the noise measurement
should reuse the existing downconvert-and-correlate datapath or get a dedicated,
minimal one. **Reuse it — allocate a channel.** Not because it is less work, but
because on the architectures FPGA correlators actually use it is also cheaper:

- On a **time-multiplexed engine** — one physical datapath across N channels at a
  high clock, per-channel state in BRAM — a noise channel is one more slot in the
  schedule. A handful of state words, essentially zero logic. A dedicated path is
  a whole second datapath.
- On a **fully parallel array** it costs *less* than one channel out of N,
  because the reference is open-loop: fixed NCO increments, no discriminator, no
  loop-filter interface, no dump-on-command synchronisation. A dedicated path
  would have to re-implement the carrier NCO, code generator and complex
  multiplier to save accumulators you were not using — so it is larger, not
  smaller.
- **No divider and no float on chip.** Report the raw integer accumulation with
  its weight and sample count; every division, the `A_c²` and `f_s` scaling, and
  all averaging happen in Tracking. That is both the cheapest option on-chip and
  the mechanism that forces the two paths to agree.

And reuse is what buys the model-freeness above: identical quantisation and
rounding, so the 1-bit/2-bit loss cancels in the ratio; identical input scaling,
so an AGC step moves numerator and denominator together; identical code
generator, so `A_c` matches exactly. A path tapping at a different point — before
a decimation filter, at a different bit width — breaks all three silently.

### What the producer must match

Only two things, because everything else is computed here:

1. **Report raw accumulations**, through one of the builders below. They ride the
   same return path as the correlator outputs, so no new wire interface is
   needed: `append_noise_observation!` sits beside
   [`append_correlator_output!`](@ref) on the host.
2. **Randomise the code phase per dump, and rotate the PRN.** Do not pin the
   phase — see above for why a stationary reference turns a chance alignment into
   a permanent one. On an FPGA this is the *easy* option: a free-running code
   generator gives an arbitrary phase, where phase 0 needs a reset. Dithering the
   carrier a few kHz either side of the nominal IF costs one more LFSR draw and
   multiplies the protection. With the phase randomised there is no need to know
   which PRNs are tracked, so the noise channel needs no input from the tracking
   channels at all.

A producer that can only dump on the code epoch matches the software default
exactly, so in the normal case there is nothing to reconcile. Where a producer
*does* differ it still does not conflict: "the two paths behave the same" means
*computed by the same code, on the same scale, with the same bias*, which the
shared window guarantees because it is `M`-weighted and bounded in time. Only the
precision differs — a code-epoch dump gives `K_n ≈ 1000` looks per second per
tap, which is indistinguishable from anything finer below 40 dB-Hz.

### The observation and its builders

```@docs
NoiseObservation
noise_observation
noise_observation_from_correlator
noise_observation_from_samples
```

All three reduce to the same `N₀` on the same white input, which is what makes
the paths interchangeable.

The weight is `M`, **the number of independent looks — not the sample count.**
For `Σ_m |B_m|²` over `M` sub-integrations of `N` samples each the density needs
only `M·N`, but its relative variance is `1/M`: a producer handing back one 1 ms
coherent despread and one handing back 64 16-chip despreads report the same total
samples and wildly different precision. Weighting the window by `M` is what makes
observations from different producers combinable at all.

A worked ingest, per fold:

```julia
# taps, per signal, exactly as before
append_correlator_output!(track_state, CorrelatorOutput(correlator, n, sample_index), prn)

# the noise channel: Σ|B|² pooled over its taps, one signal at a time. A noise
# channel is configured with a code, so this is the granularity a producer
# already has.
# `accumulated_power` is the raw integer sum off the wire; `num_taps` is how many
# accumulators it pooled, and the sample count is per tap times that.
append_noise_observation!(
    track_state,
    noise_observation_from_correlator(
        accumulated_power,
        num_taps,
        num_taps * samples_per_dump,
        sampling_frequency;
        prn = noise_prn,
        # the taps are simultaneous, not consecutive — they all integrate the
        # same samples, so the window must count that span once
        duration = samples_per_dump / sampling_frequency,
    ),
    :GPSL1CA,
)

estimate_dopplers_and_filter_prompt!(track_state, (L1 = sampling_frequency,))
```

If a signal has a noise estimator and nothing ever fills it, the fold **warns
once per signal** and every C/N₀ on it reads `-Inf dB-Hz`. That is the likeliest
integration mistake and the warning names the fix; it never fires on the
sample-driven path.

### How `append_noise_observation!` differs from `append_correlator_output!`

They are deliberately not the same mechanism:

|            | [`append_correlator_output!`](@ref)                                    | [`append_noise_observation!`](@ref)                    |
|:---------- |:---------------------------------------------------------------------- |:------------------------------------------------------ |
| scope      | one **signal** of one satellite                                        | one **signal**, across every satellite tracking it     |
| payload    | the whole correlator: every tap kept **separate**, complex, per antenna | the taps **pooled** — into one power sum, or into one spatial covariance for an array — plus `M`, the span and the PRN |
| drives     | discriminators, both loop filters, bit sync, the NCO update, the C/N₀ prompt | the noise floor only                              |
| timing     | `sample_index` is load-bearing — vector tracking needs a common grid    | none: `N₀` is slowly varying, so only "recent" matters |
| lifetime   | **drained** by the fold each chunk, buffer reused                      | **sliding window** bounded in time, never drained      |
| count      | one per completed coherent integration per signal                      | any granularity; `M` carries the weight                |

The payload row is the sharpest difference. Both traverse a multi-tap correlator,
for opposite reasons: the prompt path needs the taps kept apart because their
*differences* are the code and carrier discriminants, while the reference pools
them, because at ≥1 chip spacing they are three independent looks at the same
noise and nothing about their relative values means anything.

The **antenna** dimension is the one thing the reference does keep. With
`num_ants > 1` the pooled payload is the array's `M×M` spatial covariance
`R̂ = Σ_k b_k·b_kᴴ` rather than a scalar `Σ_k |b_k|²` — its diagonal is each
antenna's own `N₀` and its off-diagonals the antennas' noise correlation. That is
what a beamformer's C/N₀ needs: the floor its weights actually see is `wᴴR̂w`,
which no single number can answer. Every builder below takes the per-antenna form
of the input it already takes, so a hardware source reports all `M` elements
rather than collapsing them.

## Writing your own source

[`AbstractNoiseEstimator`](@ref) and its three methods are public, so a source
that is neither of the shipped paths — a front-end power monitor read over a
sideband, say — is a subtype away:

```julia
struct MyPowerMonitor <: Tracking.AbstractNoiseEstimator
    densities::Vector{typeof(1.0 / 1.0Hz)}
end

# fed from outside, so the sample-driven path is a no-op (that is the default,
# and it is why `update_noise!` need not be implemented at all here)
Tracking.append_noise_observation!(e::MyPowerMonitor, obs) =
    (push!(e.densities, obs.noise_density); e)
Tracking.get_noise_density(e::MyPowerMonitor) =
    isempty(e.densities) ? nothing : sum(e.densities) / length(e.densities)
```

Pass it in through `TrackState`'s `noise_estimators` keyword:

```julia
TrackState(; signal = GPSL1CA(), noise_estimators = (GPSL1CA = MyPowerMonitor(...),))
```

Two contracts to keep. The window must be mutated **in place** and the struct
returned unchanged — `TrackState` is immutable and never rebuilt for a noise
update. And [`Tracking.noise_density_type`](@ref) must name the concrete type
`get_noise_density` returns, so the fold can split off the `nothing` once per
signal per chunk and keep everything below it monomorphic; it defaults to
`typeof(1.0/1.0Hz)`, which every shipped builder produces.

```@docs
Tracking.noise_density_type
Tracking.NoiseUpdateContext
requires_noise_density
```
