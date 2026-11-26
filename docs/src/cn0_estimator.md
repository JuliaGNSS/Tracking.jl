# CN0 Estimator

The CN0 (carrier-to-noise density ratio) estimator measures signal quality in dB-Hz.
Each [`TrackedSignal`](@ref) holds its own estimator, fed its own prompt once per
completed integration, so a multi-signal satellite reports one CN0 per signal; the
[driver signal](tracking_state.md#Driver-signal) has no special role (on GPS L1C the
pilot carries ~75 % of the power and the data component ~25 %, 4.8 dB apart). Nothing
inside Tracking reads
[`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{NoiseReferencedCN0Estimator, Any}),
and the update costs ~40 ns against ~9 µs of correlation for a two-signal satellite
over 4000 samples, so switch a signal off with
[`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator) for clarity, not speed.

## Choosing an estimator

 1. **Sample-driven [`track!`](@ref)** →
    [`NoiseReferencedCN0Estimator`](@extref TrackingLoops.NoiseReferencedCN0Estimator), the default:
    the signal gets a
    [`CorrelatorNoiseDensityEstimator`](@extref TrackingLoops.CorrelatorNoiseDensityEstimator) that
    `track!` fills.
 2. **Correlator outputs from hardware, *with* a noise observation** → the same, plus
    one [`append_noise_observation!`](@ref) per signal per loop update (see
    [Noise Estimator](noise_estimator.md)).
 3. **Correlator outputs from hardware, *without* one** → `NWPRCN0Estimator(signal)`,
    accepting its fallback on the signals it cannot serve.
 4. **A signal whose C/N₀ nobody reads** →
    [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator); see
    [Not measuring a signal at all](#Not-measuring-a-signal-at-all).
 5. [`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator) is NWPR's
    fallback, not a standalone choice; see
    [Why not the moment method](#Why-not-the-moment-method).

|                                     | NoiseReferenced (default)                 | NWPR                                                  | Moments                  |
|:----------------------------------- |:----------------------------------------- |:----------------------------------------------------- |:------------------------ |
| needs                               | the band's samples or a noise observation | a coherent window, so bit sync                        | nothing                  |
| GPS L1C-D, Galileo E1B, pre-sync secondary codes | yes                          | no, falls back                                        | yes, reading noise as signal |
| bias                                | ≈0.13 dB at 45, ≈0.40 at 50 dB-Hz (self-leakage) | ≈+0.05 dB                                       | 20–30 dB-Hz on noise half the time |
| high-C/N₀ σ                         | keeps improving                           | saturates at `√(M/((M−1)K))` ≈0.5 dB                  | —                        |
| degenerate outputs                  | none                                      | `-Inf`/`Inf` outside `1 < μ̂ < M`: 5.6 % at 20 dB-Hz  | none                     |
| residual Doppler                    | immune (noncoherent)                     | **−12 dB at 120 Hz** (see [below](#Residual-Doppler-at-acquisition-handoff)) | immune |
| extra cost                          | one despread per signal per NCO update span | none                                                | none                     |

## The default: C/N₀ against a measured noise floor

[`NoiseReferencedCN0Estimator`](@extref TrackingLoops.NoiseReferencedCN0Estimator)
averages `|P|²/N̂₀ − 1/T` over the integrations, with `N̂₀` the measured noise density
(see [Noise Estimator](noise_estimator.md)). It averages the prompt powers
noncoherently: no bit sync, no window, no fallback, the same on every signal.

Below about 28 dB-Hz a noncoherent estimator is the noisier one at matched coherent
integration time (NWPR sums `M` prompts coherently first). A coherent sum of five
1 ms prompts *is* one 5 ms prompt, though, so lengthen the coherent integration with
[`set_preferred_num_code_periods_to_integrate!`](@ref) and the noise reference beats
NWPR at every C/N₀. Relative σ in dB over 100 ms, 4000 trials × 9 seeds:

| C/N₀ | NWPR (M = 5) | NoiseRef, 1 ms coherent | NoiseRef, 5 ms coherent |
|:---- | ------------:| ----------------------:| ----------------------:|
| 20   | 2.96         | 4.79                   | **2.76**               |
| 25   | 1.50         | 1.72                   | **1.24**               |
| 30   | 0.90         | 0.75                   | **0.65**               |
| 40   | 0.54         | 0.20                   | **0.20**               |
| 50   | 0.50         | 0.06                   | **0.06**               |

A longer coherent integration also narrows the loop (the carrier bandwidth is capped
at `0.09 / T`), and before bit or secondary sync integrations are one code period
whatever you configure.

The one bias: despreading with a wrong PRN also collects the satellite's own power,
`C/f_chip`, so on GPS L1 C/A C/N₀ reads low by 0.013 dB at 35 dB-Hz, 0.13 at 45 and
0.40 at 50, ten times less at 10.23 Mcps. Correcting it would couple every satellite's
C/N₀ to the others'.

## `NWPRCN0Estimator`: the coherent-window estimator

Van Dierendonck's narrowband/wideband power ratio over windows of `M` integrations, for
a correlator-ingest path that cannot report a noise observation. Its formulas, its
window per sync state and its coherence constraint are documented under
[`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator).

### Why not the moment method

The moments method ([`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator))
manufactures signal power out of noise at a finite window: at its default 100
integrations pure noise reads "no signal" only ~40 % of the time and otherwise mostly
20–30 dB-Hz, so ~48 % of noise-only estimates clear 25 dB-Hz, and a true 20 dB-Hz reads
the same. A code-lock detector cannot be built on that
([issue #217](https://github.com/JuliaGNSS/Tracking.jl/issues/217)).

NWPR through `track` on a data-modulated GPS L1 C/A signal at the default windows,
600 ms per run, medians over 16 seeds:

| true C/N₀  | NWPR (default)            | p10  | p90  |
|:---------- | -------------------------:| ----:| ----:|
| noise only | **−Inf** in 12 of 16 runs |      |      |
| 20 dB-Hz   | 18.3                      | 6.0  | 22.2 |
| 25 dB-Hz   | 23.9                      | 16.3 | 26.3 |
| 30 dB-Hz   | 29.0                      | 27.4 | 30.3 |
| 35 dB-Hz   | 34.1                      | 33.3 | 35.0 |
| 45 dB-Hz   | 44.3                      | 43.7 | 45.2 |

On pure noise before sync it clears 20 dB-Hz in 4 % of updates and 25 dB-Hz in
**0 %** (200 runs, highest value 22.4).

### Where NWPR has no window

Before sync a data signal without a secondary code uses an unaligned
`num_presync_narrowband_code_periods` window (5), as bit sync takes seconds and fails
below ~30 dB-Hz; it costs ~0.6 dB below 30 dB-Hz on GPS L1 C/A. Where no window is
admissible, NWPR reports its moments fallback: a secondary code before sync, an
integration as long as its window, and for good on GPS L1C-D and Galileo E1B, which
carry one symbol per code period.

!!! warning "On data-only signals NWPR reads noise as signal"

    On GPS L1C-D and Galileo E1B pure noise reads 20–30 dB-Hz about half the time, as
    on any signal with a secondary code (GPS L5, L2C, Galileo E1C, E5a) before its sync.
    Use the default estimator, which has no window; or on NWPR take the lock decision
    from the pilot of the pair (L1C-P, E1C) and give the data component a
    [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator), or use
    `NWPRCN0Estimator(signal; fallback = NoCN0Estimator())`, which reports
    `-Inf dB-Hz` exactly when no coherent window is available.

The windows (`num_narrowband_code_periods`, ~5 ms) stay shorter than a bit: decoherence
over a longer one costs a bias averaging cannot remove. The default staged loop has
narrowed to 18 Hz by the time a locked satellite fades (see [Loop Filter](@ref)), so
for it the length hardly matters; raise it for a pilot, a narrow carrier loop or a
signal that is never weak.

## Configuring the estimator

The default estimator averages 100 integrations (~100 ms for L1 C/A). Its
`num_averaged_integrations` sets another span:

```jldoctest cn0_buffer
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signal = GPSL1CA(),
           cn0_estimator = NoiseReferencedCN0Estimator(; num_averaged_integrations = 200),
       );

julia> sat = TrackedSat(
           GPSL1CA(),
           1,
           50.0,
           1000.0Hz;
           cn0_estimator = NoiseReferencedCN0Estimator(; num_averaged_integrations = 200),
       );

julia> track_state = add_satellite!(track_state, :default, sat);

julia> get_prn(track_state, 1)
1
```

Pass a `cn0_estimator` to choose another estimator or configure NWPR's windows; any
[`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator) is stored as is
(it is a type parameter of [`TrackedSignal`](@ref)):

```jldoctest cn0_estimator_kwarg
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; cn0_estimator = MomentsCN0Estimator(; num_averaged_integrations = 100));

julia> Tracking.get_cn0_estimator(sat) isa MomentsCN0Estimator
true

julia> sat = TrackedSat(
           GPSL1CA(),
           1,
           50.0,
           1000.0Hz;
           cn0_estimator = NWPRCN0Estimator(
               GPSL1CA();
               num_averaged_integrations = 400,          # ~400 ms of memory
               num_presync_narrowband_code_periods = 0,  # M2M4 until bit sync
           ),
       );

julia> Tracking.get_cn0_estimator(sat).num_averaged_integrations
400
```

`NWPRCN0Estimator(signal)` sizes its window for the code period: 5 code periods for a
1 ms code, 2 for GPS L1C-P's 10 ms one.

Each signal needs its **own** estimator instance (a shared one corrupts both); the
multi-signal constructor takes one per signal:

```julia
TrackedSat(
    (GPSL1C_P(), GPSL1CA()),
    prn,
    code_phase,
    carrier_doppler;
    cn0_estimator = (NoiseReferencedCN0Estimator(), MomentsCN0Estimator(; num_averaged_integrations = 100)),
)
```

Only a signal whose estimator reads a density gets a noise source, so this runs one
despread per NCO update span, for the noise-referenced signal.

### Not measuring a signal at all

[`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator) keeps no state, does no work
and reports `-Inf dB-Hz` (not `NaN`, which compares `>=` true against every threshold):

```julia
# the non-driver signal of a pair whose C/N₀ nobody reads
TrackedSat(
    (GPSL1C_P(), GPSL1C_D()),
    prn,
    code_phase,
    carrier_doppler;
    cn0_estimator = (NoiseReferencedCN0Estimator(), NoCN0Estimator()),
)

# NWPR reporting "no estimate" until a coherent window exists, instead of its
# moments fallback's reading of noise
NWPRCN0Estimator(GPSL1C_D(); fallback = NoCN0Estimator())
```

## What the estimator returns

  - **Before the first integration** (e.g. right after `add_satellite!`): `-Inf dB-Hz`,
    so a lock detector cannot clear on it.
  - **A noise-free signal**: `Inf dB-Hz` from NWPR and the moments estimator (`μ̂ = M`),
    a very large value from the default; the [Quick start](index.md#Quick-start) shows
    this.
  - **Noise only**: the default reports `-Inf dB-Hz` once the averaged terms fail to
    clear zero, most of the time, NWPR when `μ̂ ≤ 1`. Code that **averages** the
    estimate must handle the infinities.
  - **While the loops struggle**: NWPR measures *coherent* C/N₀, so phase noise counts
    against it, ~0.9 dB low at a true 35 dB-Hz, ~1.7 dB at 20 at 1 ms integrations and
    the default windows. The default is immune: at a true 45 dB-Hz it reads 44.9 at 1
    and 44.7 at 20 code periods per integration, where NWPR falls to 32.9.

Both NWPR and the noise reference see the other satellites' cross-correlation power
`ε`, so both report `C/(N₀+ε)`, as other receivers do.

With an antenna array the prompt is post-[`AbstractPostCorrFilter`](@extref TrackingLoops.AbstractPostCorrFilter),
so the noise reference holds the `M×M` spatial noise covariance `R̂` and each satellite
reduces it through its own current weights, `N̂₀ = wᴴR̂w`, exact for any fixed `w`. A
post-correlation filter must therefore provide
[`get_weights`](@extref TrackingLoops.get_weights) and be linear in the antennas;
[`DefaultPostCorrFilter`](@extref TrackingLoops.DefaultPostCorrFilter) selects the
`last` antenna.

### Residual Doppler at acquisition handoff

Acquisition hands over a Doppler good to a search bin (commonly ±100 Hz), which the loop
takes tens of milliseconds to pull in, just when a lock detector decides whether to
keep the satellite. Through `track!` at a true **40 dB-Hz**, 1 ms integrations, median
of 5 seeds over the first 60 ms:

| Δf     | within-integration loss | NWPR     | Moments | `NoiseRef` |
|:------ | -----------------------:| --------:| -------:| ----------:|
| 0 Hz   | 0.00 dB                 | 39.8     | 40.5    | **40.1**   |
| 30 Hz  | −0.01 dB                | 39.8     | 40.4    | **40.1**   |
| 60 Hz  | −0.05 dB                | 39.2     | 40.4    | **40.1**   |
| 120 Hz | −0.21 dB                | **27.6** | 40.5    | **39.9**   |
| 250 Hz | −0.91 dB                | **−Inf** | 40.1    | **39.1**   |

Across a window the residual turns the prompt by `Δf·M·T` cycles, 0.6 at 120 Hz over
the default 5 ms, which costs `NBP` ~5.9 dB; past 1.25 cycles `μ̂ ≤ 1` gives `-Inf`. The
noncoherent estimators see only the within-integration loss `sinc²(Δf·T)`. After
300 ms NWPR recovers to 39.2 dB-Hz up to 120 Hz; at 250 Hz the loop never pulls in.
A shorter pre-sync window helps NWPR, at more noise:

| `num_presync_narrowband_code_periods` | `Δf·M·T` at 120 Hz | reported |
|:------------------------------------ | ------------------:| --------:|
| 5 (default)                          | 0.60 cycles        | 27.6     |
| 3                                    | 0.36 cycles        | 32.4     |
| 2                                    | 0.24 cycles        | 34.9     |

A noisy L1 C/A signal at a known 45 dB-Hz, through 25 1-ms tracking cycles:

```jldoctest cn0_example; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals, Random

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_chip_rate

julia> function run_cn0_demo()
           signal = GPSL1CA()
           fs = 4e6Hz
           num_samples = 4000  # 1 ms at 4 MHz
           prn = 1
           cn0_db_hz = 45.0
           sigma = sqrt((fs/Hz) / 10^(cn0_db_hz/10) / 2)
           code_freq = get_chip_rate(signal)
           rng = MersenneTwister(0)
           track_state = TrackState(; signal)
           track_state =
               add_satellite!(track_state; prn, code_phase = 0.0, carrier_doppler = 0.0Hz)
           for _ = 1:25
               clean = gen_code(num_samples, signal, prn, fs, code_freq, 0.0)
               noise = sigma .* randn(rng, ComplexF64, num_samples)
               track_state = track(clean .+ noise, track_state, fs)
           end
           estimate_cn0(track_state, 1)
       end;

julia> run_cn0_demo()
47.9 dB-Hz
```

With unit-amplitude code samples `C/N₀ = Fs / σ²`, split across real and imaginary
parts; 25 integrations fill a quarter of the default window, so a few dB of spread is
expected.

## Custom CN0 Estimators

Subtype [`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator) and
implement:

  - `TrackingLoops.update!(cn0_estimator::MyCN0Estimator, prompt)`, returning a new
    estimator with the prompt added;
  - `estimate_cn0(cn0_estimator::MyCN0Estimator, integration_time)`, returning
    `dB-Hz`; `integration_time` is the last integration's, i.e.
    [`get_last_coherent_integration_time`](@ref).

To use the navigation-bit grid, implement
`TrackingLoops.update!(cn0_estimator::MyCN0Estimator, prompt, context::CN0UpdateContext)`
instead (see [`CN0UpdateContext`](@extref TrackingLoops.CN0UpdateContext));
post-correlation processing calls this form, whose default drops the context.

Plug it in with the `cn0_estimator` keyword of [`TrackedSignal`](@ref) or
[`TrackedSat`](@ref):

```jldoctest custom_cn0
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> struct MyCN0Estimator <: TrackingLoops.AbstractCN0Estimator end

julia> TrackingLoops.update!(estimator::MyCN0Estimator, prompt) = estimator;

julia> TrackingLoops.estimate_cn0(::MyCN0Estimator, integration_time) = 42.0 * Tracking.dBHz;

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; cn0_estimator = MyCN0Estimator());

julia> estimate_cn0(sat)
42.0 dB-Hz
```
