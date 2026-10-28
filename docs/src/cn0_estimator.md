# CN0 Estimator

The CN0 (carrier-to-noise density ratio) estimator measures signal quality in dB-Hz.
Each [`TrackedSignal`](@ref) on a [`TrackedSat`](@ref) holds its own estimator, so a
multi-signal satellite reports one CN0 per signal.

!!! note "Every signal of a satellite is estimated, not just the driver"

    Each signal's estimator is fed its own prompt once per completed record; the
    [estimator-driver signal](tracking_state.md#Estimator-driver-signal) has no
    special role, and the values differ: on GPS L1C the pilot carries ~75 % of the
    power and the data component ~25 %, 4.8 dB apart.

    Nothing inside Tracking reads
    [`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{NoiseRefCN0Estimator, Any});
    the per-record update costs well under a per cent of the record's correlation
    work (~40 ns against ~9 µs for a two-signal satellite over 4000 samples). Switch
    a signal off with [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator) for
    clarity, not speed.

## Choosing an estimator

|                                                            | [`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator) (default)                | [`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator)                                                                                                                                                 | [`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator) | [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator) |
|:---------------------------------------------------------- |:------------------------------------------------------- |:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |:----------------------------- |:------------------------ |
| Needs the band's samples (or a noise observation)?         | **yes**                                                 | no                                                                                                                                                                         | no                            | no                       |
| Needs bit sync / a coherent window?                        | no                                                      | **yes**                                                                                                                                                                    | no                            | no                       |
| Needs to know the window length `M`?                       | no                                                      | **yes**                                                                                                                                                                    | no                            | no                       |
| Works on GPS L1C-D, Galileo E1B, pre-sync secondary codes? | **yes**                                                 | no — falls back                                                                                                                                                            | yes, at its floor             | n/a                      |
| Bias, thermal only                                         | <0.03 dB                                                | ≈+0.05 dB everywhere                                                                                                                                                       | ≈27.6 dB-Hz *floor* on noise  | n/a                      |
| High-C/N₀ bias floor                                       | ≈0.13 dB at 45, ≈0.40 at 50 (self-leakage)              | none                                                                                                                                                                       | none                          | n/a                      |
| High-C/N₀ σ                                                | keeps improving, until the row above dominates          | saturates at `√(M/((M−1)K))` ≈0.5 dB                                                                                                                                       | —                             | n/a                      |
| Degenerate outputs                                         | none                                                    | `-Inf`/`Inf` outside `1 < μ̂ < M` — 5.6 % of *estimates* at 20 dB-Hz                                                                                                       | none                          | always `-Inf`            |
| Phase-noise / residual-Doppler sensitivity                 | immune (non-coherent)                                   | sensitive — the coherent sum loses power; **−12 dB at 120 Hz of handoff Doppler**, see [Residual Doppler at acquisition handoff](#Residual-Doppler-at-acquisition-handoff) | immune                        | n/a                      |
| Cost                                                       | per-record arithmetic + one despread per signal per chunk | per-record arithmetic                                                                                                                                                      | per-record arithmetic         | none                     |

The "degenerate" row counts unusable *outputs*:
[`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{NoiseRefCN0Estimator, Any})
pools every buffered window and range-checks the pooled ratio once.

### Which one do I use?

 1. **Sample-driven [`track!`](@ref)** →
    [`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator), the default:
    the signal gets a
    [`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator) that
    `track!` fills.
 2. **Correlator outputs from hardware, *with* a noise observation** → the same, plus
    one [`append_noise_observation!`](@ref) per signal per fold (see
    [Noise Estimator](noise_estimator.md)).
 3. **Correlator outputs from hardware, *without* one** → `NWPRCN0Estimator(signal)`,
    accepting its fallback on the signals it cannot serve.
 4. **A signal whose C/N₀ nobody reads** →
    [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator); see
    [Not measuring a signal at all](#Not-measuring-a-signal-at-all).
 5. [`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator) is a fallback
    component and prompt-stream tool, not a standalone default; see
    [Why not the moment method](#Why-not-the-moment-method).

## The default: C/N₀ against a measured noise floor

[`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator) computes, per record,

```
Ĉ/N₀ = ⟨|P|²⟩ / N̂₀ − 1/T
```

with `N̂₀` the measured noise **density** (see [Noise Estimator](noise_estimator.md))
and `T` the record's integration time, and averages the per-record terms. It is
**non-coherent**: no bit sync, no window, no fallback, the same behaviour on every
signal.

In the TrackingLoops manual:

- [`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator)
- [`default_cn0_estimator`](@extref TrackingLoops.default_cn0_estimator-Tuple{AbstractGNSSSignal, Int64})

### Getting the low-C/N₀ variance

Below about 28 dB-Hz a non-coherent estimator is the noisier one at matched record
length (squaring loss: NWPR sums `M` prompts coherently first, better by `√(M−1)`).
But a coherent sum of five 1 ms prompts *is* one 5 ms prompt, so take that
integration in the correlator instead, with
[`set_preferred_num_code_blocks_to_integrate!`](@ref): the noise reference then beats
NWPR at every C/N₀. Modelled over 100 ms, 4000 trials × 9 seeds, relative σ in dB:

| C/N₀ | NWPR (M = 5) | NoiseRef, 1 ms records | NoiseRef, 5 ms records |
|:---- | ------------:| ----------------------:| ----------------------:|
| 20   | 2.96         | 4.79                   | **2.76**               |
| 25   | 1.50         | 1.72                   | **1.24**               |
| 30   | 0.90         | 0.75                   | **0.65**               |
| 40   | 0.54         | 0.20                   | **0.20**               |
| 50   | 0.50         | 0.06                   | **0.06**               |

This is a trade: a longer record also narrows the loop (the carrier bandwidth is
capped at `0.09 / T_int`), and before bit/secondary sync records are one block
whatever you configure.

### The one bias it carries

Despreading with a wrong PRN also collects the tracked satellite's own power,
`ε_self/N₀ = C/f_chip`, so on GPS L1 C/A C/N₀ reads low by 0.013 dB at 35 dB-Hz,
0.042 at 40, **0.13 at 45** and **0.40 at 50**; ten times less at 10.23 Mcps (GPS L5,
Galileo E5a). It is deliberately not corrected, as that would couple every
satellite's C/N₀ to the others' (see
[`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator)).

## `NWPRCN0Estimator`: the coherent-window estimator

Van Dierendonck's **narrowband/wideband power ratio**: per window of `M` records the
coherent power `NBP = |Σ prompt|²` against the incoherent `WBP = Σ |prompt|²`, summed
over the windows in `num_prompts_for_cn0_estimation` records. It is the estimator for
a correlator-ingest path that cannot report a noise observation.

In the TrackingLoops manual:

- [`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator)
- [`NWPRCN0Estimator(::AbstractGNSSSignal)`](@extref TrackingLoops.NWPRCN0Estimator-Tuple{AbstractGNSSSignal})

### Why not the moment method

The [Moments Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1)
(M2M4, [`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator)) is a
moment ratio, and at a finite window its sample moments manufacture signal power out
of noise. On synthetic prompts (amplitude `√(C/N₀·T)` in unit-variance complex
noise), 100 one-millisecond records, 400 realizations:

| true C/N₀  | M2M4 median | p10  | p90  | NWPR (`M = 20`) median |
|:---------- | -----------:| ----:| ----:| ----------------------:|
| noise only | **27.6**    | 22.6 | 31.1 | **−Inf** ("no signal") |
| 20 dB-Hz   | 27.7        | 22.5 | 30.9 | 19.9                   |
| 25 dB-Hz   | 27.9        | 23.4 | 30.9 | 25.0                   |
| 30 dB-Hz   | 30.5        | 27.3 | 32.6 | 30.0                   |
| 35 dB-Hz   | 35.3        | 33.9 | 36.4 | 35.0                   |
| 45 dB-Hz   | 45.1        | 44.3 | 45.9 | 45.0                   |

M2M4 is blind below ~30 dB-Hz: pure noise clears a 30 dB-Hz threshold 19 % of the
time, 28 dB-Hz 45 % and 25 dB-Hz 76 %, so the usual code-lock hysteresis (lock at 34,
drop at 25) cannot be built on it. NWPR at `M = 20` has 0 % false alarms from 20 to
32 dB-Hz (Falletti, Pini & Lo Presti, *IEEE T-AES* 47(1):420–437, 2011;
[issue #217](https://github.com/JuliaGNSS/Tracking.jl/issues/217)).

End to end through `track` on a data-modulated GPS L1 C/A signal at the (shorter)
default windows, 600 ms per run, medians over 16 seeds:

| true C/N₀  | NWPR (default)            | p10  | p90  | M2M4 |
|:---------- | -------------------------:| ----:| ----:| ----:|
| noise only | **−Inf** in 12 of 16 runs |      |      | 27.9 |
| 20 dB-Hz   | 18.3                      | 6.0  | 22.2 | 27.6 |
| 25 dB-Hz   | 23.9                      | 16.3 | 26.3 | 28.3 |
| 30 dB-Hz   | 29.0                      | 27.4 | 30.3 | 30.8 |
| 35 dB-Hz   | 34.1                      | 33.3 | 35.0 | 35.3 |
| 45 dB-Hz   | 44.3                      | 43.7 | 45.2 | 45.4 |

On pure noise pre-sync, where a code-lock detector decides, the default clears a
20 dB-Hz threshold in 4 % of updates and 25 dB-Hz in **0 %** (200 runs, highest value
22.4); M2M4 on the same runs clears 25 dB-Hz in 67 %. The price is the coherent sum's
sensitivity to anything that decorrelates the prompts (see
[What the estimator returns](#What-the-estimator-returns)).

### The window follows the navigation bits

A window must not straddle a bit flip (~7 dB), so each prompt comes with the
navigation-bit state in a
[`CN0UpdateContext`](@extref TrackingLoops.CN0UpdateContext) and the window follows
it; the cases are tabulated under
[`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator). Before sync, a
data-bearing signal without a secondary code uses an unaligned
`num_presync_narrowband_code_blocks` window (5), since bit sync takes seconds and
fails below ~30 dB-Hz; it costs ~0.6 dB below 30 dB-Hz on GPS L1 C/A. Where no window
is admissible (a secondary code before sync, one symbol per code block, or a record
as long as its window, e.g. with
[`set_preferred_num_code_blocks_to_integrate!`](@ref) at a whole bit) the
[`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator) fallback, and its
floor, is reported.

!!! warning "The moment ratio's noise floor is there on data-only signals"

    GPS L1C-D and Galileo E1B carry **one navigation symbol per code block**, so the
    longest coherent window is one record, where `NBP == WBP`. They keep the moment
    ratio for good: **pure noise reads ~27.6 dB-Hz on L1C-D and E1B**, as does any
    signal with a secondary code (GPS L5, L2C, Galileo E1C, E5a) before its sync.

    No ~100-record window fixes this: without coherence only the prompt magnitudes
    carry information. (A phase-blind `(E|z|)²/E|z|²` estimator cuts false alarms,
    42 % against 76 % at 25 dB-Hz on pure noise, but its median on noise is still
    ~27.5 dB-Hz.)

    The default [`NoiseRefCN0Estimator`](@extref TrackingLoops.NoiseRefCN0Estimator)
    has no window, so these are ordinary records to it
    ([issue #217](https://github.com/JuliaGNSS/Tracking.jl/issues/217)). On NWPR,
    take the lock decision from the **pilot** of the pair (L1C-D with L1C-P, E1B
    with E1C), give the data component a
    [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator), or use
    `NWPRCN0Estimator(signal; fallback = NoCN0Estimator())`, which reports
    `-Inf dB-Hz` exactly when no coherent window is available.

### The window is capped by the loop's coherence time

The windows tile the bit from its start but stay shorter than a bit
(`num_narrowband_code_blocks`, ~5 ms): a longer window buys little, while
decoherence over it costs a bias averaging cannot remove. The default staged loop
has narrowed to 18 Hz by the time a locked satellite fades (see
[Loop Filter](@ref)) and holds phase over a whole bit, so for it the window length
hardly matters; the short default protects a loop that has not narrowed or runs
wider. Raise `num_narrowband_code_blocks` for a pilot, a narrow carrier loop or a
signal that is never weak. Measurements are under
[`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator).

## Configuring the estimator

`num_prompts_for_cn0_estimation` (default `100`, ~100 ms for L1 C/A) sets how many
records the estimate averages over (also NWPR's fallback):

```jldoctest cn0_buffer
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; num_prompts_for_cn0_estimation = 200);

julia> track_state = add_satellite!(track_state, :default, sat);

julia> get_prn(track_state, 1)
1
```

To choose a different estimator, or configure NWPR's windows, pass a
`cn0_estimator`; any
[`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator) is stored as is
(it is a type parameter of [`TrackedSignal`](@ref)):

```jldoctest cn0_estimator_kwarg
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; cn0_estimator = MomentsCN0Estimator(100));

julia> Tracking.get_cn0_estimator(sat) isa MomentsCN0Estimator
true

julia> sat = TrackedSat(
           GPSL1CA(),
           1,
           50.0,
           1000.0Hz;
           cn0_estimator = NWPRCN0Estimator(
               GPSL1CA();
               num_records = 400,                       # ~400 ms of memory
               num_presync_narrowband_code_blocks = 0,  # M2M4 until bit sync
           ),
       );

julia> Tracking.get_cn0_estimator(sat).num_records
400
```

Passing the **signal** to `NWPRCN0Estimator` sizes its window for the code period:
5 blocks for a 1 ms code, 2 for GPS L1C-P's 10 ms one (the bare constructor's `5`
would be 50 ms there).

Each signal needs its **own** estimator instance (they buffer into a vector, so a
shared one corrupts both); the multi-signal constructor takes one per signal:

```julia
TrackedSat(
    (GPSL1C_P(), GPSL1CA()),
    prn,
    code_phase,
    carrier_doppler;
    cn0_estimator = (NoiseRefCN0Estimator(), MomentsCN0Estimator(100)),
)
```

Only a signal whose estimator reads a density is provisioned a noise source, so the
case above runs one despread per chunk, for the noise-referenced signal; a signal on
NWPR runs none. The floor is per signal, not per band (see
[Why per signal and not per RF band](noise_estimator.md#Why-per-signal-and-not-per-RF-band)),
so the cost is one despread per noise-referenced signal per chunk, against `n_sats`
tracking correlations apiece.

### Not measuring a signal at all

[`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator) keeps no state, does no
per-record work and reports `-Inf dB-Hz` (not `NaN`, which compares `>=` true against
every threshold). Two uses:

```julia
# the non-driver signal of a pair whose C/N₀ nobody reads
TrackedSat(
    (GPSL1C_P(), GPSL1C_D()),
    prn,
    code_phase,
    carrier_doppler;
    cn0_estimator = (NoiseRefCN0Estimator(), NoCN0Estimator()),
)

# NWPR with no moment-ratio floor underneath it: "no estimate" until a coherent
# window exists, instead of the moment ratio's ~27.6 dB-Hz on noise
NWPRCN0Estimator(GPSL1C_D(); fallback = NoCN0Estimator())
```

The second only matters on NWPR (the default has no fallback): it turns the
fallback's noise floor into `-Inf` wherever no coherent window exists (see the
warning above).

In the TrackingLoops manual:

- [`MomentsCN0Estimator`](@extref TrackingLoops.MomentsCN0Estimator)
- [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator)

## What the estimator returns

  - **Before the first record** (e.g. right after `add_satellite!`): `-Inf dB-Hz`, so
    a lock detector cannot clear on it.
  - **A noise-free signal**: nothing to divide by, so a very large value; NWPR and
    the moment estimator report `Inf dB-Hz` (`μ̂ = M`). The
    [Quick start](index.md#Quick-start) shows this.
  - **A noisy signal**: converges to the true CN0 as the buffers fill; real signals
    typically land in 30–50 dB-Hz.
  - **Noise only**: the default reports `-Inf dB-Hz` once the averaged terms fail to
    clear zero, which is most of the time, and has no positive floor; NWPR reports
    `-Inf` when `μ̂ ≤ 1`. That is the median, not a guarantee: code that
    **averages** the estimate must handle the infinities.
  - **While the loops struggle**: NWPR measures *coherent* C/N₀, so phase noise
    counts against it: at 1 ms records and default windows ~0.9 dB low at a true
    35 dB-Hz, ~1.1 dB at 25 and ~1.7 dB at 20, much worse with
    `num_narrowband_code_blocks` past the loop's coherence time. The default is
    **immune**: through `track!` at a true 45 dB-Hz it reads 44.9 at 1-block and 44.7
    at 20-block records, where NWPR falls to 32.9.

All estimators report C/N₀ against the **effective** noise density: the other
satellites' cross-correlation power `ε` is seen by NWPR and the noise reference
alike, so both report `C/(N₀+ε)`, as other receivers do. Only the reference's
self-leakage differs (see [The one bias it carries](#The-one-bias-it-carries)).

!!! note "An antenna array's floor follows each satellite's own beamformer"

    The prompt is post-[`AbstractPostCorrFilter`](@extref TrackingLoops.AbstractPostCorrFilter),
    so the floor must go through the same weights. With `num_ants > 1` the noise
    reference's window (one per signal) holds the `M×M` spatial noise covariance
    `R̂`, and each satellite reduces it at C/N₀ time through its own current weights,

    ```
    N̂₀ = wᴴ R̂ w
    ```

    exact for any fixed `w` since `E[|wᴴn|²] = wᴴRw`. Hence
    [`get_weights`](@extref TrackingLoops.get_weights) is required of a post-correlation
    filter, and the contract is **linear in the antennas**.
    [`DefaultPostCorrFilter`](@extref TrackingLoops.DefaultPostCorrFilter) selects
    the `last` antenna, whose diagonal entry `wᴴR̂w` then is, as in a single-antenna
    run on that column.

### Residual Doppler at acquisition handoff

Acquisition hands over a Doppler good to a search bin (commonly ±100 Hz), which the
loop takes tens of milliseconds to pull in, just when a lock detector decides
whether to keep the satellite. Through `track!` at a true **40 dB-Hz**, 1 ms records,
seeded with a carrier-Doppler error, median of 5 seeds over the first 60 ms:

| Δf     | real within-record loss | NWPR     | Moments | `NoiseRef` |
|:------ | -----------------------:| --------:| -------:| ----------:|
| 0 Hz   | 0.00 dB                 | 39.8     | 40.5    | **40.1**   |
| 30 Hz  | −0.01 dB                | 39.8     | 40.4    | **40.1**   |
| 60 Hz  | −0.05 dB                | 39.2     | 40.4    | **40.1**   |
| 120 Hz | −0.21 dB                | **27.6** | 40.5    | **39.9**   |
| 250 Hz | −0.91 dB                | **−Inf** | 40.1    | **39.1**   |

At 120 Hz NWPR under-reports by **12 dB**, and at 250 Hz reports `-Inf dB-Hz` on a
satellite being tracked perfectly well. Across the window a residual `Δf` turns the
prompt by `Δf·M·T` **cycles**: 0.6 at 120 Hz over the default 5 ms, spreading the
five phasors over 216° and costing `NBP` ~5.9 dB; `μ̂` collapses toward 1, and past
1.25 cycles `μ̂ ≤ 1` gives `-Inf`. This is the constraint `M·T ≪ 1/(2·Δf)` under
[`NWPRCN0Estimator`](@extref TrackingLoops.NWPRCN0Estimator), which at the default window means `Δf ≪ 100 Hz`.

The non-coherent estimators see only the real **within-record** loss,
`sinc²(Δf·T)`: 0.21 dB at 120 Hz, 0.91 dB at 250 Hz (40 − 0.91 ≈ the measured 39.1).
The open-loop noise reference cannot be perturbed by any satellite's Doppler error.

After 300 ms NWPR recovers to 39.2 dB-Hz at every offset up to 120 Hz. At 250 Hz it
stays at `-Inf` and the noise reference reads 38.9: the loop never pulled in, so the
residual and its 0.9 dB loss remain.

On NWPR, a shorter pre-sync window helps, though a two-record window is also the
noisiest:

| `num_presync_narrowband_code_blocks` | `Δf·M·T` at 120 Hz | reported |
|:------------------------------------ | ------------------:| --------:|
| 5 (default)                          | 0.60 cycles        | 27.6     |
| 3                                    | 0.36 cycles        | 32.4     |
| 2                                    | 0.24 cycles        | 34.9     |

```julia
NWPRCN0Estimator(GPSL1CA(); num_presync_narrowband_code_blocks = 2)
```

A noisy L1 C/A signal at a known 45 dB-Hz, through 25 1-ms tracking cycles:

```jldoctest cn0_example; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals, Random

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_chip_rate

julia> function run_cn0_demo()
           sys = GPSL1CA()
           fs = 4e6Hz
           num_samples = 4000  # 1 ms at 4 MHz
           prn = 1
           cn0_db_hz = 45.0
           sigma = sqrt((fs/Hz) / 10^(cn0_db_hz/10) / 2)
           code_freq = get_chip_rate(sys)
           rng = MersenneTwister(0)
           track_state = TrackState(; signal = GPSL1CA())
           track_state =
               add_satellite!(track_state; prn, code_phase = 0.0, carrier_doppler = 0.0Hz)
           for _ = 1:25
               clean = gen_code(num_samples, sys, prn, fs, code_freq, 0.0)
               noise = sigma .* randn(rng, ComplexF64, num_samples)
               track_state = track(clean .+ noise, track_state, fs)
           end
           estimate_cn0(track_state, 1)
       end;

julia> run_cn0_demo()
47.9 dB-Hz
```

With unit-amplitude code samples, `C/N₀ = 1 / (σ² / Fs)` gives
`σ² = Fs / 10^(CN0_dB/10)`, split across real and imaginary parts. 25 records fill a
quarter of the default 100-record window, so a few dB of spread is expected.

In the TrackingLoops manual:

- [`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{NoiseRefCN0Estimator, Any})

## Custom CN0 Estimators

In the TrackingLoops manual:

- [`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator)
- [`CN0UpdateContext`](@extref TrackingLoops.CN0UpdateContext)

Subtype [`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator) and
implement:

  - `TrackingLoops.update(cn0_estimator::MyCN0Estimator, prompt)`, returning a new
    estimator with the prompt added;
  - `estimate_cn0(cn0_estimator::MyCN0Estimator, integration_time)`, returning
    `dB-Hz`; `integration_time` is the last record's, i.e.
    [`get_last_fully_integrated_integration_time`](@ref).

To use the navigation-bit grid (boundaries, sync, soft bits), implement
`TrackingLoops.update(cn0_estimator::MyCN0Estimator, prompt, context::CN0UpdateContext)`
instead; the tracking loop calls this form, whose default drops the context.

In the TrackingLoops manual:

- [`update(::MomentsCN0Estimator, ::Any)`](@extref TrackingLoops.update-Tuple{MomentsCN0Estimator, Any})
- [`update(::AbstractCN0Estimator, ::Any, ::CN0UpdateContext)`](@extref TrackingLoops.update-Tuple{TrackingLoops.AbstractCN0Estimator, Any, TrackingLoops.CN0UpdateContext})
- [`update(::NWPRCN0Estimator, ::Any, ::CN0UpdateContext)`](@extref TrackingLoops.update-Tuple{NWPRCN0Estimator, Any, TrackingLoops.CN0UpdateContext})
- [`update(::NWPRCN0Estimator, ::Any)`](@extref TrackingLoops.update-Tuple{NWPRCN0Estimator, Any})
- [`update(::NoiseRefCN0Estimator, ::Any, ::CN0UpdateContext)`](@extref TrackingLoops.update-Tuple{NoiseRefCN0Estimator, Any, TrackingLoops.CN0UpdateContext})
- [`update(::NoiseRefCN0Estimator, ::Any)`](@extref TrackingLoops.update-Tuple{NoiseRefCN0Estimator, Any})
- [`update(::NoCN0Estimator, ::Any)`](@extref TrackingLoops.update-Tuple{NoCN0Estimator, Any})

Plug it in with the `cn0_estimator` keyword of [`TrackedSignal`](@ref) or
[`TrackedSat`](@ref):

```jldoctest custom_cn0
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> struct MyCN0Estimator <: TrackingLoops.AbstractCN0Estimator end

julia> TrackingLoops.update(estimator::MyCN0Estimator, prompt) = estimator;

julia> TrackingLoops.estimate_cn0(::MyCN0Estimator, integration_time) = 42.0 * Tracking.dBHz;

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; cn0_estimator = MyCN0Estimator());

julia> estimate_cn0(sat)
42.0 dB-Hz
```
