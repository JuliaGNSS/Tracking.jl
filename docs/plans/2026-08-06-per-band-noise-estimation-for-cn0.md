# Per-band noise estimation for C/N₀

**Date:** 2026-08-06
**Status:** **Implemented** 2026-08-12 — Phases 0 through 3 all landed. See
"Where the implementation diverged" below for the four places it did not follow
this document, and why.

Earlier: Planned (Phase 0 implemented; single noise estimator as of 2026-08-06;
review pass 2026-08-09 — self-leakage promoted from a parenthetical to a quantified
risk, the "no density" skip respecified, anchors corrected against `b9709191`)

## Where the implementation diverged

Four places, each because building it surfaced something this document had not.

1. **`NoiseObservation` stores a `duration`, not `total_samples`.** The window is
   bounded in time, and the plan's own trim rule (`Σ total_samplesᵢ >
   window_duration · f_s`) needs a sampling frequency the observation does not
   carry. Worse, the field is not even a proxy for span in the case the software
   source actually produces: its `M` looks are the correlator's *taps*, which all
   integrate the same `N` samples, so `M·N` over-counts the span threefold.
   The builders still take `total_samples` and `sampling_frequency` exactly as
   specified — only the stored field changed — and
   `noise_observation_from_correlator` gained a `duration` keyword for the
   simultaneous-looks case.

2. **`NoiseUpdateContext` carries the backend.** The plan's interface signature
   omits it, and the plan's own "Backends" table requires it: the one- and
   two-bit accumulators are popcount *counts* rather than sample sums, so a
   float-kernel reference would compare `|P|²` and `N̂₀` on incompatible scales.
   The context was explicitly designed as "one struct so adding a field later is
   not a signature change", which is what that escape hatch was for. Each backend
   supplies one `_correlate_noise_reference!` method.

3. **`requires_noise_density` is a trait on the *type*, and provisioning reads a
   group's slot type rather than a satellite value.** The plan's value-based form
   works only while the default requires no density: once it does, the
   `isempty(satellites)` branch leaves the whole `noise_estimators` NamedTuple —
   keys included — inferring as a union of "provisioned" and "not", which infects
   `TrackState`'s own type. The slot type is fixed at construction and is already
   the right answer either way.

4. **`NWPRCN0Estimator(signal)` is new.** `default_cn0_estimator` used to size
   NWPR's coherent window from the signal's code period (5 blocks for a 1 ms
   code, 2 for L1C-P's 10 ms one). Phase 2 replaced that function's body, which
   would have silently dropped the sizing for the one path NWPR is still
   recommended on. It moved onto the estimator.

One bug the plan did not anticipate, found by the cadence test:
`_chunk_last_sample(cd, index − 1, …) + 1` clamps to `num_samples + 1` without
chunking, so an unchunked whole-buffer pass measured nothing at all. Hence
`_chunk_first_sample`.

Deferred exactly as written: the optional coherent pre-sum, the acquisition CFAR
hook, `noise_gain(::PostCorrFilter)`, and the `Σ|x|²` diagnostic source.

> **Anchors are as of `b9709191` (6.0.1).** Phase 0.5 deletes
> `src/cn0_estimation.jl` and moves every symbol in it, so all
> `src/cn0_estimation.jl:*` line references die at that point. After Phase 0.5,
> locate by symbol name rather than by line.

## Goal

Estimate the noise floor **once per band** and hand it to the C/N₀ estimator, so
C/N₀ follows from a measured noise reference instead of being inferred from the
prompt's own statistics. The noise source is pluggable behind
`AbstractNoiseEstimator`, with two classes of implementation:

- **software** — derived from the band's samples inside `track!` (a shared per-band
  noise correlator, one despread per band per integration),
- **hardware** — the statistic is supplied from outside (FPGA / ASIC correlator,
  or a front-end power monitor) via `append_noise_observation!`, exactly parallel
  to the existing [`append_correlator_output!`](@ref) producer path.

The consumer is one new estimator, `NoiseRefCN0Estimator`, that works on **every**
signal with **one** set of characteristics and **no fallback**.

The final deliverable is documentation that **compares all the estimators and says
when to pick which** — with four CN0 estimators and one noise estimator that can be
filled two ways, that becomes the main thing a user needs from the docs.

## Motivation

The current default, `NWPRCN0Estimator`, is accurate where it applies but does
not apply uniformly:

- It needs a coherent narrowband window of `M ≥ 2` records tiling a navigation
  bit. Where no such window exists it produces nothing and defers to
  `fallback` — a *different estimator with different bias*. That is permanent,
  not a warm-up, for `blocks_per_bit == 1` (GPS L1C-D, Galileo E1B) and for any
  secondary-coded signal pre-sync (`_narrowband_window`,
  `src/cn0_estimation.jl:446-477`). The fallback is `MomentsCN0Estimator`, whose
  noise floor on pure noise is ≈27.6 dB-Hz (`src/cn0_estimation.jl:249-260`,
  `docs/src/cn0_estimator.md:40-90`). This is issue #217.
- `μ̂ → M` saturates, so its relative error asymptotes to `√(M/((M−1)K))`
  regardless of signal strength.
- The ratio-of-sums inversion carries a small positive bias at every C/N₀
  (`src/cn0_estimation.jl:219-247`).
- At low C/N₀ a fraction of *estimates* land outside `1 < μ̂ < M` and report
  `-Inf`/`Inf` dB-Hz (`src/cn0_estimation.jl:705-706`). Note where that bound is
  applied: `estimate_cn0` pools every buffered window — `μ̂ = Σ NBP / Σ WBP` over all
  `num_records ÷ M` of them — and range-checks the *pooled* ratio once
  (`:700-706`). Individual windows are never screened or discarded, so this is a
  rate of unusable outputs, not a per-window discard rate.

Measured by `test/cn0_estimator_comparison.jl` (Phase 0): 100 ms of observation at
1 ms records, NWPR at its default `M = 5`, 4000 trials × 9 seeds per point (median),
**NWPR given its best case** (perfect carrier phase, bit-aligned windows, no
transitions):

| C/N₀ (dB-Hz) | λ = (C/N₀)·T | NWPR σ | NWPR bias | NWPR degenerate | NoiseRef σ (1 ms) | bias | NoiseRef σ (5 ms) | bias |
|---|---|---|---|---|---|---|---|---|
| 20 | 0.1  | 2.96 | **+0.45** | **5.6 %** | 4.77 | +0.03 | **2.75** | +0.01 |
| 25 | 0.32 | 1.51 | +0.08 | 0 | 1.76 | −0.01 | **1.25** | −0.01 |
| 30 | 1    | 0.90 | +0.07 | 0 | 0.76 | +0.01 | **0.65** | 0.00 |
| 35 | 3.2  | 0.64 | +0.05 | 0 | 0.37 | −0.01 | **0.35** | −0.01 |
| 40 | 10   | 0.54 | +0.05 | 0 | 0.20 | −0.00 | **0.20** | −0.00 |
| 45 | 32   | 0.51 | +0.05 | 0 | 0.11 | 0.00 | **0.11** | 0.00 |
| 50 | 100  | 0.50 | +0.05 | 0 | 0.06 | −0.00 | **0.06** | −0.00 |

σ is the relative standard deviation in dB; bias is of the linear-domain mean. The
NWPR column is the *shipped* `NWPRCN0Estimator` folded through its public
`update`/`estimate_cn0`; the two NoiseRef columns are reference implementations of
what `NoiseRefCN0Estimator` must reproduce.

**Two things the NoiseRef columns assume, both of which cost something real:**

- a **variance-free** noise reference. "How much averaging the reference needs"
  derives what that costs and sets the defaults that get within a few percent of it.
- a noise floor that is exactly the thermal one — i.e. **an empty sky**. A despread
  reference also collects every satellite's cross-correlation, and unlike NWPR it
  collects the tracked satellite's *own*. That is a bias, it grows with the
  satellite's own C/N₀, and above ≈45 dB-Hz it is the dominant error term rather
  than σ. See "Multi-access interference and the reference's self-leakage".

NWPR's σ and bias are computed over its **non-degenerate** draws only, which is the
generous reading — the discarded draws are its worst. At 20 dB-Hz that exclusion is
what turns a wide low tail into the visible +0.45 dB bias: 5.6 % of estimates fall
outside `1 < μ̂ < M` and are reported as -Inf or Inf dB-Hz, and the rest sit high.

Reading of the table:

- At **matched record length** the non-coherent noise-reference estimator is
  worse below ≈28 dB-Hz (squaring loss: NWPR sums `M = 5` prompts coherently
  before squaring, so the low-SNR ratio is exactly `√(M−1) = 2`) and better
  above, without bound.
- A coherent sum of five 1 ms prompts **is** one 5 ms prompt. So NWPR's
  narrowband window is not extra information — it is a longer coherent
  integration reconstructed post-correlation. Taking that integration in the
  correlator instead (`preferred_num_code_blocks_to_integrate`) makes the
  noise-reference estimator better than NWPR at **every** C/N₀ *in σ*, and it is the
  better place to spend it because a longer coherent record also helps the
  discriminators. (In *total* error the picture flips above ≈45 dB-Hz once
  self-leakage is counted — see below.)
- NWPR is unbiased nowhere, saturates at 0.50 dB, and throws away 1 output in 18
  at 20 dB-Hz. The noise-reference estimator is unbiased to <0.03 dB throughout
  against a thermal-only floor, and its σ keeps improving with C/N₀.

**Where this actually wins, stated honestly.** The σ advantage is unbounded on
paper, but self-leakage puts a bias floor under the top of the range: ≈0.13 dB at
45 dB-Hz and ≈0.40 dB at 50. So the case for the change is *not* "arbitrarily
accurate strong signals". It is:

1. **It works at all on GPS L1C-D, Galileo E1B and pre-sync secondary-coded
   signals**, where NWPR produces nothing and hands over to a fallback with a
   ≈27.6 dB-Hz floor. This is issue #217 and it is the reason for the change.
2. **20–40 dB-Hz**, where lock and loss decisions are made: no degenerate outputs
   and no +0.45 dB low-end bias. On σ, at matched 1 ms records it wins only above
   the ≈28 dB-Hz crossover; with records lengthened to NWPR's window it wins across
   the whole span. Self-leakage is ≤0.04 dB there (0.013 dB at 35, 0.042 dB at 40),
   i.e. an order of magnitude under NWPR's +0.05 dB.
3. **Uniformity** — one estimator, one set of characteristics, no fallback with a
   different bias, no `M`, no bit sync.

Above 45 dB-Hz the two are comparable in total error and it stops mattering, since
nothing consuming C/N₀ distinguishes 0.4 dB there.

Secondary benefits: a per-band noise floor is also what acquisition needs for a
real CFAR threshold, and it is a natural interference/AGC monitor.

## Design

### The invariant: a noise *density*, not a noise power

For any correlator normalised as we already normalise
(`normalize`, `src/correlators/correlator.jl:181`, dividing by
`integrated_samples * code_amplitude`), white input noise of per-sample variance
`σ²` gives

```
E|P|² = σ²/N = N₀/T
```

`code_amplitude` is defined as the replica's RMS amplitude, so it cancels — CBOC
included. `N₀ = σ²/f_s` is therefore **independent of the integration time** and
is the right thing to store and share:

```
N̂₀  = |B|² / (N_ref · A_c² · f_s)                (from a despread accumulation)
N̂₀  = σ̂²/f_s                                     (from a front-end power monitor)
C/N₀ = ⟨|P|²⟩/N̂₀ − 1/T                          (consumer, per record)
```

`N₀` carries Unitful dimensions of `1/Hz`, so `⟨|P|²⟩/N₀` and `1/T` are both in
`Hz` and the arithmetic is dimension-checked. Consequences:

- The noise reference needs **no alignment to correlator records**. It can run on
  whatever grid is convenient (per chunk) and be averaged over a much longer
  window than the per-satellite prompt buffer, because `N₀` is a slowly-varying
  front-end property. Making the reference's own variance negligible is what
  buys the numbers in the table above.
- Records with different `T` (1 block pre-sync vs
  `preferred_num_code_blocks_to_integrate` after) are handled correctly by
  construction, because `T` enters per record. This retires the class of bug
  behind `dc817828` and `9c9e49d8` structurally rather than by guarding it.

### `AbstractNoiseEstimator`

New file `src/noise_estimation.jl`. Interface, all methods on the abstract type:

```julia
abstract type AbstractNoiseEstimator end

# Software sources: measure one band's chunk of samples and append the observation.
# Mutates the buffer only; the struct is unchanged. No-op for external sources.
update_noise!(est, measurement::BandMeasurement, first_sample, last_sample, ctx) -> est

# Consumer side: the band's current noise density, or `nothing` if none yet.
# A read, not a drain -- the window keeps sliding.
get_noise_density(est) -> Union{Nothing,typeof(1.0/1Hz)}

# Producer side (hardware). Appends to the same buffer.
append_noise_observation!(est, obs::NoiseObservation) -> est
```

`ctx` carries what a source may need but a `BandMeasurement` does not have: a
signal instance for the band, a PRN to use, and the chunk index. Kept as a
single struct so adding a field later is not a signature change.

Implementations:

```julia
struct CorrelatorNoiseEstimator{D} <: AbstractNoiseEstimator
    window_duration::typeof(1.0s)              # default 1 s
    buffered::Vector{NoiseObservation{D}}       # FIFO, written in place
end
```

**One concrete type, and no parameter for who produced the observations.** An earlier
revision carried a `NoiseSource` type parameter (`SoftwareMeasured` /
`ExternallyProvided`) to make `update_noise!` a no-op on the hardware path. It is
unnecessary, because the two paths live in **disjoint call graphs**:

| | fills the window via | enters `downconvert_and_correlate!`? |
|---|---|---|
| software | `track!` → `downconvert_and_correlate!` → `update_noise!` | yes |
| hardware | `append_noise_observation!`, beside `append_correlator_output!` | **no** |

`update_noise!` has exactly one call site, inside `downconvert_and_correlate!`. A hardware
producer never calls that function — it injects correlator outputs and then calls
`estimate_dopplers_and_filter_prompt!(track_state, (L1 = fs,))` with a bare per-band
sampling-frequency NamedTuple and no sample buffers at all (`_band_sampling_frequency`,
`src/conventional_pll_and_dll.jl:694-696`; the path is exercised by
`test/external_correlator_producer.jl`). So `update_noise!` is simply never reached on that
path, and there is nothing for a type parameter to disambiguate.

The only case it would have served is software correlation *plus* a hardware power
monitor — and that combination is not worth a type parameter, because a software despread
is the better reference by construction (matched spectral weighting, same tap point, known
scaling) and is already available whenever the samples are. Anyone who wants it can define
their own `AbstractNoiseEstimator` with a no-op `update_noise!`; the abstract type and its
interface are public.

A band with no entry at all (`nothing` in the NamedTuple) is the "do not measure" case,
which `CN0UpdateContext`'s `Nothing` parametrisation already encodes — so there is no
`NoNoiseEstimator` either.

An earlier draft also shipped a `SampleVarianceNoiseEstimator` (`Σ|x|²` over the chunk).
It is not a separate algorithm: at a sub-integration of **one sample** the replica is a
single ±1 value, the despread degenerates to `|x·c|² = |x|²`, and the window mean is
exactly `σ̂²`. The two "sources" are the two ends of one continuum in sub-integration length —
minimum variance and flat weighting at one end, full spectral fidelity at the other. This
plan sits at the spectral-fidelity end and does **not** expose the length as a knob (see "How
much averaging the reference needs"), so that limit is explanatory rather than reachable; the
diagnostic follow-up would implement it separately.

The window is bounded in **time**, not in observation count: entries are dropped from the
front while `Σᵢ total_samplesᵢ` exceeds `window_duration · f_s`. That is what lets one
producer report 16-chip accumulations and another a single pre-averaged 0.2 s figure, with
the same configuration and no scalar state — the sum is computable from the FIFO.

### The reference is open-loop

Worth stating explicitly, because it is what makes the hardware side cheap: the noise
correlator has **no feedback of any kind**. It runs at the band's nominal IF with zero
Doppler and an off-peak code phase from a rotating PRN, so there is no discriminator, no
loop filter, and no NCO update ever written to it. Consequences:

- It is correct even with **zero satellites tracked**, which is what makes it usable as an
  acquisition CFAR floor.
- Nothing about it interacts with where the PLL/DLL run. In a split system where the FPGA
  does only downconversion and correlation and the loops run on the host, the noise
  channel needs no host round-trip at all — it free-runs and dumps, and its accumulations
  ride the same return path as the correlator outputs. No new wire interface.
- On an FPGA a noise channel is therefore a **strict subset** of a tracking channel: same
  datapath, feedback disconnected, fixed NCO increments, one accumulator instead of three,
  and no dump-on-command synchronisation. Cheaper than the "one channel out of N" estimate
  below, not more expensive.

### `NoiseObservation`

The producer-side value, parallel to
[`CorrelatorOutput`](@ref) (`src/correlators/correlator.jl:42-46`). It stores the
already-normalised density plus a weight, so the internal contract is one scalar
and the unit conversion happens at the boundary where the caller knows their
hardware:

```julia
struct NoiseObservation{D}
    noise_density::D             # N₀, dimension 1/Hz
    num_sub_integrations::Int    # M — statistical weight
    total_samples::Int           # M · N — window-length bookkeeping
    prn::Int16                   # which PRN measured it; carries the rotation position
end
```

**The weight is `M`, the number of independent sub-integrations — not the sample count.**
An earlier draft used samples, which is wrong: for `Σ_m |B_m|²` over `M` sub-integrations
of `N` samples each,

```
N̂₀ = Σ_m |B_m|² / (M · N · A_c² · f_s) = Σ_m |B_m|² / (total_samples · A_c² · f_s)
```

so the *density* needs only the total sample count, but its **relative variance is `1/M`**.
A producer that hands back one 1 ms coherent despread and one that hands back 64 16-chip
despreads report the same total samples and wildly different precision. Weighting the
window by `M` is what makes observations from different producers combinable.

The natural granularity is **one observation per dump, `M = 1`** — a producer that dumps
every 16 chips pushes 64 observations per 1 ms chunk, exactly as the software path emits
`num_sub` per chunk. `M > 1` exists only so a high-rate producer can pre-sum on-chip to cut
host traffic; pre-summing must be **incoherent** (`Σ|B_m|²`, never `|ΣB_m|²`), which is why
the pre-summed builder takes a power rather than a complex value.

Three builders:

```julia
# One dump, as a raw complex accumulation Σ x·c. The most faithful form: Tracking
# squares it, so nothing but the accumulate happens off-host.
noise_observation(accumulation::Complex, integrated_samples, sampling_frequency;
                  code_amplitude = 1)

# M dumps pre-summed on-chip: Σ_m |B_m|², each over total_samples ÷ M samples.
noise_observation_from_correlator(accumulated_power, num_sub_integrations,
                                  total_samples, sampling_frequency;
                                  code_amplitude = 1)

# A front-end / AGC power monitor: Σ|x|² over num_samples. Equivalent to the above
# with one-sample sub-integrations, so M == num_samples.
noise_observation_from_samples(accumulated_power, num_samples, sampling_frequency)
```

All three reduce to the same `N₀` for white noise, which is the check that the paths are
interchangeable and belongs in the tests.

#### `append_noise_observation!` vs `append_correlator_output!`

They are deliberately not the same mechanism, and the docs should say so plainly:

| | `append_correlator_output!` | `append_noise_observation!` |
|---|---|---|
| scope | one **signal** of one satellite (`src/sat_state.jl:273`) | one **band** |
| payload | the whole correlator: every tap kept **separate**, complex, per antenna, because E, P and L drive different discriminators | the taps **pooled** into one power sum, plus `M`, the sample count and the PRN |
| drives | code and carrier discriminators, both loop filters, bit sync, the NCO update, the C/N₀ prompt | the noise floor only — no loop, no discriminator, no bit sync |
| timing | `sample_index` is load-bearing; vector tracking needs outputs on a common grid | none — `N₀` is slowly varying, so only "recent" matters |
| lifetime | **drained** by the fold each chunk, buffer reused | **sliding window** bounded in time, never drained |
| count | one per completed coherent integration per signal | any granularity the producer likes; `M` carries the weight, so a 1-dump and a 64-dump observation combine correctly |

The payload row is the sharpest difference. Both traverse the same multi-tap correlator — see
"Use every tap, spaced wide" — but for opposite reasons. The prompt path needs the taps kept
apart because their *differences* are the code and carrier discriminants. The reference pools
them, because at ≥1 chip spacing they are three independent looks at the same scalar and
nothing about their relative values means anything.

So a noise observation is not a degenerate `CorrelatorOutput` and does not travel through
`TrackedSignal.correlator_outputs`. Reusing that path would mean giving the noise channel a
PRN key in `SignalGroup.satellites`, a bit buffer, a post-corr filter and a C/N₀ estimator
it must never use — and would surface it in `estimate_cn0` and every result accessor.

### Where the state lives

`TrackState` (`src/Tracking.jl:232-235`) gains a third field:

```julia
struct TrackState{G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    groups::G
    doppler_estimator::DE
    noise_estimators::NE      # NamedTuple keyed by band id
end
```

`NoiseEstimators = NamedTuple{<:Any,<:Tuple{Vararg{AbstractNoiseEstimator}}}`,
keyed by `GNSSSignals.get_band_id` symbols — the same idiom as `BandMeasurements`
(`src/band_measurement.jl:115`), so lookups fold to compile-time constants and
`_band_sampling_frequency` (`src/conventional_pll_and_dll.jl:694-696`) can be
copied verbatim for `_band_noise_density`.

Rejected alternatives, and why:

- **`SignalGroup`** — a group is not a band. Grouping is by signal-tuple shape;
  two groups can share one band (`src/sat_state.jl:936-947`), so per-group state
  would duplicate the estimator and split the averaging.
- **The `dc` backend struct** — wrong lifetime. `downconvert_and_correlator`
  defaults to a freshly constructed value on every `track!` call
  (`src/track.jl:176`), so the average would reset unless the user hoists it. It
  is also not a backend concern: the same noise estimator must work across
  backends.

Because the state is mutated in place, `TrackState` itself never needs rebuilding
for a noise update, and the per-band keying means one estimator instance is never
shared between bands (the same reasoning as `_per_signal_cn0_estimators`,
`src/sat_state.jl:617-659`).

### Where it is updated

At the three `downconvert_and_correlate!` entry points, which are the only places
that see `measurements` whole and therefore the only per-band-exactly-once sites:
`src/downconvert_and_correlate_cpu.jl:856` (CPU, CPUThreaded, and Int16 by
inheritance), `src/downconvert_and_correlate_onebit.jl:1349`,
`src/downconvert_and_correlate_twobit.jl:1838`.

One shared helper defined on `AbstractDownconvertAndCorrelator`, called
immediately before `_foreach_group!` at each site:

```julia
_update_band_noise!(track_state, measurements, chunk_index, chunk_duration, samples_unchanged)
```

It walks the `noise_estimators` NamedTuple by **tuple recursion**, in the style of
`_foreach_group!` (`src/sat_state.jl:1147-1154`) — not with a runtime
`for (k, v) in pairs(...)` loop. Different bands may hold different estimator
types, so a runtime loop over the NamedTuple's values is type-unstable and would
allocate on every chunk. Recursion unrolls it and each `update_noise!` call
devirtualises. For each band it recomputes the chunk sample range with the
existing `_chunk_last_sample` (`src/downconvert_and_correlate_cpu.jl:885-895`) and
calls `update_noise!`.

`_dc_one_group!` was rejected as the site: it fires once per *group*, so a band
with two groups would be measured twice.

#### Same kernel, not the same channel

In hardware the noise reference *is* a normal downconvert-and-correlate on a PRN nobody is
tracking, with the feedback disconnected. The software side should be the same thing — but
"the same thing" means reusing the **kernel**, not modelling a fake satellite.

Reuse the single-satellite convenience form already in the codebase,
`downconvert_and_correlate!(signal_type, signal, correlator, code_replica, code_phase,
carrier_phase, code_frequency, carrier_frequency, sampling_frequency, signal_start_sample,
num_samples_left, prn)` (`src/downconvert_and_correlate_cpu.jl:974-1012`). It is
`gen_code_replica!` plus the fused kernel and nothing else — bit-identical arithmetic to the
prompt path, which is exactly the model-freeness the whole approach rests on, with none of
the channel machinery. Call it once per sub-integration.

Do **not** plumb it as a `TrackedSat`/`TrackedSignal`, for three concrete reasons:

- **The dump cadence is wrong for a channel.** Every completion path in the correlate step
  is built around code-block boundaries (`calc_num_code_blocks_to_integrate`,
  `_calc_min_samples_and_completed`). The noise reference dumps on a *time* grid derived from
  the primary code period of the band's reference signal, or the chunk where that is shorter.
  It is not re-anchored to absolute code phase, and it carries no partial across calls.
- **It would drag in everything it must not use**: a bit buffer, a post-corr filter, a C/N₀
  estimator, a doppler-estimator state, and a PRN key in `SignalGroup.satellites` that would
  then surface in `estimate_cn0` and every result accessor.
- **The estimator would try to steer it.** `_update_tracked_sat_doppler` would run
  discriminators and write an NCO update to a channel that must stay at nominal.

So: same kernel, same arithmetic, same replica generator; no channel, no feedback, no
satellites-dictionary entry.

#### Integration boundaries, and how many observations per call

The prompt correlators integrate boundary-to-boundary — code phase 0 to code length —
because the despread signal only accumulates coherently over a whole code period. **The
noise correlator has no such constraint.** It despreads with a deliberately wrong PRN at
a deliberately wrong phase, so there is no signal to accumulate coherently; the estimate

```
N̂₀ = |B|² / (N_ref · A_c² · f_s)
```

is unbiased for *any* `N_ref` and *any* starting phase, because `A_c` is the replica's RMS
amplitude by definition and the noise is white. Nothing about a code-period boundary
enters.

That matters, because following code-period boundaries would be actively harmful here: a
chunk generally does not contain a whole number of code periods, so a boundary-aligned
noise correlator would have to **carry a partial accumulator across calls** — a complex
sum plus a sample count, i.e. scalar state, i.e. exactly what has no anchor per band (see
"Immutability and allocation"). Integrating slices of the chunk carries nothing.

So each sub-integration simply starts its replica at code phase 0 with the current
rotating PRN and runs for its slice. Successive observations are independent because their
*sample* ranges are disjoint — the replica repeating is irrelevant.

**But the count per call is not one.** Correcting what this plan said a revision ago:
one-per-call would weld the noise window's time span to `doppler_update_interval`. The
window is `K ≈ 200` observations, so at a 1 ms chunk it spans 200 ms, but a user who sets
`doppler_update_interval = 20ms` would get a 4 s window, over which `N₀` is no longer
stationary. The fix is to make the sub-integration length a property of the *source*, not
of the chunk:

```julia
num_sub = max(1, round(Int, chunk_duration / target_integration_time))
```

where `target_integration_time` **is** the primary code period of the band's reference signal — see "How much
averaging the reference needs" for why it is derived rather than configurable, and why the
variance that costs is bought back with the window instead. At the default
`doppler_update_interval` the two coincide and `num_sub` is 1; a longer update interval yields
several. The chunk is split into `num_sub` **equal**
slices — equal rather than
fixed-length so no remainder is wasted and every window entry is statistically identical.
The `max(1, …)` matters: a band whose signal has a 10 ms code sampled with 1 ms chunks
would otherwise yield zero observations forever.

Per call, per band: `num_sub` observations, pushed **individually**, not pre-averaged.
Pre-averaging would put one entry per chunk in the window, which re-welds its time span to
the chunk length and unequally weights a buffer's trailing partial chunk. Pushing
individually decouples C/N₀ quality from the Doppler update rate entirely — worth having,
since that coupling is already a documented wart elsewhere in this design.

Per `track!` over a buffer of `C` chunks: `C · num_sub` observations per band.

**`samples_unchanged` must not gate the measurement.** An earlier draft said a re-passed
buffer "is not re-measured", by analogy with the one-bit backend's shared band pack
(`src/downconvert_and_correlate_onebit.jl:1267-1290`). The analogy is wrong: that pack is a
**whole-buffer** transform, so `samples_unchanged` tells it the cache is still valid while
the chunk loop keeps advancing through `chunk_last_sample`. Since `track!` passes
`samples_unchanged = chunk_index > 0`, gating on it would measure chunk 0 and then freeze
the window for the rest of the buffer.

The call that *does* need skipping is the **final drain pass** (`src/track.jl:224-229`), which
passes neither `chunk_index` nor `chunk_duration`. With `chunk_duration === nothing` the
chunk range defaults to the whole buffer (`_chunk_last_sample`,
`src/downconvert_and_correlate_cpu.jl:885-895`), so measuring there would re-measure
everything already covered. The satellites avoid this through their per-signal
`signal_start_sample`; a per-band source has no such pointer, and a pointer would be the
scalar state that has no anchor.

The discriminating condition needs no pointer:

| call | `chunk_duration` | `samples_unchanged` | measure? |
|---|---|---|---|
| chunked call `k` | set | `k > 0` | yes, `num_sub` slices of chunk `k` |
| final drain pass | `nothing` | `true` | **no** |
| unchunked single call | `nothing` | `false` | yes, whole buffer in `num_sub` slices |

So: skip iff `chunk_duration === nothing && samples_unchanged`. The cost is that the
trailing partial's samples never reach the window — at most one chunk per `track!` call,
which changes nothing, since omitting an entry biases `σ̂²` not at all.

Verified against every shape a call can take. Note first that **inside `track!` the
unchunked case never arises**: `_resolve_doppler_update_interval` (`src/track.jl:271-275`)
falls back to `_smallest_code_period`, so `chunk_duration` is always a concrete time there.
`track!` therefore has exactly two call shapes, and the third is reachable only by a direct
caller:

| call | `chunk_duration` | `samples_unchanged` | gate | measures |
|---|---|---|---|---|
| chunked loop, iteration `k` (`:206-218`) | set | `k > 0` | no skip | chunk `k` in `num_sub` slices |
| drain (`:224-229`) | `nothing` | `chunk_index > 0` | **skip** | — |
| direct unchunked `downconvert_and_correlate!` | `nothing` | `false` | no skip | whole buffer |

A buffer shorter than one chunk is covered by the **loop**, not the drain: `_chunks_left`
(`:258-265`) seeds from `_chunk_last_sample(cd, -1, …) = 0 < n`, so the loop always runs at
least once for a non-empty buffer, and `_chunk_last_sample` clamps that single chunk to the
buffer end. The drain then sees `chunk_index = 1` and is skipped, which is right — the loop
already measured everything. The drain is reached with `samples_unchanged = false` only when
the loop ran zero times, i.e. an empty buffer, where there is nothing to measure anyway.

**This does overload a public kwarg, which is worth stating rather than discovering.**
`samples_unchanged` is documented as a *cache-validity hint* — "the measurement buffers are
fixed for the whole call, so sample-derived backend caches … are built on the very first pass
and reused ever after" (`src/downconvert_and_correlate_cpu.jl:849-854`) — and
`downconvert_and_correlate!` is public API. After this change it also means "do not measure",
so a direct caller who passes `samples_unchanged = true` on an unchunked call, truthfully,
gets their noise measurement silently skipped. The gate is right for `track!` and a
pointer-free discriminator is worth having, but the kwarg's docstring must gain a sentence
saying so, and the Verification cadence test should cover the direct-caller shape as well as
`track!`'s.

### How much averaging the reference needs

The Motivation table's noise-reference columns assume the reference contributes no
variance of its own. That has to be earned. Writing `K_rec` for the C/N₀ ring depth (100
records) and `K_n` for the number of noise observations averaged:

```
sd_dB ≈ 4.343 · sqrt( (1+2λ)/K_rec + (1+λ)²/K_n ) / λ
```

The second term is the reference's contribution. Verified against Monte Carlo to two
decimals at every C/N₀ in the table. The requirement for it to be negligible,

```
K_n  ≳  10 · K_rec · (1+λ)²/(1+2λ)
```

**grows with λ**, because the prompt term's relative error keeps falling as `√(2/(λK_rec))`
while the reference's is floored at `1/√K_n`. So there is no fixed `K_n` that is negligible
everywhere, and a too-small one quietly eats exactly the high-C/N₀ advantage that motivates
this whole change:

| C/N₀ | `K_n` = 200 | 2 000 | 20 000 | variance-free floor |
|---|---|---|---|---|
| 20 | 5.85 | 4.88 | 4.78 | 4.76 |
| 30 | 0.98 | 0.78 | 0.76 | 0.75 |
| 40 | 0.40 | 0.23 | 0.20 | 0.20 |
| 45 | 0.34 | 0.15 | 0.11 | 0.11 |
| 50 | 0.32 | 0.12 | 0.07 | 0.06 |

(σ in dB, `K_rec` = 100, simulated and analytic agreeing to ±0.01.) At `K_n` = 200 the
50 dB-Hz figure is **0.32 dB against a 0.06 dB floor** — five times worse, and barely
better than NWPR's 0.50 dB saturation. An earlier revision of this plan defaulted to
≈200 observations of one code period each. That was asserted, not derived, and it was
wrong.

#### The sub-integration length is the knob — but the default is a full code period

**Coherent integration buys nothing for noise-power estimation.** This is why a shorter dump
raises `K_n` at all. For a complex-Gaussian accumulation, `E|B|² = σ²·N·⟨c²⟩` *and*
`Var(|B|²) = (E|B|²)²` — so one dump carries 100 % relative error **regardless of how long it
is**. A 1 ms coherent dump and a 16-chip coherent dump are each worth exactly one look. All
the information about `σ²` lives in the *number of independent looks*, and coherent
integration destroys looks by collapsing `N` samples into one complex number.

(That is the opposite of the prompt correlator, where coherent integration is the whole point
because there *is* a signal to accumulate. It is why the noise reference cannot simply
inherit the prompt's integration boundaries.)

So `K_n` per second is `1/T_sub`. But three things push back on shortening `T_sub`, and
together they win:

- **SIMD.** The fused kernel is vectorised, and a short integration amortises per-call setup,
  tail handling and NCO initialisation over far fewer vector iterations — 64 kernel calls per
  1 ms chunk instead of one. It *could* be avoided by snapshotting and resetting the
  accumulator inside a single kernel pass (a horizontal reduction every ~47 vector iterations
  is a few percent), but that needs a **new kernel variant**, forfeiting the "reuse the
  existing kernel, bit-identical arithmetic" property the whole approach rests on.
- **DC balance.** A full-period Gold code is balanced (`E[(Σc)²]` = 1.0, measured); any
  sub-period window behaves like iid signs (16.3 at 16 chips, 57.0 at 64), so an ADC DC
  offset biases a short despread and not a full-period one. See Risks.
- **Cadence parity with hardware.** A hardware correlator naturally dumps on the code epoch.
  Defaulting software to the same makes the two paths agree in cadence as well as in code.

So `T_sub` defaults to **one full primary code period**, a derived quantity, not a
parameter. The variance is bought back with the
**window** instead, which is nearly free — a 1 s window is 1000 `Float64` densities, 7.8 kB
per band:

| C/N₀ | floor | 0.2 s (`K_n`=200) | **1 s (`K_n`=1000)** | 2 s (`K_n`=2000) |
|---|---|---|---|---|
| 20 | 4.76 | 5.83 (+1.08) | **4.99 (+0.23)** | 4.88 (+0.12) |
| 25 | 1.755 | 2.17 (+0.42) | **1.845 (+0.09)** | 1.80 (+0.05) |
| 30 | 0.752 | 0.97 (+0.22) | **0.801 (+0.05)** | 0.78 (+0.02) |
| 40 | 0.199 | 0.39 (+0.19) | **0.250 (+0.05)** | 0.226 (+0.03) |
| 50 | 0.062 | 0.32 (+0.25) | **0.152 (+0.09)** | 0.116 (+0.05) |

(σ in dB; brackets are the absolute penalty over the variance-free floor.)

**Read this in absolute dB, not as a ratio.** The percentage penalty looks alarming at
50 dB-Hz (+414 % at `K_n` = 200) and harmless at 20 dB-Hz (+23 %), but the *absolute* cost
runs the other way: **+1.08 dB at 20 dB-Hz** against +0.25 dB at 50. Low C/N₀ is where lock
and loss decisions are made, so that is the end that actually needed fixing — and the 1 s
window fixes it, bringing every point to ≤0.23 dB and everything above 25 dB-Hz to ≤0.09 dB.
At that level it is genuinely negligible and no further tuning is warranted.

**And there is no knob for it.** An earlier revision kept `sub_integration_chips` configurable
"for the one case that wants it". Priced out, that case does not exist: shortening the
sub-integration buys ≈0.08 dB at 50 dB-Hz, ≈0.06 dB at 45 and nothing measurable below 40 —
against a public exported parameter, a second SIMD regime to benchmark, a DC-offset failure
mode to document and test, and cadence divergence from every hardware producer. Nothing that
consumes C/N₀ can tell 0.07 dB from 0.15 dB on a 50 dB-Hz signal. So the field is gone.

`window_duration` remains the one tunable, because it trades genuinely: longer means lower
variance (table above) but a longer smear across AGC changes (see Risks).

The chunk still decides the length in exactly one configuration — when it is **shorter** than
the reference signal's code period, which another band can force via a smaller
`doppler_update_interval` (a band carrying only GPS L1C-P at 10 ms alongside one carrying
L1 C/A at 1 ms). Then `max(1, …)` gives `num_sub = 1` and the sub-integration is a *fraction*
of a code period. That is the precise and only trigger for the DC-offset risk.

#### Use every tap, spaced wide

The reference correlator should use **all** of its accumulators, not just the prompt. An
earlier revision implicitly used one, which threw away a free factor of three.

The catch is that taps are only worth having if they are statistically **independent**, and at
the default E/P/L spacing they are not. For jointly circular Gaussian accumulations,
`corr(|B_i|², |B_j|²) = |ρ_ij|²` where `ρ_ij` is the sampled-code correlation between the two
replica shifts. Measured over GPS L1 C/A PRNs 1-32 at 4 samples/chip:

| tap spacing | mean `ρ` | max `\|ρ\|` | max power correlation |
|---|---|---|---|
| 0.5 chip | 0.502 | 0.531 | 0.282 |
| 1.0 chip | 0.005 | 0.064 | 0.004 |
| 1.5 chips | 0.000 | 0.032 | 0.001 |
| 2.0 chips | −0.005 | 0.064 | 0.004 |

So at the tracking default of ±0.5 chip, three taps are worth **2.25** independent looks, not
3. At **≥1 chip** they are worth **2.98** — effectively fully independent, the residual being
the Gold code's off-peak autocorrelation floor (±64/1023).

The reference therefore uses a **wide tap spacing, ≥1 chip** (1.5 chips as the default: it
sits in the autocorrelation null and stays clear of both the 1- and 2-chip sidelobe values).
This costs nothing extra — the taps are already in the correlator type and the kernel, the
downconversion is shared across them, and it is still one channel per band. `M` per dump
becomes the number of taps, which `NoiseObservation`'s weight field already carries.

What it buys, at code-period dumps and a 1 s window:

| C/N₀ | 1 tap (`K_n`=1000) | **3 taps (`K_n`=3000)** | 5 taps (`K_n`=5000) |
|---|---|---|---|
| 20 | +0.234 | **+0.079** | +0.048 |
| 30 | +0.049 | **+0.017** | +0.010 |
| 40 | +0.051 | **+0.018** | +0.011 |
| 50 | +0.090 | **+0.039** | +0.026 |

(absolute penalty in dB over the variance-free floor.) Three taps bring every point to
≤0.08 dB. Equivalently — and this is the more useful way to spend it — the window can drop to
≈0.35 s for the same accuracy as one tap at 1 s, cutting the AGC-smear exposure in Risks by
three. `VeryEarlyPromptLateCorrelator`'s five taps extend this further at no structural cost.

This also strengthens the hardware argument: an FPGA tracking channel **already computes three
accumulators**. Reporting all three costs nothing beyond programming the spacing wide, and if a
design's spacing is fixed at ±0.5 chip it degrades gracefully to 2.25 looks rather than
failing. Where the spacing is fixed, the producer should say so, since it changes `M`.

#### This is also what became of the sample-variance source

At `T_sub` = **one sample** the replica is a single ±1 value, the despread degenerates to
`|x·c|² = |x|²`, and the window mean is exactly `σ̂²`. So the "two sources" an earlier
revision shipped were never alternatives — they are the two ends of this one knob, with
minimum variance and flat weighting at one end and full spectral fidelity at the other.
Dropping `SampleVarianceNoiseEstimator` as a *type* therefore costs nothing but the
diagnostic value of running both weightings at once, which is what the Deferred follow-up
records.

### Where it is read

`_est_one_group!` (`src/conventional_pll_and_dll.jl:704-715`) already holds
`get_band_id(g.band)` and already does a band-keyed lookup for the sampling
frequency. Add a second band-keyed lookup there, read that band's
`get_noise_density`, and thread it down to
`_apply_correlator_output` (`src/conventional_pll_and_dll.jl:413`), where it joins
the `CN0UpdateContext` built at `:459-463`.

The fold reads one scalar per band, never the samples — `measurements` stays in this call
purely for `_band_sampling_frequency`, as today. All sample processing is in the
correlate step.

`_foreach_group!` currently receives only `track_state.groups` and
`sampling_frequencies` (`src/conventional_pll_and_dll.jl:738`), so it also needs
`track_state.noise_estimators` passed through — one extra argument at one site.

Functions gaining **two** arguments — the density and a `Bool` saying whether it is
meaningful this chunk (see "Two different 'no density' conditions" for why the
availability flag cannot be folded into the density's type, and why it cannot be
resolved higher up) — all underscore-private:
`_update_tracked_sat_doppler` (`:284`), `_process_estimator_driver_signal`
(`:499`, **two** dispatch families — the vector-tracking one is
`src/vector_pll_and_dll.jl:264`), `_process_passenger_signals` /
`_process_one_passenger_signal` (`:600`, `:620`), `_apply_correlator_output`
(`:413`). Three `_apply_correlator_output` call sites:
`src/conventional_pll_and_dll.jl:534`, `:636`, `src/vector_pll_and_dll.jl:295`.

Both arguments are scalars of a type fixed per call site, so every one of these stays
monomorphic; the flag is what keeps the density from ever needing to be
`Union{Nothing,D}`.

`CN0UpdateContext` (`src/cn0_estimation.jl:48-54`) gains two fields:

```julia
noise_density::N       # N === Nothing when no reference is available
integration_time::T    # this record's own T
```

`N` is a **type parameter, not a `NaN` sentinel**: a given setup consistently has
or has not got a reference, so every call site monomorphises and the
allocation-free contract at `src/cn0_estimation.jl:629-633` holds. It also gives
the misconfiguration check for free — see below.

### `NoiseRefCN0Estimator`

```julia
struct NoiseRefCN0Estimator <: AbstractCN0Estimator
    num_records::Int
    buffered_cn0::Vector{Float64}   # linear Hz, written in place
    current_index::Int
    filled_length::Int
end
```

It requires a noise density on its band, supplied by **any** `AbstractNoiseEstimator` —
`CorrelatorNoiseEstimator` being the only shipped one. On a hardware path you still configure
that type; you simply fill it with `append_noise_observation!` instead of letting
`update_noise!` fill it (see "One concrete type, and no parameter for who produced the
observations"). What happens when no density is available is set out under "Two different 'no
density' conditions" below.

One ring, advancing once per correlator record. Noise averaging is the band's job, so
nothing about `N₀` is stored here: the estimator is handed a density per record and
uses it immediately.

The scalar index is legal here, unlike in the per-band source, because this estimator
lives in `TrackedSignal`, which `_build_new_signals` rebuilds and writes back through
`SignalGroup.satellites`.

`update` computes, per record, `(|P|²/N₀ − 1/T)/1` in linear Hz and rings it in.
`estimate_cn0` averages and converts once with `dBHz`.

Three deliberate properties:

- **No `fallback` field.** The three reasons NWPR needs one — warm-up, no
  admissible window pre-sync, and no window ever at `blocks_per_bit == 1` — all
  vanish. One record plus a density is a valid (noisy) estimate, on every signal.
- **`estimate_cn0`'s `integration_time` argument is ignored**, because `T` is
  applied per record at update time. The signature stays for interface uniformity;
  this gets an explicit comment so it does not read as an oversight.
- **Clamp only the average, never the per-record term.** Individual terms go
  negative at low C/N₀ and that is exactly what makes the mean unbiased.
  Clamping per record would reintroduce a floor. The final average is floored the
  way PocketSDR floors its `0.01·sumN`.

**Two different "no density" conditions, deliberately kept apart.** They were conflated in an
earlier revision, which would have made a legitimate warm-up throw:

- **No source configured** — a *static* property of the setup. `CN0UpdateContext`'s `N` type
  parameter is `Nothing`, and `update` throws (below). Compile-time, type-stable.
- **Window not filled yet** — a *runtime* condition, where `get_noise_density` returns
  `nothing`. This must not throw and must not reach the context, or the context field would
  have to become a `Union` and lose type stability.

The second is resolved one level up, in `_est_one_group!`, which reads the band's density
*before* building any context. But **where the resulting branch lands, and how wide it is,
both need pinning** — an earlier revision said only "does not update the C/N₀ estimators for
that band", which is not implementable as written and would have been wrong if it were.

**The branch cannot live in `_est_one_group!`.** That function calls
`_update_tracked_sat_doppler`, which reaches `_apply_correlator_output`
(`src/conventional_pll_and_dll.jl:413-489`) — and that function folds a record as one
indivisible unit: it normalises the correlator (`:422-423`), advances the post-corr filter
(`:424-426`), pushes the filtered prompt (`:428`), computes the bit-block count
(`:429-434`), updates the C/N₀ estimator (`:459-463`) and feeds the bit buffer
(`:464-470`) — everything a record contributes except the discriminators themselves, which
run one level up in `_process_estimator_driver_signal` (`:499`) off the state this function
produces. Skipping it for a band would stall bit sync and the prompt history the loops read.
The only thing that may be skipped is the `update(get_cn0_estimator(…), …)` call itself, so
the "have we got a density" answer has to reach `_apply_correlator_output`.

Thread it as a `Bool` (or `Val`), not as a `Union{Nothing,D}` density: the density argument
stays a plain scalar of the band's type, and a separate flag says whether it is meaningful.
Both branches inside `_apply_correlator_output` are then type-stable and the argument type
never varies. The cost is honest and worth stating: this is a **second** value threaded
through the same five private functions the density already goes through
(`_update_tracked_sat_doppler`, `_process_estimator_driver_signal` in both dispatch families,
`_process_passenger_signals` / `_process_one_passenger_signal`, `_apply_correlator_output`),
not a free branch at the top.

**The skip must be per estimator, not per band.** A band is one-to-many with signals, and
only signals whose estimator returns `requires_noise_density() == true` care. Skipping the
whole band would stall an `NWPRCN0Estimator` sharing it — and for NWPR that is not "one
fewer sample", it is **state corruption**: `_update_nwpr` drops the open narrowband window
whenever `window_block_index != num_accumulated_code_blocks` — the test is
`src/cn0_estimation.jl:554-555` and the drop itself `:556-557` — so records missing from the
bit grid silently break window
alignment and NWPR degrades to its fallback. Gate the skip on the same
`requires_noise_density` trait that drives provisioning, evaluated per signal — it is a
compile-time constant on the estimator's type, so this costs nothing at run time.

The context's density field is then unconditionally a plain scalar whenever
`N !== Nothing`, and both branches are type-stable — one calls `update`, the other does not.

In the software path the second condition is at most a formality: `downconvert_and_correlate!`
measures before the fold reads, so the window holds an observation from the very first chunk.
It is reachable on the external path, where a producer may fold before appending, and on a
`track!` call whose buffer is shorter than one sub-integration. An entirely empty C/N₀ ring
reports `-Inf dB-Hz`, matching the house convention that a missing estimate is `-Inf` and
never `NaN` (`src/cn0_estimation.jl:168-172`).

**A configured-but-never-fed estimator must not stay silent.** The two conditions above get
loud treatment at different strengths, and the gap between them is the likeliest hardware
integration mistake: a producer that configures a `CorrelatorNoiseEstimator` per the docs and
then forgets `append_noise_observation!` has a *configured* source whose window is empty
forever. The static check does not fire (a source exists), and the runtime skip repeats
silently, so every C/N₀ reads `-Inf dB-Hz` indefinitely.

That is a loud enough *symptom* — nobody mistakes "every satellite reports `-Inf`" for correct
output — but it is not a loud enough *diagnosis*. So `_est_one_group!` emits a single
`@warn … maxlog = 1` at the skip site, naming the fix:

> band `:L1` has a noise estimator but no noise density, so `NoiseRefCN0Estimator` will not
> update. Call `append_noise_observation!` before
> `estimate_dopplers_and_filter_prompt!`, or use a source that measures from samples.

The skip site is the right place because it is the only point that knows a fold actually ran
*and* the density was unavailable. It needs no state (`maxlog` is handled by the logging
macro), it never fires on the software path, and on the short-buffer edge it fires once with
information that is arguably worth having anyway.

**`maxlog` is keyed per callsite, not per band**, so with two misconfigured bands only the
first to reach the skip ever warns — and since the message interpolates the band id, the
other band's name never appears. Two options, and the second is preferred: either pass an
explicit per-band `_id` (`@warn … _id = Symbol(:no_noise_density_, band_id), maxlog = 1`),
or drop the band id from the message and name the general fix. Take the `_id` route — the
band id is the single most useful thing in the message for a multi-band setup, which is
exactly where the mistake is most likely.

**Warn rather than throw, deliberately**, and the asymmetry with the static case is the point:
"no source configured" is unambiguously a mistake, whereas "window empty at this instant" has
a legitimate transient reading, so making it fatal would break a caller streaming buffers
shorter than one sub-integration.

Where `noise_density` is `Nothing`, `update(::NoiseRefCN0Estimator, prompt,
::CN0UpdateContext{<:Any,<:Any,Nothing})` throws an `ArgumentError` naming the
fix ("no band noise estimator configured; pass `cn0_estimator =
NWPRCN0Estimator()` for externally supplied correlator outputs, or set
`noise_estimators`"). This fires on the first record and is loud, which is right
for a wiring mistake — a silent substitution would hide a backend that fails to
populate the reference.

### Immutability and allocation

**No `mutable struct` anywhere in this work, and every steady-state path
allocation-free.** Both match existing contracts: the CN0 estimators
(`src/cn0_estimation.jl:629-633`) and the backend scratch buffers
(`src/downconvert_and_correlate_cpu.jl:10-11`). Three paths run per chunk or per
record and are in scope: `_update_band_noise!` and `update_noise!` (once per band
per chunk), `update(::NoiseRefCN0Estimator, …)` (once per record per signal), and
`append_noise_observation!`.

What this constrains:

- **Immutable structs only, and no mutable containers for scalar state** — no
  `mutable struct`, and no `Ref` either (`Base.RefValue` *is* a mutable struct).
- **Append and ring buffers are `Vector` fields written in place**, as
  `TrackedSignal.correlator_outputs` and `NWPRCN0Estimator`'s power rings already are.
  A buffer is not scalar state; that is the existing house pattern
  (`src/cn0_estimation.jl:629-633`, `src/downconvert_and_correlate_cpu.jl:10-11`) and
  it stays.
- **Scalar state only where a rebuilt struct has somewhere to be written back.**
  `NWPRCN0Estimator` gets away with `ratio_current_index::Int` because it lives in
  `TrackedSignal`, which `_build_new_signals` rebuilds and writes back through the
  `Dictionary` at `SignalGroup.satellites`. Per-band state has no such anchor, which
  is what shapes the split below.

#### Where per-band accumulation is anchored — resolved

A band has exactly **one** noise density, so it must be averaged **once per band**.
The obstacle was that `TrackState` is immutable and the only reason `track!` works in
place is that `SignalGroup.satellites` is a `Dictionary` for rebuilt
`TrackedSat`/`TrackedSignal` values to land in — there is no equivalent anchor for a
per-band scalar, and a ring index is a scalar.

A ring index is not actually needed. A **length-managed FIFO** carries its own state:

```julia
push!(buffered_densities, density)
length(buffered_densities) > num_observations && popfirst!(buffered_densities)
```

The `Vector`'s own length *is* the position, so there is no scalar to write back, the
struct is never rebuilt, and `get_noise_density` is the mean of the window (`nothing`
while empty). Measured allocation, with the buffer `sizehint!`-ed to
`4 * num_observations + 1`:

| window | pushes | allocated |
|---|---|---|
| 200 | 10 000 000 | **0 bytes** |
| 1000 | 10 000 000 | **0 bytes** |

The four-times headroom is what makes it exactly zero: at 1× or 2× Julia periodically
shifts the front offset back and reallocates (≈0.3 bytes/call at a 1000-deep window).
That is an implementation detail of `Vector`, so the Verification requirement is an
`@allocated == 0` test at the shipped default, not a claim about Julia.

Consequences worth stating:

- The measurement stays in `downconvert_and_correlate!`, where the
  thousands-of-samples work belongs, alongside the despreading it already does.
- `estimate_dopplers_and_filter_prompt!` **keeps receiving `measurements` only for the
  sampling frequency**, exactly as today (`_band_sampling_frequency`,
  `src/conventional_pll_and_dll.jl:694-696`). It must never read `samples`.
- The fold **reads** the band's `get_noise_density`; it drains nothing. The window keeps
  sliding across chunks and across `track!` calls.
- `NoiseRefCN0Estimator` therefore needs **no noise ring of its own** — it is handed a
  density per record and keeps only its C/N₀ ring.

- **`noise_estimators` stays a NamedTuple**, not a `Dictionary`: a dictionary would
  need an abstract value type as soon as two bands hold different estimator types,
  which is type-unstable and allocates on every chunk.
- **`CN0UpdateContext`'s noise density is a type parameter, not `Union{Nothing,…}`
  and not a `NaN` sentinel.** A small union in a struct field risks boxing and
  kills the devirtualisation of `update`; the parametric form monomorphises per
  call site. This is the same choice already argued for on type-stability grounds —
  allocation is the other half of the reason.
- **NamedTuple walks are tuple recursion**, per the section above.
- **`NoiseObservation` is `isbits`** (a Unitful density plus an `Int`), so it is
  stack-allocated and `append_noise_observation!` costs nothing.
- **`N₀` is a Unitful quantity**, which is `isbits` — the dimension checking is
  free at run time.
- No closures capturing mutated locals in the per-sample accumulation loops, and no
  intermediate arrays: accumulate scalars over indices (or `@views` over dense
  ranges), never `sum(abs2, x[a:b])`.

This is a *preference*, not a hard gate, in one place: if the correlator source's
replica handling on a quantised backend cannot be made allocation-free without
distorting the design, take the allocation and record it in the benchmark suite
rather than contorting the abstraction. Everything else should measure zero.

### Multi-access interference and the reference's self-leakage

**No multi-access-interference correction**, and the two estimators are *nearly* but
not exactly alike in what they pick up. This distinction was a parenthetical in an
earlier revision and it should not have been: it is the design's only bias term, and
it lands squarely on the high-C/N₀ end the σ tables advertise.

A satellite of carrier power `C` spreads over the code band when despread with the
wrong code, contributing `ε = C/f_chip` to a measured noise density. At
`C/N₀ = 45 dB-Hz` on L1 C/A that is `10^4.5 / 1.023e6 = 0.031` of `N₀`, i.e.
**0.13 dB**. It is Doppler-independent: the product of two different codes is spread
over ±`f_chip`, so a ±5 kHz offset does not null it.

**The shared term.** For the 11 *other* satellites, both estimators see the same
thing. NWPR sees `E[NBP] = M²P_s + M(P_n+ε)` and `E[WBP] = M(P_s+P_n+ε)`, giving
`λ̂ = P_s/(P_n+ε)`; the reference measures `N₀+ε` directly and yields the same. Both
report C/N₀ against the effective noise density including MAI, which is what other
receivers report and the more useful number for a lock detector. Summed over a
typical sky this is ≈0.9 dB for both — which is the figure quoted under "The
rotation policy", where it is the *base* the PRN-choice spread is measured against.
It is not a difference between the estimators and needs no correction. This is also
what removes the cross-satellite feedback loop an earlier sketch of this design had.

**The term that is *not* shared**, and the one that matters: the reference despreads
with a wrong PRN, so it also collects **the tracked satellite's own power**. NWPR
does not — it despreads with the correct PRN, where the satellite's power appears
properly in both `NBP` and `WBP` and cancels in the ratio. So the reference carries
a bias NWPR does not, and it scales with the satellite's own C/N₀:

| own C/N₀ | ε_self/N₀ | bias (C/N₀ reads low by) |
|---|---|---|
| 35 | 0.0031 | 0.013 dB |
| 40 | 0.0098 | 0.042 dB |
| 45 | 0.031 | **0.13 dB** |
| 50 | 0.098 | **0.40 dB** |

(L1 C/A, `f_chip` = 1.023 MHz; `10log10(1 + ε_self/N₀)`. A wider spreading bandwidth
scales it down proportionally, so GPS L5 and E5a at 10.23 Mcps are ≈10× smaller —
0.013 dB even at 50 dB-Hz. For BOC/CBOC signals the split spectrum is wider than the
chip rate suggests, so `C/f_chip` is an over-estimate there and E1B is somewhat
better than the L1 C/A figures above, not equal to them. Measure it rather than
assume it — that is what the Verification MAI sweep is for.)

**Consequences, and why this is still accepted:**

- Above ≈45 dB-Hz the reference is **bias-limited, not variance-limited**: 0.40 dB of
  bias at 50 dB-Hz against 0.06 dB of σ. The "keeps improving without bound" claim is
  true of σ only, and the Motivation table now says so.
- At and below 40 dB-Hz it is ≤0.042 dB and irrelevant — under NWPR's own +0.05 dB
  bias there. That is the range that motivated the change, and where NWPR's
  degenerate outputs and +0.45 dB low-end bias live.
- **Correcting it is worse than carrying it.** The correction is
  `N̂₀ ← N̂₀ − Σᵢ Ĉᵢ/f_chip` over the tracked satellites, which reintroduces exactly
  the cross-satellite feedback loop this design removed, and makes each satellite's
  C/N₀ a function of every other satellite's estimate. A stable, documented,
  monotone bias at the top of the range is the better trade.
- It **must not go undocumented**, because it is the one axis on which NWPR is
  better. The Phase 3 comparison table carries a "high-C/N₀ bias floor" row for
  this, and the Verification MAI item measures the NoiseRef−NWPR delta rather than
  assuming the two shift together.

**No in-estimator coherent pre-sum**, initially. Deferred — see below.

## Backends

| Backend | Sample-variance source | Correlator source |
|---|---|---|
| CPU / CPUThreaded | `Σ\|x\|²` over `Complex{Float32/64}` | `downconvert_and_correlate_fused!` (`src/downconvert_and_correlate_fused.jl:51`) |
| Int16 | same, integer accumulate | `_int16_hybrid_blocked!` (`src/downconvert_and_correlate_int16.jl:740`) |
| OneBit | `σ_q² ≡ 1` **analytically** — no pass needed | `_onebit_hybrid_blocked!` (`src/downconvert_and_correlate_onebit.jl:643`) |
| TwoBit | level count against `threshold` | `_twobit_hybrid_blocked!` (`src/downconvert_and_correlate_twobit.jl:1042`) |

The correlator source needs a real kernel call on every backend, which is the bulk of
step 1b's work. In exchange it is **model-free**: it traverses the identical
quantise→despread→accumulate path as the prompt, so the measured `N₀` already includes the
quantisation loss and the quantiser's operating point under load, with no per-backend
analysis and no closed-form correction. A `Σ|x|²` source would have needed one line per
backend and still missed the second-order coupling where a strong signal shifts that
operating point.

## Hardware producers (FPGA/ASIC)

The question a hardware implementer will ask is whether the noise measurement should reuse
the existing downconvert-and-correlate datapath or get a dedicated, minimal one. **Reuse
it — allocate a channel.** Not because it is less work, but because on the architectures
FPGA GNSS correlators actually use it is *also* cheaper, and because a separate path
forfeits the reason for measuring noise through a correlator at all.

### Why reuse is cheaper, not just better

- **Time-multiplexed engines**: most FPGA correlators run one physical datapath (carrier
  NCO, complex multiplier, code generator, accumulators) across N channels at a high clock,
  with per-channel state in BRAM. Adding a noise channel is **one more slot in the
  schedule** — a handful of state words and one time slot, essentially zero logic. A
  dedicated path is a whole second datapath, strictly more area.
- **Fully parallel arrays**: reuse costs *less* than one channel out of N (≈3 % at 32
  channels), because the reference is open-loop — see "The reference is open-loop". A noise
  channel is a tracking channel with the feedback disconnected: fixed NCO increments, one
  accumulator instead of three, no discriminator, no loop-filter interface, no
  dump-on-command sync. A dedicated path would have to re-implement the carrier NCO, code
  generator and complex multiplier to save the accumulators you already were not using, so
  it is *larger*, not smaller.
- **No divider, no float.** The producer reports the raw integer accumulation `Σ_m|B_m|²`
  with `M` and the sample count; every division, the `A_c²` and `f_s` scaling, and all
  averaging happen in Tracking. That is both the cheapest option on-chip and the mechanism
  that forces the two paths to agree.

### What reuse buys that a dedicated path cannot

The entire argument for a despread noise reference over a power meter is that it traverses
the **identical** chain as the prompt, so the measured `N₀` already contains every
imperfection of that chain and needs no model:

- identical quantisation and rounding, so 1-bit/2-bit loss cancels in the ratio;
- identical input scaling, so an AGC step moves numerator and denominator together — a
  dedicated path tapping at a different point (before a decimation filter, at a different
  bit width) breaks this silently;
- identical code generator, so `A_c` matches exactly, CBOC scaling included;
- identical saturation headroom under interference;
- one datapath to verify rather than two.

A separate path re-introduces exactly the per-implementation modelling this design exists to
avoid.

### What the producer must match

Only two things, because everything else is computed in Tracking:

1. **Report raw accumulations**, per the `noise_observation_from_correlator` builder:
   `Σ_m|B_m|²`, `M`, total samples. Not a density. They ride the same return path as the
   correlator outputs, so no new wire interface is needed — `append_noise_observation!` sits
   beside `append_correlator_output!` on the host.
2. **Use an off-peak code phase and rotate the PRN** — a PRN not currently tracked is
   simplest, since then any phase is safe. LFSR-generated Gold codes make this free; memory
   codes already hold every PRN in ROM.

**A producer that can only dump on the code epoch matches the software default exactly**, so
there is nothing to reconcile in the normal case. And where a producer does differ, it still
does not conflict. It is worth
being explicit, because "the two paths must behave the same" could be read as requiring
identical cadence. It does not: the producer's dump length and `window_duration` do **not**
have to match the software default. Behaving the same means *computed by the same code, on the
same scale, with the same bias* — which the shared window guarantees, because it is
`M`-weighted and time-bounded. What differs is only precision, and deliberately so:
degrading the software path to match the coarsest hardware would throw away a free 4×
variance improvement at high C/N₀ for nothing. The window is bounded in time and weighted by `M`, so a producer that can only dump
on the code epoch (1 ms) is handled correctly — it just delivers a smaller `M` per second.
The consequence is quantified in "How much averaging the reference needs": at one
sub-integration per code period and the default 1 s window, `K_n ≈ 1000`, versus 200 at a
0.2 s window — 0.32 dB versus 0.07 dB at 50 dB-Hz, and indistinguishable below 40 dB-Hz.
Three ways to recover the top end, if it matters:

1. **Several noise channels in parallel** at different PRNs, summed — multiplies `M`
   directly, and each is cheaper than a tracking channel because it is open-loop.
2. **A programmable dump counter**, which is a comparator on the code-phase accumulator.
3. **Accept it.** Below 40 dB-Hz there is no measurable penalty, and that is where lock and
   loss decisions are actually made.

A code-epoch dump is also DC-immune, a full-period Gold code being balanced by construction —
which is one of the three reasons the software default matches it.

## External APIs

New exports (**11**): `AbstractNoiseEstimator`, `CorrelatorNoiseEstimator`,
`requires_noise_density`, `NoiseObservation`, `noise_observation_from_correlator`,
`noise_observation`, `noise_observation_from_samples`, `append_noise_observation!`,
`update_noise!`, `get_noise_density`, `NoiseRefCN0Estimator`.

`update_noise!` is on that list deliberately, and an earlier revision left it off by
mistake. It is one of the three methods that *define* `AbstractNoiseEstimator`, and
the design leans on users implementing it — "anyone who wants [a hardware power
monitor alongside software correlation] can define their own `AbstractNoiseEstimator`
with a no-op `update_noise!`; the abstract type and its interface are public". An
unexported method cannot be extended by name, so leaving it out would have made that
escape hatch inaccessible and the interface only two-thirds public.

`append_noise_observation!` mirrors `append_correlator_output!`
(`src/sat_state.jl:273`, `src/tracking_state.jl:834-842`) with band-id selection:

```julia
append_noise_observation!(track_state, obs)              # single band
append_noise_observation!(track_state, obs, :L1)         # explicit band
```

The hardware path is then symmetric with what an FPGA producer already does:
`append_correlator_output!` for the taps, `append_noise_observation!` for the
noise, `estimate_dopplers_and_filter_prompt!` to fold.

`TrackState`'s kwarg constructor (`src/tracking_state.jl:43`) gains
`noise_estimators = nothing`, which provisions a `CorrelatorNoiseEstimator` **only for bands
that need one** — see below.

### Provision only where a C/N₀ estimator asks for it

A noise estimator is pure cost if nothing consumes its density. `MomentsCN0Estimator`,
`NWPRCN0Estimator` and `NoCN0Estimator` never read one; only `NoiseRefCN0Estimator` does. An
earlier draft defaulted to one per band unconditionally, which would have run a despread per
band per chunk — O(N), ~1-3 % of correlation work — for nobody, and shown up as an
unexplained regression in `SUITE["track"]` for anyone who kept NWPR.

So the requirement becomes a trait on the C/N₀ estimator:

```julia
requires_noise_density(::AbstractCN0Estimator) = false
requires_noise_density(::NoiseRefCN0Estimator) = true
```

`TrackState`'s constructor walks its groups' signals and provisions a `CorrelatorNoiseEstimator`
for a band only where some signal on it returns `true`. Bands with no requiring signal get **no
entry** in the `noise_estimators` NamedTuple, so `_update_band_noise!`'s tuple recursion does no
work for them and the cost is exactly zero — which is the answer for anyone who deliberately
stays on `MomentsCN0Estimator` or `NWPRCN0Estimator`.

The trait belongs on the public extension point rather than being a hard-coded type check,
because `update(::AbstractCN0Estimator, prompt, ::CN0UpdateContext)`
(`src/cn0_estimation.jl:92`) is already documented as user-extensible: a third-party estimator
that wants a density must be able to say so, and one that does not must not pay for it.

Two consequences worth stating:

- With the trait in place, the "no source configured" throw becomes reachable **only** by
  explicitly passing `noise_estimators` that omit a band whose signals need one — i.e. exactly
  the deliberate mistake the throw exists to catch, rather than something the defaults can
  produce.
- `add_satellite!` can in principle add a requiring signal to a band that had none (a band that
  was empty, or held only non-requiring estimators). Provisioning on add would mean rebuilding
  `TrackState`'s type; the warn-once at the skip site already diagnoses it, so leave it to that
  and document the limitation rather than building machinery for an edge case. Post-Phase 2 it
  is nearly unreachable anyway, since `default_cn0_estimator` returns a requiring estimator and
  any band with satellites will already have been provisioned.

### Breaking changes

Per `AGENTS.md`, judged against 6.0.1. Both land as `feat(cn0)!` with a
`BREAKING CHANGE:` footer:

- **`TrackState` gains a field and a type parameter.** Positional construction
  and the copy constructor (`src/tracking_state.jl:254`, which pins type params)
  change shape. Internal positional `TrackState(groups, doppler_estimator)` call
  sites to fix: `src/tracking_state.jl:82`, `:196`, `:205`, `:247`, and the two
  type-parameterised forms `:259` (inside the copy constructor) and `:413`.
  `src/tracking_state.jl:295` is `TrackState(track_state; groups = new_groups)` — the
  **kwarg** copy constructor, so the new field flows through it untouched and it needs
  no third argument, but it is the call that makes the copy constructor at `:254` the
  place the field must be threaded. (An earlier revision listed
  `src/tracking_state.jl:101`, `:131`, `:282` and `src/sat_state.jl:1088` here — those
  are all `SignalGroup` constructors and are **not** affected.)
- **`CN0UpdateContext` gains two fields.** It is exported
  (`src/Tracking.jl:63`), its 5-field positional constructor is called from
  `test/cn0_estimation.jl:184`, `:347`, `:381`, and
  `update(::AbstractCN0Estimator, prompt, ::CN0UpdateContext)`
  (`src/cn0_estimation.jl:92`) is a documented extension point. The convenience
  constructor (`:63-77`) keeps its current call shape by defaulting the new
  fields, so most callers are unaffected.

Explicitly **not** breaking: `estimate_cn0(estimator, integration_time)` keeps its
signature (`src/cn0_estimation.jl:137`, `:698`, `src/sat_state.jl:1195-1201`,
`src/tracking_state.jl:789-816`), and `NWPRCN0Estimator` / `MomentsCN0Estimator`
/ `NoCN0Estimator` are untouched.

## Implementation phases

### Phase 0 — baseline and harness (no `src/` changes) — **done**

`test/cn0_estimator_comparison.jl`: a seeded Monte-Carlo comparison of the shipped
`NWPRCN0Estimator` against reference non-coherent and coherent-`M` noise-reference
estimators, at known λ with the noise density known exactly. 4000 trials × 9 seeds
per point, medianed — the house convention at `test/cn0_estimation.jl:625-634`, and
necessary rather than decorative here: at one seed the sampling error on `bias_db`
is ≈0.014 dB at 30 dB-Hz against a 0.02 dB bound, so a single-seed version passed on
its own seed and failed on ≈13 % of others. At 9 seeds every assertion holds across
40 alternative seed bases with ≈40 % headroom; the measured worst cases are recorded
in the comment above `_sweep`. Runs in ≈4 s and prints the Motivation table above,
so it stays reproducible. Five testsets pin the
structural claims the plan rests on: NWPR's degenerate estimates, its bias at every
C/N₀, its σ saturation, the coherent noise reference beating it everywhere, and —
the honest half — the non-coherent squaring loss below the crossover.

Claims 1-3 are assertions about Tracking's *real* behaviour (the estimator is driven
through its public `update`/`estimate_cn0`), not about a model of it. When
`NoiseRefCN0Estimator` lands in step 1a, its measured bias and σ must match the
reference columns this file already pins.

**No manual benchmark baseline is needed.** An earlier draft of this phase called for
recording `SUITE["track"]` numbers by hand; `.github/workflows/benchmark_pr.yml`
already benchmarks the PR head against its base branch on both x86 and ARM, reporting
minimum times. A hand-captured baseline would be redundant and, being
machine-specific, misleading if a later phase runs elsewhere. Phase 1 just needs its
own `SUITE["noise estimation"]` group so the per-band pass is visible to that
comparison.

### Phase 0.5 — split the estimators into one file each (AFTER the harness, BEFORE any new code)

`src/cn0_estimation.jl` is 737 lines holding four estimators plus the shared
abstraction, and `NoiseRefCN0Estimator` would make it five. Split it first, as
pure motion, so the new estimator arrives as a new file and its diff is
reviewable instead of being buried in a reorganisation.

Mirror the existing `src/correlators/` convention exactly — a folder whose shared
file is the singular of the folder name, then one file per concrete type:

```
src/cn0_estimators/
  cn0_estimator.jl   # AbstractCN0Estimator, CN0UpdateContext, the 3-arg `update`
                     # extension point, default_cn0_estimator
  moments.jl         # MomentsCN0Estimator                         (:95-146)
  no_cn0.jl          # NoCN0Estimator                             (:148-190)
  nwpr.jl            # NWPRCN0Estimator, _num_ratios, _narrowband_window,
                     # _update_nwpr, _with_window_state, _with_open_window  (:192-709)
```

**The ranges include each type's docstring**, which is the whole point of stating
them: `NWPRCN0Estimator`'s alone is `:192-351`, 160 lines — more than the struct and
all its methods put together — and an earlier revision's `:352-709` would have left
it stranded in the shared file, splitting a type from its documentation and quietly
making the commit not-pure-motion. Line accounting for the 737-line original:
`nwpr.jl` takes 518 lines, `moments.jl` 52, `no_cn0.jl` 43, leaving ~124 for
`cn0_estimator.jl`.

Constraints and consequences:

- **Pure motion only.** No behaviour, name, signature or docstring changes in this
  commit, so it lands as `refactor(cn0):` — no release per `AGENTS.md`. Verify with
  a full test run, not by inspection.
- **Include order.** The shared file must still follow `bit_buffer.jl`, because
  `CN0UpdateContext` has a `BitBuffer` *field* (`src/cn0_estimation.jl:53`) and
  that is a compile-time dependency. The explanatory comment at
  `src/Tracking.jl:201-203` moves with it and needs rewording for the new paths
  (the two `include`s it explains are `:204-205`).
  `default_cn0_estimator` only *references* `NWPRCN0Estimator` and
  `MomentsCN0Estimator` inside a function body, so it is resolved at call time and
  can stay in the shared file regardless of include order.
- **Docs are unaffected.** `docs/src/cn0_estimator.md` uses `@docs` blocks keyed on
  symbols, not paths.
- Do the same for the noise sources when they arrive in Phase 1, rather than
  creating `src/noise_estimation.jl` and splitting it later:

```
src/noise_estimators/
  noise_estimator.jl   # AbstractNoiseEstimator, NoiseObservation + its three builders,
                       # the update_noise!/get_noise_density/
                       # append_noise_observation! interface
  correlator.jl        # CorrelatorNoiseEstimator
```

Tests follow the same split (`test/cn0_estimators/`, `test/noise_estimators/`),
each still its own `module XyzTest` with explicit import lists, added to the flat
`include` list in `test/runtests.jl`.

### Phase 1 — architecture and both fill paths

Two landable steps, so there is a working, tested intermediate state at each commit. The
whole of Phase 1 is opt-in, because the default `cn0_estimator` still comes from `default_cn0_estimator`
(post-0.5: `src/cn0_estimators/cn0_estimator.jl`) until Phase 2.

#### Step 1a — architecture and the externally-provided path

`src/noise_estimators/{noise_estimator,correlator}.jl`, the `TrackState` field,
`_update_band_noise!` at the three entry points, the `CN0UpdateContext` fields and the
threading through the six private functions, `src/cn0_estimators/noise_ref.jl`,
`append_noise_observation!`.

`append_noise_observation!` is what keeps this step end-to-end testable
before any kernel work exists: a test pushes synthetic `NoiseObservation`s at a known `N₀`
and checks that `NoiseRefCN0Estimator` reproduces the bias and σ that
`test/cn0_estimator_comparison.jl` already pins. That is also exactly the hardware/FPGA
path, so it lands first rather than as an afterthought. Step 1b then adds `update_noise!` —
the only new code is the measurement; the window it feeds is already tested.

#### Step 1b — `CorrelatorNoiseEstimator`

One extra despread per band per chunk, using **all** the correlator's taps at ≥1 chip spacing
(see "Use every tap, spaced wide"): rotating PRN, code phase offset from the
tracked phase of that PRN (preferring PRNs not currently tracked, where any phase
is safe), nominal IF as the carrier frequency and zero Doppler — ±5 kHz is 0.5 %
of a 1 MHz code band, so no satellite's state is needed and the reference runs
even with nothing tracked.

**Keep the downconversion.** Skipping it is tempting (the noise power is
invariant under the carrier multiply, so no rescaling would be needed) but it
moves the reference's passband to DC ±f_chip — into ADC offset, 1/f and clock
spurs, and for real IF into a different part of the spectrum than the signal
occupies, where any filter tilt biases `N₀` silently. It also costs √2 in the
reference's own std, since a real output's `Var(B²) = 2v²` against a complex
`Var(|A|²) = v²`. And it saves nothing: this is *one* correlation per band, so
~1/(n_sats·n_signals) of correlation cost either way, and reusing the existing
single-satellite template (`src/downconvert_and_correlate_cpu.jl:974-1012`) is
less new code than a bespoke no-downconvert path.

**The rotation policy, and how it avoids a counter.** Advancing a PRN needs a position, and a
position is scalar state with no per-band anchor (see "Immutability and allocation"). The
solution is the same trick as the FIFO: carry it in the buffer. `NoiseObservation` gains a
`prn` field — legitimate metadata, since the observation really was measured with that PRN,
and useful to the diagnostic follow-up — and the next PRN is derived from `last(buffered).prn`,
advancing to the next PRN of the chosen code family that is **not currently tracked on that
band** (derivable each chunk from `SignalGroup.satellites`, no state). An empty buffer starts
at the lowest untracked PRN. If every PRN is tracked, reuse any and offset the code phase by
half a code period from that satellite's tracked phase.

Rotation is worth this much and no more. It is not load-bearing: the leakage term is a sum
over ~12 satellites, so it is already averaged over 12 cross-correlation draws. What a fixed
PRN costs is a *slowly drifting* bias — with Gold cross-correlation power roughly
exponentially distributed, the sum's relative spread across PRN choices is ≈1/√12 ≈ 29 %,
which on the ≈0.9 dB whole-sky leakage bias derived in "Multi-access interference and the
reference's self-leakage" is ≈±0.28 dB that moves as the constellation geometry does.
Rotation replaces that drift with a stable mean. Since the leakage bias is deliberately
uncorrected, a stable bias is the whole objective.

Note which part of that bias rotation can and cannot touch. The ≈0.9 dB is dominated by the
term NWPR carries identically — the other satellites' interference — and rotation only
stabilises *which* draw of it you get. The ≈0.13 dB self-leakage term is not affected at all:
every wrong PRN collects the tracked satellite's power to the same expected degree, so no
choice of reference code removes it. That is why it is treated as a documented bias rather
than as something the rotation policy is expected to solve.

Rotating the PRN randomises the deterministic sidelobe/cross-correlation term
rather than reducing its mean — it buys tail-risk reduction (Gold sidelobes span
≈6 dB between RMS and peak), not bias and not variance. The estimate stays
**strictly per band**; only the borrowed *code* rotates.

**Which code family, when two groups share a band?** This needs pinning, because a band is
one-to-many with groups (`:legacy_gps` = `(GPSL1CA(),)` and `:galileo` = `(GalileoE1B(),)`
are both on `L1()`, `src/sat_state.jl:936-947`) and `_validate_signal_group` only guarantees
a single chip rate *within* a group. Rule: **the first signal of the first group on that
band, in `groups` order** — deterministic, needs no new configuration, and its chip rate is
whose primary code period sets the sub-integration length. For white noise the choice is
immaterial
(`N̂₀ = σ²/f_s` is code-independent once `A_c` is divided out), so this only selects *which*
spectral weighting a tilted front end is measured with. Do **not** rotate across code
families within a band: that would mix two weightings in one window, and a BOC replica
changes the spectral
weighting (harmless for white noise, not for a tilted front end).

Scratch: a dedicated `Int8` replica slot rather than borrowing
`_with_code_replica_buffer` (`src/downconvert_and_correlate_cpu.jl:155-162`), to
avoid aliasing the per-signal path. Mind that `gen_code_replica!` writes *at*
`start_sample` while the kernels read the replica offset by `start_sample - 1`
(`src/downconvert_and_correlate_fused.jl:470-476`).

### Phase 2 — validate, then flip the default

Change `default_cn0_estimator` (post-0.5: `src/cn0_estimators/cn0_estimator.jl`) to return
`NoiseRefCN0Estimator` when a band reference exists. `NWPRCN0Estimator` stays
exported and documented as the estimator to pass explicitly for
externally-supplied correlator outputs without a noise observation — the one
place it remains necessary.

### Phase 3 — comparison and selection docs (final step)

With four C/N₀ estimators and two ways to feed the noise reference, "which one do I use"
becomes the
main documentation question. `docs/src/cn0_estimator.md` is currently written as an
extended justification of *one* default (`## Default Estimator`, `### Why not the
moment method`, `### The window follows the navigation bits`, …). Restructure it
around a comparison and a decision, and add a companion page for the noise sources.

**Hard gate first:** `docs/make.jl:10` sets `checkdocs = :exports`, so every new
exported symbol needs a docstring or the docs build fails. That is 11 new exports.
Add `noise_estimator.md` to the `pages` list at `docs/make.jl:11-21`.

`docs/src/cn0_estimator.md` gains a comparison table up front, one row per
estimator, with these columns — chosen because each one is a reason somebody picks
a different estimator:

| Axis | Why it decides |
|---|---|
| Needs band samples? | Rules out `NoiseRef` on a pure correlator-ingest path with no noise observation |
| Needs bit sync / a coherent window? | Rules out `NWPR` on L1C-D, E1B and pre-sync secondary-coded signals |
| Needs to know `M`? | `NWPR` only; couples the estimator to the bit grid |
| Bias (thermal only) | `NWPR` ≈ +0.05 dB everywhere; `Moments` ≈27.6 dB-Hz floor on noise; `NoiseRef` <0.03 dB |
| High-C/N₀ bias floor | **`NoiseRef` only**: self-leakage, ≈0.13 dB at 45 dB-Hz and ≈0.40 dB at 50. The one axis on which `NWPR` is better — see "Multi-access interference and the reference's self-leakage" |
| High-C/N₀ σ | `NWPR` saturates at `√(M/((M−1)K))`; `NoiseRef` keeps improving (until the row above dominates) |
| Degenerate outputs | `NWPR` reports -Inf/Inf dB-Hz outside `1 < μ̂ < M` — 5.6 % of *estimates* at 20 dB-Hz, the pooled ratio over all windows, not a per-window discard rate |
| Phase-noise sensitivity | Non-coherent `NoiseRef` and `Moments` are immune; `NWPR` and long coherent records are not |
| Cost | Per-record arithmetic, plus O(N) per band per chunk for a software noise source |

Then a selection guide as an explicit decision list, not prose:

1. Sample-driven `track!` → `NoiseRefCN0Estimator`. This is the default and needs no
   configuration.
2. Correlator outputs from hardware **with** a noise observation → the same, plus
   `append_noise_observation!` per band. The docs must carry the "Hardware producers"
   guidance: reuse the tracking datapath with its feedback disconnected, report raw
   accumulations.
3. Correlator outputs from hardware **without** one → `NWPRCN0Estimator`, accepting
   its fallback on the signals it cannot serve.
4. A signal whose C/N₀ nobody reads (the passenger of a co-tracked pair) →
   `NoCN0Estimator`, which already has its own section at
   `docs/src/cn0_estimator.md:235`.
5. `MomentsCN0Estimator` — document it as a fallback component and a prompt-stream
   tool, not a recommended standalone default, and say plainly why: the ≈27.6 dB-Hz
   floor on noise.

Two further sections, both carrying material that has no home today:

- **"Getting the low-C/N₀ variance"** — the coherent-gain point, with the
  measured table from the Motivation section of this plan. A coherent sum of five
  1 ms prompts *is* one 5 ms prompt, so the gain belongs in
  `preferred_num_code_blocks_to_integrate`, not in a post-correlation window; and
  spending it there also helps the discriminators. Must state the coupling honestly:
  record length also sets discriminator behaviour and loop bandwidth, so this is a
  trade, not a free win. Include the pre-sync caveat (1-block cap,
  `src/sat_state.jl:36-40`).
- **"What C/N₀ means here"** — extend `## What the estimator returns`
  (`docs/src/cn0_estimator.md:267`) to say that every estimator reports against the
  *effective* noise density including multi-access interference, that `NWPR` and
  `NoiseRef` shift together for that reason, and that the post-corr filter is assumed
  to have unity noise gain (`src/post_corr_filter.jl:27`) — so a beamformer
  invalidates `NoiseRef` until `noise_gain` exists.

New `docs/src/noise_estimator.md`: the `N₀` invariant and why a density rather than a
power, the two fill paths and why they need no type distinction, the derived
sub-integration length and why it is not tunable, what a hardware producer must supply, and
the `NoiseObservation` builders as the
hardware-integration recipe with a worked FPGA example.

Lands as `docs(cn0):` — no release per `AGENTS.md`.

## Files

New:

| Path | Contents |
|---|---|
| `src/cn0_estimators/cn0_estimator.jl` | moved: `AbstractCN0Estimator`, `CN0UpdateContext`, 3-arg `update`, `default_cn0_estimator` |
| `src/cn0_estimators/moments.jl` | moved: `MomentsCN0Estimator` |
| `src/cn0_estimators/no_cn0.jl` | moved: `NoCN0Estimator` |
| `src/cn0_estimators/nwpr.jl` | moved: `NWPRCN0Estimator` + window helpers |
| `src/cn0_estimators/noise_ref.jl` | new: `NoiseRefCN0Estimator` |
| `src/noise_estimators/noise_estimator.jl` | new: `AbstractNoiseEstimator`, `NoiseObservation`, interface |
| `src/noise_estimators/correlator.jl` | new: `CorrelatorNoiseEstimator` |

Deleted: `src/cn0_estimation.jl`.

Modified:

| Path | Change |
|---|---|
| `src/Tracking.jl` | `include`s (`:205`) and the order comment (`:201-205`); exports (`:59-63`, `:106`); `TrackState` struct (`:232-235`) |
| `src/tracking_state.jl` | `TrackState` constructors (`:43`, `:184`, `:199`, `:238`, `:254`) and internal call sites (`:101`, `:131`, `:282`); `append_noise_observation!` alongside `:834-842` |
| `src/sat_state.jl` | `TrackState` call site (`:1088`) |
| `src/downconvert_and_correlate_cpu.jl` | `_update_band_noise!` + call at `:856`; shared helper on `AbstractDownconvertAndCorrelator` |
| `src/downconvert_and_correlate_onebit.jl` | call at `:1349` |
| `src/downconvert_and_correlate_twobit.jl` | call at `:1838` |
| `src/conventional_pll_and_dll.jl` | `_band_noise_density` beside `:694-696`; threading through `:284`, `:413`, `:499`, `:600`, `:620`; `CN0UpdateContext` build at `:459-463` |
| `src/vector_pll_and_dll.jl` | threading through `:264`, `:295` |
| `docs/src/cn0_estimator.md` | restructured around the comparison table + selection guide (Phase 3); `@docs` entries for the new types |
| `docs/src/noise_estimator.md` | new page: the `N₀` invariant, source comparison, hardware-integration recipe |
| `docs/make.jl` | `noise_estimator.md` in `pages` (`:11-21`); `checkdocs = :exports` (`:10`) means all 11 new exports need docstrings |
| `test/runtests.jl`, `test/cn0_estimators/`, `test/noise_estimators/` | split + new test modules |
| `benchmark/benchmarks.jl` | `SUITE["noise estimation"]`, feature-gated |

## Verification

- **The Phase 0.5 split is pure motion**, so the gate is the full test suite plus
  Aqua passing byte-identically before and after, with no test edited other than
  its module path. If a test needs changing, the commit is not pure motion and the
  behaviour change belongs in its own commit.
- **Harness reproduction** (Phase 0): bias and σ vs C/N₀ for every estimator at
  known `N₀`, seeded, median-over-seeds with the measured spread written into the
  comment — the convention at `test/cn0_estimation.jl:625-634`.
- **Source equivalence**: `noise_observation_from_correlator` and
  `noise_observation_from_samples` agree on the same white input; the software
  and external paths give bit-identical C/N₀ when fed the same observations.
- **Closed-loop sweep** using `track_noisy` (`test/cn0_estimation.jl:579-615`):
  15→45 dB-Hz in 1 dB steps, `T` = 1 / 4 / 20 ms, bias and σ, on CPU, Int16,
  OneBit and TwoBit. This is where the **phase-noise wash** is measured — a 5 ms
  coherent record loses as much to residual phase noise as a 5-record coherent
  sum does, so the "5 ms" column inherits NWPR's coherence-loss bias
  (`src/cn0_estimation.jl:290-310`); only the 1 ms non-coherent form is immune.
  That axis is deliberately absent from the Phase 0 model and is the one number
  that must come from a real loop before the default flips.
- **Signals NWPR cannot serve**: GPS L1C-D and Galileo E1B, plus a
  secondary-coded signal pre-sync, must land inside 1 dB with no floor — the #217
  acceptance criterion.
- **MAI, and the self-leakage delta**: a 12-satellite constellation. The shared
  term is confirmed as before — both read against `P_n + ε`, so they shift together
  on the *other* satellites' interference. But the test must **measure the
  NoiseRef−NWPR difference** rather than assert it is zero, because it is not: the
  reference also collects the tracked satellite's own power, which NWPR does not.
  Sweep the satellite under test across 35/40/45/50 dB-Hz and pin the delta against
  the predicted `10log10(1 + (C/N₀)/f_chip)` — 0.013 / 0.042 / 0.13 / 0.40 dB — to
  within a few hundredths. That both validates the model in "Multi-access
  interference and the reference's self-leakage" and turns the design's one bias
  term into a regression rather than a footnote. It is also the gate on the Phase 3
  comparison table's "high-C/N₀ bias floor" row being accurate.
- **Multi-band**: L1 + L5 with different `f_s`, confirming the densities stay
  separate.
- **Immutability**: `grep -rn "^mutable struct" src/` stays empty, and
  `ismutabletype` is `false` for every new type — a one-line test over the new
  types, since nothing else enforces it.
- **Tap independence**: the measured power correlation between the reference's taps stays
  below ~0.01, i.e. the spacing really is ≥1 chip after `get_correlator_sample_shifts`
  rounds it to whole samples at the band's `f_s`. A regression on that rounding, not on the
  code statistics.
- **No provisioning without a consumer**: a `TrackState` whose signals all use
  `MomentsCN0Estimator`/`NWPRCN0Estimator`/`NoCN0Estimator` gets an empty `noise_estimators`,
  performs no despread, and shows no change against `SUITE["track"]`. The mixed case — one band
  requiring, another not — provisions exactly one.
- **Configured-but-never-fed**: a band with a `CorrelatorNoiseEstimator` that is never
  appended to and never measured warns exactly once and leaves C/N₀ at `-Inf dB-Hz` without
  throwing; and the software path never emits that warning. This is the likeliest hardware
  integration mistake, so it gets its own regression. With **two** misconfigured bands, both
  must warn — that is what pins the per-band `_id` on the `@warn`, since a bare
  `maxlog = 1` is keyed per callsite and would silence the second band.
- **Mixed estimators on one band**: a band carrying one `NoiseRefCN0Estimator` signal and one
  `NWPRCN0Estimator` signal, driven through a warm-up where the density is unavailable. The
  NWPR signal's `num_records_per_ratio` and `ratios_are_bit_aligned` must be unaffected — the
  regression against skipping the whole band's C/N₀ fold, which would drop NWPR's open window
  on every skipped record (`src/cn0_estimation.jl:554-557`) and silently demote it to its
  fallback.
- **Observation cadence**: over a buffer of `C` chunks, exactly `C · num_sub` observations
  reach each band's window, none from the final drain pass, and `num_sub` tracks
  `chunk_duration / code_period` across at least 1 ms, 4 ms and 20 ms Doppler update
  intervals and across two sampling frequencies, and collapses to 1 on a band whose reference
  code period exceeds the chunk. This is the regression that catches a wrong `samples_unchanged` gate,
  which would otherwise freeze the window silently after chunk 0, and a `num_sub` that
  collapses to zero on a band whose code period exceeds the chunk. Cover the **direct
  `downconvert_and_correlate!` caller** as well as `track!`, since the gate reads a public
  kwarg whose documented meaning is a cache hint (see "Integration boundaries"): all four
  `track!` shapes, plus an explicit `samples_unchanged = true` unchunked call, whose skip is
  intended but surprising.
- **Per-band window**: `@allocated == 0` over ≥1 M pushes at each source's default
  window length, which is what pins the `sizehint!` headroom (see "Where per-band
  accumulation is anchored"). Also that a band with two `SignalGroup`s is measured once,
  and that two bands at different `f_s` keep separate windows.
- **Allocation**: `@allocated == 0` in steady state for `update_noise!`,
  `update(::NoiseRefCN0Estimator, …)` and `append_noise_observation!`, following the
  existing idioms — the direct form at `test/cn0_estimation.jl:104`, and the
  wrap-the-measurement-in-a-function form at `test/track_in_place.jl:84-118` (needed
  because `@allocated` at module scope picks up boxing from untyped global lookups,
  `test/bit_buffer.jl:43-47`). Extend the existing per-stage assertions in
  `test/track_in_place.jl:119` so a regression shows up against
  `downconvert_and_correlate!` and `estimate_dopplers_and_filter_prompt!` as wholes,
  not only in isolation. The threaded backends already carry a documented small
  residual (`test/track_in_place.jl:139-142`) — match that tolerance rather than
  tightening it.
- `@inferred` on `update` / `estimate_cn0` / `get_noise_density` for every source
  and both `CN0UpdateContext` parametrisations; Aqua; `format(["src", "test"])`
  before committing.
- `benchmark/benchmarks.jl`: a new `SUITE["noise estimation"]` group, feature-gated
  in the established style (`benchmark/benchmarks.jl:41-63`), plus confirmation
  that `SUITE["track"]` is unchanged within noise.

## Risks / watch-items

- **Post-corr filter noise gain.** `N̂₀` is measured at the band, but the prompt
  reaching the estimator is post-`PostCorrFilter`
  (`src/conventional_pll_and_dll.jl:424-427`). `DefaultPostCorrFilter` is
  `last(x)` (`src/post_corr_filter.jl:27`), gain 1, fine. A real beamformer
  changes the prompt's noise scale and silently invalidates the ratio. The same item covers
  **which antenna the reference despreads** when `num_ants > 1`: it must be the one
  `DefaultPostCorrFilter` selects (`last`), or the ratio compares different channels. Phase 1
  documents unity as a requirement; a `noise_gain(filter)` hook is the follow-up.
- **Pre-sync record length.** `preferred_num_code_blocks_to_integrate` is capped
  at 1 block pre-sync (`src/sat_state.jl:36-40`), so the coherent gain in the
  table's third column is unavailable there and pre-sync sits on the 1 ms column.
  NWPR's pre-sync path buys its lower σ with a bias the harness does not charge
  it: a free-running 5 ms window on L1 C/A straddles a bit transition ≈20 % of
  the time, and a straddle collapses NBP. Unbiased-but-noisy is the better trade
  for a sync decision, but this should be stated in the docs, not assumed.
- **Tuning coupling.** Recommending a longer record length to recover the
  coherent gain couples C/N₀ quality to the loop's integration time, which also
  sets discriminator behaviour and loop bandwidth. It is a real coupling, not a
  capability gap — but the docs must say so rather than presenting the 5 ms
  column as free.
- **AGC contemporaneity.** `N̂₀` must reflect the samples a record actually
  integrated. Records span chunks, so if a front end re-gains mid-buffer the
  EMA lags. Acceptable at ≈1 s time constants; revisit if a fast-AGC front end
  shows up.
- **DC-offset bias on a partial code period** — a non-default concern, since `T_sub` now
  is a full, DC-balanced code period. It applies only in the one configuration named
  above — a band whose reference code period exceeds `doppler_update_interval`, so the
  sub-integration is a *fraction* of a period. A despread over `W` chips accumulates the
  code sum `Σc`, and for `W ≪ 1023` a Gold code behaves like iid signs, so
  `E[(Σc)²] ≈ W` (measured: 16.3 at `W` = 16, 57.0 at `W` = 64, versus **1.0** over the full
  1023-chip period, which is balanced by construction). An ADC DC offset `d` therefore adds
  a *positive* term to every `|B|²` — it does not average away, because the sign of `Σc`
  is squared — inflating `N̂₀` and making C/N₀ read **low** by
  `10log10(1 + (d/σ)²·f_s/f_chip)`, independent of `W` for any `W ≪ 1023`:

  | `d/σ` | short sub-integration | full code period |
  |---|---|---|
  | 0.01 | +0.010 dB | +0.000 dB |
  | 0.022 | +0.049 dB | +0.000 dB |
  | 0.05 | +0.248 dB | +0.000 dB |
  | 0.10 | +0.915 dB | +0.001 dB |

  Mitigation is mostly free: the reference **downconverts**, so an ADC DC offset lands at
  `−f_IF` and the baseband code rejects it whenever `f_IF` is outside the code band. This is
  therefore a **zero-IF / complex-baseband concern only**. For such a front end either keep
  `d/σ ≲ 0.02`, subtract a running DC estimate before despreading, or set
  `doppler_update_interval` to at least the band's reference code period, so whole periods are
  integrated again. Worth a test with an injected offset either way.
- **`Int16` has no own entry point** (`src/downconvert_and_correlate_int16.jl:1249-1256`),
  it inherits CPU's. Confirm `_update_band_noise!` is reached exactly once there
  and that the sample type check does not reject the noise pass.

## Deferred follow-ups

- **Optional coherent pre-sum of `M` records** inside
  `NoiseRefCN0Estimator`, for the low-C/N₀ variance without lengthening records.
  The reason this would not recreate NWPR's problem: `M = 1` is a valid operating
  point — same formula, same code path, unbiased, only noisier — whereas NWPR at
  `M = 1` is undefined (`NBP ≡ WBP`, `μ̂ ≡ 1`). Graceful degradation instead of a
  cliff. Add only if the closed-loop sweep shows the 20 dB-Hz gap matters.
- **Expose `N̂₀` to acquisition** for a real CFAR threshold, replacing the
  search-grid-mean noise reference.
- **`noise_gain(::PostCorrFilter)`** so beamforming setups are supported rather
  than documented-as-unsupported.
- **A `Σ|x|²` second source, purely as a diagnostic.** Dropped from the shipped design to
  keep one algorithm, but running it alongside the correlator source makes their ratio a
  direct read-out of multi-access leakage plus front-end spectral tilt, for one extra O(N)
  pass. Predicted ratio for 12 satellites at 45 dB-Hz on L1 C/A is ≈1.30, because a
  correlator reference concentrates every satellite's power in the code band
  (`Σᵢ(C/N₀)ᵢ/f_chip`) while `σ̂²` spreads it over the sampled band
  (`Σᵢ(C/N₀)ᵢ/f_s`) — smaller by `f_s/f_chip`. Residue beyond that prediction is tilt or an
  in-band spur, which has no analytic form. Worth adding if the correlator source's
  calibration is ever in doubt; not worth a second supported code path before then.
