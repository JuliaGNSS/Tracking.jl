# CN0 Estimator

The CN0 (Carrier-to-Noise density ratio) estimator provides a measure of
signal quality in dB-Hz. Each [`TrackedSignal`](@ref) on a
[`TrackedSat`](@ref) holds its own CN0 estimator, so a multi-signal
satellite produces one CN0 value per signal.

!!! note "Every signal of a satellite is estimated, not just the driver"

    The estimator of *each* signal is fed its own prompt once per completed
    record — the [estimator-driver signal](tracking_state.md#Estimator-driver-signal)
    has no privileged role here, and the per-signal values are not
    interchangeable: on GPS L1C the pilot carries ~75 % of the power against the
    data component's ~25 %, a 4.8 dB difference.

    Nothing inside Tracking reads [`estimate_cn0`](@ref) — it is a reporting API,
    so this costs only the per-record estimator update, which is well under a per
    cent of the correlation work for that record (~40 ns against ~9 µs for a
    two-signal satellite over 4000 samples). `estimate_cn0` itself is of the same
    order and runs only when you call it, so a signal whose C/N₀ you never read
    costs almost nothing either way.

    The reason to switch a signal off is therefore not speed but honesty.
    [`NoCN0Estimator`](@ref) is the per-signal opt-out.

## Choosing an estimator

Four estimators ship, and each column below is a reason somebody picks a
different one.

|                                                            | [`NoiseRefCN0Estimator`](@ref) (default)                | [`NWPRCN0Estimator`](@ref)                                                                                                                                                 | [`MomentsCN0Estimator`](@ref) | [`NoCN0Estimator`](@ref) |
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
| Cost                                                       | per-record arithmetic + one despread per band per chunk | per-record arithmetic                                                                                                                                                      | per-record arithmetic         | none                     |

The "degenerate" row is a rate of unusable *outputs*, not a per-window discard
rate: [`estimate_cn0`](@ref) pools every buffered window and range-checks the
pooled ratio once.

### Which one do I use?

 1. **Sample-driven [`track!`](@ref)** → [`NoiseRefCN0Estimator`](@ref). This is
    the default and needs no configuration: the band is provisioned a
    [`CorrelatorNoiseEstimator`](@ref) automatically and `track!` fills it.
 2. **Correlator outputs from hardware, *with* a noise observation** → the same,
    plus one [`append_noise_observation!`](@ref) per band per fold. See
    [Noise Estimator](noise_estimator.md) for what a producer has to report; the
    short version is *raw accumulations of an untracked PRN*, with every division
    and all averaging done here.
 3. **Correlator outputs from hardware, *without* one** →
    `NWPRCN0Estimator(signal)`, accepting its fallback on the signals it cannot
    serve. This is the one place NWPR is still the right choice.
 4. **A signal whose C/N₀ nobody reads** (the passenger of a co-tracked pair) →
    [`NoCN0Estimator`](@ref); see [Not measuring a signal at all](#Not-measuring-a-signal-at-all).
 5. [`MomentsCN0Estimator`](@ref) is a fallback component and a prompt-stream
    tool, not a recommended standalone default — see
    [Why not the moment method](#Why-not-the-moment-method) for the ≈27.6 dB-Hz
    floor that rules it out.

## The default: C/N₀ against a measured noise floor

[`NoiseRefCN0Estimator`](@ref) computes, per record,

```
Ĉ/N₀ = ⟨|P|²⟩ / N̂₀ − 1/T
```

with `N̂₀` the band's measured noise **density** and `T` that record's own
integration time, then averages the per-record terms. It is **non-coherent**: no
bit sync, no window, no `M`, no fallback, no saturation ceiling — one set of
characteristics on every signal.

Where the density comes from is a separate, pluggable thing: see
[Noise Estimator](noise_estimator.md). On the sample-driven path it is automatic.

```@docs
NoiseRefCN0Estimator
Tracking.default_cn0_estimator
```

### Getting the low-C/N₀ variance

Below about 28 dB-Hz a *non-coherent* estimator is the noisier of the two at
matched record length, and the reason is squaring loss: NWPR sums `M` prompts
coherently before squaring, so at low SNR its σ is better by exactly `√(M−1)`.

The way to buy that back is not a post-correlation window. **A coherent sum of
five 1 ms prompts *is* one 5 ms prompt** — NWPR's narrowband window is a longer
coherent integration reconstructed after the fact, not extra information. Taking
the same integration in the *correlator* instead, with
[`set_preferred_num_code_blocks_to_integrate!`](@ref), makes the noise reference
better than NWPR at every C/N₀ and helps the discriminators as well. Modelled
over 100 ms of observation, 4000 trials × 9 seeds, relative σ in dB:

| C/N₀ | NWPR (M = 5) | NoiseRef, 1 ms records | NoiseRef, 5 ms records |
|:---- | ------------:| ----------------------:| ----------------------:|
| 20   | 2.96         | 4.79                   | **2.76**               |
| 25   | 1.50         | 1.72                   | **1.24**               |
| 30   | 0.90         | 0.75                   | **0.65**               |
| 40   | 0.54         | 0.20                   | **0.20**               |
| 50   | 0.50         | 0.06                   | **0.06**               |

Two honest caveats. The record length also sets the discriminator behaviour and
the loop's update interval, so this is a **trade, not a free win** — the
conventional estimator rescales the carrier bandwidth by `1/N` to keep the loop
stable, but a longer record still means a slower loop. And the length is capped
at one block until bit/secondary sync is found, so pre-sync you are on the 1 ms
column whatever you configure.

### The one bias it carries

The reference despreads with a *wrong* PRN, so besides the thermal floor it also
collects **the tracked satellite's own power** — `ε_self/N₀ = C/f_chip`. NWPR
does not: it despreads with the correct code, where the satellite's power appears
in both `NBP` and `WBP` and cancels in the ratio. On GPS L1 C/A:

| own C/N₀ | C/N₀ reads low by |
|:-------- | -----------------:|
| 35 dB-Hz | 0.013 dB          |
| 40 dB-Hz | 0.042 dB          |
| 45 dB-Hz | **0.13 dB**       |
| 50 dB-Hz | **0.40 dB**       |

A wider spreading bandwidth scales it down in proportion, so GPS L5 and Galileo
E5a at 10.23 Mcps are ≈10× better — 0.013 dB even at 50 dB-Hz.

It is carried rather than corrected on purpose: the correction
`N̂₀ ← N̂₀ − Σᵢ Ĉᵢ/f_chip` would make every satellite's C/N₀ a function of every
*other* satellite's estimate, which is a cross-satellite feedback loop this
design exists without. Below 40 dB-Hz the term sits under NWPR's own +0.05 dB
bias, and that is the range where lock and loss decisions are made.

## `NWPRCN0Estimator`: the coherent-window estimator

Van Dierendonck's **narrowband/wideband power ratio**. Per narrowband window of
`M` records it forms the coherent power `NBP = |Σ prompt|²` and the incoherent
power `WBP = Σ |prompt|²`, divides the sums of both over as many windows as fit
in `num_prompts_for_cn0_estimation` records, and turns that mean power ratio into
a C/N₀. It was the default until the per-band noise reference landed, and it is
what to configure on a correlator-ingest path that cannot report a noise
observation.

```@docs
NWPRCN0Estimator
NWPRCN0Estimator(::GNSSSignals.AbstractGNSSSignal)
```

### Why not the moment method

The [Moments Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1)
(M2M4, [`MomentsCN0Estimator`](@ref)) was the default up to and including
Tracking 5.1. It is a *moment ratio*, and at a finite window the sample moments
fluctuate enough to manufacture signal power out of noise. Measured on
synthetic prompts (amplitude `√(C/N₀·T)` in unit-variance complex noise), 100
one-millisecond records, 400 realizations:

| true C/N₀  | M2M4 median | p10  | p90  | NWPR (`M = 20`) median |
|:---------- | -----------:| ----:| ----:| ----------------------:|
| noise only | **27.6**    | 22.6 | 31.1 | **−Inf** ("no signal") |
| 20 dB-Hz   | 27.7        | 22.5 | 30.9 | 19.9                   |
| 25 dB-Hz   | 27.9        | 23.4 | 30.9 | 25.0                   |
| 30 dB-Hz   | 30.5        | 27.3 | 32.6 | 30.0                   |
| 35 dB-Hz   | 35.3        | 33.9 | 36.4 | 35.0                   |
| 45 dB-Hz   | 45.1        | 44.3 | 45.9 | 45.0                   |

M2M4 is accurate from ~35 dB-Hz up and blind below ~30: a true 20 dB-Hz signal
and pure noise both read ~27.7. As a per-update detection statistic that makes
pure noise clear a 30 dB-Hz threshold 19 % of the time, a 28 dB-Hz threshold
45 %, and a 25 dB-Hz threshold 76 % — so the usual code-lock hysteresis (lock at
34 dB-Hz, drop at 25) cannot be built on it. At `M = 20` NWPR's false-alarm rate
against any threshold from 20 to 32 dB-Hz is 0 %, and its spread at low C/N₀ is
much tighter (p10–p90 of 15.8–22.2 at a true 20 dB-Hz). This mirrors the
published comparisons (Falletti, Pini & Lo Presti, *IEEE T-AES* 47(1):420–437,
2011). See [issue #217](https://github.com/JuliaGNSS/Tracking.jl/issues/217) for
the full measurements.

Those are the numbers for one 20-record window; the shipped default runs shorter
windows (see below), which trades some of that margin for immunity against the
loop's phase noise. End to end through `track` on a data-modulated GPS L1 C/A
signal at the defaults, 600 ms per run, medians over 16 seeds:

| true C/N₀  | NWPR (default)            | p10  | p90  | M2M4 |
|:---------- | -------------------------:| ----:| ----:| ----:|
| noise only | **−Inf** in 12 of 16 runs |      |      | 27.9 |
| 20 dB-Hz   | 18.3                      | 6.0  | 22.2 | 27.6 |
| 25 dB-Hz   | 23.9                      | 16.3 | 26.3 | 28.3 |
| 30 dB-Hz   | 29.0                      | 27.4 | 30.3 | 30.8 |
| 35 dB-Hz   | 34.1                      | 33.3 | 35.0 | 35.3 |
| 45 dB-Hz   | 44.3                      | 43.7 | 45.2 | 45.4 |

On pure noise, pre-sync — where a code-lock detector has to make its call — the
default reports a median of `-Inf`, clears a 20 dB-Hz threshold in 4 % of updates
and a 25 dB-Hz threshold in **0 %** (200 runs, highest value seen 22.4). M2M4 on
the same runs clears 25 dB-Hz in 67 % of updates.

The trade-off is that `NBP` is a coherent sum, so it is sensitive to anything
that decorrelates the prompts inside a window — see
[What the estimator returns](#What-the-estimator-returns).

### The window follows the navigation bits

A narrowband window must not straddle a navigation-bit flip (a mid-window flip
costs ~7 dB). Where the bit boundaries are is something the tracking loop knows
and a downstream consumer of `get_filtered_prompts` does not, so each
record's prompt is handed to the estimator together with the navigation-bit
state in a [`CN0UpdateContext`](@ref), and the window is derived from it:

| signal state                                         | narrowband window                                                          |
|:---------------------------------------------------- |:-------------------------------------------------------------------------- |
| sync found, data-bearing signal (GPS L1 C/A, L5I, …) | `num_narrowband_code_blocks` (5), tiling the navigation bit from its start |
| sync found, pilot (GPS L1C-P, L2CL, Galileo E1C, …)  | `num_narrowband_code_blocks` (no bit grid to respect)                      |
| sync not found, data-bearing without secondary code  | `num_presync_narrowband_code_blocks` (5), unaligned                        |
| sync not found, signal with a secondary code         | none — the unknown overlay flips every code block                          |
| one symbol per code block (GPS L1C-D, Galileo E1B)   | none — no window longer than one record is coherent                        |
| record at least as long as its own window            | none — a one-record window has `NBP == WBP`                                |

Where no window is admissible the estimator falls back to a
[`MomentsCN0Estimator`](@ref), which needs no coherence at all — and with it to
that estimator's floor. Note the last row: with
[`set_preferred_num_code_blocks_to_integrate!`](@ref) at a whole navigation bit,
a window closes on a single record and NWPR reports its fallback for good.

!!! warning "The moment ratio's noise floor is there on data-only signals"

    GPS L1C-D and Galileo E1B carry **one navigation symbol per code block** (a
    10 ms and a 4 ms code period against a 100 sym/s and 250 sym/s data rate), so
    no two consecutive records are guaranteed to share a sign and the longest
    coherent window is a single record — where `NBP == WBP` identically. Those two
    signals therefore keep the moment ratio for good, and with it its floor:
    **pure noise still reads ~27.6 dB-Hz on L1C-D and E1B**, and a true 20 dB-Hz
    signal is indistinguishable from it. The same applies to any signal *before*
    its sync is found if it carries a secondary code (GPS L5, L2C, Galileo E1C,
    E5a), where the unknown overlay flips the sign every code block.

    There is no fix inside a ~100-record window: without coherence only the
    prompt *magnitudes* carry information, and at a per-record SNR of 0.1–0.3
    their distribution is barely distinguishable from the noise-only one at that
    sample count. (A phase-blind `(E|z|)²/E|z|²` estimator does better than M2M4
    on false alarms — 42 % against 76 % at a 25 dB-Hz threshold on pure noise, 6 %
    against 19 % at 30 dB-Hz — but its median on noise is still ~27.5 dB-Hz, so it
    moves the floor rather than removing it. Its floor shrinks as `N^{-1/2}`
    against the moment ratio's `N^{-1/4}`, so it would need seconds of averaging,
    not tens of seconds.)

    **This is what the default estimator exists to retire.**
    [`NoiseRefCN0Estimator`](@ref) has no window at all, so those signals and
    phases are ordinary records to it — that is
    [issue #217](https://github.com/JuliaGNSS/Tracking.jl/issues/217), and the
    reason the default moved. If you are on NWPR anyway (a correlator-ingest path
    with no noise observation), the mitigations are: take the lock decision from
    the **pilot** of the pair, which is the component the loops track anyway
    (L1C-D with L1C-P, E1B with E1C, both on one [`TrackedSat`](@ref)); give the
    data component a [`NoCN0Estimator`](@ref); or give it an
    `NWPRCN0Estimator(signal; fallback = NoCN0Estimator())`, which reports
    `-Inf dB-Hz` exactly when no coherent window is available.

The pre-sync window matters: the bit-edge detector needs seconds to lock at
35 dB-Hz and does not lock below ~30 dB-Hz, so a strictly bit-aligned window
would leave the low-C/N₀ regime — the one this estimator exists for — on the
moment ratio. An unaligned window of `M` records inside an `L`-block bit
straddles a flip with probability `(M−1)/L`, so the default five blocks of a
20-block GPS L1 C/A bit costs ~0.6 dB below 30 dB-Hz (and up to ~6 dB at a strong
signal, for the fraction of a second until sync is found) while removing the
moment ratio's noise floor completely.

### The window is capped by the loop's coherence time

A bit-aligned window does not have to span the whole bit, and it should not: the
`num_narrowband_code_blocks` windows tile the bit from its start, so none
straddles a flip *and* none runs past the coherence time of the loop that
produced the prompts. The longest flip-free window is not the best one — a longer
coherent window buys only a little (the same `num_records` records are combined
either way, so only their partition changes: at a true 25 dB-Hz, going from a
5-record to a 20-record window buys 0.8 dB of spread), while decoherence over a
long window costs a *bias*, the one error more averaging cannot remove.

The case that shows it is a satellite that fades **after** bit sync was found,
which is the normal way a receiver reaches low C/N₀. Locked in at 45 dB-Hz, then
faded to a true 25 dB-Hz, medians over 96 runs:

| coherent window       | reported C/N₀ | p10  | runs reading `-Inf` |
|:--------------------- | -------------:| ----:| -------------------:|
| 5 blocks (default)    | **23.6**      | 19.5 | 1 / 96              |
| one full 20-block bit | 21.0          | 12.0 | 17 / 96             |

A whole-bit coherent sum reports "no signal" on a satellite that is being
tracked. Raise `num_narrowband_code_blocks` for a pilot, a narrow carrier loop,
or a signal that is never weak; lower it if the loop is noisier than the default.

## Configuring the estimator

`num_prompts_for_cn0_estimation` (default `100`, ~100 ms for L1 C/A) sets how
many records the estimate averages over — for the default estimator, and for
NWPR's fallback as well:

```jldoctest cn0_buffer
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; num_prompts_for_cn0_estimation = 200);

julia> track_state = add_satellite!(track_state, :default, sat);

julia> get_prn(track_state, 1)
1
```

To choose a different estimator — or to configure NWPR's windows explicitly —
pass a `cn0_estimator` instance. It is a type parameter of
[`TrackedSignal`](@ref), so any [`AbstractCN0Estimator`](@ref) is stored as is:

```jldoctest cn0_estimator_kwarg
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; cn0_estimator = MomentsCN0Estimator(100));   # the default up to 5.1

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

Passing the **signal** to `NWPRCN0Estimator` is what sizes its coherent window
for that signal's code period — 5 blocks for a 1 ms code, 2 for GPS L1C-P's
10 ms one. The window is a count of blocks, so leaving it at the bare
constructor's `5` makes it 50 ms long on a 10 ms code.

Each signal needs its **own** estimator instance — they buffer into a vector, so
sharing one between two signals of a satellite corrupts both. The multi-signal
`TrackedSat` constructor therefore takes a tuple, one estimator per signal:

```julia
TrackedSat(
    (GPSL1C_P(), GPSL1CA()),
    prn,
    code_phase,
    carrier_doppler;
    cn0_estimator = (NoiseRefCN0Estimator(), MomentsCN0Estimator(100)),
)
```

A band is provisioned a noise source if **any** signal on it uses an estimator
that reads one, and the mixed case above works: the noise-referenced signal gets
its density, the other ignores it. A band where *no* signal asks for one is not
provisioned at all, so it runs no per-band despread and costs exactly nothing —
which is the answer for anyone who deliberately stays on NWPR.

### Not measuring a signal at all

[`NoCN0Estimator`](@ref) keeps no state, does no per-record work and reports
`-Inf dB-Hz`. It is the way to say *"this signal's C/N₀ is not measured"* — which
is worth saying in two situations:

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

The second only matters if you are on NWPR at all — the default estimator has no
fallback and no floor. Where you are, it turns the fallback's noise floor into an
honest `-Inf` for every signal and phase that admits no coherent window (see the
warning above). The per-record saving of the first is real but tiny — one
estimator update is well under a per cent of the correlation work for that
record — so choose it for clarity, not for speed.

It reports `-Inf dB-Hz` rather than `NaN dB-Hz` deliberately: with Unitful's
`Level` comparison, `NaN dB-Hz >= threshold` is `true` for *every* threshold, so a
`NaN` would clear every lock detector it met. `-Inf` compares `false` against any
finite threshold, which is the safe answer to "is this signal locked?".

```@docs
MomentsCN0Estimator
NoCN0Estimator
```

## What the estimator returns

  - **Before the first record has been folded in** — nothing has been measured and
    `estimate_cn0` returns `-Inf dB-Hz`. This is the value seen immediately after
    `add_satellite!`, and it is `-Inf` rather than a finite number precisely so a
    lock detector cannot clear on it.
  - **Driving a noise-free signal through `track`** — there is no noise to divide
    by, so the result runs off to a very large value; NWPR's coherent and
    incoherent powers agree exactly (`μ̂ = M`) and it reports `Inf dB-Hz`, as does
    the moment estimator. The [Quick start](index.md#Quick-start) shows this.
  - **Driving a noisy signal** — the estimate converges to the underlying CN0 as
    the buffers fill. Real-world signals typically land in 30–50 dB-Hz.
  - **Driving noise only** — the noise-referenced default reports whatever the
    averaged per-record terms come to, floored at `-Inf dB-Hz` once the mean fails
    to clear zero — which on pure noise is most of the time, and never a positive
    floor. NWPR reports `-Inf dB-Hz` whenever the mean ratio sits at the incoherent
    floor (`μ̂ ≤ 1`). Either way that is the *median* outcome and not a guarantee:
    code that **thresholds** the estimate needs no special case, code that
    **averages** it must handle the infinities.
  - **While the loops are struggling** — this is where the two differ most. NWPR
    measures *coherent* C/N₀, so residual phase noise counts against it: at 1 ms
    records and the default windows it reads ~0.9 dB low at a true 35 dB-Hz,
    ~1.1 dB low at 25 dB-Hz and ~1.7 dB low at 20 dB-Hz, and raising
    `num_narrowband_code_blocks` past the loop's coherence time makes it much
    worse. The non-coherent default is **immune** — measured through `track!` at a
    true 45 dB-Hz it reads 44.9 at 1-block records and 44.7 at 20-block ones, where
    NWPR falls to 32.9.

Every estimator here reports C/N₀ against the **effective** noise density,
including multi-access interference: the other satellites' cross-correlation adds
incoherent power `ε` that both NWPR and the noise reference see, so both report
`C/(N₀+ε)`. That is what other receivers report and the more useful number for a
lock detector. The one term they do *not* share is the reference's self-leakage —
see [The one bias it carries](#The-one-bias-it-carries).

!!! warning "The post-correlation filter is assumed to have unity noise gain"

    `N̂₀` is measured at the band, but the prompt the estimator sees is
    post-[`AbstractPostCorrFilter`](@ref). [`DefaultPostCorrFilter`](@ref) selects
    one antenna (`last`) and has gain 1, which is fine — and the noise reference
    despreads that same antenna, so the two sides of the ratio match. A real
    beamformer changes the prompt's noise scale and silently invalidates the
    ratio; there is no `noise_gain` hook yet.

### Residual Doppler at acquisition handoff

Acquisition hands over a Doppler estimate good to a search-bin width — commonly
±100 Hz or so — and the loop takes tens of milliseconds to pull it in. That
transient is where the coherent and non-coherent estimators differ most, and it
is also exactly when a lock detector is being asked whether to keep the
satellite.

Measured through `track!` at a true **40 dB-Hz**, 1 ms records, the satellite
seeded with a deliberate carrier-Doppler error, median of 5 seeds, read over the
first 60 ms while the loop is still converging:

| Δf     | real within-record loss | NWPR     | Moments | `NoiseRef` |
|:------ | -----------------------:| --------:| -------:| ----------:|
| 0 Hz   | 0.00 dB                 | 39.8     | 40.5    | **40.1**   |
| 30 Hz  | −0.01 dB                | 39.8     | 40.4    | **40.1**   |
| 60 Hz  | −0.05 dB                | 39.2     | 40.4    | **40.1**   |
| 120 Hz | −0.21 dB                | **27.6** | 40.5    | **39.9**   |
| 250 Hz | −0.91 dB                | **−Inf** | 40.1    | **39.1**   |

At 120 Hz NWPR under-reports by **12 dB**; at 250 Hz it reports `-Inf dB-Hz` —
"no signal" — on a satellite that is being tracked perfectly well.

The reason is the coherent sum. A residual `Δf` ramps the prompt's phase by
`2π·Δf·T` per record, so across the window the ramp is `Δf·M·T` **cycles**: 0.6
of a cycle at 120 Hz over the default 5 ms window, which spreads the five phasors
over 216° and costs `NBP` ~5.9 dB. `μ̂` collapses toward 1 and the inversion
`(μ̂−1)/(M−μ̂)` falls off a cliff; past 1.25 cycles `μ̂ ≤ 1` and the estimate is
`-Inf`. This is the constraint `M·T ≪ 1/(2·Δf)` stated under
[`NWPRCN0Estimator`](@ref), which at the default window means `Δf ≪ 100 Hz`.

The non-coherent estimators see only the **within-record** loss, `sinc²(Δf·T)` —
0.21 dB at 120 Hz, 0.91 dB at 250 Hz. That is real signal loss the correlator
genuinely suffers, so reporting it is correct rather than a defect; the measured
39.1 dB-Hz at 250 Hz is 40 minus exactly that.

What makes the default robust here is that the **noise reference is open-loop**:
it runs at the band's nominal IF with zero Doppler and measures *noise*, which is
white. No amount of Doppler error on any satellite can perturb the denominator,
so `N̂₀` is trustworthy precisely when nothing derived from the loop is.

Once the loop converges, NWPR recovers — 39.2 dB-Hz at every offset up to 120 Hz
after 300 ms. At 250 Hz it stays at `-Inf`, and the noise reference's 38.9 says
why: the loop never pulled in from there, so the residual — and its real 0.9 dB
loss — is still present. The two agree about the physics and disagree about
whether to report it.

If you are on NWPR anyway (a correlator-ingest path with no noise observation),
shorten the window that applies before bit sync — but it only goes so far,
because a two-record window is also the noisiest one:

| `num_presync_narrowband_code_blocks` | `Δf·M·T` at 120 Hz | reported |
|:------------------------------------ | ------------------:| --------:|
| 5 (default)                          | 0.60 cycles        | 27.6     |
| 3                                    | 0.36 cycles        | 32.4     |
| 2                                    | 0.24 cycles        | 34.9     |

```julia
NWPRCN0Estimator(GPSL1CA(); num_presync_narrowband_code_blocks = 2)
```

The doctest below builds a noisy L1 C/A signal at a known 45 dB-Hz CN0
and drives 25 1-ms tracking cycles through it:

```jldoctest cn0_example; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals, Random

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_code_frequency

julia> function run_cn0_demo()
           sys = GPSL1CA()
           fs = 4e6Hz
           num_samples = 4000  # 1 ms at 4 MHz
           prn = 1
           cn0_db_hz = 45.0
           sigma = sqrt((fs/Hz) / 10^(cn0_db_hz/10) / 2)
           code_freq = get_code_frequency(sys)
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

julia> run_cn0_demo()  # 25 ms is pre-sync, so the unaligned five-block window is in use
47.9 dB-Hz
```

The 1-σ noise amplitude is derived from the target CN0 by inverting
`C/N₀ = signal_power / (σ² / Fs)`. With unit-amplitude code samples this
collapses to `σ² = Fs / 10^(CN0_dB/10)`, split evenly across the complex
sample's real and imaginary parts. The seeded `MersenneTwister(0)` keeps
the doctest deterministic; 25 records only fill a quarter of the default
100-record window, so a few dB of spread is expected. The value is also a
pre-sync one: five of those records land in a window that straddles a data-bit
flip once the signal carries data, which this noise-free-code demo does not.

```@docs
estimate_cn0
```

## Custom CN0 Estimators

```@docs
AbstractCN0Estimator
CN0UpdateContext
```

You can implement your own estimator by creating a subtype of
[`AbstractCN0Estimator`](@ref) and implementing:

  - `Tracking.update(cn0_estimator::MyCN0Estimator, prompt)` — return a new
    estimator with the latest prompt added (immutable update).
  - `estimate_cn0(cn0_estimator::MyCN0Estimator, integration_time)` —
    return the CN0 estimate as a `dB-Hz` quantity. `integration_time` is the
    integration time of the *record* the last prompt came from, i.e.
    [`get_last_fully_integrated_integration_time`](@ref).

If the estimator can profit from the navigation-bit grid — bit boundaries,
whether sync has been found, the decoded soft bits — implement the
three-argument form instead, which is what the tracking loop calls:

  - `Tracking.update(cn0_estimator::MyCN0Estimator, prompt, context::CN0UpdateContext)`

The default three-argument method drops the context and calls the two-argument
one, so implementing only the latter is fine.

```@docs
Tracking.update(::MomentsCN0Estimator, ::Any)
Tracking.update(::Tracking.AbstractCN0Estimator, ::Any, ::CN0UpdateContext)
Tracking.update(::NWPRCN0Estimator, ::Any, ::CN0UpdateContext)
Tracking.update(::NWPRCN0Estimator, ::Any)
```

Plug it in with the `cn0_estimator` keyword of [`TrackedSignal`](@ref) or
[`TrackedSat`](@ref) — the estimator is a type parameter of `TrackedSignal`, so
your type is stored as is:

```jldoctest custom_cn0
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> struct MyCN0Estimator <: Tracking.AbstractCN0Estimator end

julia> Tracking.update(estimator::MyCN0Estimator, prompt) = estimator;

julia> Tracking.estimate_cn0(::MyCN0Estimator, integration_time) = 42.0 * Tracking.dBHz;

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz; cn0_estimator = MyCN0Estimator());

julia> estimate_cn0(sat)
42.0 dB-Hz
```
