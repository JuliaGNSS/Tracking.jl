# CN0 Estimator

The CN0 (Carrier-to-Noise density ratio) estimator provides a measure of
signal quality in dB-Hz. Each [`TrackedSignal`](@ref) on a
[`TrackedSat`](@ref) holds its own CN0 estimator, so a multi-signal
satellite produces one CN0 value per signal.

## Default Estimator

The default estimator is Van Dierendonck's **narrowband/wideband power ratio**
(NWPR), implemented as [`NWPRCN0Estimator`](@ref). Per narrowband window of `M`
records it forms the coherent power `NBP = |Σ prompt|²` and the incoherent power
`WBP = Σ |prompt|²`, divides the sums of both over as many windows as fit in
`num_prompts_for_cn0_estimation` records, and turns that mean power ratio into a
C/N₀.

```@docs
NWPRCN0Estimator
Tracking.default_cn0_estimator
```

### Why not the moment method

The [Moments Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1)
(M2M4, [`MomentsCN0Estimator`](@ref)) was the default up to and including
Tracking 5.1. It is a *moment ratio*, and at a finite window the sample moments
fluctuate enough to manufacture signal power out of noise. Measured on
synthetic prompts (amplitude `√(C/N₀·T)` in unit-variance complex noise), 100
one-millisecond records, 400 realizations:

| true C/N₀   | M2M4 median | p10  | p90  | NWPR (`M = 20`) median |
|:----------- | -----------:| ----:| ----:| ----------------------:|
| noise only  | **27.6**    | 22.6 | 31.1 | **−Inf** ("no signal") |
| 20 dB-Hz    | 27.7        | 22.5 | 30.9 | 19.9                   |
| 25 dB-Hz    | 27.9        | 23.4 | 30.9 | 25.0                   |
| 30 dB-Hz    | 30.5        | 27.3 | 32.6 | 30.0                   |
| 35 dB-Hz    | 35.3        | 33.9 | 36.4 | 35.0                   |
| 45 dB-Hz    | 45.1        | 44.3 | 45.9 | 45.0                   |

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

| true C/N₀   | NWPR (default) | p10  | p90  | M2M4 |
|:----------- | --------------:| ----:| ----:| ----:|
| noise only  | **−Inf** in 12 of 16 runs | | | 27.9 |
| 20 dB-Hz    | 18.3           | 6.0  | 22.2 | 27.6 |
| 25 dB-Hz    | 23.9           | 16.3 | 26.3 | 28.3 |
| 30 dB-Hz    | 29.0           | 27.4 | 30.3 | 30.8 |
| 35 dB-Hz    | 34.1           | 33.3 | 35.0 | 35.3 |
| 45 dB-Hz    | 44.3           | 43.7 | 45.2 | 45.4 |

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

| signal state                                            | narrowband window                                      |
|:------------------------------------------------------- |:------------------------------------------------------ |
| sync found, data-bearing signal (GPS L1 C/A, L5I, …)     | `num_narrowband_code_blocks` (5), tiling the navigation bit from its start |
| sync found, pilot (GPS L1C-P, L2CL, Galileo E1C, …)      | `num_narrowband_code_blocks` (no bit grid to respect)  |
| sync not found, data-bearing without secondary code      | `num_presync_narrowband_code_blocks` (5), unaligned    |
| sync not found, signal with a secondary code             | none — the unknown overlay flips every code block      |
| one symbol per code block (GPS L1C-D, Galileo E1B)       | none — no window longer than one record is coherent    |
| record at least as long as its own window                | none — a one-record window has `NBP == WBP`            |

Where no window is admissible the estimator falls back to a
[`MomentsCN0Estimator`](@ref), which needs no coherence at all — so those
signals and phases behave exactly as they did before NWPR became the default.
Note the last row: with
[`set_preferred_num_code_blocks_to_integrate!`](@ref) at a whole navigation bit,
a window closes on a single record and NWPR reports its fallback for good.

!!! warning "The moment ratio's noise floor is still there on data-only signals"
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

    In practice both of these are the data component of a pair whose **pilot**
    does get NWPR: L1C-D with L1C-P, E1B with E1C. Put both on one
    [`TrackedSat`](@ref) and take the lock decision from the pilot's estimate,
    which is the component the loops track anyway.

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

| coherent window | reported C/N₀ | p10  | runs reading `-Inf` |
|:--------------- | -------------:| ----:| -------------------:|
| 5 blocks (default) | **23.6**   | 19.5 | 1 / 96              |
| one full 20-block bit | 21.0    | 12.0 | 17 / 96             |

A whole-bit coherent sum reports "no signal" on a satellite that is being
tracked. Raise `num_narrowband_code_blocks` for a pilot, a narrow carrier loop,
or a signal that is never weak; lower it if the loop is noisier than the default.

## Configuring the estimator

`num_prompts_for_cn0_estimation` (default `100`, ~100 ms for L1 C/A) sets how
many records the estimate averages over, for the default estimator and its
fallback alike:

```jldoctest cn0_buffer
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        num_prompts_for_cn0_estimation = 200);

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

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        cn0_estimator = MomentsCN0Estimator(100));   # the default up to 5.1

julia> Tracking.get_cn0_estimator(sat) isa MomentsCN0Estimator
true

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        cn0_estimator = NWPRCN0Estimator(;
                            num_records = 400,                     # ~400 ms of memory
                            num_presync_narrowband_code_blocks = 0,  # M2M4 until bit sync
                        ));

julia> Tracking.get_cn0_estimator(sat).num_records
400
```

Each signal needs its **own** estimator instance — they buffer into a vector, so
sharing one between two signals of a satellite corrupts both. The multi-signal
`TrackedSat` constructor therefore takes a tuple, one estimator per signal:

```julia
TrackedSat((GPSL1C_P(), GPSL1CA()), prn, code_phase, carrier_doppler;
           cn0_estimator = (NWPRCN0Estimator(), MomentsCN0Estimator(100)))
```

```@docs
MomentsCN0Estimator
```

## What the estimator returns

- **Before the first record has been folded in** — nothing has been measured
  and `estimate_cn0` returns `0.0 dB-Hz`. This is the value seen immediately
  after `add_satellite!`.
- **Driving a noise-free signal through `track`** — the coherent and incoherent
  powers agree exactly (`μ̂ = M`), so the result is `Inf dB-Hz`. The moment
  estimator's noise moment goes to zero and it reports `Inf dB-Hz` too. The
  [Quick start](index.md#Quick-start) shows this behavior.
- **Driving a noisy signal** — the estimate converges to the underlying CN0 as
  the power buffers fill. Real-world signals typically land in 30–50 dB-Hz.
- **Driving noise only** — NWPR reports `-Inf dB-Hz` whenever the mean ratio sits
  at the incoherent floor (`μ̂ ≤ 1`), i.e. "no detectable signal". That is the
  median outcome, not a guarantee: a short window has an upper tail on noise, and
  at the default the estimate clears 20 dB-Hz in ~4 % of updates (0 % at
  25 dB-Hz), so a lock threshold wants some margin or a dwell. Code that
  *thresholds* the estimate needs no special case; code that *averages* it must
  handle the infinities.
- **While the loops are struggling** — NWPR measures *coherent* C/N₀, so residual
  phase noise counts against it. With the conventional PLL at 1 ms records and
  the default windows it reads ~0.9 dB low at a true 35 dB-Hz, ~1.1 dB low at
  25 dB-Hz, and ~1.7 dB low at 20 dB-Hz, with a wide spread and occasional `-Inf`
  there (see the table above). Raising `num_narrowband_code_blocks` past the
  loop's coherence time makes this much worse, not better.

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
           track_state = add_satellite!(track_state; prn, code_phase = 0.0, carrier_doppler = 0.0Hz)
           for _ in 1:25
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
