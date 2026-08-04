# CN0 Estimator

The CN0 (Carrier-to-Noise density ratio) estimator provides a measure of
signal quality in dB-Hz. Each [`TrackedSignal`](@ref) on a
[`TrackedSat`](@ref) holds its own CN0 estimator, so a multi-signal
satellite produces one CN0 value per signal.

## Default Estimator

The default estimator is Van Dierendonck's **narrowband/wideband power ratio**
(NWPR), implemented as [`NWPRCN0Estimator`](@ref). Per narrowband window of `M`
records it forms the coherent power `NBP = |Σ prompt|²` and the incoherent power
`WBP = Σ |prompt|²`, averages `NBP / WBP` over as many windows as fit in
`num_prompts_for_cn0_estimation` records, and turns the mean ratio into a C/N₀.

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
34 dB-Hz, drop at 25) cannot be built on it. NWPR's false-alarm rate against any
threshold from 20 to 32 dB-Hz is 0 %, and its spread at low C/N₀ is much tighter
(p10–p90 of 15.8–22.2 at a true 20 dB-Hz). This mirrors the published
comparisons (Falletti, Pini & Lo Presti, *IEEE T-AES* 47(1):420–437, 2011). See
[issue #217](https://github.com/JuliaGNSS/Tracking.jl/issues/217) for the full
measurements.

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
| sync found, data-bearing signal (GPS L1 C/A, L5I, …)     | one navigation bit, started on a bit boundary          |
| sync found, pilot (GPS L1C-P, L2CL, Galileo E1C, …)      | `num_narrowband_code_blocks` (no bit grid to respect)  |
| sync not found, data-bearing without secondary code      | `num_presync_narrowband_code_blocks` (2), unaligned    |
| sync not found, signal with a secondary code             | none — the unknown overlay flips every code block      |
| one symbol per code block (GPS L1C-D, Galileo E1B)       | none — no window longer than one record is coherent    |

Where no window is admissible the estimator falls back to a
[`MomentsCN0Estimator`](@ref), which needs no coherence at all — so those
signals and phases behave exactly as they did before NWPR became the default.

The pre-sync window matters: the bit-edge detector needs seconds to lock at
35 dB-Hz and does not lock below ~30 dB-Hz, so a strictly bit-aligned window
would leave the low-C/N₀ regime — the one this estimator exists for — on the
moment ratio. An unaligned window of `M` records inside an `L`-block bit
straddles a flip with probability `(M−1)/L`, so the default two blocks of a
20-block GPS L1 C/A bit costs ~1 dB of bias while still removing the moment
ratio's noise floor completely.

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
  the ratio buffer fills. Real-world signals typically land in 30–50 dB-Hz.
- **Driving noise only** — NWPR reports `-Inf dB-Hz` once the mean ratio drops
  to the incoherent floor (`μ̂ ≤ 1`), i.e. "no detectable signal". Code that
  *thresholds* the estimate needs no special case; code that *averages* it must
  handle the infinities.
- **While the loops are struggling** — NWPR measures *coherent* C/N₀, so
  residual phase noise counts against it: with the conventional PLL at 1 ms
  records it reads ~1.5 dB low at a true 35 dB-Hz, ~5 dB low at 25 dB-Hz, and
  bottoms out at `-Inf` at 20 dB-Hz, where the loop no longer holds phase. That
  is the property that makes it a lock detector, and the reason a receiver
  reporting C/N₀ for weak, barely-tracked satellites may prefer
  `MomentsCN0Estimator`.

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

julia> run_cn0_demo()  # 25 ms is pre-sync, so the two-block window is in use
48.8 dB-Hz
```

The 1-σ noise amplitude is derived from the target CN0 by inverting
`C/N₀ = signal_power / (σ² / Fs)`. With unit-amplitude code samples this
collapses to `σ² = Fs / 10^(CN0_dB/10)`, split evenly across the complex
sample's real and imaginary parts. The seeded `MersenneTwister(0)` keeps
the doctest deterministic; 25 records only fill a quarter of the default
100-record window, so a few dB of spread is expected.

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
