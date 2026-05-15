# Multi-Signal Tracking Design (v2.0.0)

## Goal

Support tracking multiple signals per satellite — e.g. one satellite carrying GPS L1 C/A, L1C-D, and L1C-P — sharing a single carrier downconvert per satellite for performance, while keeping the track-state structure fully type-stable.

## Motivation

GNSSSignals v2 introduces signal types beyond the legacy one-signal-per-band view: GPS L1 now exposes `GPSL1CA`, `GPSL1C_D`, and `GPSL1C_P` as distinct `AbstractGNSSSignal` types. A modern satellite transmits all three on the same carrier; tracking them independently wastes work because the carrier Doppler is shared. A single downconvert per satellite, fanned out to a correlator per signal, is the natural primitive.

The current `SatState` holds exactly one correlator, one CN0 estimator, one bit buffer — i.e. one signal per satellite. The restructure lifts those per-signal fields into a new `TrackedSignal` and generalises the per-satellite state to hold a heterogeneous tuple of signals.

## Non-goals

- **Multi-band tracking.** `track` still consumes a single signal input. GPS L1 + GPS L5 is out of scope; the user constructs separate `TrackState`s per band and runs them independently.
- **Per-signal Doppler estimation.** All signals on a satellite share one carrier Doppler and one code Doppler, estimated from the first signal in the tuple.
- **Backward compatibility shims.** v2.0.0 is a clean break; the previous `TrackState(system, sat_states)` / `merge_sats` API is removed, not deprecated.

## Final types

```julia
struct TrackedSignal{Sig<:AbstractGNSSSignal, C<:AbstractCorrelator, PCF<:AbstractPostCorrFilter}
    signal::Sig
    integrated_samples::Int
    is_integration_completed::Bool
    correlator::C
    last_fully_integrated_correlator::C
    last_fully_integrated_filtered_prompt::ComplexF64
    cn0_estimator::MomentsCN0Estimator
    bit_buffer::BitBuffer
    post_corr_filter::PCF
    filtered_prompts::Vector{ComplexF64}
end

struct TrackedSat{Signals<:Tuple{Vararg{TrackedSignal}}, D}
    prn::Int
    code_phase::Float64              # mod max_code_length(signals) — incl. secondary
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    carrier_doppler::typeof(1.0Hz)
    signal_start_sample::Int
    signals::Signals                 # signals[1] = Doppler source
    doppler_estimator_state::D       # e.g. SatConventionalPLLAndDLL
end

struct TrackState{Sats<:NamedTuple, DE<:AbstractDopplerEstimator}
    satellites::Sats                 # (modern_gps = Dictionary{Int, TrackedSat{...}}, ...)
    doppler_estimator::DE
end
```

The signal-types tuple per group is recoverable from `Dictionary{Int, TrackedSat{Signals, ...}}`'s type parameter; there is no separate `signal_types` field on `TrackState`.

The current branch already introduced a `TrackedSat` wrapper that held `sat_state::SatState` + `estimator_state::E` as a refactor step. That wrapper is collapsed in this restructure: `SatState` goes away, and the per-sat estimator state moves into the new flat `TrackedSat` as `doppler_estimator_state`.

## Public API

### Construction

Multi-signal-set form:

```julia
track_state = TrackState(;
    signals = (
        legacy_gps = (GPSL1CA(),),
        modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
        galileo    = (GalileoE1B(),),
    ),
)
```

Each named entry is a *capability* — a group of satellites all tracked on the same signal-tuple shape. Different capabilities exist as different concrete `Dictionary` types side-by-side in the `satellites` NamedTuple, so the iteration over capabilities is type-stable.

Single-signal-set shortcut:

```julia
track_state = TrackState(; signal = GPSL1CA())
# desugars internally to (default = (GPSL1CA(),),)
```

### Adding satellites

```julia
# easy mode — library builds default TrackedSignals
add_satellite!(track_state;
    prn = 11,
    capability = :modern_gps,        # may be omitted iff one capability exists
    code_phase = 0.0,
    code_doppler = 0.0Hz,
    carrier_phase = 0.0,
    carrier_doppler = 1234.0Hz,
)

# immutable variant
track_state = add_satellite(track_state; prn=11, capability=:modern_gps, ...)

# escape hatch — pre-built TrackedSat for non-default chip spacing etc.
add_satellite!(track_state, :modern_gps, tracked_sat::TrackedSat)
```

`add_satellite!` cannot change `Sats` type parameters — correlator and post-corr-filter *types* are frozen at `TrackState` construction. Runtime knobs (initial Doppler, etc.) vary per call. Users who need non-default correlator types or post-corr-filter types build a `TrackedSat` themselves and hand it to the escape-hatch overload.

## Doppler source

The first signal in the tuple drives both the PLL and the DLL. The user controls this through tuple order at capability declaration:

```julia
modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA())
#             ^^^^^^^^^  Doppler source — pilot signal, recommended
```

**Putting a pilot signal first is encouraged** when one is available. Pilot signals carry no data-bit modulation, which lets the PLL run longer coherent integrations and reach lower phase-noise floors than a data-bearing signal can.

The data-bearing signals (L1C-D, L1 C/A) still recover their navigation bits independently — each `TrackedSignal` carries its own `bit_buffer` and processes its own prompt regardless of which signal drives the loop filter. The Doppler-source role and the data-recovery role are independent.

## Integration cadence

Each outer iteration of the tracking loop integrates to the **shortest** signal's next primary-code boundary:

```
samples_to_integrate = min(
    samples_to_next_boundary_for_each_signal...,
    samples_left_in_buffer
)
```

Walkthrough with `modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA())` over 20 ms of input:

- Each outer iteration ≈ 1 ms (L1 C/A primary period = 1023 chips / 1.023 Mcps, the shortest of the three).
- Every iteration: one shared downconvert over 1 ms; three `correlate!` calls.
- L1 C/A's correlator completes every iteration → `is_integration_completed = true`, its `last_fully_integrated_correlator` updates, working correlator resets to zero, its `bit_buffer` advances.
- L1C-D and L1C-P have 10 ms primary periods (10230 chips at 1.023 Mcps). Their correlators **accumulate across 10 iterations without resetting**; `is_integration_completed` stays false during those 10 iterations and flips to true on the 10th.
- After 20 ms of input: 20 L1 C/A integrations have completed, 2 L1C-D and 2 L1C-P integrations have completed.

### Doppler updates within long integrations

Carrier and code Doppler update at the rate of `signals[1]`'s completed integrations. When the Doppler source is a long-cadence signal like L1C-P, the PLL/DLL fires every 10 ms. Between firings, all signals' correlators accumulate against the *previous* iteration's Doppler.

When the Doppler source is a short-cadence signal (e.g. C/A first), the PLL/DLL fires every 1 ms, and long-cadence signals see ~10 piecewise Doppler updates across their integration window. This is the natural per-iteration-Doppler-correction behaviour of a real receiver: long integrations are not done against a single stale Doppler.

Choosing the Doppler source is therefore a tradeoff:
- **Pilot first** (e.g. L1C-P): cleanest phase reference, longest coherent integration, but Doppler updates only at the pilot's cadence (10 ms).
- **Short data-signal first** (e.g. C/A): Doppler updates at 1 ms, better dynamic response, but the PLL discriminator must contend with bit transitions.

Either is valid; the right choice depends on dynamics and C/N₀.

### Why MIN, not MAX

An alternative was to integrate to the **longest** signal's boundary per iteration. That would have aligned all signals' completion to the same iteration boundary, but at three costs:

1. The shared downconvert would only be valid for the full long window if no Doppler updates occurred mid-window — forcing Doppler updates to be deferred to the long cadence even when the Doppler source is a short signal.
2. The short signal's `last_fully_integrated_correlator` field would have to become a `Vector{C}` to hold the burst of completed correlators per iteration.
3. Downstream consumers would see bursty output instead of a steady per-primary-period stream.

MIN-boundary iteration keeps the data layout simple (scalar `last_fully_integrated_correlator`), gives the user a free choice of Doppler-source cadence, and shares the downconvert at the iteration grain.

## Code-phase wrapping

The shared `sat.code_phase` must wrap at the longest code period across all signals — including the secondary code. For L1C-P this is 10230 primary chips × 1800 secondary chips ≈ 18.4 M chips ≈ 18 s at 1.023 Mcps. Wrapping earlier would discard secondary-code phase tracking.

The wrap constant is recoverable from the signal-tuple type at compile time via tuple recursion:

```julia
@inline _max_code_length(::Tuple{}) = 0
@inline _max_code_length(t::Tuple) = max(
    get_code_length(first(t).signal) * get_secondary_code_length(first(t).signal),
    _max_code_length(Base.tail(t)),
)
@inline max_code_length(signals::Tuple{Vararg{TrackedSignal}}) = _max_code_length(signals)
```

A `@generated` function would also work, but tuple recursion is idiomatic Julia and produces the same constant-folded IR — verified to compile to `return 18414000` for the `modern_gps` tuple type, against `~20 ns` runtime cost for a plain `for` loop over the heterogeneous tuple. Plain `for` loops on heterogeneous tuples do not unroll at compile time and were rejected for this reason.

Per-signal replica-relative phase is derived from the shared phase:

```julia
get_code_phase(sat::TrackedSat, ::Type{Sig}) where {Sig} =
    mod(sat.code_phase, get_code_length(Sig()))
```

The secondary-code position within the long signal's overlay is `floor(sat.code_phase / get_code_length(Sig())) mod get_secondary_code_length(Sig())`.

## Downconvert dispatch — fused for N=1, tile-share for N≥2

The original design assumed that for multi-signal sats, a shared-downconvert path would beat calling the fused kernel N times. After two probes (the first one unfair, see git history of this file), the **empirical answer is yes** — but only if the shared-downconvert path reuses the existing SIMD downconvert infrastructure, not a naive Julia-loop replacement.

### Two existing fused kernels

`downconvert_and_correlate_fused!` already ships in two overloads:

- **Static-shifts** (`sample_shifts::SVector{NC}`): `@generated` unroll over antennas × taps; downconverted samples stay in SIMD registers and are immediately accumulated against each tap's code replica.
- **Dynamic-shifts** (`sample_shifts::AbstractVector`, requires caller-supplied SoA tile buffers): downconverts samples into a tile (one per antenna × num_samples), then runs a per-tap `@simd` correlate loop reading the tile back. Structurally already a "downconvert-once + per-tap correlate" implementation, just for one correlator.

### Measurements (5000 samples, 1 antenna, 3-tap EPL)

From `claude_scratch/dynamic_vs_static_microbench.jl`:

| Kernel | Time |
|---|---|
| Static-shifts (in-register accumulate) | 3.89 µs |
| Dynamic-shifts (downconvert-to-tile + correlate) | 4.82 µs |
| Pure correlate-from-tile (no downconvert) | 1.21 µs |

So for the dynamic-shifts path: ~3.61 µs of downconvert + ~1.21 µs of correlate. The tile materialization costs ~24% at N=1 vs keeping samples in registers — that's why single-signal sats should keep using the static-shifts kernel.

### A "tuple-of-correlators" kernel — what it costs

Extending the dynamic-shifts kernel to share its tile across N correlators (each with its own code replica and sample_shifts) adds one extra correlate-from-tile per added signal. Predicted costs vs fused-N-times:

| N | fused-N-times (static) | tile-share (downconvert + N × correlate) | savings |
|---|---|---|---|
| 1 | 3.89 µs | 4.82 µs | **−24% (loss)** |
| 2 | 7.77 µs | 6.03 µs | **+22%** |
| 3 | 11.66 µs | 7.25 µs | **+38%** |

Crossover at N > 1.35, so the tile-share kernel wins from N=2 upward. The savings grow with N because the downconvert (~3.6 µs) is amortised across more correlate work.

### Resulting design — compile-time dispatch on `length(sat.signals)`

```julia
if length(typeof(sat.signals).parameters) == 1
    # N=1: in-register fused kernel (keeps the static-shifts hot path).
    downconvert_and_correlate_fused!(sat.signals[1].correlator, ...)
else
    # N≥2: tile-share kernel — one downconvert into the shared tile,
    # then N correlate-from-tile passes against each signal's code
    # replica.
    downconvert_and_correlate_fused_tuple!(sat.signals, ...)
end
```

The branch resolves at type-inference time, so single-signal sats see no regression. Multi-signal sats pay one downconvert plus N×correlate-from-tile, beating fused-N-times by ≥20% from N=2 onward.

### Implementation arrives in two commits

Step 3 ships fused-N-times — the conservative path that always works. Step 4 then adds the tile-share kernel and the compile-time dispatch, keeping the change focused on a single optimization with its own benchmark to validate.

## Type-stability story

The hot loop walks two heterogeneous containers:

1. **`track_state.satellites`** — a `NamedTuple` of `Dictionary{Int, TrackedSat{...}}`s. Each named entry has its own concrete value type. Iteration via tuple recursion (`values(satellites)` returns a tuple; a recursive walk traverses it).
2. **`sat.signals`** — a `Tuple{TrackedSignal{...}, TrackedSignal{...}, ...}`. Each element has its own concrete type. Iteration via the same recursive pattern.

Both walks are type-stable and allocation-free. The existing test infrastructure for `@allocated == 0` in steady-state extends directly to the multi-signal case.

## Commit staging

The restructure lands on branch `ss/tracked-sat-wrapper` (PR #113) as a sequence of commits:

0. `docs: design plan for v2 multi-signal tracking` (this file).
1. `build(deps)!: bump to GNSSSignals v2 and update for renamed API` — pure rename pass: `AbstractGNSS` → `AbstractGNSSSignal`, `GPSL1` → `GPSL1CA`, `GPSL5` → `GPSL5I`, etc. No Tracking.jl semantics change.
2. `refactor!: introduce TrackedSignal wrapper around per-sat signal state` — moves per-signal fields out of `SatState` into `TrackedSignal`. Renames `SatState` → `TrackedSat` and folds the existing wrapper's `estimator_state` into `doppler_estimator_state`. Still single-signal-per-sat — `signals::Tuple{TrackedSignal{...}}` (one-tuple).
3. `refactor!: generalize TrackedSat.signals to Tuple{Vararg{TrackedSignal}}` — teaches the tracking loop to walk the tuple via recursive tuple-walks. Calls the existing fused kernel once per signal in the tuple. Multi-signal test fixtures added. Establishes the baseline against which step 4 will measure savings.
4. `perf: tile-share fused kernel for multi-signal sats` — adds a `downconvert_and_correlate_fused_tuple!` kernel that reuses the existing dynamic-shifts kernel's downconvert-to-tile path and a `@generated` sample-outer fused correlate that accumulates into all N×NC accumulators in one tile-streaming pass. Compile-time dispatch on `length(sat.signals)` keeps single-signal sats on the in-register static-shifts path. End-to-end measurement: **36% faster at N=2, 46% faster at N=3** vs step 3's fused-N-times baseline. Zero per-call allocations across N=1..4.
5. `feat!: TrackState carries satellites NamedTuple; add_satellite!/add_satellite API` — new constructor, new add-satellite API, band parameter dropped. Replaces `TrackState(system, sat_states)` / `merge_sats`.
6. `docs: rewrite for v2 multi-signal API` — README, docs/src/*, docstrings on new types/methods.

Each `!` commit carries a `BREAKING CHANGE:` footer in conventional-commit format; the CI lifts these into the generated v2.0.0 CHANGELOG entry.

## Open questions

None at design-time. Step 3's measurements confirmed that the original "shared-downconvert + per-signal correlate" design is correct, with a refinement: the shared part is the existing dynamic-shifts kernel's tile path, and the dispatch is by compile-time `length(sat.signals)`.

### Post-step-4 investigation: option C (no tile, all in registers)

A natural follow-up question: could a fully in-register multi-signal kernel (no tile materialization, downconverted samples never leave SIMD registers, accumulated into N×NC named-variable accumulators) beat option B's tile-share approach?

Two probes addressed this:

1. **Lower-bound probe** (`claude_scratch/option_c_inregister_probe.jl`): a cheating measurement that reuses a stale tile and skips the tile-write pass. Showed ≤5% headroom at N=2-3.

2. **Hand-coded SIMD prototype** (`claude_scratch/option_c_real_kernel_probe.jl`): an actual fused in-register kernel for N=2 with the same SIMD-intrinsic style as the production static-shifts fused kernel (fast `_to_vec` Int16→Float32 convert, `vload` SIMD loads, all 12 SIMD accumulators in named locals).

Against option B at N=2:

| Variant | Time |
|---|---|
| Option B (production) | 5.32 µs |
| Option C hand-coded SIMD prototype | 5.75 µs |
| Δ | **−8%** (option C is slower) |

**Option C is actually slower in practice.** The two-pass structure of option B has hidden advantages:

1. **Better pipelining.** Each of option B's passes has a shorter dependency chain, letting the CPU's out-of-order execution overlap multiple sample iterations in flight. Option C's longer per-iteration body (sincos → downconvert → N×NC code loads → 2N×NC muladds) limits ILP.
2. **Lower register pressure.** Option B's pass 1 only needs carrier-generation registers; pass 2 only needs the N×NC accumulators. Option C needs both live across the entire loop, pushing toward the AVX-2 16-register limit and risking spills.
3. **Better cache prefetching.** Option B's pass-2 tile reads are pure linear streams that the hardware prefetcher loves. Option C's per-iteration scattered code-replica loads interleaved with carrier generation don't prefetch as cleanly.

**Conclusion: option C is not worth implementing.** Option B is the right design — and not just by a small margin.
