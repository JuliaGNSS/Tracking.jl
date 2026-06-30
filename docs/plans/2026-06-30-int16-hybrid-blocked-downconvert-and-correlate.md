# Integer (Int16) hybrid-blocked DownconvertAndCorrelator

## Context

The current Tracking.jl downconvert→correlate path (`src/downconvert_and_correlate_fused.jl`)
runs a **Float32** pipeline: `Complex{Float32}` downconvert, carrier from FastSinCos
(`fast_sincos_u100k`, polynomial), Float32 accumulate with Float64 flush, and code from
GNSSSignals 2.2.2's fixed-point `gen_code!`.

Two upstream optimizations together replace this with a much faster **integer** pipeline, and an
integration benchmark exercising it lives in GNSSSignals:

- **GNSSSignals.jl PR #90** (stacked on #69, branch `feat/signal-embedded-lut`): per-signal
  embedded SIMD LUT code generator. `gen_code!` is now **Int8-only**, ~2× faster (~3.8× for
  L1C-P), plus a fused value engine `code_engine(…, Val(K))` and a value-threaded continuing fill
  engine `CodeFillEngine`. (`CodeReplicaLUT` and `CodeGeneratorLUT` were both removed.)
- **SinCosLUT.jl**: register-resident (`vpermb`) Int8 carrier LUT — sin/cos is a single hardware
  permute, no multiplies, no range limit.

The GNSSSignals `benchmark/benchmarks.jl` compares four correlation kernels (fused / unfused /
hybrid / **hybrid-blocked**) on a 12-bit-ADC `Complex{Int16}` measurement, Int8 carrier, Int8
code, integer `vpmaddwd` accumulate. **hybrid-blocked** (strip-mine into L1-resident blocks,
regenerate code+carrier per block) is the all-round winner — fused-class LLC footprint with
hybrid-class per-sample cost — and is also the most general (uses the array/continuing code fill,
which unlike the fused register engine handles non-baked secondaries like L5I NH10).

**Goal:** add a new `Int16`-input `DownconvertAndCorrelator` backend (separate from the Float32
one) using the new code+carrier generation with the hybrid-blocked kernel. It is **opt-in**
(selected upfront by the user, hoisted) and operates on `Complex{Int16}` measurements. Build it
up incrementally — EPL → parameterized correlator count → dynamic count → multiple antennas —
confirming the kernel choice by benchmark and never regressing tracking metrics.

## Decisions

- **Separate, opt-in backend.** `track!`'s default stays `CPUThreadedDownconvertAndCorrelator`
  (chosen upfront, hoisted); the user selects the integer backend manually. **No** auto-selection
  from sample eltype — the backend type can't wait until the first measurement is seen, and
  per-call selection defeats hoisting.
- Provide a **threaded** integer backend (`Int16ThreadedDownconvertAndCorrelator`) so it can be
  benchmarked fairly against the threaded Float32 default. A single-threaded
  `Int16DownconvertAndCorrelator` comes along for parity/debugging; the comparison of record is
  threaded-vs-threaded.
- Errors (clear message) if given a non-`Complex{Int16}` (e.g. Float) measurement. Constructor must
  be **hoistable** (build once, reuse — zero per-call alloc in steady state), like the existing CPU
  backends.
- **Accuracy bar: tracking metrics within tolerance** (not bit-identical). The Int8 carrier is
  deliberately low-precision; downstream consumers are all ratio/normalized so the carrier
  amplitude cancels.
- **Dependencies via `[sources]`** pinned to the PR branches; revert to registered versions once
  merged.

## Key facts

- Accumulators must stay `ComplexF64` / `SVector{M,ComplexF64}` — the correlator structs are
  `T<:Complex` and `sat_state.jl` / `cn0_estimation.jl` hardcode `ComplexF64`. The integer kernel
  computes Int sums (scaled by carrier amplitude) and **converts to `ComplexF64`** at finalize.
  No downstream type changes: discriminators (`atan(im/re)`, `(E−L)/(E+L)`), C/N0 (power ratios),
  bit buffer (sign) are all amplitude-invariant after the per-`integrated_samples` normalize.
- `Complex{Int16}` vectors/matrix-columns are stored as interleaved `[re,im,…]` Int16,
  bit-identical to the benchmark's `meas::Vector{Int16}` — the kernel can
  `vload(Vec{2W,Int16}, ptr, 2n+1)` directly. `BandMeasurement` already enforces dense layout.
- The backend is dispatched purely on the `downconvert_and_correlator` object today; there is **no**
  eltype hook, and we are **not** adding one. `track!`'s default kwarg stays
  `CPUThreadedDownconvertAndCorrelator()` (`src/track.jl:149`).
- Existing correlators are static `SVector{3}` (EPL) / `SVector{5}` (VEPL);
  `get_correlator_sample_shifts` returns an `SVector`. Dynamic-count infra exists
  (`get_initial_accumulator(num_ants, ::Integer)` returns a `Vector`) but no concrete type uses it.
- `code_engine(…, Val(K))` (the **fused register engine** for `code_lookup`/`code_advance`)
  **errors on non-baked secondaries** (L5I NH10, L1C-P TMBOC). The continuing `CodeFillEngine`
  (`code_engine(…)` without `Val(K)`) and the one-shot array `gen_code!` **do** apply them
  per-period. hybrid-blocked uses the `CodeFillEngine` fill path → works for all signals.

## External APIs

- **Carrier (SinCosLUT):** `tbl = SinCosTable(Int8; steps=64, amplitude=amp)`;
  `eng = carrier_engine(tbl; frequency, sampling_frequency)`;
  `st = carrier_state(eng, start; phase=<cycles::Real>)` (fractional initial phase in cycles);
  `carrier_lookup(eng, st) -> (sin,cos)::Tuple{Vec{W,Int8},Vec{W,Int8}}`;
  `carrier_advance(eng, st, nchunks)`; `carrier_width(eng)`.
- **Code (primary, continuing):** value-threaded fill engine —
  `eng = code_engine(signal, prn, fs, fc; start_phase, start_index_shift)` (no `Val(K)`) →
  immutable, isbits `CodeFillEngine` (one-time ~40 ns DDA setup); `st = code_state(eng)` → isbits
  `CodeFillState`; `gen_code!(out::Vector{Int8}, eng, st) -> st'` fills the next `length(out)`
  samples and returns the advanced state. **Allocation-free, thread-safe (caller threads `st`),
  exact across block seams, applies non-baked secondaries (L5I NH10) across boundaries.** Per-block
  code path for hybrid-blocked: build `eng`/`st` once per integration (per sat), thread `st` across
  blocks.
- **Code (one-shot option):** array `gen_code!(out::Vector{Int8}, signal, prn, fs, fc, start_phase,
  start_index_shift)` — same embedded LUT, Int8, 0-alloc, thread-safe; simpler when a single fill
  suffices.
- Constraint everywhere: `sampling_frequency ≥ code_frequency · subchip_factor`.

## Benchmark matrix (used in Phases 0a / 0b / 1+)

Sweep **signals × sampling rates**, threaded, single antenna, EPL, 1 ms, so the kernel choice and
speedups are representative across modulations and oversampling ratios (not just high-rate L1 C/A):

- **Signals:** GPS L1 C/A (BPSK, subchip 1), GPS L5I (BPSK + non-baked NH10 secondary, subchip 1),
  Galileo E1B (CBOC, subchip 12), GPS L1C-D (BOC), GPS L1C-P (TMBOC + 1800-chip overlay).
- **Rates per signal:** a **low-oversampling** point just above `fs_min = code_frequency ·
  subchip_factor`, plus a higher-OSR point. Concretely include **GPS L1 C/A at 2 MHz** (≈2× OSR),
  5 MHz and 40 MHz; for E1B / L1C use rates ≥ ~12.276 MHz (subchip 12); for L5I use ≥ 10.23 MHz.
  Skip any (signal, rate) with `fs < fc · subchip_factor` (error early in the kernel).

## Implementation phases

### Phase 0a — capture the CURRENT baseline first (no changes at all)
On the current registered deps (GNSSSignals 2.2.2), benchmark the **current implementation as-is**:
the Float32 fused path with its **Int16** code replicas and `CPUThreadedDownconvertAndCorrelator`,
across the benchmark matrix (all signals; rates incl. low-OSR like 2 MHz L1 C/A), threaded. Record
as the authoritative "before" baseline. GNSSSignals #90 removes the old Int16 `gen_code!`, so the
current implementation **cannot coexist** with the new deps in one environment — the baseline must
be taken now, and the current-vs-new head-to-head is **cross-environment** (current registered deps
vs PR-branch deps) on the same machine and scenario.

### Phase 0b — add deps + confirm the kernel choice (no Tracking source changes yet)
- Add `[sources]` to `Project.toml`: `GNSSSignals` → branch `feat/signal-embedded-lut`,
  `SinCosLUT` → `main`; add `SinCosLUT` to `[deps]`/`[compat]`; relax `GNSSSignals` compat to the
  branch's declared version. `Pkg.instantiate`.
- Reproduce the GNSSSignals integration benchmark on this host and **confirm hybrid-blocked is
  fastest** for the parallel multi-channel regimes; record numbers. If a different variant wins
  locally, surface it before committing (default stays hybrid-blocked per the PR's broad results).

### Phase 0.5 — migrate the existing Float32 path to Int8 (AFTER the baseline is recorded)
The new `gen_code!` is Int8-only and the legacy generator is gone, so the Float32 path must adopt it
to build on the new dep.
- Migrate `gen_code_replica!` (`src/code_replica.jl`) and the `ScratchBuffers.code_replica` /
  `_with_code_replica_buffer` element type (`src/downconvert_and_correlate_cpu.jl`) from
  `get_code_type(signal)` to **`Int8`**. The fused Float32 kernel reads code as generic
  `CT = eltype(code_replica)` and widens to Float32, so Int8 code works unchanged.
- Run `Pkg.test()`; **re-benchmark the migrated Float32 path** to isolate the Int8-code effect from
  the full integer-path effect (compare to the Phase-0a Int16 baseline).

### Phase 1 — minimal integer backend: EPL, single antenna, single signal/sat
- New file `src/downconvert_and_correlate_int16.jl` defining `Int16DownconvertAndCorrelator` and
  `Int16ThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator`, each holding
  hoistable per-(thread) scratch: block code buffer `extb` (`BLK + 2·max_shift` Int8), `csb`/`ccb`
  (BLK Int8 sin/cos), and a cached `SinCosTable`. Mirror the existing `ScratchBuffers` lazy-grow /
  `Threads.threadid()` pattern.
- Port the benchmark's integer helpers, generalized: `_wipe` (Int16/Int32 paths), `_madd*`/`_maddacc`
  (x86 `vpmaddwd`, gated to `Sys.ARCH`), and a `correlate_hybrid_blocked!` kernel generalized from
  E/P/L to **arbitrary `sample_shifts`** (NC offset reads from `extb`), producing `ComplexF64`
  accumulators. **Per-block code via the value-threaded `CodeFillEngine`**:
  `eng = code_engine(signal, prn, fs, fc; start_phase = signal_code_phase)`, `st = code_state(eng)`,
  then `st = gen_code!(view(extb, …), eng, st)` per block; carrier via the value engine filled into
  `csb`/`ccb` (4-way unrolled, as in `_epl_fill_carrier!`). Add **per-block Int64 flush** of the
  Int32 lane accumulators sized so worst-case `code_amp · |DI| · samples/lane` cannot overflow
  Int32 (mirrors the #152 blocked-accumulation reasoning, in integers).
- Wire the per-sat flow with methods dispatched on the new types, paralleling the CPU backend:
  `downconvert_and_correlate(!)`, `_dc_one_group!` reuse, `_dc_group_loop!` (serial + `@batch`),
  `_update_tracked_sat_correlator`, integer `_correlate_signals(::Tuple{TrackedSignal}, …)`. Reuse
  `_signal_replica_params`, boundary/`is_integration_completed` logic, and `update_accumulator`.
- Choose carrier amplitude via the benchmark's `choose_carrier(max_meas)`; expose it as a
  constructor argument with a 12-bit default. Leave integer sums unscaled when converting to
  `ComplexF64` (amplitude-invariant downstream); document this.
- **No auto-selection.** The user opts into `Int16ThreadedDownconvertAndCorrelator()` explicitly. An
  `Int16*` backend handed a non-`Complex{Int16}` measurement errors with a clear message at the
  group dispatch boundary.
- Validate: existing test signals track to convergence; compare discriminator outputs / C/N0 /
  final Doppler vs the Float32 path within tolerance. Benchmark the threaded integer backend vs the
  threaded Float32 default — against both the Phase-0a Int16 baseline and the Phase-0.5 Int8 path.

### Phase 2 — parameterized (static) number of correlators
- Generalize the kernel's tap unrolling over `NC = length(sample_shifts)` (`@generated` on NC).
  Cover `VeryEarlyPromptLateCorrelator` (NC=5) and any `NumAccumulators{NC}`-built EPL. Re-validate +
  benchmark.

### Phase 3 — dynamic number of correlators
- Support correlators whose `accumulators` is a `Vector` (the `get_initial_accumulator(_, ::Integer)`
  path). Add a runtime-NC kernel overload (loop over taps), analogous to the Float32 path's
  `AbstractVector`-shifts fallback. Re-validate + benchmark; document the static-vs-dynamic cost.

### Phase 4 — multiple antenna channels
- Matrix `samples` (columns = antennas). Per-antenna carrier wipe reusing the shared code block;
  accumulators become `SVector{M,ComplexF64}`. Generalize the kernel's accumulator set to `M·NC`
  (antenna-outer passes to bound live registers, mirroring the Float32 tile-share kernel). Re-validate
  (M=1 unchanged) + benchmark. Add the tile-share multi-signal analog last if needed.

## Files

- `Project.toml` — `[sources]`, `[deps]`/`[compat]` for `SinCosLUT`, relaxed `GNSSSignals`.
- `src/code_replica.jl`, `src/downconvert_and_correlate_cpu.jl` — Int8 code buffers (Phase 0.5).
- `src/downconvert_and_correlate_int16.jl` — **new**: backend structs, integer kernels, per-sat flow.
- `src/Tracking.jl` — `using SinCosLUT`, `include` the new file, export the new backend type(s).
  (`src/track.jl` unchanged — default backend stays `CPUThreadedDownconvertAndCorrelator`.)
- `test/` — new tests for the Int16 backend (correctness vs Float32 within tolerance; opt-in usage
  via the `downconvert_and_correlator` kwarg; clear error on non-`Complex{Int16}` input).

## Risks / watch-items

- **Int32 accumulator overflow** for multi-level CBOC code and long (≥20 ms) integrations — handled
  by per-block Int64 flush with an overflow-safe block size; verify against worst case.
- **`CodeFillEngine` build cost per integration** (~40 ns DDA setup; engine/state are isbits, so the
  steady-state fill is allocation-free). Confirm construction stays off the per-sample hot path and
  doesn't allocate.
- **GNSSSignals bump is breaking** — Phase 0.5 must land with green tests before Phase 1.
- **Non-x86 hosts** — `vpmaddwd` helpers are x86-gated; ensure a correct (if slower) NEON/portable
  path so tests pass everywhere.
- **`subchip_factor` minimum-rate constraint** — error early if a sat's `fs < fc·subchip_factor`.

## Verification

- `Pkg.test()` green at each phase (Float32 path unchanged after Phase 0.5; new Int16 tests added).
- Tracking-metric comparison: run both backends on the same synthetic `Complex{Int16}` capture;
  assert PLL/DLL discriminator outputs, C/N0, and converged carrier/code Doppler agree with the
  Float32 path within tolerance, and that lock is maintained over the full capture.
- Benchmarks at each phase vs the Phase-0a baseline across the benchmark matrix (all signals;
  low-OSR + high-OSR rates incl. 2 MHz L1 C/A; single-channel and Polyester multi-channel) showing
  speedup and no regression; record in the benchmark suite.
