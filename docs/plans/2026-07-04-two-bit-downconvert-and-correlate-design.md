# Two-bit (sign+magnitude) hybrid-blocked DownconvertAndCorrelator

## Context

Two quantized downconvert→correlate backends already exist as opt-in alternatives to the Float32
default (`downconvert_and_correlate_fused.jl` / `CPUThreadedDownconvertAndCorrelator`):

  - **`Int16DownconvertAndCorrelator`** (`src/downconvert_and_correlate_int16.jl`, PR #160) — a full
    integer pipeline: `Complex{Int16}` measurement, Int8 carrier (SinCosLUT register-resident LUT),
    Int8 ±1 code (GNSSSignals embedded SIMD LUT), integer carrier-wipe + `vpmaddwd` accumulate.
    Essentially no quantization SNR loss; 1.3–2.9× faster than Float32.
  - **`OneBitDownconvertAndCorrelator`** (`src/downconvert_and_correlate_onebit.jl`, PR #162) —
    hard-limits **everything** to one sign bit, so the carrier wipe-off collapses to **XOR** and the
    tap accumulate to **popcount**. Fastest at wideband (5.6× Float32 at 20 MHz), at the cost of
    ≈3 dB of correlation SNR (≈2 dB input + ≈1 dB carrier).

This adds a **third point on that curve**: a two-bit **sign+magnitude** backend that sits between
one-bit and Int16 — clearly better SNR than one-bit, in the same speed class as Int16 (and far
faster than Float32), so a receiver that finds one-bit too lossy but Int16 too slow has a middle
option. It is a
direct extension of the one-bit backend (it reuses the bit-plane packing, funnel-shifted code
planes, band-shared measurement, and XOR+popcount correlate machinery), so it is built on
`feat/one-bit-on-160`.

**Goal:** add an opt-in `TwoBit(Threaded)DownconvertAndCorrelator` on `Complex{Int16}`
measurements. The measurement is always two-bit; the carrier precision is a **type parameter**
`CB ∈ {1, 2}` (1-bit sign-only carrier, or 2-bit sign+magnitude carrier); the code stays 1-bit
(±1), which is exact for binary modulations (BPSK/BOC/TMBOC). Never regress tracking metrics.

## The two-bit representation

Each measurement component is stored as **sign + magnitude**:

```
value = s·(1 + 2·b)          s = sign bit (as today);  b = magnitude bit (b=1 ⇔ |sample| ≥ threshold)
      ∈ {±1, ±3}
```

`{±1, ±3}` with the threshold near 1σ of the input is the classic near-optimal 4-level quantizer
(≈0.55 dB loss vs float, against ≈1.96 dB for 1-bit). Because the code and (when `CB=1`) the
carrier stay ±1, a product's **sign is still an XOR** of sign bits — only the **weight** changes,
from 1 to a product of `(1+2·b)` factors ∈ {1, 3, 9}. That weight expands into masked popcounts.

## The math (masked-popcount weight expansion)

For the `mᵣ·cos` term of tap *x* (measurement real `mr`, carrier `cos`), let the combined sign
plane be `S = codeₓ ⊻ smr ⊻ scos` (bit=1 ⇔ product is −1, exactly the one-bit "A" plane). With
measurement magnitude plane `Gmr` (bit=1 ⇔ `|mr| ≥ thr_m`) and, for `CB=2`, carrier magnitude
plane `Hcos` (bit=1 ⇔ `|cos| ≥ thr_c`):

```
Σ codeₓ·mᵣ·cos = Σ S·(1+2·Gmr)(1+2·Hcos)
   = (N − 2·pc S)                    ← one-bit result
   + 2·(pc Gmr − 2·pc(Gmr & S))      ← measurement magnitude correction
   + 2·(pc Hcos − 2·pc(Hcos & S))    ← carrier magnitude correction        (CB=2 only)
   + 4·(pc(Gmr&Hcos) − 2·pc(Gmr&Hcos&S))                                    (CB=2 only)
```

`pc Gmr`, `pc Hcos`, `pc(Gmr&Hcos)` are per-block constants (computed once). Every per-tap term is
a `count_ones` of `S` optionally ANDed with the per-sample magnitude planes — the same
`Vec{VW,UInt64}` VPOPCNTQ primitive the one-bit kernel already uses. The four sign planes of the
complex correlation are the one-bit `A/B/C/E` planes (`mr·cos`, `mi·sin`, `mi·cos`, `mr·sin`);
`Iₓ = (mr·cos)+(mi·sin)`, `Qₓ = (mi·cos)−(mr·sin)`.

**Validated** (standalone, no GNSS deps): the bit-plane kernel is bit-exact against the direct
elementwise quantized product sum across 2000 random cases (both `CB`), and in a weak-signal
phase-noise Monte-Carlo the two-bit (2-bit meas + 1-bit carrier) correlation loses **1.52 dB** vs
float where one-bit loses **3.08 dB** — a **1.56 dB** recovery, matching theory.

## Cost — popcounts per (antenna, tap)

| backend             | popcounts/tap | ANDs/tap | notes                              |
|:------------------- |:------------- |:-------- |:---------------------------------- |
| one-bit             | 4 (`A,B,C,E`) | 0        | baseline                           |
| **two-bit, `CB=1`** | **8**         | 4        | +1 masked popcount per sign plane  |
| **two-bit, `CB=2`** | **16**        | 12       | +3 masked popcounts per sign plane |

Only the correlate inner loop grows; code gen, carrier gen, and packing are shared. Because one
`count_ones(Vec{8,UInt64})` consumes 512 samples/instruction (vs Int16's `vpmaddwd` at 32
samples/instruction), the bit-plane loop keeps a throughput edge, but `CB=1`'s extra popcount +
the second (magnitude) measurement plane in packing put it **in the same performance class as
Int16** rather than decisively ahead.

**Measured** (this host, threaded, random 12-bit capture, GPS L1 C/A; microbenchmarks noisy on a
busy 24-core box, so read these as ranges): `CB=1` is ≈2.2–3.3× faster than Float32 and trades
blows with Int16 — ≈par-to-slightly-faster at 20 MHz (≈2.4 vs 2.5 µs / 1 sat, ≈7.6 vs 8.4 µs /
8 sats), a bit slower at 5 MHz. `CB=2` (2-bit carrier, Int8-LUT) recovers a further ≈0.6 dB and
runs ≈0.7–0.9× of Int16 — close to par, not the blowout the SNR-max option might suggest. The
wired-up CI benchmark (`track! Int16 vs Float32`, with `TwoBit`/`TwoBitC2` leaves) gives the
authoritative apples-to-apples numbers.

The clean, robust win is **SNR over one-bit** (≈1.4–1.8 dB, confirmed three ways: theory, a
standalone weak-signal Monte-Carlo (1.56 dB), and an end-to-end noisy tracking test (1.76 dB))
while staying ≈2–5× faster than Float32 and competitive with Int16.

## Decisions

  - **Separate, opt-in backend**, chosen upfront and hoisted (like Int16/one-bit); no auto-selection
    from sample eltype. Requires `Complex{Int16}` samples; errors clearly on Float measurements and
    on CBOC (non-binary, amplitude-carrying) codes — same gates as one-bit.
  - **Measurement always two-bit; carrier bits `CB ∈ {1,2}` a type parameter** on the backend
    (`TwoBitDownconvertAndCorrelator{CB}`), so the `@generated` kernel emits exactly the popcount
    terms needed and `CB=1` pays nothing for the carrier-magnitude path. Default `CB=1` (the
    speed/SNR sweet spot).
  - **Carrier planes.** `CB=1` reuses `SinCosLUT.generate_carrier_signs!` (the 1-bit NCO, sign
    plane only — its fast square-wave-run path). `CB=2` reuses SinCosLUT's **Int8 carrier LUT** —
    the exact generator the Int16 backend drives via `_int16_fill_carrier!` — built at **amplitude
    3** so a lane holds `round(3·sin) ∈ {0,±1,±2,±3}`; the sign plane is `signbit(lane)` and the
    magnitude plane is `|lane| ≥ 2` (since `|round(3·sin)| ≥ 2 ⇔ |sin| ≥ 2/3`, exactly the
    `{±1,±3}` crossover). Both planes come from one proven, register-resident (`vpermb`) source and
    are self-consistent — no bespoke phase math. The table's phase resolution is the `steps` knob.
    (An earlier draft hand-rolled the `CB=2` magnitude plane from a UInt32 phase-arc test; the Int8
    LUT is cleaner and measurably faster.)
  - **Measurement magnitude threshold** is a constructor knob `threshold` (ADC counts), analogous to
    Int16's `max_meas`. Default ≈1σ for a properly-AGC'd 12-bit front end; documented as
    per-front-end tunable (like Int16's amplitude choice). Under/over-declaring only shifts the
    quantizer operating point, never overflows (accumulators are `Int64`).
  - **Provide a threaded variant** (`TwoBitThreadedDownconvertAndCorrelator{CB}`), one scratch per
    thread, for a fair comparison against the threaded Float32 default and the other backends.
  - **Scope** mirrors one-bit: static tap counts (EPL NC=3, VEPL NC=5), any antenna count M,
    band-shared measurement packing for >1 sat, per-sat direct packing for a single sat. A
    multi-signal-per-sat correlates its signals in turn (no carrier-share tile — a future
    optimisation, as in one-bit).

## Implementation plan

New file `src/downconvert_and_correlate_twobit.jl`, closely paralleling
`downconvert_and_correlate_onebit.jl`:

 1. **Scratch/band buffers** — the one-bit planes plus magnitude planes: measurement `gmr`/`gmi`
    (and band `gmrband`/`gmiband`), and for `CB=2` carrier `hsin`/`hcos`. Grown lazily, reused.
 2. **Packing** — extend `_ob_pack_meas!` / `_ob_pack_band!` / `_ob_realign_meas!` to also emit the
    magnitude plane per component (a second `movemask` of `|re|≥thr` / `|im|≥thr`). Code packing
    (`_ob_pack_code!` + funnel `_ob_shift_plane!`) is unchanged (code is 1-bit).
 3. **Carrier NCO** — `_twobit_gen_carrier!(sinw, cosw[, hsinw, hcosw], len, cps; phase, ::Val{CB})`
    generating sign planes (both `CB`) and magnitude planes (`CB=2`) from the phase accumulator.
 4. **Kernel** — `_twobit_hybrid_blocked!`, `@generated` over `(NC, M, CB)`. Same strip-mine and
    plane layout as one-bit; the per-(antenna,tap) accumulate emits the masked-popcount terms above
    (2 or 4 popcounts per sign plane by `CB`). Finalize to `ComplexF64`/`SVector{M,ComplexF64}`.
 5. **Plumbing** — `_correlate_signals` (single + multi-signal), `_update_tracked_sat_correlator`,
    `_dc_one_group!`, `_dc_group_loop!` (threaded + single), `downconvert_and_correlate[!]` —
    mechanical copies of the one-bit versions with the `_TwoBitDC` union and threshold plumbed in.
 6. **Wire-up** — `include` + `export` in `src/Tracking.jl`.

Tests `test/downconvert_and_correlate_twobit.jl` mirror the one-bit suite (ratios vs Float32, VEPL,
multi-antenna, single-vs-threaded agreement, multi-signal, band-shared vs per-sat, error gates)
**plus** an SNR test asserting the two-bit correlation tracks Float32 better than one-bit, and a
`CB=1`-vs-`CB=2` agreement-of-shape test. Benchmark cases and `docs/src/track.md` updated to list
the third backend.
