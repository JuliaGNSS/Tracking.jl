# Bit-Sync / Secondary-Code Detection Redesign

Part of the Tracking.jl 1.x → 2.0 bump. No backward-compat shims;
`is_upcoming_integration_new_bit`'s signature, `BitBuffer`'s field layout,
and `TrackedSignal`'s type parameters all change.

**Upstream dependency.** This redesign requires GNSSSignals.jl ≥ v2.2.1,
which fixes `get_data_frequency(::GPSL1C_D)` to return the broadcast
100 sps symbol rate instead of the post-LDPC 50 bps info-bit rate
([GNSSSignals.jl PR #61](https://github.com/JuliaGNSS/GNSSSignals.jl/pull/61),
merged). Without it, `_calc_num_code_blocks_that_form_a_bit` computes 2
instead of 1 for L1C-D and the L1C-D row of the design below is invalid.

*Current state of this branch:* v2.2.1 has not yet propagated through
the General registry ([JuliaRegistries/General#155956](https://github.com/JuliaRegistries/General/pull/155956)),
so Tracking.jl is `Pkg.develop`'d against the local
`~/Code/GNSSSignals.jl` checkout (`Manifest.toml` reflects this; the
`[compat]` entry stays at `"2"` until the registry catches up, at which
point it gets pinned to `"2.2.1"` and the dev-dep is replaced).

## Goal

Replace the current "exact match against an unshifted template" sync detector
with one that **(a)** works for arbitrary template lengths up to ~1800 bits,
**(b)** runs an explicit phase search rather than waiting for self-alignment,
and **(c)** tolerates bit errors. This unlocks secondary-code lock on GPS
L1C-P and tightens the existing detectors on GPS L1 C/A and GPS L5I.
GPS L1C-D and Galileo E1B drop to one-block-per-symbol trivial cases (no
edge search needed — each broadcasts one channel symbol per primary code
period).

## Motivation

Today every signal funnels through `is_upcoming_integration_new_bit(signal,
code_block_bits, num_code_blocks)`. Each implementation compares the buffered
primary-code-block signs against a single fixed bit pattern at the buffer's
*current* shift:

- GPS L1 C/A — match against `0xfffff` / `0xfffff00000` over 40 blocks.
- Galileo E1B — match against `0xf` / `0xf0` over 8 blocks.
- GPS L5I — XOR against `0x35` (NH10) over 10 blocks, accept all-zero or
  all-one result.
- GPS L1C-D — `return false` (deferred).
- GPS L1C-P — `return false` (deferred); a `UInt128` buffer can't hold the
  1800-chip overlay.

Three problems:

1. **Worst-case lock time is double the code length.** The detector only
   fires when the receive buffer happens to align with the template. For
   L5I that's 10 ms × 20 ≈ 200 ms; for an 1800-chip L1C-P overlay it would
   be 18 s × 2 = 36 s.
2. **No error tolerance.** A single bit-flip from noise blocks lock even at
   the right phase.
3. **`code_block_buffer::UInt128` is too narrow for L1C-P.** Its 1800-chip
   overlay does not fit.

## Non-goals

- **Preamble search.** L1C-D's CNAV-2 preamble (and the corresponding
  decoder-layer phase resolution) lives in GNSSDecoder.jl. This redesign
  stops at bit-edge / overlay-phase detection.
- **Coherent integration past the secondary-code period.** Once L1C-P locks
  the overlay phase, the existing `max_code_length` / `code_phase` wrap
  machinery handles indexing. Adopting >18 s coherent integration is a
  follow-up.

## Background — signal parameters

| Signal     | Primary period | Symbol rate | Primary blocks per symbol | Secondary code | Sync detector input length |
|------------|----------------|-------------|---------------------------|----------------|----------------------------|
| GPS L1 C/A | 1 ms           | 50 Hz       | 20                        | none           | 40 bits                    |
| Galileo E1B| 4 ms           | 250 Hz      | 1                         | none           | n/a (no edge to find)      |
| GPS L5I    | 1 ms           | 100 Hz      | 10 (= secondary code)     | NH10 (shared)  | 20 bits                    |
| GPS L1C-D  | 10 ms          | 100 Hz      | 1                         | none           | n/a (no edge to find)      |
| GPS L1C-P  | 10 ms          | 0 Hz        | n/a (pilot)               | 1800 (per-PRN) | 1800 bits                  |

The "sync detector input length" is `2 × num_blocks_per_symbol` (so the
buffer spans at least one full symbol-edge transition) for data-bearing
signals with multiple primary periods per symbol, and the full secondary-
code length for pilot signals. L1C-D and Galileo E1B both broadcast one
channel symbol per primary code period, so they have no sub-symbol
boundary for the tracker to find — `is_upcoming_integration_new_bit`
returns `true` immediately and the tracker hands one prompt-correlator
sum per primary period to the downstream decoder.

> **Note** — `get_data_frequency` in GNSSSignals.jl reports the channel
> symbol rate (post-FEC, on-air), not the post-decode information bit
> rate. For L1C-D that's 100 sps, with 50 bps information after rate-½
> LDPC decoding handled by GNSSDecoder.jl. See GNSSSignals.jl PR #61.

## Benchmark summary (length 1800)

Measured on the dev machine, `--startup-file=no`:

| Approach                              | Single Hamming | Single shifted Hamming | Full 1800-phase search |
|---------------------------------------|----------------|------------------------|------------------------|
| `BitIntegers.@define_integers 1800`   | 16.6 ns        | 49.9 ns                | **70.6 μs**            |
| `BitIntegers.@define_integers 1856` (padded) | 24.1 ns | 70.5 ns                | 106.2 μs               |
| `Vector{UInt64}` (29 limbs)           | 7.7 ns         | 59.3 ns                | 104.1 μs               |
| FFT cross-correlation (length ~3600)  | n/a            | n/a                    | 113.2 μs (+170 KiB)    |

**Exact-width `UInt1800` wins.** The padded variant pays a mask-everywhere
tax; the limbs variant matches it on a single op but loses on the shift because
multi-limb rotation in scalar Julia is awkward; FFT is no faster and
allocates. Bench scripts live in
[`claude_scratch/bench_bitintegers_1800.jl`](../../claude_scratch/bench_bitintegers_1800.jl)
and [`claude_scratch/bench_bitintegers_1800_exact.jl`](../../claude_scratch/bench_bitintegers_1800_exact.jl).

## Design

### Sync detector contract

Replace the `Bool`-returning detector with a typed result:

```julia
struct SyncResult
    found::Bool
    phase::Int       # 0 for no-secondary signals; secondary-chip offset otherwise
    polarity::Int8   # +1 or -1; which match orientation locked
end

is_upcoming_integration_new_bit(
    signal::AbstractGNSSSignal,
    code_block_bits::B,        # B<:Unsigned, per-signal width
    num_code_blocks::Int,
)::SyncResult
```

Existing one-shot match detectors collapse to this with `phase = 0` and the
matched polarity passed through.

### Per-signal buffer widths

The current `code_block_buffer::UInt128` becomes a parameter of `BitBuffer`:

```julia
struct BitBuffer{B<:Unsigned}
    code_block_buffer::B
    code_block_buffer_length::Int
    found::Bool
    secondary_phase::Int       # 0 until found; secondary-chip offset after
    buffer::UInt128            # decoded navigation bits (post-sync)
    length::Int
    prompt_accumulator::ComplexF64
    prompt_accumulator_integrated_code_blocks::Int
end
```

Per-signal width is set by a new trait:

```julia
get_code_block_buffer_type(::AbstractGNSSSignal) = UInt64
get_code_block_buffer_type(::GPSL1CA)            = UInt64   # need 40 bits
get_code_block_buffer_type(::GalileoE1B)         = UInt8    # 1 block per symbol — buffer unused; type still required
get_code_block_buffer_type(::GPSL5I)             = UInt32   # need 20 bits
get_code_block_buffer_type(::GPSL1C_D)           = UInt8    # 1 block per symbol — buffer unused; type still required
get_code_block_buffer_type(::GPSL1C_P)           = UInt1800 # need 1800 bits
```

`UInt1800` is introduced via
`BitIntegers.@define_integers 1800` once at module load. The post-sync output
`buffer::UInt128` stays a `UInt128` — it holds decoded navigation bits, not
sync-search state.

`BitBuffer{B}` is parametric, so the concrete type of `code_block_buffer`
flows into `TrackedSignal{Sig, C, PCF}` through the existing
`bit_buffer::BitBuffer` field. The `B` parameter is determined at
`TrackedSignal` construction from `get_code_block_buffer_type(signal)`, so
the parameter chain stays type-stable. No `Union{...}` needed.

### Detector internals

A generic Hamming-window matcher handles every signal:

```julia
@inline function _try_match(
    masked_buffer::B,
    template::B,
    mask::B,
    max_errors::Int,
) where {B<:Unsigned}
    masked = masked_buffer & mask
    dist_pos = count_ones(masked ⊻ template)
    dist_pos ≤ max_errors && return (found=true,  polarity=Int8(+1))
    dist_neg = count_ones(masked ⊻ (template ⊻ mask))
    dist_neg ≤ max_errors && return (found=true,  polarity=Int8(-1))
    return (found=false, polarity=Int8(0))
end
```

Per-signal `is_upcoming_integration_new_bit`:

- **GPS L1 C/A** — `_try_match(code_block_bits, 0xfffff, 0xffffffffff, max_errors)`,
  where `max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * 40)`.
  Wait until `num_code_blocks ≥ 40`. Note: the *signed* `±1` representation of a
  bit-edge is "20 same followed by 20 opposite", so the canonical template is
  `0xfffff` (20 ones followed by 20 zeros, i.e. `0x0000_000f_ffff`) and the
  negated-polarity branch picks up `0xfffff00000`. Default tolerance 2.5 %
  discretizes to 1 bit-flip allowed over 40 blocks.
- **Galileo E1B** — 1 channel symbol per primary code period (250 sym/s,
  4 ms primary period; Galileo OS SIS ICD Tables 11 & 15). No sub-symbol
  boundary to find, so:
  ```julia
  is_upcoming_integration_new_bit(::GalileoE1B, _, _, _) = SyncResult(true, 0, Int8(+1))
  ```
  Polarity ambiguity is resolved by GNSSDecoder.jl via the I/NAV preamble.
- **GPS L5I** — template `0x000ff` for NH10 (10-bit secondary code `0000110101`
  expanded into a 20-bit window as `10-bits ++ negated-10-bits`),
  `max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * 10)`.
  Default tolerance 2.5 % discretizes to 0 errors (exact match) — the 10-block
  window doesn't carry enough chips to tolerate even one flip at 2.5 %.
  Equivalent to the pre-redesign XOR-with-`0x35` check, generalised across
  all 10 possible alignments of the buffer rather than only the aligned case.
- **GPS L1C-D** — 1 channel symbol per primary code period (100 sps,
  10 ms primary period). No sub-symbol boundary to find, so:
  ```julia
  is_upcoming_integration_new_bit(::GPSL1C_D, _, _, _) = SyncResult(true, 0, Int8(+1))
  ```
  Identical contract to Galileo E1B (also 1 block per symbol). Polarity
  ambiguity (whose "+1" is the data symbol's "0" bit?) is GNSSDecoder.jl's
  problem — it runs its CNAV-2 preamble search at both polarities and locks
  the one whose CRC checks out.
- **GPS L1C-P** — full 1800-phase shifted search against the per-PRN overlay
  matrix from `get_secondary_code(signal)`. Two-stage:
  1. **Fill** (1800 primary periods, 18 s): just shift bits into the buffer.
     `_try_match` is *not* called.
  2. **Search** (one shot, ~70 μs): once `num_code_blocks == 1800`, run
     `_full_phase_search` with
     `max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * 1800)`,
     store the best phase, set `found = true`. Default tolerance 2.5 %
     discretizes to 45 bit-flips allowed; user-overridable via dispatch on
     the trait.

The full-phase search reuses the rotate-XOR-popcount inner loop from the
bench:

```julia
@inline function _full_phase_search(received::B, overlay::B, max_errors::Int) where {B<:Unsigned}
    best_k, best_dist, best_pol = -1, typemax(Int), Int8(0)
    for k in 0:1799
        shifted = (received << k) | (received >> (1800 - k))
        d_pos = count_ones(shifted ⊻ overlay)
        d_neg = 1800 - d_pos
        if d_pos < best_dist; best_dist, best_k, best_pol = d_pos, k, Int8(+1); end
        if d_neg < best_dist; best_dist, best_k, best_pol = d_neg, k, Int8(-1); end
    end
    best_dist ≤ max_errors ?
        SyncResult(true, best_k, best_pol) :
        SyncResult(false, 0, Int8(0))
end
```

(`UInt1800`'s exact width means no mask is needed — the rotation is the
natural 1800-bit rotate.)

### Phase plumbing

When `SyncResult.found == true && phase != 0`, the recovered secondary-chip
offset must seed the shared `sat.code_phase` so the wrap-at-`max_code_length`
machinery in [`sat_state.jl:138`](../../src/sat_state.jl#L138) indexes the
overlay starting from the correct chip.

For pilot signals the offset is meaningful immediately (it tells you "we
just integrated primary block N of the secondary code, where N = phase").
For data signals it's always 0 — no plumbing change needed beyond the
record-keeping.

The seeding happens in `buffer()` ([`bit_buffer.jl:81-132`](../../src/bit_buffer.jl#L81)),
where today the post-sync transition stores `found = true` and a one-shot
decoded-bit conversion. We add:

- Carry `SyncResult.phase` and `polarity` through into the post-sync
  `BitBuffer`.
- In the same `track` step that flips `found`, return the recovered phase
  alongside the new `BitBuffer` and let `_process_doppler_source_signal`
  (or the per-passenger equivalent) snap `sat.code_phase` to the correct
  offset.

This is the only piece that touches outside `bit_buffer.jl` and the
per-signal detectors. It's small: one extra field on `BitBuffer`, one
branch in the post-sync handover.

## Error tolerance defaults

A single package-wide trait
`get_bit_edge_or_secondary_code_tolerance(::AbstractGNSSSignal) = 0.025`
expresses the tolerance as a fraction of the search window. Each detector
converts to an integer error budget at its call site via
`floor(Int, tolerance × window_size)`. The 2.5 % choice covers every signal
without per-signal overrides:

| Signal    | Window | Effective `max_errors` at 2.5 % | Justification |
|-----------|--------|----------------------------------|---------------|
| L1 C/A    | 40     | 1                                | Uniform 2.5 % ceiling across detectors |
| E1B       | n/a    | n/a                              | Trivial — `SyncResult(true, ...)` returned unconditionally (1 block per symbol) |
| L5I       | 10     | 0 (exact match)                  | 10-block window can't discretize 2.5 % into a non-zero integer |
| L1C-D     | n/a    | n/a                              | Trivial — `SyncResult(true, ...)` returned unconditionally |
| L1C-P     | 1800   | 45                               | 18 s coherent gain — still well within the original 2 % envelope |

The trait is `@inline`'d and folds at each detector's call site, so user
overrides via dispatch (e.g.
`Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL1CA) = 0.05`) take
effect at the next call to `is_upcoming_integration_new_bit` without
rebuilding any `TrackState`. Out of scope for this redesign: tying
tolerance to the live CN0 estimator.

## Implementation order

1. **Parameterize `BitBuffer` on `B<:Unsigned`.** Mechanical: thread the
   parameter through `TrackedSignal{Sig, B, C, PCF}`. Existing detectors
   keep their `Bool` return; the `SyncResult` refactor is step 3.
   No behavior change. Touched files: `src/bit_buffer.jl`,
   `src/sat_state.jl`, `src/tracking_state.jl`,
   `src/conventional_pll_and_dll.jl`, `src/downconvert_and_correlate_cpu.jl`
   (read-only access to the bit-buffer field).
   - **Tests:** update `test/bit_buffer.jl` — the existing `BitBuffer()`
     `@inferred` constructor test stays but the default width changes;
     add a `BitBuffer{UInt8}()` case to confirm the parameter flows. Also
     run `test/sat_state.jl`, `test/tracking_state.jl`, `test/multi_signal.jl`
     to confirm `TrackedSignal`'s new type parameter doesn't break
     inference.
2. **Add `get_code_block_buffer_type` trait** with widths from the table
   above. Wire it into `TrackedSignal`'s constructor so each signal picks
   its width.
   - **Tests:** new `@testset` per signal in `test/gps_l1ca.jl`,
     `test/galileo_e1b.jl`, `test/gps_l5i.jl`, `test/gps_l1c_d.jl`,
     `test/gps_l1c_p.jl` calling `get_code_block_buffer_type(signal)` and
     asserting the exact width. (`UInt1800` arrives in step 4 — until
     then L1C-P's trait can return `UInt64` as a placeholder, with the
     test marked `@test_skip`.)
3. **Introduce `SyncResult` and refactor the three existing detectors**
   (L1 C/A, E1B, L5I) to return it. Pure refactor, behavior preserved
   bit-for-bit. Add `_try_match` as a shared helper but keep `max_errors = 0`
   so the existing exact-match tests pass unchanged.
   - **Tests:** update the existing per-signal
     `is_upcoming_integration_new_bit` tests to assert against
     `SyncResult` fields (`found`, `phase`, `polarity`) instead of `Bool`.
     Add `@inferred` calls. Add `_try_match` unit tests in
     `test/bit_buffer.jl` covering: exact match, single-bit error rejected
     at tolerance 0, polarity flip, mask-out-of-band bits.
4. **Add BitIntegers dependency, define `UInt1800`,** and implement L1C-P
   `is_upcoming_integration_new_bit` with the full-phase search. Wire the
   phase offset into `sat.code_phase` seeding.
   - **Tests:** in `test/gps_l1c_p.jl`, drop a synthetic 1800-bit buffer
     with a known phase and polarity, confirm `_full_phase_search` returns
     the expected phase / polarity / found state. Tolerance: inject up to
     36 random bit-flips, still lock; 37+ flips, reject. End-to-end:
     extend `test/track.jl`'s L1C-P case to assert that `found` becomes
     true after 1800 primary periods of clean simulated data and that
     `sat.code_phase` snaps to the injected offset.
5. **Replace L1C-D's `return false` with `SyncResult(true, 0, +1)`.**
   One-line change in `src/gpsl1c_d.jl`. Requires the
   [GNSSSignals.jl PR #61](https://github.com/JuliaGNSS/GNSSSignals.jl/pull/61)
   bump to have landed in `Project.toml`'s `[compat]`; otherwise
   `_calc_num_code_blocks_that_form_a_bit` returns 2 (from the stale
   50 Hz `get_data_frequency`) and the post-sync prompt-accumulator
   integrates pairs of symbols, halving the symbol rate to the decoder.
   - **Tests:** in `test/gps_l1c_d.jl`, assert
     `is_upcoming_integration_new_bit(GPSL1C_D(), 0x0, 1)` returns
     `SyncResult(found=true, phase=0, polarity=+1)`. Integration test in
     `test/track.jl`: a few hundred primary periods of simulated L1C-D
     yield one bit per primary period in `get_bits(tracked_signal)`.
6. **Raise tolerances on the three existing detectors** to the table
   values. This is the only behavioral change to existing signals; gate
   it behind its own commit and run the existing track tests to confirm
   no regressions.
   - **Tests:** update `test/gps_l1ca.jl`, `test/galileo_e1b.jl`,
     `test/gps_l5i.jl` to add an "N-bit error at lock window" case per
     signal: inject the per-signal tolerance number of errors into a
     correctly-phased buffer and assert `found = true`; inject one more
     and assert `found = false`.
7. **Documentation.** Three doc touchpoints:
   - `docs/src/tracking_state.md` — the existing
     `has_bit_or_secondary_code_been_found` row in the accessors table
     ([line 135](../../docs/src/tracking_state.md#L135)) and the
     secondary-code wrap paragraph ([line 191](../../docs/src/tracking_state.md#L191))
     mention sync. Add a short subsection right after explaining the
     detector contract (per-signal sync window, Hamming tolerance, L1C-P
     two-stage fill+search behavior, where `phase` lands in
     `sat.code_phase`).
   - **Docstrings.** New `SyncResult` and `_try_match` get full
     docstrings. The five per-signal `is_upcoming_integration_new_bit`
     methods grow a one-line docstring stating their window / template /
     tolerance / phase semantics — replacing the long "deferred" prose
     in [gpsl1c_d.jl:1-14](../../src/gpsl1c_d.jl#L1-L14) and
     [gpsl1c_p.jl:1-24](../../src/gpsl1c_p.jl#L1-L24). Add a docstring
     to `get_code_block_buffer_type` listing the per-signal widths.
   - `CHANGELOG.md` — under the 2.0.0 entry, add a "Bit-sync /
     secondary-code detection" subsection summarizing: new `SyncResult`
     return type, Hamming-tolerance matching, L1C-D one-symbol-per-block
     sync (no edge search), L1C-P 1800-chip overlay lock, `BitBuffer` is
     now parameterized on `B<:Unsigned`. Also call out the GNSSSignals.jl
     compat bump and the symbol-rate semantics of `get_data_frequency`.

Each step is independently revertible. Step 1 is the biggest mechanical
change; 3–5 are localized; 7 is doc-only and can land alongside or after
step 6.

## Open questions

- **Should L1C-P run an incremental search** at e.g. `num_code_blocks =
  225, 450, ...` for early-lock under good conditions? 8 partial searches
  cost ~560 μs total, run once per multi-second integration boundary. I'd
  defer until someone complains about 18 s acquisition.
- **What's the right BitIntegers version pin?** Last public release is
  0.3.x; the API we use (`@define_integers`, `count_ones`, `xor`, shifts)
  is stable since 0.2. Pin `[6e96ef4d] BitIntegers = "0.3"`.
- **Does `count_ones` on `UInt1800` codegen to a popcnt loop or
  byte-table?** Bench shows ~17 ns for a 1800-bit popcount, which is
  consistent with ~28 × 64-bit popcnt instructions (Zen 4: ~1 cycle each).
  Compiler is doing the right thing; no special handling needed.
