# Bit and Secondary-Code Sync

Once a satellite is being tracked, the next milestone is locking onto its
**bit/symbol boundary** (and, for pilots like GPS L1C-P, the
**secondary-code phase**). These are what let the receiver:

- Switch from 1-ms coherent integrations to full-symbol coherent
  integrations (e.g. 20 ms for L1 C/A), driving the PLL phase-noise floor
  down.
- Hand over a decoded bit stream to a navigation-message parser.
- Reconstruct the absolute code phase modulo the longest secondary-code
  cycle, so downstream consumers (e.g. PositionVelocityTime.jl) can
  pseudorange.

This page is a deep dive into the mechanism. The high-level user only
needs to know: call `has_bit_or_secondary_code_been_found` per signal
(see the [per-signal accessor table](tracking_state.md#What-you-can-read))
and `get_bits` once it returns `true`. The rest of the machinery is
internal.

## Code-phase wrap period

The shared `TrackedSat.code_phase` wraps at two distinct values depending on the per-signal sync state:

- **Before any signal has synced** — wrap is the largest primary code length across the sat's signals. For a sat tracking only L1 C/A this is 1023 chips. The wrap is intentionally narrow at this stage because we don't yet know which data bit or secondary-code chip the current primary period belongs to.

- **After a signal syncs** — its contribution to the wrap widens to one full *symbol period*: `primary × secondary_code_length` for a pilot, or `primary × blocks_per_data_bit` for a data-bearing signal. The shared wrap is then the `max` across signals, so the longest synced signal pins it.

Concrete values:

| Sat tracks         | Before sync | After full sync                          |
|--------------------|-------------|------------------------------------------|
| GPS L1 C/A only    | 1023        | 1023 × 20 = **20460** (one 50 Hz data bit) |
| GPS L5I only       | 10230       | 10230 × 10 = **102300** (one NH10 cycle / one data bit) |
| GPS L1C-P only     | 10230       | 10230 × 1800 ≈ **18.4 M** (one overlay-code cycle, ≈ 18 s) |
| L1C-P + L1C-D + L1CA | 10230 (longest primary) | 18.4 M (L1C-P dominates) |

The post-sync widening is what lets downstream consumers (e.g. PositionVelocityTime.jl) distinguish which primary-code period within the symbol they're currently in — `mod(code_phase, primary)` gives the per-signal replica phase, while `div(code_phase, primary)` gives the symbol-internal position.

[`max_code_length`](@ref) returns the *upper bound* (the post-full-sync value) at compile time. [`current_code_wrap`](@ref) returns the *runtime* value honoring the current per-signal sync state — this is what the inner loop's `mod` actually uses.

```@docs
max_code_length
current_code_wrap
```

## Bit-sync and secondary-code-sync detection

Each per-signal `BitBuffer` runs an `detect_bit_or_secondary_code_sync` detector against the running buffer of primary-code-block signs. The detector returns a `SyncResult` containing whether sync was found, the secondary-code phase (chip offset within the secondary code, used for code-phase seeding — see below), and the locked polarity (±1).

There are three detector families:

- **Soft, maximum-energy CFAR** — GPS L1 C/A (bit-edge, `_detect_bit_edge_cfar`) and the short-secondary-code signals GPS L5I/L5Q and Galileo E1C/E5a-I/E5a-Q (overlay-rotation, `_detect_secondary_code_cfar`). These accumulate per-hypothesis, coherently-summed bin energy in `PhaseAccumulators` and lock when the peak hypothesis beats its runner-up with a Student-t confidence (`get_bit_edge_detection_confidence`, default `0.999`). They self-pace with C/N₀ and fire only at the winning hypothesis's own boundary — so they report **`phase = 0`** (the upcoming integration starts a fresh data bit / secondary-code period). Selected by `uses_soft_bit_edge_detection` / `uses_soft_secondary_code_detection`.
- **Hard-decision rotation/Hamming sweep** (`_secondary_code_search`) — among the currently implemented signals, GPS L1C-P alone. It matches packed prompt *signs* against the known overlay at every rotation, accepting the best within `get_bit_edge_or_secondary_code_tolerance`; it locks in **one** full period worst case (matches at any alignment) and reports the recovered secondary-chip `phase`.
- **Trivial / none** — Galileo E1B, GPS L1C-D and GPS L2CM broadcast one channel symbol per primary period, so their detector returns `SyncResult(true, 0, +1)` immediately; GPS L2CL is a dataless pilot with no sync.

*Min-to-fire* is the smallest `num_code_blocks` the detector accepts before it can lock. The soft CFAR detectors need `2 × period` blocks before a runner-up (and hence a decision) exists, then fire at the next winning-hypothesis boundary; the hard sweep needs `1 × period`. *Buffer width* is `get_code_block_buffer_type(signal)`; note it sizes the packed prompt-sign buffer that **only the hard sweep consults** — the soft detectors read `PhaseAccumulators` instead, so for them the packed buffer is vestigial.

| Signal | Detector | Min-to-fire | Buffer width | Accept rule | Phase | Blocks per symbol |
|--------|----------|-------------|--------------|-------------|-------|---------------------|
| GPS L1 C/A | soft bit-edge CFAR | 40 blocks (2 × 20) | `UInt64` (vestigial) | confidence 0.999 | 0 | 20 |
| Galileo E1B | trivial | n/a | `UInt8` (unused) | always | 0 | 1 |
| GPS L5I | soft secondary CFAR | 20 blocks (2 × NH10) | `UInt32` (vestigial) | confidence 0.999 | 0 | 10 |
| GPS L5Q | soft secondary CFAR | 40 blocks (2 × NH20) | `UInt32` (vestigial) | confidence 0.999 | 0 | 20 (pilot) |
| GPS L1C-D | trivial | n/a | `UInt8` (unused) | always | 0 | 1 |
| GPS L1C-P | **hard** rotation sweep | 1800 blocks | `UInt1800` (exact width) | 45 errors (2.5 %) | `0..1799` | n/a (pilot) |
| GPS L2CM | trivial | n/a | `UInt8` (unused) | always | 0 | 1 |
| GPS L2CL | never fires | `UInt8` (unused) | none (dataless pilot) | n/a | 0 | n/a (pilot) |
| Galileo E1C | soft secondary CFAR | 50 blocks (2 × CS25) | `UInt32` (vestigial) | confidence 0.999 | 0 | 25 (pilot) |
| Galileo E5a-I | soft secondary CFAR | 40 blocks (2 × CS20) | `UInt32` (vestigial) | confidence 0.999 | 0 | 20 |
| Galileo E5a-Q | soft secondary CFAR | 200 blocks (2 × CS100) | `UInt128` (vestigial) | confidence 0.999 | 0 | 100 (pilot) |

The buffer-width type threads through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF}` as a type parameter. The L1C-P case uses an exact-width `UInt1800` defined via `BitIntegers.@define_integers 1800`. The soft secondary CFAR detectors take the ±1 overlay chip directly from `get_secondary_code(signal)` per rotation, so a new short-secondary-code signal needs no bespoke template — only for `uses_soft_secondary_code_detection` to return `true` (secondary length `1 < N ≤ 100`).

The soft CFAR confidence is package-wide and adjustable per signal via `get_bit_edge_detection_confidence`; the hard-sweep Hamming tolerance is adjustable per (hard-path) signal via `get_bit_edge_or_secondary_code_tolerance`:

```jldoctest tolerance_override
julia> using Tracking, GNSSSignals

julia> Tracking.get_bit_edge_or_secondary_code_tolerance(GPSL1C_P())  # default
0.025

julia> # Loosen the L1C-P ceiling to 5 % (= 90 errors over 1800 blocks) for low-C/N₀ work.
       Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL1C_P) = 0.05;

julia> Tracking.get_bit_edge_or_secondary_code_tolerance(GPSL1C_P())  # after override
0.05
```

The hard sweep reads the tolerance at its call site and converts to an integer error budget via `floor(Int, tolerance × window_size)`, so the override picks up the next time `detect_bit_or_secondary_code_sync` runs — no `TrackState` rebuild needed. It has no effect on the soft-detector signals (tune their confidence instead); the trivial detectors ignore both.

### Lifecycle of a `BitBuffer`

Two distinct phases, separated by the `found::Bool` flag:

1. **Pre-sync search** (`found = false`). Each completed integration folds the prompt into the running state and the detector is called; while it returns `SyncResult(false, ...)` the loop keeps integrating one primary code period at a time. The **soft CFAR** detectors (GPS L1 C/A bit-edge; GPS L5I/L5Q and Galileo E1C/E5a-I/E5a-Q overlay) fold each prompt into `PhaseAccumulators` and lock only at the winning hypothesis's own boundary — so the call that flips `found` lands exactly on a data-bit / secondary-code-period boundary and reports `phase = 0` (the upcoming integration starts a fresh period). The **hard rotation-sweep** detector (GPS L1C-P) instead shifts each prompt sign into `code_block_buffer::B` and locks as soon as one full overlay period has been buffered — at *any* alignment — reporting the recovered `phase` (the upcoming integration's secondary chip); the integration cadence then re-aligns to that boundary (see below).

2. **Post-sync accumulation** (`found = true`). The post-sync branch in `Tracking.buffer` ignores `code_block_buffer` and instead accumulates the complex prompt into `prompt_accumulator`. (Pilot signals such as GPS L1C-P carry no data bits, so for them this branch is a no-op — the buffer just retains `found` / `secondary_phase` / `polarity`.) Because each post-sync integration spans to the next secondary-code (data-bit) boundary, bit decoding stays aligned even when sync fired mid-period. Each integration also bumps `prompt_accumulator_integrated_code_blocks`. Once that counter reaches the per-signal "blocks per symbol" value above — `_calc_num_code_blocks_that_form_a_bit(signal) = get_code_frequency(signal) / (get_code_length(signal) * get_data_frequency(signal))` — one decoded bit is committed as a **soft bit** — the polarity-corrected coherent prompt sum, sign = hard decision, magnitude = confidence — pushed onto the unbounded `soft_bits` vector (the packed hard bits of `get_bits` are derived from the signs on demand, newest 128 only), and the accumulator resets to zero. For Galileo E1B and GPS L1C-D the counter is `1`, so one symbol commits per integration; for GPS L1 C/A it's `20`, so the loop counts 20 primary-code periods (≈ 20 × 1023 chips = 20460 chips of `code_phase` advance, modulo wrap) per data bit; for GPS L5I it's `10`. The polarity flag flips the accumulator's sign at commit time when the detector locked at negative polarity, so downstream consumers always see `1 = data symbol 0`.

Pilot signals (`get_data_frequency = 0 Hz`, e.g. GPS L1C-P) never enter the post-sync accumulation branch — their `bit_buffer` carries the recovered secondary-code phase but no decoded bits, and the post-sync work is purely the [code-phase seeding](#Code-phase-seeding-from-the-secondary-code-phase) described next.

### Code-phase seeding from the secondary-code phase

When a signal with a secondary code (`secondary_code_length > 1`) syncs — GPS L5I/L5Q, Galileo E1C/E5a-I/E5a-Q, or GPS L1C-P — its `secondary_phase` seeds `TrackedSat.code_phase` so subsequent wrap-mod-[`current_code_wrap`](@ref) arithmetic gives the absolute position in the longest secondary-code cycle. The seeding follows a fallback chain: the synced signal with the largest `(primary × secondary)` code length wins. For the soft-CFAR signals `secondary_phase` is always `0` (they fire on a period boundary, so the upcoming integration is at chip 0); for the hard-sweep L1C-P it is the recovered chip offset. Either way it is a multiple-of-primary snap into the correct secondary window.

Signals with `secondary_code_length == 1` (bit-edge only, e.g. GPS L1 C/A) do **not** carry an explicit `secondary_phase` to snap — there's no per-PRN overlay to recover a chip offset from. Instead, their post-sync bit-edge alignment is captured by **two complementary mechanisms**:

1. The `BitBuffer.prompt_accumulator_integrated_code_blocks` counter tracks "how many of the next N primary periods have I integrated since the last bit commit." `reset(bit_buffer)` preserves it, so the bit cadence survives intra-call resets without re-syncing.

2. The wrap returned by [`current_code_wrap`](@ref) widens from `primary` to `primary × blocks_per_data_bit` (e.g. 1023 → 20460 for L1 C/A) the moment `bit_buffer.found` flips to `true`. From that point on `mod(code_phase, primary)` continues to give the replica phase, while `div(code_phase, primary)` reads off which primary period within the data bit we're in. The transition is one-shot: on the call that flips `found`, `code_phase` is implicitly the start of a fresh data bit because the detector only matches at a true bit boundary (the `[0, primary)` range of `code_phase` then represents primary period 0 of the new bit), so no explicit snap is needed.
