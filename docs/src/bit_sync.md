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
and `get_soft_bits` once it returns `true` — each entry is the coherently
accumulated complex prompt of one navigation bit, from which the receiver
forms its own hard decision. The rest of the machinery is internal.

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

Per-signal contract. *Min-to-fire* is the smallest `num_code_blocks` the detector accepts before it tries a match (smaller windows return `SyncResult(false, …)` without searching). *Buffer width* is the sliding window's container type, picked by `get_code_block_buffer_type(signal)`. The GPS L1 C/A bit-edge detector needs `2 × symbol` (40 blocks) so a full bit can slide past either polarity; the secondary-code signals (GPS L5I, GPS L1C-P) only need **one** full secondary period, because their rotation search (`_secondary_code_search`) matches at any alignment and recovers the phase — so they lock in `1 ×` the period worst case, not `2 ×`.

| Signal | Min-to-fire | Buffer width | Template | Tolerance | Phase | Blocks per symbol |
|--------|-------------|--------------|----------|-----------|-------|---------------------|
| GPS L1 C/A | 40 blocks | `UInt64` (40 bits used) | bit-edge `0xfffff` (20+20) | 1 error (2.5 %) | 0 | 20 |
| Galileo E1B | n/a | `UInt8` (unused) | trivial (1 block per symbol) | n/a | 0 | 1 |
| GPS L5I | 10 blocks | `UInt32` (low 10 bits used) | NH10 `0x035` (rotation search) | 0 errors (2.5 % → exact match) | `0..9` | 10 (secondary code) |
| GPS L5Q | 20 blocks | `UInt32` (low 20 bits used) | NH20 (rotation search) | 0 errors (2.5 % → exact match) | `0..19` | 20 (secondary code, pilot) |
| GPS L1C-D | n/a | `UInt8` (unused) | trivial (1 block per symbol) | n/a | 0 | 1 |
| GPS L1C-P | 1800 blocks | `UInt1800` (exact width) | per-PRN overlay (rotation search) | 45 errors (2.5 %) | `0..1799` | n/a (pilot) |
| GPS L2CM | n/a | `UInt8` (unused) | trivial (1 block per symbol) | n/a | 0 | 1 |
| GPS L2CL | never fires | `UInt8` (unused) | none (dataless pilot, no secondary code) | n/a | 0 | n/a (pilot) |
| Galileo E1C | 25 blocks | `UInt32` (low 25 bits used) | CS25 (rotation search) | 0 errors (2.5 % → exact match) | `0..24` | 25 (secondary code, pilot) |
| Galileo E5a-I | 20 blocks | `UInt32` (low 20 bits used) | CS20 (rotation search) | 0 errors (2.5 % → exact match) | `0..19` | 20 (secondary code) |
| Galileo E5a-Q | 100 blocks | `UInt128` (low 100 bits used) | per-PRN CS100 (rotation search) | 2 errors (2.5 %) | `0..99` | 100 (secondary code, pilot) |

The buffer-width type threads through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF}` as a type parameter. The L1C-P case uses an exact-width `UInt1800` defined via `BitIntegers.@define_integers 1800`; the other signals use built-in `UInt8` / `UInt32` / `UInt64` / `UInt128`. The secondary-code signals other than GPS L5I / L1C-P get their packed reference from the generic `_packed_secondary_code(::Type{B}, ::AbstractGNSSSignal, prn)`, which derives it directly from `get_secondary_code(signal)` — so a new secondary-coded signal needs no bespoke template, only a wide enough `get_code_block_buffer_type`.

The 2.5 % tolerance is a package-wide default and adjustable per-signal via dispatch:

```jldoctest tolerance_override
julia> using Tracking, GNSSSignals

julia> Tracking.get_bit_edge_or_secondary_code_tolerance(GPSL1CA())  # default
0.025

julia> # Loosen the L1 C/A ceiling to 5 % (= 2 errors over 40 blocks) for low-C/N₀ work.
       Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL1CA) = 0.05;

julia> Tracking.get_bit_edge_or_secondary_code_tolerance(GPSL1CA())  # after override
0.05
```

Each detector reads the trait at its call site and converts to an integer error budget via `floor(Int, tolerance × window_size)`, so the override picks up the next time `detect_bit_or_secondary_code_sync` runs — no `TrackState` rebuild needed. Galileo E1B and GPS L1C-D have trivially-true detectors and ignore the trait.

### Lifecycle of a `BitBuffer`

Two distinct phases, separated by the `found::Bool` flag:

1. **Pre-sync search** (`found = false`). Each completed integration shifts one bit (the sign of the prompt's real part) into `code_block_buffer::B`. `detect_bit_or_secondary_code_sync` is called on every shift; while it returns `SyncResult(false, ...)` the loop keeps integrating one primary code period at a time. The **bit-edge** detector (GPS L1 C/A) only matches at a true bit boundary, so the call that flips `found` lands exactly on the edge between two data symbols and reports `phase = 0`. The **secondary-code** detectors (GPS L5I, GPS L1C-P) instead lock as soon as one full secondary period has been buffered — at *any* alignment — and report the recovered `phase` (the upcoming integration's secondary chip); the integration cadence then re-aligns to the secondary-code boundary (see below).

2. **Post-sync accumulation** (`found = true`). The post-sync branch in `Tracking.buffer` ignores `code_block_buffer` and instead accumulates the complex prompt into `prompt_accumulator`. (Pilot signals such as GPS L1C-P carry no data bits, so for them this branch is a no-op — the buffer just retains `found` / `secondary_phase` / `polarity`.) Because each post-sync integration spans to the next secondary-code (data-bit) boundary, bit decoding stays aligned even when sync fired mid-period. Each integration also bumps `prompt_accumulator_integrated_code_blocks`. Once that counter reaches the per-signal "blocks per symbol" value above — `_calc_num_code_blocks_that_form_a_bit(signal) = get_code_frequency(signal) / (get_code_length(signal) * get_data_frequency(signal))` — one decoded bit is committed to `buffer::UInt128` and the accumulator resets to zero. For Galileo E1B and GPS L1C-D the counter is `1`, so one symbol commits per integration; for GPS L1 C/A it's `20`, so the loop counts 20 primary-code periods (≈ 20 × 1023 chips = 20460 chips of `code_phase` advance, modulo wrap) per data bit; for GPS L5I it's `10`. The polarity flag flips the accumulator's sign at commit time when the detector locked at negative polarity, so downstream consumers always see `1 = data symbol 0`.

Pilot signals (`get_data_frequency = 0 Hz`, e.g. GPS L1C-P) never enter the post-sync accumulation branch — their `bit_buffer` carries the recovered secondary-code phase but no decoded bits, and the post-sync work is purely the [code-phase seeding](#Code-phase-seeding-from-the-secondary-code-phase) described next.

### Code-phase seeding from the secondary-code phase

When a signal whose detector exposes a `secondary_phase` syncs — GPS L5I and GPS L1C-P — that phase is used to seed `TrackedSat.code_phase` so subsequent wrap-mod-[`current_code_wrap`](@ref) arithmetic gives the absolute position in the longest secondary-code cycle. The seeding follows a fallback chain: the synced signal with the largest `(primary × secondary)` code length wins.

Signals with `secondary_code_length == 1` (bit-edge only, e.g. GPS L1 C/A) do **not** carry an explicit `secondary_phase` to snap — there's no per-PRN overlay to recover a chip offset from. Instead, their post-sync bit-edge alignment is captured by **two complementary mechanisms**:

1. The `BitBuffer.prompt_accumulator_integrated_code_blocks` counter tracks "how many of the next N primary periods have I integrated since the last bit commit." `reset(bit_buffer)` preserves it, so the bit cadence survives intra-call resets without re-syncing.

2. The wrap returned by [`current_code_wrap`](@ref) widens from `primary` to `primary × blocks_per_data_bit` (e.g. 1023 → 20460 for L1 C/A) the moment `bit_buffer.found` flips to `true`. From that point on `mod(code_phase, primary)` continues to give the replica phase, while `div(code_phase, primary)` reads off which primary period within the data bit we're in. The transition is one-shot: on the call that flips `found`, `code_phase` is implicitly the start of a fresh data bit because the detector only matches at a true bit boundary (the `[0, primary)` range of `code_phase` then represents primary period 0 of the new bit), so no explicit snap is needed.
