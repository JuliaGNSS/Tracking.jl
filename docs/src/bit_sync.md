# Bit and Secondary-Code Sync

Until a signal's **bit boundary** (or, for a pilot, its **secondary-code
phase**) is found it integrates one primary code period at a time. Sync unlocks
longer coherent integrations (e.g. 20 ms for GPS L1 C/A), decoded bits, and the
code phase modulo a whole symbol, which a pseudorange needs.

Users only need `has_bit_or_secondary_code_been_found` per signal (see the
[per-signal accessor table](tracking_state.md#What-you-can-read)) to gate
`get_soft_bits`. The detectors and the bit buffer are TrackingLoops'; see
[`BitBuffer`](@extref TrackingLoops.BitBuffer).

## Code-phase wrap period

The shared `TrackedSat.code_phase` wraps at the least common multiple of the
signals' wraps (#129), each of which depends on the signal's sync state:

- **Before sync**: its primary code length, since the bit or secondary chip of
  the current period is unknown.
- **After sync**: one *symbol period*, `primary × secondary_code_length` for a
  pilot and `primary × blocks_per_data_bit` for a data signal.

| Sat tracks         | Before sync | After full sync                          |
|--------------------|-------------|------------------------------------------|
| GPS L1 C/A only    | 1023        | 1023 × 20 = **20460** (one 50 Hz data bit) |
| GPS L5I only       | 10230       | 10230 × 10 = **102300** (one NH10 cycle / one data bit) |
| GPS L1C-P only     | 10230       | 10230 × 1800 ≈ **18.4 M** (one overlay-code cycle, ≈ 18 s) |
| L1C-P + L1C-D + L1CA | 10230 | 18.4 M (L1C-P dominates) |

`mod(code_phase, primary)` is a signal's replica phase, `div(code_phase, primary)`
its position within the symbol. [`max_code_wrap`](@ref) is the post-sync bound,
[`current_code_wrap`](@ref) the runtime wrap.

```@docs
max_code_wrap
current_code_wrap
```

## Sync detection

Each signal's detector (`detect_bit_or_secondary_code_sync`) comes in one of three
families:

- **Soft, maximum-energy CFAR**: GPS L1 C/A's bit edge and every short secondary
  code (`1 < N ≤ 100`). It locks when the peak hypothesis beats its runner-up at
  the confidence `get_bit_edge_detection_confidence` (default `0.999`), needs
  two periods first, and fires on a period boundary, so it reports **`phase = 0`**.
- **Hard rotation sweep**: the 1800-chip overlay pilots GPS L1C-P and BeiDou
  B1C-P. It matches the prompt signs against the overlay at every rotation within
  `get_bit_edge_or_secondary_code_tolerance`, locks after one period at any
  alignment and reports the recovered chip `phase`.
- **Trivial**: one symbol per primary period (Galileo E1B, E6-B, GPS L1C-D, L2CM,
  BeiDou B2b-I, B1C-D), and Galileo E5a-QP (no data, no overlay), sync on the first
  block. GPS L2CL, a pilot without secondary code, never syncs and stays pre-sync.

*Min-to-fire* is the number of blocks before a detector can lock.

| Signal | Detector | Min-to-fire | Accept rule | Phase | Blocks per symbol |
|--------|----------|-------------|-------------|-------|---------------------|
| GPS L1 C/A | soft bit-edge CFAR | 40 blocks (2 × 20) | confidence 0.999 | 0 | 20 |
| Galileo E1B | trivial | n/a | always | 0 | 1 |
| GPS L5I | soft secondary CFAR | 20 blocks (2 × NH10) | confidence 0.999 | 0 | 10 |
| GPS L5Q | soft secondary CFAR | 40 blocks (2 × NH20) | confidence 0.999 | 0 | 20 (pilot) |
| GPS L1C-D | trivial | n/a | always | 0 | 1 |
| GPS L1C-P | **hard** rotation sweep | 1800 blocks | 45 errors (2.5 %) | `0..1799` | n/a (pilot) |
| GPS L2CM | trivial | n/a | always | 0 | 1 |
| GPS L2CL | none (dataless pilot) | n/a | n/a | 0 | n/a (pilot) |
| Galileo E1C | soft secondary CFAR | 50 blocks (2 × CS25) | confidence 0.999 | 0 | 25 (pilot) |
| Galileo E5a-I | soft secondary CFAR | 40 blocks (2 × CS20) | confidence 0.999 | 0 | 20 |
| Galileo E5a-Q | soft secondary CFAR | 200 blocks (2 × CS100) | confidence 0.999 | 0 | 100 (pilot) |
| Galileo E5b-I | soft secondary CFAR | 8 blocks (2 × CS4) | confidence 0.999 | 0 | 4 |
| Galileo E5b-Q | soft secondary CFAR | 200 blocks (2 × CS100) | confidence 0.999 | 0 | 100 (pilot) |
| Galileo E5a-QP | trivial | n/a | always | 0 | n/a (pilot) |
| Galileo E6-B | trivial | n/a | always | 0 | 1 |
| Galileo E6-C | soft secondary CFAR | 200 blocks (2 × CS100) | confidence 0.999 | 0 | 100 (pilot) |
| BeiDou B1I | soft secondary CFAR | 40 blocks (2 × NH20) | confidence 0.999 | 0 | 20 |
| BeiDou B3I | soft secondary CFAR | 40 blocks (2 × NH20) | confidence 0.999 | 0 | 20 |
| BeiDou B2b-I | trivial | n/a | always | 0 | 1 |
| BeiDou B2a data | soft secondary CFAR | 10 blocks (2 × 5) | confidence 0.999 | 0 | 5 |
| BeiDou B2a pilot | soft secondary CFAR | 200 blocks (2 × 100) | confidence 0.999 | 0 | 100 (pilot) |
| BeiDou B1C data | trivial | n/a | always | 0 | 1 |
| BeiDou B1C pilot | **hard** rotation sweep | 1800 blocks | 45 errors (2.5 %) | `0..1799` | n/a (pilot) |

BeiDou B1I/B3I GEO satellites (PRN 1-5, 59-63) carry no NH20 overlay and
broadcast D2, so the detector never locks on them: they track and range but stay
pre-sync (see [Known limitations](signals.md#Known-limitations)).

The confidence (soft detectors) and the tolerance (hard sweep only) can be
overridden per signal; the hard sweep reads its tolerance on every call, so no
`TrackState` rebuild is needed:

```jldoctest tolerance_override
julia> using Tracking, TrackingLoops, GNSSSignals

julia> TrackingLoops.get_bit_edge_or_secondary_code_tolerance(GPSL1C_P())  # default
0.025

julia> # Loosen the L1C-P ceiling to 5 % (= 90 errors over 1800 blocks) for low-C/N₀ work.
       TrackingLoops.get_bit_edge_or_secondary_code_tolerance(::GPSL1C_P) = 0.05;

julia> TrackingLoops.get_bit_edge_or_secondary_code_tolerance(GPSL1C_P())  # after override
0.05
```

### Code-phase seeding

When a signal with a secondary code syncs, its `secondary_phase` snaps
`TrackedSat.code_phase` into the right secondary-code window, by a multiple of
the primary length; the within-period phase is kept (#117). Of several synced
signals, the one with the longest `primary × secondary` length wins.

A signal without secondary code (e.g. GPS L1 C/A) needs no snap: its detector
locks at a true bit boundary, so the widened wrap starts at the bit's first
period, and the bit buffer's block counter, which `reset` preserves, keeps the
bit cadence across `track!` calls.
