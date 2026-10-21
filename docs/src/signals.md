# Signal support matrix

Every concrete signal type [GNSSSignals.jl](https://github.com/JuliaGNSS/GNSSSignals.jl)
defines, and what Tracking.jl does with it.
[`test/signal_coverage.jl`](https://github.com/JuliaGNSS/Tracking.jl/blob/master/test/signal_coverage.jl)
checks the table against the live type tree, so a GNSSSignals release that adds
a signal fails the test suite rather than user code.

## What "supported" means here

**Tracking is not navigation decoding, and neither is PVT.**

- *Sync feature* is what the tracker locks onto to leave its one-block pre-sync
  integration: a data-bit edge, a secondary-code phase, or nothing where the
  primary code period already is the integration unit.
- *Navigation data* is what the signal carries. A pilot or acquisition aid
  carries none: it yields code phase, carrier phase, Doppler and C/N₀, and
  `get_soft_bits` stays empty. Decoding is
  [GNSSDecoder.jl](https://github.com/JuliaGNSS/GNSSDecoder.jl)'s job, and a
  satellite can be ranged on without it.

**Evidence level.** Every row is verified in this repository's test suite
against clean synthetic GNSSSignals replicas: `test/<constellation>_<signal>.jl`
pins defaults and sync behaviour, and `test/signal_coverage.jl` runs each signal
through [`track`](@ref). Simulated-FPGA, hardware-replay and live-RF validation
are tracked in
[GNSSReceiver.jl#130](https://github.com/JuliaGNSS/GNSSReceiver.jl/issues/130).

## The matrix

*Blocks per integration* is the coherent integration length in primary code
blocks a fresh [`TrackedSignal`](@ref) starts at
([`default_num_code_blocks_to_integrate`](@extref TrackingLoops.default_num_code_blocks_to_integrate)),
and in brackets the ceiling
([`max_num_code_blocks_to_integrate`](@extref TrackingLoops.max_num_code_blocks_to_integrate))
[`set_preferred_num_code_blocks_to_integrate!`](@ref) may raise it to: one data
bit, overlay period or code cycle. Until sync, one block is integrated regardless.

| Signal | Band | Primary code | Blocks per integration | Sync feature | Navigation data | PRNs |
|:--|:--|:--|:--|:--|:--|:--|
| `GPSL1CA` | `L1` | 1023 chips, 1 ms | 1 (max 20) | soft CFAR bit-edge search | LNAV, 50 sym/s | 1-37 |
| `GPSL1C_D` | `L1` | 10230 chips, 10 ms | 1 (max 1) | none needed — 1 symbol per code period | CNAV-2, 100 sym/s | 1-63 |
| `GPSL1C_P` | `L1` | 10230 chips, 10 ms | 1 (max 1800) | hard rotation sweep, 1800-chip overlay | none (pilot) | 1-63 |
| `GPSL2CM` | `L2` | 10230 chips, 20 ms | 1 (max 1) | none needed — 1 symbol per code period | CNAV, 50 sym/s | 1-63 |
| `GPSL2CL` | `L2` | 767250 chips, 1.5 s | 1 (max 1) | none — the code period is the integration | none (pilot) | 1-63 |
| `GPSL5I` | `L5` | 10230 chips, 1 ms | 1 (max 10) | soft CFAR overlay search, NH10 | CNAV, 100 sym/s | 1-37 |
| `GPSL5Q` | `L5` | 10230 chips, 1 ms | 1 (max 20) | soft CFAR overlay search, NH20 | none (pilot) | 1-37 |
| `GalileoE1B` | `L1` | 4092 chips, 4 ms | 1 (max 1) | none needed — 1 symbol per code period | I/NAV, 250 sym/s | 1-50 |
| `GalileoE1B_BOC11` | `L1` | 4092 chips, 4 ms | 1 (max 1) | none needed — 1 symbol per code period | I/NAV, 250 sym/s | 1-50 |
| `GalileoE1C` | `L1` | 4092 chips, 4 ms | 1 (max 25) | soft CFAR overlay search, CS25 | none (pilot) | 1-50 |
| `GalileoE1C_BOC11` | `L1` | 4092 chips, 4 ms | 1 (max 25) | soft CFAR overlay search, CS25 | none (pilot) | 1-50 |
| `GalileoE5aI` | `L5` | 10230 chips, 1 ms | 1 (max 20) | soft CFAR overlay search, CS20 | F/NAV, 50 sym/s | 1-50 |
| `GalileoE5aQ` | `L5` | 10230 chips, 1 ms | 1 (max 100) | soft CFAR overlay search, per-PRN CS100 | none (pilot) | 1-50 |
| `GalileoE5aQP` | `L5` | 330 chips, 64.5 µs | 31 (max 31) | none — no data, no overlay | none (acquisition aid) | 1-40 |
| `GalileoE5bI` | `E5b` | 10230 chips, 1 ms | 1 (max 4) | soft CFAR overlay search, CS4 | I/NAV, 250 sym/s | 1-50 |
| `GalileoE5bQ` | `E5b` | 10230 chips, 1 ms | 1 (max 100) | soft CFAR overlay search, per-SVID CS100 | none (pilot) | 1-50 |
| `GalileoE6B` | `E6` | 5115 chips, 1 ms | 1 (max 1) | none needed — 1 symbol per code period | C/NAV, 1000 sym/s | 1-50 |
| `GalileoE6C` | `E6` | 5115 chips, 1 ms | 1 (max 100) | soft CFAR overlay search, per-PRN CS100 | none (pilot) | 1-50 |
| `BeiDouB1I` | `B1I` | 2046 chips, 1 ms | 1 (max 20) | soft CFAR overlay search, NH20 | D1, 50 sym/s (see GEO note) | 1-63 |
| `BeiDouB3I` | `B3I` | 10230 chips, 1 ms | 1 (max 20) | soft CFAR overlay search, NH20 | D1, 50 sym/s (see GEO note) | 1-63 |
| `BeiDouB1C_D` | `L1` | 10230 chips, 10 ms | 1 (max 1) | none needed — 1 symbol per code period | B-CNAV1, 100 sym/s | 1-63 |
| `BeiDouB1C_P` | `L1` | 10230 chips, 10 ms | 1 (max 1800) | hard rotation sweep, 1800-chip overlay | none (pilot) | 1-63 |
| `BeiDouB2aI` | `L5` | 10230 chips, 1 ms | 1 (max 5) | soft CFAR overlay search, 5-chip | B-CNAV2, 200 sym/s | 1-63 |
| `BeiDouB2aQ` | `L5` | 10230 chips, 1 ms | 1 (max 100) | soft CFAR overlay search, per-PRN 100-chip | none (pilot) | 1-63 |
| `BeiDouB2bI` | `E5b` | 10230 chips, 1 ms | 1 (max 1) | none needed — 1 symbol per code period | B-CNAV3, 1000 sym/s (see below) | 6-58 |

`PRNs` is the range the upstream code table defines; this package does not check
it, and other PRNs give an all-zero or out-of-bounds code. BeiDou B2b-I's table
has 63 columns, but the ICD defines codes for PRN 6-58 only.

Signals sharing a band can be tracked in one [`BandMeasurement`](@ref)
regardless of constellation: Galileo E5b and BeiDou B2b share the 1207.14 MHz
carrier (`E5b` for both), BeiDou B2a shares `L5` with GPS L5 and Galileo E5a,
and BeiDou B1C shares `L1` with GPS L1 C/A, GPS L1C and Galileo E1.

## Galileo E5a-QP, the acquisition aid

E5a-QP (OS SIS ICD v2.2 §2.3.1.4) is a dataless BPSK(5) quasi-pilot at
5.115 Mcps with a 330-chip, 64.5 µs primary code, repeated 31 times per 2 ms with
no overlay, broadcast so a receiver can find Galileo with a search 31× smaller
than E5a-I's.

- **As an acquisition aid.** Its carrier Doppler transfers directly to E5a-I/Q
  (same carrier). Its code phase pins E5a-I only modulo 64.5 µs, leaving 31
  hypotheses, still a 330× reduction. At *half* the E5a-I/Q chip rate it cannot
  share their [`SignalGroup`](@ref)
  ([#151](https://github.com/JuliaGNSS/Tracking.jl/issues/151)); track it in its
  own group and hand over.
- **As a tracked signal** it yields the usual observables. One 64.5 µs block
  would run the loop at 15.5 kHz, so it integrates a whole 31-block, 2 ms code
  cycle by default: a 500 Hz loop, where the carrier loop pulls in at 45 Hz and
  narrows to 18 Hz.
- **It never provides navigation data**, for itself or for E5a.

Only part of the constellation transmits it (see the
[GSC list](https://www.gsc-europa.eu/galileo/services/galileo-open-service/quasipilot)),
and the ICD says it will eventually be discontinued.

## Known limitations

- **BeiDou GEO satellites (PRN 1-5, 59-63) never sync on B1I/B3I.** They carry
  no NH20 overlay and broadcast D2 at 500 sym/s instead of D1 at 50 sym/s, so
  the secondary-code detector has nothing to lock onto. They track and can be
  *ranged* on, but are not decoded: a D2 bit-edge search would need a per-PRN
  data rate, while `get_data_frequency` reports D1's for every PRN. See
  TrackingLoops' `src/beidou/b1i.jl`.
- **BeiDou B2b-I** correlates to an all-zero prompt on a clean replica, which
  takes the loop down through a NaN Doppler
  ([#234](https://github.com/JuliaGNSS/Tracking.jl/issues/234)); its defaults
  and sync detector are checked, but it is not usable end to end.
- **Long-overlay pilots** (`GPSL1C_P`, `BeiDouB1C_P`) need a full 18 s overlay
  period before the hard rotation sweep can lock
  ([#215](https://github.com/JuliaGNSS/Tracking.jl/issues/215)).
- **Secondary-code phase is re-searched at handoff** even when the acquisition
  stage already recovered it
  ([#120](https://github.com/JuliaGNSS/Tracking.jl/issues/120)).
- **One chip rate per [`SignalGroup`](@ref)**, which keeps E5a-QP out of an E5a
  group
  ([#151](https://github.com/JuliaGNSS/Tracking.jl/issues/151)).
- **Pilots without a secondary code** (GPS L2 CL, Galileo E5a-QP) keep the
  two-quadrant Costas discriminators, at half the four-quadrant range, as no
  sync gives them the prompt sign; secondary-coded pilots switch to
  four-quadrant once synced
  ([#230](https://github.com/JuliaGNSS/Tracking.jl/issues/230)).
