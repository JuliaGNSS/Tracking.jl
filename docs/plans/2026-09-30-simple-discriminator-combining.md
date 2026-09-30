# Simple Multi-Signal Discriminator Combining

A low-complexity variant of the design worked out on branch
`sc/tracking-multi-signal` (`docs/plans/2026-08-17-multi-signal-discriminator-combining.md`
there). The goal, the weights, the group-delay referral and the reasons for
combining at the discriminator level are that design's and are condensed here.
What differs is the timing: this variant combines only records that **end on the
same sample**, and keeps no state between driver records or chunks.

## Goal

Make multi-signal tracking *use* its extra signals. The
[2026-05-15 multi-signal design](2026-05-15-multi-signal-tracking.md) built the
structure — one `TrackedSat` carrying a tuple of `TrackedSignal`s, one shared
carrier downconvert, a correlator per signal — but `signals[1]` (the **driver**)
drove the loops and `signals[2:end]` (the **passengers**) only recovered their
own bits and C/N₀. A satellite tracked on a pilot/data pair paid for two
correlators and closed its loops on one.

With `discriminator_combining = true` on a `SignalGroup`, a passenger record that
ends on the same sample as a driver record is folded into that driver record's
PLL, FLL and DLL update. Off by default; a single-signal satellite, and any
satellite of a group with the flag off, is bit-identical to before.

## Why the code loop needs a group delay

Combining the *carrier* discriminators needs nothing from the caller. All signals
of a group sit on one carrier (a `SignalGroup` is single-band and single-chip-rate),
so after de-rotating onto the driver's phase frame each signal's `pll_disc` /
`fll_disc` measures the same carrier phase / frequency error.

Combining the *code* discriminators does not have that property. A combined DLL
drives the satellite-shared `code_phase` towards a weighted average of the
signals' code phases, and those differ by the satellite's payload group delay
differences. Uncorrected, the bias **moves** as the weights move, and a
downstream consumer applying the group-delay correction of one ranging signal
(`PositionVelocityTime.jl`) applies it to a phase that is no longer that
signal's. An L1C-P passenger can outweigh an L1 C/A driver by roughly 80:1, so a
1–3 ns difference lands as 0.3–0.9 m.

**Resolution:** each `TrackedSignal` carries a `group_delay` (a time, on a datum
the caller chooses), and a passenger's code discriminator has
`(delay(signals[1]) − delay(passenger)) · f_code` chips subtracted at the loop
update, at the code frequency in effect. `sat.code_phase` keeps meaning the
driver's code phase. The sign is fixed by the loop itself: a positive `dll_disc`
raises the code frequency and advances the replica, so `dll_disc` carries the
sign of `(true − replica)`, and a passenger less delayed than the driver (larger
code phase) reads `e_driver + δ`.

## Unknown vs. zero, and who decides

Two missing states, and conflating them is wrong in both directions:

  - **Known to be zero.** Every Galileo pair (E1B/E1C, E5aI/E5aQ, E5bI/E5bQ,
    E6B/E6C) is one composite out of one payload chain, and the broadcast group
    delays are cross-band only. The caller says `0.0s` and the code loop combines
    from the first integration.
  - **Not yet known.** GPS broadcasts a per-component ISC, BeiDou does for B1C and
    B2a — a statement that the components do not share a group delay. The value
    has to be decoded.

Telling them apart needs constellation knowledge and a parsed navigation message,
so **Tracking does neither**. Every slot starts at `nothing`; `nothing` on a
passenger withholds that passenger from the code loop only, `nothing` on the
driver withholds every passenger's code contribution, and the carrier loops
combine throughout. Guessing zero for an unknown difference is exactly the moving
bias above, and indistinguishable from a value the caller forgot. The two ICD
families reference their corrections differently (GPS: write `−ISC_x` on each
signal; BeiDou: `0.0s` on the pilot, `+ISC` on the data component); the manual's
**Group delay** section spells that out, and the BeiDou sign has still to be
checked against real data.

## Weighting

`w = P · N / G` per record: the nominal post-integration SNR over the
discriminator's noise gain. `P` is `GNSSSignals.get_relative_power(signal)` —
the ICD power split (L1C-P `0.75` / L1C-D `0.25`; E1B/E1C, E5aI/E5aQ, L5I/L5Q
`0.5` each) — and `N` the record's sample count. Every factor common to a
satellite's signals (noise density, sampling frequency, code amplitude) cancels
and is never formed.

The power is nominal rather than a measured `|P|²` because within one
constellation and band the ICD ratio is the stable quantity — elevation, block,
path loss and front-end gain cancel out of it — while a measured `|P|²` brings
the estimator's own noise into the loop gain. It is also a compile-time constant
and needs no C/N₀ estimate. The price: nominal weights cannot notice a component
that is not transmitted, so declaring one is a configuration error to fix at the
source (a pre-Block-III GPS satellite belongs in a C/A-only group). The one
record excluded automatically is a passenger's correlated with a pre-sync replica
on a secondary-coded signal, whose coherent sum partially cancels.

| Discriminator | `G` | Notes |
|:--- |:--- |:--- |
| Code, early-late (spacing `d` chips) | `d / 4` | Kaplan & Hegarty Table 5.6; the `d` dependence is the E/L noise correlation. |
| Code, VEML on BOC(1,1) | `1 / (2 (s_i + s_o)²)` | Same envelope model as the discriminator's chip calibration; tap-pair noise correlation not modelled, so it over-states VEML noise (the conservative direction). `1/32` for the default taps. |
| Carrier phase | signal-independent | Cancels in the normalized mean; weight is the bare `P · N`. |
| Carrier frequency (FLL) | `∝ 1/T²` | `fll_disc` divides a phase difference by `2π·T`, so the weight is `P · N³`. A record with no previous prompt has weight zero. |

A wrong gain costs combining efficiency, never bias — each discriminator is
calibrated in its own units. There is no `AbstractCorrelator` fallback for
`dll_disc_noise_gain`: a correlator with its own `dll_disc` must state its gain.

## Discriminator-level, not prompt-level

Coherently summing de-rotated prompts would be the ML phase estimator for a
common carrier, and is wrong here: `pll_disc` and `fll_disc` are Costas
discriminators, blind to the ±1 data sign, which is what lets a data component be
averaged with a pilot; their prompts would cancel on every bit flip. De-rotation
is needed for the PLL only — `dll_disc` is noncoherent and `fll_disc` forms
`conj(p_prev) · p_cur`, in which a constant rotation cancels.

The code loop has a stronger reason of its own: `dll_disc` is built from tap
magnitudes, and each method's chip calibration is tied to its own signal's
correlation shape. Taps summed across a BOC and a BPSK signal have no single
slope to divide by, so the result would not be calibrated in chips — the property
that makes two signals' code errors averageable. The ICD power split does not
force the choice either way: coherent maximal-ratio combining would weight by
`sqrt(P · N)`, the same input one square root away.

## Weighted **mean**, not sum

The combined value per loop is `Σ wᵢ dᵢ / Σ wᵢ` over the driver's record and
the coincident passenger records, so the loop gain does not change with the
number of contributors and only the measurement noise does. Where no passenger
contributed to a loop, the loop reads the driver's own discriminator as it is —
not `(w·d)/w`, which is not `d` bit for bit — so a driver record with nothing
coincident closes the loops exactly as the driver-only path does. Zero total
weight reads as zero.

## Timing: coincident records only, and the equal-length assumption

The loop filters fire per driver record. For each one, a passenger record
contributes **iff its `sample_index` equals the driver record's**. The fold walks
each passenger's records in order with one cursor: records ending before the
driver record are applied to their signal (prompt, post-correlation filter, C/N₀,
bit buffer, `last_fully_integrated_*`) and dropped from the loops; a record ending
on the same sample is applied and combined; after the chunk's last driver record
the rest are drained, applied only. The combination for one driver record is a
local value: nothing is parked on the per-satellite state, no sample-index
windows are formed, and nothing crosses a chunk boundary.

This rests on one assumption, documented and not enforced here: **the driver
integrates at least as long as every passenger, and for full aiding equally
long.** GNSSReceiver.jl checks it before enabling combining. It holds by default:
every signal starts at one primary code block (Galileo E5a-QP, at 31, is tracked
alone), and every intra-band pilot/data pair GNSSSignals defines shares one
primary code period between its components, so their records end on the same
sample. All signals of a satellite share one sub-step grid in the correlate
phase, so equal boundaries give equal `sample_index` values; an external
producer has to keep a satellite's signals on one per-chunk origin for the same
reason, and a mismatch means only that a passenger record is not combined.

What is lost, stated plainly:

  - **A passenger integrating shorter than the driver contributes one record in
    `k`** — the one that ends with the driver's — at its own (`1/k`) weight. A
    10 ms L1C-P driver with a 1 ms C/A passenger uses one C/A record in ten.
    The window-based design recovered the other `k − 1` by accumulating every
    passenger record inside the driver's window, and carried them across chunk
    boundaries to do so; this variant gives that up.
  - **A caller lengthening the driver alone** (`set_preferred_num_code_blocks_to_integrate!`)
    puts its passengers in exactly that position, and loses `(k − 1)/k` of
    their aiding until the passengers are lengthened to match.
  - **A passenger integrating longer than the driver** breaks the assumption and
    aids only the driver records it happens to end with.

In exchange there is no accumulator to thread, park, drop on a `vt_on`
transition, drop on a group-delay change or drop on `reset_loop_filters!`, and no
question of which update a record belongs to.

## Vector tracking

`VectorPLLAndDLL` runs through the same fold. Its per-satellite accumulators
become one `(count, sum)` pair per signal and per loop, filled under `vt_on` by
every record of every signal with the **raw** `dll_disc` / `fll_disc` — no
group-delay referral (under `vt_on` the code loop is the navigation filter's),
whatever the combining flag says, coincident or not. The pre-sync-invalidated
passenger record skips both of its slots together. The loops still closed here
follow `vt_on`: in the scalar fallback all three combine exactly as under the
conventional estimator; under vector closure only the carrier phase loop does,
and the combined code and frequency values the fold forms are ignored.

## API

```julia
# Opt-in, per group.
TrackState(; signals = (galileo_e1 = (GalileoE1C(), GalileoE1B()),), discriminator_combining = true)
SignalGroup((GalileoE1C(), GalileoE1B()); discriminator_combining = true)

# Per-signal group delay, a time; every slot starts at `nothing`.
set_group_delay!(ts, :galileo_e1, 11, 1, 0.0s)
set_group_delay!(ts, :galileo_e1, 11, GalileoE1B, 0.0s)
get_group_delay(ts, :galileo_e1, 11, GalileoE1B)

# Vector tracking: one raw measurement per signal.
mean_code_discr(ts, :galileo_e1, 11, GalileoE1B)
```

`group_delay` is the last positional field of `TrackedSignal`
(`Maybe{typeof(1.0s)}`), so the kwarg-update constructor takes it wrapped in
`Some` — `nothing` is a legal value. The value is supplied rather than derived
because deriving it means knowing which constellation broadcasts which term, on
what datum and with which sign, and then parsing a navigation message; putting
that in the receiver also puts the bias and the downstream group-delay
correction in one place.

## Scope

  - **Cross-band combining** (L1 + L5 on one satellite) is out: a `SignalGroup`
    is single-band, and a cross-band combination needs an ionospheric divergence
    model, not a constant.
  - **Replica offsets** are out: the group delay is applied to the
    discriminator, not by shifting each passenger's replica, which at ISC
    magnitudes (≤ ~0.01 chips) costs negligible correlation power.
  - **Unequal integration lengths** are supported only as described above; the
    window-based accumulation is the documented alternative if a receiver needs
    it.

## Validation on real data

The measurements below were taken with the window-based implementation on
`sc/tracking-multi-signal`, not re-run for this variant. They carry over **for
the equal-length pairs only**, where this variant forms the same combination by
construction: with driver and passenger integrating equally long, every passenger
record ends on a driver record's sample and completes in the same chunk, so the
window-based fold never had a record pending or a window holding more than one
passenger record, and both folds combine the same two records into every update.

Two static captures (TEX-CUP's first ~600 s; the Fraunhofer Spirent scenario),
both configurations over the same samples from identical handoffs, differing only
in the combining flag. Noise is the std of the second difference of the Doppler
series over √6, reported as the off/on ratio (`> 1`: combining reduced it). Every
satellite carries a lock check — the two configurations must agree on the mean
Doppler. All three pairs split power evenly and use one correlator type, so the
combined discriminator is a plain average of two measurements: **prediction
√2 ≈ 1.41.**

| capture | group | n | carrier | code |
|:--- |:--- |:--- |:--- |:--- |
| TEX-CUP | Galileo E1 (E1C+E1B) | 7 | 1.39 | 1.52 |
| TEX-CUP | Galileo E5a (E5aQ+E5aI) | 5 | 1.39 | 1.40 |
| TEX-CUP | GPS L5 (L5Q+L5I) | 5 | 1.40 | 1.04 → 1.41 |
| Fraunhofer | Galileo E1 (E1C+E1B) | 8 | 1.23 | 1.48 |
| Fraunhofer | Galileo E5a (E5aQ+E5aI) | 8 | 1.34 | 1.43 |
| Fraunhofer | GPS L5 (L5Q+L5I) | 8 | 1.27 | 1.02 → 1.40 |

41 locked satellites, every mean-Doppler difference within 0.11 Hz. The GPS L5
code column is `1.04`/`1.02` with the passenger held out of the code loop and
`1.41`/`1.40` once it was let in, while the carrier ratio did not move — the gate
is visible in the data. Those two L5 code figures were measured under an earlier
revision in which Tracking itself withheld the code loop until an ISC was
decoded; the mechanism is the one an unknown `group_delay` triggers here, but
the numbers were not re-taken with the caller supplying the value. The carrier figures sit a little below
√2, consistent with a common noise floor (reference oscillator, residual
dynamics) the combination cannot reduce.

**Not applicable to this variant:**

  - The reference plan's **GPS L1 C/A + L1C** results (a 1 ms C/A driver with
    10 ms L1C passengers, and the three-component sweeps on TEX-CUP PRN 4). They
    mix code periods, so the window-based fold combined records this variant
    drops; the numbers say nothing about this fold.
  - The **unequal-integration-length** runs (a pilot driver at 2× its data
    passenger). There the window-based fold combined both of the passenger's
    records per driver window; this one combines one of them, so the reported
    "no change in the ratios" does not transfer.

Re-measuring this variant on those configurations, and on the equal-length pairs
to confirm the by-construction argument end to end, is follow-up work.
