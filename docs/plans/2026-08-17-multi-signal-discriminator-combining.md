# Multi-Signal Discriminator Combining

## Goal

Make multi-signal tracking *use* its extra signals. The
[2026-05-15 multi-signal design](2026-05-15-multi-signal-tracking.md) built the
structure — one `TrackedSat` carrying a tuple of `TrackedSignal`s, one shared
carrier downconvert, a correlator per signal — but explicitly listed
"per-signal Doppler estimation" as a non-goal: `signals[1]` drove the loops and
`signals[2:end]` only recovered their own bits and C/N₀. So a satellite tracked
on L1C-P + L1C-D + L1 C/A paid for three correlators and closed its loops on
one, discarding the pilot's 75 % power share, its 10× longer coherent
integration, and its ~3× steeper BOC S-curve unless it happened to be the signal
placed first.

This step folds **every** signal's discriminator output into the loop update.

## Why this needs broadcast correction parameters

Combining the *carrier* discriminators is free. All signals in a group sit on one
carrier at one frequency (a `SignalGroup` enforces same band and same chip rate),
so after de-rotating onto the driver's phase frame each signal's `pll_disc` /
`fll_disc` measures the *same* carrier phase / frequency error. Nothing has to be
decoded.

Combining the *code* discriminators does not have that property. A combined DLL
drives the satellite-shared `code_phase` towards a weighted average of the
signals' code phases, and those differ by the satellite's **differential payload
group delay** — which is exactly what the broadcast inter-signal corrections
(ISCs) describe.

Two things make this more than a rounding error:

1. **The bias moves.** The weights track C/N₀, so an uncorrected differential
   makes the shared `code_phase` wander between the signals' code phases as
   conditions change — worse than a constant bias, which a receiver could
   calibrate.
2. **The downstream correction is applied to a phase that is no longer that
   signal's.** `PositionVelocityTime.jl` reads the shared
   `TrackedSat.code_phase` and applies the ISC of whichever ranging signal it is
   handed (`correct_by_group_delay` dispatches on the ranging signal). If the
   shared phase is a mixture, that per-signal correction is simply wrong.

The magnitudes rule out ignoring it. An L1C-P passenger outweighs an L1 C/A
driver by roughly 80:1 — comparable received power (−158.25 vs −158.5 dBW), 10×
the integration, ~8× the discriminator gain — so the shared code phase becomes
essentially the *passenger's*, and a 1–3 ns differential lands as 0.3–0.9 m of
bias. That is larger than the jitter the combining just bought.

**Resolution:** refer every passenger's code discriminator to the *driver's* code
phase by subtracting the differential group delay. `sat.code_phase` then keeps
meaning "the driver signal's code phase", which is exactly what the existing
downstream group-delay correction already assumes — no change needed in
`PositionVelocityTime.jl`.

### Sign convention

Fixed by what the pseudorange does with it. A receiver reduces a signal's
transmit-time estimate by `−T_GD + ISC_signal` (IS-GPS-200 §20.3.3.3.3.2, and
IS-GPS-705/800 for the ISCs), so two signals of one satellite measuring the same
geometric range satisfy `t_raw(a) − t_raw(b) = ISC(a) − ISC(b)`, with `t_raw`
growing with measured code phase. A signal whose ISC exceeds the driver's
therefore sits at a **larger** code phase, by `(ISC − ISC_driver) · f_code`
chips — the amount to subtract from its discriminator.

Independently confirmed by the loop's own stability: a positive `dll_disc`
raises the code frequency, which advances the replica phase, so `dll_disc` must
carry the sign of `(true phase − replica phase)`. The passenger's raw
discriminator therefore reads `e_driver + δ`, and `δ` comes off.

## Unknown vs. zero, and who decides

A passenger's differential group delay has two very different missing states, and
conflating them is wrong in both directions:

  - **Known to be zero.** Every Galileo pair (E1B/E1C, E5aI/E5aQ, E5bI/E5bQ,
    E6B/E6C) is one CBOC or QPSK composite out of one payload chain; the
    broadcast BGDs are cross-band (E1↔E5a/E5b) and no intra-band term exists.
    Waiting for one would disable code combining **forever**.
  - **Not yet known.** IS-GPS broadcasts a *per-component* ISC (`ISC_L1CA`,
    `ISC_L1CD`, `ISC_L1CP`, `ISC_L2C`, `ISC_L5I5`, `ISC_L5Q5`), which is
    precisely a statement that the components are not assumed to share a group
    delay; BeiDou does the same for its B1C and B2a pairs (`ISC_B1Cd`,
    `ISC_B2ad`, in B-CNAV1 / B-CNAV2). The real value has to be decoded.

Telling those apart needs constellation knowledge and, in the second case, a
parsed navigation message — so **Tracking does neither**. Every signal's
`differential_group_delay` starts at `nothing`, meaning *unknown*, and the caller supplies what
it knows. `nothing` keeps the passenger out of the **code** loop only, leaving its
carrier contribution untouched, and it joins the moment a bias is set. A caller
that knows the value is zero by construction says `0.0` and combines from the
first integration.

This is the safe direction rather than the convenient one: guessing zero for an
unknown differential is exactly the moving bias this section exists to prevent,
and it is indistinguishable from a bias the caller meant to supply and forgot.

How long the wait is depends on the message, not the signal: `ISC_L1CD` /
`ISC_L1CP` ride CNAV-2 subframe 2 next to the ephemeris (every 18 s), while
`ISC_L1CA` sits in subframe 3 page 1 and the CNAV `ISC_L5I5` / `ISC_L5Q5` /
`ISC_L2C` only in message type 30 — both far rarer, which is why
GNSSDecoder.jl's `is_decoding_completed_for_positioning` does not wait for them
either, and advises treating a missing one as zero for *positioning*. LNAV
carries no ISC at all, so a legacy-only receiver has nothing to derive an
L1 C/A + L1C bias from.

## Weighting

`w = P · N / G` per record: nominal post-integration SNR over the
discriminator's noise gain.

**Power from the ICD split.** `P` is `GNSSSignals.get_relative_power(signal)` —
a dimensionless linear power giving the component's share of its composite, with
composites on the same constellation and band scaled against one another (L1C-P
`0.75` / L1C-D `0.25`; E1B/E1C, E5aI/E5aQ, L5I/L5Q, L2CM/L2CL `0.5` each;
GPSL1CA `0.708`, scaled against the L1C composite). `P · N` is the
post-integration SNR up to factors common to a satellite's signals, and each of
those — front-end noise density, sampling frequency, code amplitude — cancels
and is never formed. A `SignalGroup` is single-band by construction and a
`TrackedSat` is one satellite, so the "same constellation and band" precondition
on those ratios always holds where they are used.

Chosen over weighting by the record's own measured `|P|²` because within one
constellation and band the ICD ratio is the *stable* quantity: elevation,
satellite block, free-space loss, antenna gain and front-end gain are common to
the components and cancel out of it. A measured `|P|²` instead carries the
estimator's own noise into the loop gain, biased upward by `σ²/N` and so
compressing the weights towards equal exactly where correct weighting matters
most. Nominal power also folds to a compile-time constant and needs no C/N₀
estimate, so it works with `NoCN0Estimator`.

The cost is that nominal weights cannot notice a component that is not there,
which a measured weight muted for free. That job moves to an explicit
presence gate — see below.

**`N` handles unequal integration lengths by construction — in one direction.** A
passenger integrating `1/k` of the driver's length completes `k` records per
driver window, each with `1/k` of the weight, so its *total* weight over the
window equals one driver-length record: the combination is energy-fair for any
passenger whose records are no *longer* than the driver's. The FLL is weighted
`P · N³` instead, because `fll_disc` divides an inter-prompt phase difference by
`2π·T`: one 20 ms measurement genuinely beats five averaged 4 ms ones, and the
extra `N²` is what says so.

The reverse case is **not** covered, and the driver must therefore be the
longest-integrating signal in the group. The accumulator is zeroed at every
driver record, so a passenger reporting once per `k` driver records reaches the
loop in one update out of `k` and carries no weight in the other `k − 1`;
its influence over the window is ~`1/k` of its energy share. This costs most of
the gain rather than a little of it: on the static Spirent L1 capture, over nine
satellites, a `(GPSL1CA(), GPSL1C_P(), GPSL1C_D())` group (1 ms driver, 10 ms
passengers) reduces carrier-Doppler noise by a median 1.05x where its own weights
predict 1.55x, and the identical group with the driver raised to 10 code blocks
reaches a median 1.40x carrier and 4.10x code. It
is dormant for every shipped intra-band pilot/data pair (GPS
L5, GPS L1C, Galileo E1, Galileo E5a all share one primary code period between
their components; GPS L2C's pilot is the longer one) and for the usual tuning of
lengthening a *pilot driver*, and it bites only when a group mixes code periods
(C/A with L1C) or a *data passenger* is integrated past its pilot driver.

Two ways to close the gap, neither implemented: weight each passenger record by
its **sample overlap** with the driver's window and let one long record
contribute to every window it spans (rate-preserving, but reuses one noise
realisation across `k` updates and raises the applied lag from ~`T/2` to ~`T`),
or decimate the driver to the slowest signal's rate (exact, but gives up the fast
update rate and rescales the `1/N` carrier bandwidth). The second is available
today as `set_preferred_num_code_blocks_to_integrate!` on the driver.

**No presence check.** A record correlated with a pre-sync replica on a
secondary-coded signal is excluded (its coherent sum is partially cancelled by the
missing wipe-off, so it measures nothing), but there is no check that a component
is being received at all. A signal is in a `TrackedSat` only because the caller
declared it, so "declared but not transmitted" is a configuration error to fix at
the source, not something to detect per record in the hot path. The consequence
worth documenting: a pre-Block-III GPS satellite carries no L1C, so it must go in a
C/A-only group rather than an `(L1C_P, L1C_D, L1CA)` one, where two noise channels
at `0.75 + 0.25` would outvote C/A at `0.708`.

**Noise gains**, all in the convention `σ² ≈ G / SNR` with `SNR ≡ |P|²/σ²` and
per-tap envelope noise `σ²/2`:

| Discriminator | `G` | Notes |
|:--- |:--- |:--- |
| Code, early-late (spacing `d` chips) | `d / 4` | Kaplan & Hegarty Table 5.6. The `d` dependence is entirely the E/L noise correlation `ρ ≈ 1 − d`. |
| Code, VEML on BOC(1,1) | `1 / (2 k² (e_i + e_o)²)` | Same envelope model as the discriminator's own chip calibration; `k` is `_veml_discriminator_slope`. Tap-pair noise correlation **not** modelled, so it over-states VEML noise — the conservative direction. |
| Carrier phase | signal-independent (`≈ 1/2`) | Cancels in the normalized mean; weight is the bare SNR. |
| Carrier frequency (FLL) | `∝ 1/T²` | `fll_disc` divides an inter-prompt phase difference by `2π·T`, so accuracy scales with `T`. Weight carries an extra `N²` — a 10 ms passenger is worth 100× a 1 ms one at equal SNR. |

For the default taps this puts VEML about 8× above a 1-chip early-late gain,
matching the ~3× steeper BOC(1,1) main peak (variance ~9×). A wrong gain costs
*combining efficiency*, never bias — each discriminator is individually
calibrated, so any positive weights give a consistent estimate. That is what
licenses first-order models and an `AbstractCorrelator` fallback instead of a
`MethodError` for custom correlators.

## Discriminator-level, not prompt-level

The obvious alternative — coherently sum the de-rotated prompts and take one
`atan` — is the ML phase estimator for a common carrier and would weight by
amplitude automatically. It is **wrong here**: `pll_disc` and `fll_disc` are
Costas discriminators, blind to the unknown ±1 navigation-data sign, which is
what lets a data component's discriminator be averaged with a pilot's. Their
*prompts* would partially cancel on every bit flip.

Combining at the discriminator level also handles differing cadences naturally,
and needs de-rotation for the PLL only: `dll_disc` is noncoherent (tap
magnitudes) and `fll_disc` forms `conj(p_prev) · p_cur`, in which a constant
unit rotation cancels. The de-rotation is not hypothetical — GPS L1 C/A's
carrier phase offset is −π/2 while L1C-D/P sit at 0, so C/A really is in
quadrature with an L1C driver.

The data sign is the argument for the two **carrier** loops. The **code** loop
has its own, and it is the stronger of the two, because it survives the data sign
being resolved: `dll_disc` is not a function of a prompt at all. It is built from
the tap *magnitudes* (`|VE|`, `|E|`, `|L|`, `|VL|`), and each method's chip
calibration is tied to its own signal's correlation shape — the `(2 - d) / 2`
factor of the early-late form assumes a triangular BPSK autocorrelation, the VEML
form divides by the S-curve slope of a sine-BOC(1,1) envelope
(`_veml_discriminator_slope`). Taps summed across a BOC signal and a BPSK one
produce an envelope with no single slope to divide by, so the combined statistic
would not be calibrated in chips — which is the one property that makes two
signals' code errors averageable in the first place. There is no prompt-level
code combination to choose against, whatever the data sign does.

Nor is the *weighting* what forces the choice, which is worth stating because it
is the natural suspicion: the ICD power split is not something only a
discriminator-level combination can use. Coherent prompt combining is
maximal-ratio combining, whose optimum weight is amplitude-over-noise —
`sqrt(P · N)` where the discriminator form uses `P · N`, i.e. the same
`get_relative_power` input one square root away. Either level can weight by the
power split. What rules prompt summing out is the data sign for the carrier
loops, and the calibration for the code loop; nothing about the weights.

## Weighted **mean**, not sum

The combined value is `(Σ wᵢ dᵢ) / (Σ wᵢ)`. This is what makes combining safe to
default to: the result stays calibrated in the discriminator's own units however
many signals contributed at a given epoch, so the loop *gain* never changes and
only the measurement *noise* varies. A plain sum would multiply the loop gain by
the contributor count and destabilise the filter every time a long-period
passenger reported.

The `iszero(passenger_weight)` short-circuit returns the driver's discriminator
**bit-identically**, so a single-signal satellite — and a multi-signal one in a
chunk where no passenger completed — is numerically untouched.

## Timing: exact windows, and why the accumulator must outlive the chunk

The loop filter fires per driver record. The fold therefore walks the driver's
and the passengers' records in `sample_index` order: before closing the loops on
driver record `k`, every passenger is advanced through the records that completed
at or before `k`'s end sample, so each loop update sees exactly the measurements
belonging to its own window.

Records completing after the last driver record of a chunk stay pending in a
`DiscriminatorAccumulator` carried in the per-satellite estimator state. This is
a **requirement, not an optimisation**: `track!`'s default chunk is the *shortest*
primary code period across all tracked signals, so with a 10 ms L1C-P driver and
a 1 ms C/A passenger, nine chunks in ten contain C/A records and no driver record
at all. A per-chunk accumulator would discard nine measurements in ten.

Two consequences worth stating:

- A passenger's `found_before_chunk` sync flag is captured once per chunk and
  threaded, not re-read per window. Re-reading it would clear
  `correlated_pre_sync` for records in later windows that were nonetheless
  correlated with pre-sync replicas.
- A long-period passenger's measurement represents the average error over its
  own (longer) window, so it enters the loop with up to `T_passenger/2` of lag.
  At a 1 Hz DLL bandwidth (~160 ms time constant) a 5 ms lag is negligible; it
  is another reason the existing advice to put the longest-integration pilot
  first is good.

## API

```julia
# Opt-in: the caller enables it, having checked the driver-ordering precondition.
ConventionalAssistedPLLAndDLL(signal_combining = true)

# Per-passenger differential group delay, relative to the driver, addressed like
# `set_preferred_num_code_blocks_to_integrate!`. The driver's own is 0.0.
set_differential_group_delay!(ts, :gps_l1, 7, GPSL1CA, -1.0e-9)
get_differential_group_delay(sat, GPSL1CA)
```

`differential_group_delay` is a `Maybe{Float64}` field on `TrackedSignal`, so the kwarg-update
constructor takes it wrapped in `Some` — `nothing` is a legal value ("unknown"),
and the usual `isnothing(x) ? keep : set` test cannot express "clear it"
otherwise.

### Why the bias is supplied rather than derived

Deriving it means knowing which constellation broadcasts which term, on what
datum, with which sign — and then parsing a navigation message to get it. None of
that is signal processing, and a table of it inside Tracking would have to grow
with every signal GNSSSignals adds and be re-verified each time. So the value
crosses the boundary already reduced to the one quantity the code loop uses: a
per-passenger differential in seconds, positive when the passenger sits at the
larger code phase.

That keeps Tracking free of GNSSDecoder.jl entirely, and puts the bias and the
*downstream* group-delay correction (`PositionVelocityTime.jl`'s, applied to the
same shared `code_phase`) in one place — the receiver — where they can be held
consistent. Under the earlier design Tracking set one end and PVT applied the
other, with the invariant asserted in docstrings and enforced nowhere.

#### Sign conventions the caller has to reconcile

The two ICD families disagree, so a caller forwarding decoded fields must not
treat them alike. Both correct SV time as `t = t_SV − Δt_SV_x`, and per-signal:

```
IS-GPS (705/800):      Δt_SV_x = Δt_SV − T_GD       + ISC_x
BeiDou (B1C/B2a §7.6): Δt_SV_x = Δt_SV − T_GD_pilot − ISC_data
```

Writing each as `Δt_SV − G_x` for a hardware delay `G_x` gives
`G = T_GD − ISC_x` for GPS but `G = T_GD_pilot (+ ISC_data)` for BeiDou. Two
signals of one satellite measuring the same range satisfy
`t_raw(a) − t_raw(b) = G_b − G_a`, and this package's bias is defined so that
difference equals `bias(a) − bias(b)` — hence `bias = −G`. For GPS that is
`ISC_x − T_GD`, and `T_GD` is common to the satellite's signals and cancels, so
**a plain difference of two GPS ISCs is already correct**. For BeiDou nothing
cancels the sign, so `bias(B1C_D) = −ISC_B1Cd` relative to the pilot (and
likewise `−ISC_B2ad` for B2a).

**The BeiDou sign is unconfirmed against real data.** It follows from the ICDs and
agrees with the loop-stability argument above, but the Doppler-noise ratios below
are structurally blind to it: a sign error there is a *code-phase bias*, so
catching it needs a differential-range or double-difference check against a known
baseline. Treat a BeiDou pair's combined code loop as provisional until such a
check exists.

## Scope

Out of scope, and why:

- **Cross-band combining (L1 + L5 on one satellite).** A `SignalGroup` is
  single-band by construction, so cross-band signals of one satellite are
  separate `TrackedSat`s today. Combining them needs more than an ISC: the
  ionospheric delay differs between bands *and* drifts, so it needs a
  divergence model or a per-band replica offset, not a constant. This step's
  machinery (per-signal group delay, weighted accumulator) is the prerequisite.
- **`VectorPLLAndDLL`.** Its discriminators are consumed by an external
  navigation filter with its own semantics, so it keeps the driver-only fold
  (its own `_process_estimator_driver_signal` method is untouched). Additive to
  enable later.
- **Replica offsets.** The differential group delay is applied to the
  discriminator, not by shifting each passenger's replica. At ISC magnitudes
  (≤ ~0.01 chips) the correlation-power loss from an off-centre passenger
  replica is negligible; a cross-band step would need real replica offsets.

## Validation on real data

Two independent captures, both **static** (TEX-CUP's first ~600 s — `ground_truth.log`
gives mean speed 0.004 m/s and max 0.007 m/s over t < 300 s, with sustained motion
only from t ≈ 1759 s; the Fraunhofer Spirent scenario is static by construction). A
static antenna is what makes the metric trustworthy: the true carrier Doppler is
smooth, so high-frequency content in the Doppler estimate is receiver noise.

Method: both configurations run over the **same** samples in one pass from identical
acquisition handoffs, differing only in `signal_combining`. Noise is the std of the
second difference of the Doppler series divided by √6 (a strong high-pass; the first
difference and a moving-average detrend agree with it within a few percent). Reported
as the off/on ratio, so `> 1` means combining reduced the noise. Every satellite
carries a **lock check**: two locked configurations must agree on the *mean* Doppler,
so a non-trivial `Δdopp` marks a mistracking run whose ratio is meaningless. This
check is what caught every false result below.

GPS and Galileo are never in one group — a satellite is one constellation.

### Headline: intra-band pilot/data pairs

Both Galileo pairs and the GPS L5 pair split power evenly
(`get_relative_power == 0.5` each) and use the same correlator, so at equal
integration lengths the weights are equal and the combined discriminator is a plain
average of two independent measurements. **Prediction: √2 ≈ 1.41.**

| capture | group | n | carrier | code |
|:--- |:--- |:--- |:--- |:--- |
| TEX-CUP | Galileo E1 (E1C+E1B) | 7 | 1.39 | 1.52 |
| TEX-CUP | Galileo E5a (E5aQ+E5aI) | 5 | 1.39 | 1.40 |
| TEX-CUP | GPS L5 (L5Q+L5I) | 5 | 1.40 | 1.04 → **1.41** |
| Fraunhofer | Galileo E1 (E1C+E1B) | 8 | 1.23 | 1.48 |
| Fraunhofer | Galileo E5a (E5aQ+E5aI) | 8 | 1.34 | 1.43 |
| Fraunhofer | GPS L5 (L5Q+L5I) | 8 | 1.27 | 1.02 → **1.40** |

41 locked satellites, every `Δdopp` within 0.11 Hz. The code figures match the √2
prediction; the carrier figures sit a little below it, consistent with a common noise
floor the combination cannot reduce (reference oscillator, and whatever dynamics the
high-pass leaves in).

### The code-loop gate is visible in the data

Measured under the earlier design, where Tracking itself withheld the code loop
until an ISC was decoded: on GPS L5 the code-Doppler ratio sat at 1.04x while the
gate held and reached 1.41x once it was lifted, with the carrier ratio unchanged
either way — on both captures independently. The *mechanism* is unchanged (an
unset `differential_group_delay` still holds the passenger's code contribution back while its
carrier contribution combines), but the numbers measured a flag this package no
longer owns, so they are recorded as history rather than carried forward as a
claim about the current API. Re-measuring needs a run with the receiver supplying
the bias.

### Unequal integration lengths

Doubling the pilot's coherent integration (E1 8/4 ms, E5a 2/1 ms) changed the ratios
by nothing — carrier 1.39, code 1.47 against 1.33 / 1.47 at equal lengths. That is
what weighting by `P · N` predicts: the passenger's two shorter records carry,
together, the weight of the driver's one long one.

**A 5× pilot integration (E1C at 20 ms) does not track at all** on this capture, and
that is a coherence limit, not a combining problem: a coherent integration of `T`
tolerates a residual frequency error of only about `1/(2T)` = 25 Hz, while the
acquisition hands over a Doppler on a ~333 Hz grid, and the `1/N`-scaled 0.9 Hz
carrier bandwidth cannot pull that in before the 20 ms integrations cancel. Symptoms:
E1C C/N₀ collapses to 12–17 dBHz and `std(pll_disc) ≈ 0.9 rad` (i.e. uniform over
the discriminator's range — no lock). To test long integrations, converge first at one
block and lengthen afterwards.

### Two anomalies, both resolved

**PRN 31's carrier "getting worse" (ratio 0.161) was an unlocked receiver, not
combining.** It appeared only in the 20 ms-pilot configuration described above. In
the static window at the default 4 ms the same satellite gives carrier 1.32 / code
1.54 with `Δdopp = 0.000 Hz`, and 1.34 / 1.51 in the all-systems run — no anomaly,
across several independent runs. The diagnostic that settles it is the mean-Doppler
check: in the 20 ms run the two configurations differed by **132 Hz**, so at least
one was mistracking and the noise ratio was meaningless. That check is now part of
every result.

**PRN 33's E1C C/N₀ collapse is a C/N₀-*estimator* artifact, not a tracking
problem.** On the same samples: NWPR (the default) reports 12.8 dBHz while
`MomentsCN0Estimator` reports **36.5 dBHz**, and its prompt magnitude (114 600)
matches its own E1B component's (118 600). Every other satellite agrees between the
two estimators within ~1 dB. NWPR sums prompts coherently over the signal's
bit / secondary-code window, so for a pilot it depends on the recovered CS25 phase;
E1C on PRN 33 also takes far longer than its band-mates to declare sync. So the fault
lies in the secondary-code phase handed to the estimator, not in the loops — and note
that under **nominal** weighting the C/N₀ estimator is not in the loop at all, so this
class of bug can no longer influence combining (it could have under a measured
`|P|²` weight). Worth a separate look at the CS25 detector; out of scope here.

### GPS L1's three components

L1C exists on far fewer satellites than one might assume, and that dominates what can
be measured:

- **TEX-CUP (2019-05-09):** L1C-D/L1C-P acquire on **PRN 4 only** — GPS III SV01, the
  single Block III satellite broadcasting L1C at the time. C/A acquires on 14.
- **Fraunhofer (2014 Spirent GSS8000):** L1C is **not simulated** — no detections with
  the correct spectral orientation. (An apparent PRN 21 detection under the *wrong*
  orientation was a false alarm.)

On TEX-CUP PRN 4, sweeping the combinations (all seeded from L1C-D, whose 10230-chip
code is the frame the shared code phase must be in when an L1C component is present):

| combination (driver first) | carrier | code |
|:--- |:--- |:--- |
| L1C-P + L1C-D | **1.19** | **1.15** |
| L1C-P + L1CA | 0.81 | 1.08 |
| L1C-D + L1CA | 0.92 | 1.23 |
| L1C-P + L1C-D + L1CA | 0.91 | 1.21 |

**Open item:** the two L1C components combine well, but adding C/A degrades the
carrier loop while still improving the code loop. Two candidate mechanisms were
tested and **ruled out**:

1. *A mismodelled C/A↔L1C carrier phase.* GNSSSignals puts C/A at −π/2 and L1C at 0,
   so combining de-rotates C/A by a quarter turn; the measured relative phase with
   L1C driving is ≈ 0 rad, suggesting they are co-phased. Overriding
   `get_carrier_phase_offset(GPSL1CA) = 0.0` moved the ratio only from 0.811 to
   0.854 — not the cause.
2. *A BOC side-peak lock on the L1C VEML discriminator* (which would put the shared
   code phase ~0.5 chip off C/A's broad BPSK peak). The converged code phase agrees
   to **0.005 chips** across L1C-driven, C/A-driven and cross-seeded runs — no side
   lock.

A confound to clear first: this group is **very** sensitive to the acquisition seed,
because an L1C driver's default carrier bandwidth is 1.8 Hz (sized from its 10 ms
primary period). A 7.65 Hz difference in the handover Doppler moves the *driver-only*
carrier noise from 0.20 Hz to 3.01 Hz, and C/A's own C/N₀ inside the group reads
33.1 dBHz under one seed and 42.2 dBHz under another at the same converged code
phase. Until the L1 handover is made robust (finer acquisition, or a wider initial
bandwidth that narrows after lock), the L1C numbers above are the least trustworthy
in this document.

### Capture notes worth keeping

- **TEX-CUP `ntlab.bin`:** 4 RF channels bit-interleaved, one byte per sampling
  instant, 2-bit sign/magnitude MSB-first, real IF, 79.5 MHz. ch0 = L1/E1 at
  IF −14.58 MHz, ch2 = L5/E5a at IF −23.55 MHz. `max_meas = 3`.
- **Fraunhofer Flexiband L125 III-1b:** one file per band, `Complex{Int8}` I/Q
  byte-interleaved, IF = 0, L1 at 20 MHz and L5 at 40 MHz, 212.5 s.
  `max_meas = 15` (4-bit values stored as the odd integers −15…+15).
  **Neither band is to be conjugated** with the plain
  `Complex{Int16}(raw[2k-1], raw[2k])` reader used here — settled by tracking a
  single-signal GPS L1 C/A group: without conjugation all 8 satellites hold
  (drift < 12 Hz over 10 s, C/N₀ 43–47.5 dBHz), with it 4 of 8 diverge (drift up to
  10 kHz, C/N₀ 15–21 dBHz). L5 likewise: 8/8 hold at C/N₀ 49.5–52.1 dBHz without
  conjugation, 8/8 diverge with it.
- **L5 and E5a acquisition must integrate coherently over a whole number of secondary-code
  periods** (NH10 → 10 for L5I, CS20 → 20 for E5aI). Ignoring this still *detects* the
  satellites but hands over a Doppler poor enough that the loops never lock — which
  surfaces much later as unexplained mistracking, not as a failed acquisition.

## Commit staging

1. `feat(discriminators): de-rotated PLL discriminator and noise-gain models`
2. `feat(track): per-signal differential group delay on TrackedSignal`
3. `feat(pll): combine every signal's discriminators into the loop update`
4. `docs: multi-signal discriminator combining and the differential group delay`

Landed as one commit in the end. With combining opt-in, no existing user's
numerics change: a single-signal satellite was always bit-identical, and a
multi-signal one now needs `signal_combining = true` to see any difference. The
commit is breaking only because `TrackedSignal` gained a field, so all-positional
calls to its default constructor need the new trailing argument.

## Follow-up: the receiver side

Not in this commit, and the reason the policy layer left Tracking in the first
place. GNSSReceiver.jl already owns everything needed:

- `CombinedSignal{pilot,data}` fixes the group ordering structurally — the pilot
  is always `signals[1]`, so the driver-ordering precondition is met by
  construction for all five shipped pairs (GPS L5 / L1C / L2C, Galileo E1 / E5a),
  where the pilot's primary code period is equal to or longer than the data
  component's.
- It holds the per-satellite decoder, so it can derive the bias and reconcile the
  two ICD sign conventions above in one place.
- It applies the downstream group-delay correction too, so both ends of the
  "shared `code_phase` means the driver's code phase" invariant sit in one
  package.

So the receiver enables `signal_combining` and supplies the bias — `0.0` for a
Galileo pair, the (correctly signed) decoded differential for a GPS or BeiDou
one — and the BeiDou sign gets its first real-data check there.
