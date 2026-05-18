# Multi-band tracking — design plan

**Date:** 2026-05-15
**Status:** Draft (design only — no implementation yet)
**Builds on:** `2026-05-15-multi-signal-tracking.md` (v2 multi-signal-per-sat)

## Goal

Let `track` / `track!` consume **multiple incoming measurements at once**, one
per RF band. Each measurement is its own complex sample vector/matrix, may use
a different sampling rate, but all share the same observation time window. A
tracked satellite is bound to exactly one band; routing the satellite's
correlators to the correct measurement happens by construction, with no
runtime band lookup in the hot loop.

The change extends — does not replace — the v2 multi-signal-per-sat design.
Within one band a satellite can still track multiple signals (L1 C/A + L1C-P +
L1C-D), as today. The new axis is *across bands*.

## Scope: in vs. deferred

**In scope (this plan):** the structural restructure that *enables*
cross-band joint tracking — `SignalGroup`s with explicit `band` fields,
per-band measurement routing, TrackState-level Doppler estimator preserved
so it can see all of a physical satellite's signals at once.

**Deferred to a follow-up:** the joint-tracking algorithm itself. Until
that follow-up lands, a sat tracked on multiple bands is registered as
multiple independent `TrackedSat`s (one per group, sharing a PRN) that
evolve independently. The cross-band fusion logic — the estimator
recognising those entries as one physical satellite and combining their
discriminators — is the next step, not this one. See "Cross-band satellite
linkage" below for the open mechanism question.

**Non-goals:**

- **Heterogeneous bit-sample formats across bands.** All measurements still
  use the same eltype family (`Complex{Int16/Int32/Float32/Float64}`).
- **Run-time changing of the signal-group set.** The groups tuple is frozen
  at TrackState construction (same as today's group tuple). Adding a
  new signal group needs a new TrackState; adding a *sat* to an existing
  group stays runtime-cheap.

## Design principles

### Single-band users write today's code

The over-arching constraint: a user with one band, one group, one signal
writes exactly the same code they write today. They never see `groups = (...)`,
never see a NamedTuple of measurements, never type the word "band":

```julia
# Single-band, single-group, single-signal — unchanged from today
track_state = TrackState(; signal = GPSL1CA())
add_satellite!(track_state; prn = 11, carrier_doppler = 1234.0Hz)
track!(buf, track_state, 4e6Hz)
```

```julia
# Single-band, multi-group — also unchanged from today
track_state = TrackState(; signals = (
    legacy_gps = (GPSL1CA(),),
    modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
    galileo    = (GalileoE1B(),),
))
add_satellite!(track_state; prn = 11, group = :modern_gps)
track!(buf, track_state, 4e6Hz)
```

Multi-band uses the same one-keyword `signals = (...)` construction —
the band is inferred from each group's first signal, so no extra keyword is
needed. `track` accepts a bare `VecOrMat` + scalar `sampling_frequency`
whenever the TrackState has one band.

### One code path, many ergonomic entry points

Internally everything is always multi-group — there is no "single-band fast
path" code duplication. The single-band constructors wrap their content in a
single auto-named group; the single-band `track(buf, state, fs)` method wraps
its inputs in a single-entry NamedTuple. One core implementation; many
shortcuts. The compiler sees a one-element NamedTuple and unrolls the group
loop to one iteration, so there's no runtime cost for the single-band case.

## Today's shape (one band only)

```julia
TrackState
├── satellites      :: NamedTuple{caps...}  # one Dictionary{Int, TrackedSat} per capability
├── signal_groups   :: NamedTuple{caps...}  # one Tuple{Vararg{AbstractGNSSSignal}} per capability
└── doppler_estimator
```

- `track(signal::AbstractVecOrMat, state, sampling_frequency)` — one signal
  buffer, one sampling rate, walks every capability and every sat against it.
- `num_ants` is per-capability today, recovered from the dict's slot type
  ([`_num_ants_for_capability`](src/tracking_state.jl#L565)).
- Every capability's signals are implicitly on the same band: there is no
  band field anywhere; the receiver code assumes all capabilities share one
  RF carrier in one measurement.

## Target shape — `SignalGroup`

The grouping layer exists for **type stability**: each group has one
signal-tuple shape (and therefore one concrete `TrackedSat` dictionary value
type). Today these groups are named `:legacy_gps` / `:modern_gps` / `:galileo`
and implicitly share one band. Multi-band makes the band an *explicit field*
on the group rather than a hidden global assumption.

Two groups may share a band (e.g. `:legacy_gps` and `:galileo` both on L1) —
the grouping is not by band, it is by signal-tuple shape. Band is metadata
each group carries so the right measurement is routed to it.

```julia
SignalGroup{B<:Band, S<:Dictionary{Int,<:TrackedSat}, Sigs<:Tuple{Vararg{AbstractGNSSSignal}}, NA<:NumAnts}
├── band         :: B            # e.g. L1(), L5() — from GNSSSignals
├── satellites   :: S            # Dictionary{Int, TrackedSat} for this group's shape
├── signals      :: Sigs         # signal-instance tuple (e.g. (GPSL1C_P(), GPSL1C_D(), GPSL1CA()))
└── num_ants     :: NA           # antenna count for this group's band
```

`TrackState` shrinks to:

```julia
TrackState{Groups<:NamedTuple{<:Any,<:Tuple{Vararg{SignalGroup}}}, DE}
├── groups            :: Groups
└── doppler_estimator :: DE
```

Notes:

- The `band` field uses GNSSSignals' [`Band`](src/bands.jl) type
  ([`L1`](src/bands.jl#L18), [`L5`](src/bands.jl#L26), …) and is
  auto-inferred from `get_band(signals[1])` when the user doesn't pass it
  explicitly. All signals in one group must share a band — invariant
  checked at construction.
- `num_ants` is per-group today (recovered from the slot type); making it
  an explicit field is a side-improvement that falls out for free and lets
  per-band antenna counts work naturally — two groups on the same band
  should agree on `num_ants` (validated at construction), groups on
  different bands may differ.
- `doppler_estimator` stays where it is. Per-sat estimator state on
  `TrackedSat` already carries band-specific gains; keeping the estimator
  TrackState-level leaves the door open for cross-band fusion of the same
  physical sat (e.g. PRN 5's L1 C/A and L5I instances sharing a Doppler
  estimate via the known L1↔L5 frequency ratio).

## API surface

### Construction

**One keyword does it all.** `signals = (...)` declares one or many
signal-tuple groups; the band of each group is inferred from
`get_band(signals[1])`. Single-band and multi-band share one entry point.

```julia
# Single signal, single group, single band — unchanged from today
TrackState(; signal = GPSL1CA())

# Multiple groups, all on one band — unchanged from today
TrackState(; signals = (
    legacy_gps = (GPSL1CA(),),
    modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
    galileo    = (GalileoE1B(),),
))

# Multiple groups, two bands — same keyword
TrackState(; signals = (
    legacy_gps_l1 = (GPSL1CA(),),
    modern_gps_l1 = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
    galileo       = (GalileoE1B(),),
    gps_l5        = (GPSL5I(),),
))
```

The fourth example is four groups across two bands (L1 for the first three,
L5 for the fourth). Multi-band is reached purely by adding a group whose
first signal lives on a different band — no new keyword, no nested config.

The group key (`:modern_gps_l1`, `:gps_l5`) is the user's choice — same
semantics as today's `capability` key, plus the freedom to encode band into
the name if the user wants.

**Per-band antenna count and other overrides.** When the user needs to set
`num_ants` (or, rarely, override the inferred `band`) per group, they
construct a `SignalGroup` directly:

```julia
TrackState(; signals = (
    legacy_gps_l1 = SignalGroup((GPSL1CA(),); num_ants = NumAnts(2)),
    gps_l5        = SignalGroup((GPSL5I(),);  num_ants = NumAnts(1)),
))
```

The `SignalGroup` constructor takes the signal tuple positionally; `band`
and `num_ants` are kwargs. `band` defaults to `get_band(signals[1])`;
`num_ants` defaults to `NumAnts(1)`. The `satellites` dictionary inside the
group starts empty and is populated by `add_satellite!`.

An entry's value in the `signals = (...)` NamedTuple may therefore be
either:

- a bare `Tuple{Vararg{AbstractGNSSSignal}}` — the common case, internally
  wrapped as `SignalGroup(tuple)`; or
- a `SignalGroup` instance — for the rare cases that need per-group
  overrides.

The two forms can be mixed inside one call. There is exactly one struct
(`SignalGroup`) that represents a group throughout the codebase; no
parallel config types.

### Adding / removing / querying satellites

The existing `capability` kwarg is renamed to `group` to match the new field
name. It hasn't been released yet, so the rename is a clean break — no
deprecation alias needed:

```julia
add_satellite!(track_state; prn = 11, group = :modern_gps_l1,
               code_phase = ..., carrier_doppler = ...)
```

- For a single-group TrackState (auto-wrapped single-band cases), `group`
  defaults to `:default`. Today's single-capability call sites that pass
  nothing keep working.
- For multi-group TrackStates that happen to be single-band (today's
  shape), the user passes `group = :modern_gps` exactly as they pass
  `capability = :modern_gps` today.

`remove_satellite!`, `remove_satellite`, `add_satellite`, `merge_sats`,
`get_sat_states`, `get_sat_state`, `get_system`, `estimate_cn0` all gain the
renamed `group` kwarg with the same defaulting rule.

### `Measurement` — one band's incoming buffer + its sample rate

Sample buffer and sample rate are inseparable in practice — a buffer at the
wrong claimed rate is a bug, not a feature. Bundling them into a single
struct removes the chance of mismatched parallel NamedTuples and shrinks the
`track` signature back to one positional argument.

```julia
Measurement{S<:AbstractVecOrMat, F}
├── samples              :: S            # complex sample buffer (Vector or Matrix-of-antennas)
├── sampling_frequency   :: F            # e.g. 4e6Hz
└── intermediate_frequency :: F          # defaults to 0.0Hz
```

Constructors:

```julia
Measurement(samples, sampling_frequency)                           # IF defaults to 0.0Hz
Measurement(samples, sampling_frequency, intermediate_frequency)
Measurement(samples; sampling_frequency, intermediate_frequency = 0.0Hz)  # kwarg form
```

`intermediate_frequency` lives on the measurement (not the TrackState)
because in a real receiver IF is a property of the front-end / digitiser
configuration for *this band*, alongside the sample rate. Moving it onto
`Measurement` co-locates everything the receiver knows about one band's
incoming chunk.

### The `track` call

**Single-measurement form (unchanged from today):**

```julia
track(buf, track_state, 4e6Hz)              # bare VecOrMat, bare scalar
track(buf, track_state, 4e6Hz; intermediate_frequency = 1.575e6Hz)
```

This dispatches when the TrackState has exactly one distinct band across all
its groups (which is also the auto-wrapped single-band case). Internally
wrapped to `Measurement(buf, 4e6Hz, …)` and routed through the multi-band
form.

**Multi-measurement form:**

```julia
track(
    (l1 = Measurement(l1_buf, 4e6Hz),
     l5 = Measurement(l5_buf, 25e6Hz)),
    track_state;
    downconvert_and_correlator = ...,
    preferred_num_code_blocks_to_integrate = 1,
)
```

The NamedTuple is keyed by **band**, not by group, since groups sharing a
band share the measurement. The L1 measurement is consumed by every group
whose `band == L1()`; the L5 measurement by every group whose `band == L5()`.

A single `Measurement` can also be passed for a single-band TrackState:

```julia
track(Measurement(buf, 4e6Hz), track_state)
```

So the dispatch ladder for the first positional argument is:

| Argument                  | Meaning                                          |
|---------------------------|--------------------------------------------------|
| `AbstractVecOrMat`        | Today's bare buffer. Single-band TS only. `sampling_frequency` passed as second positional arg, `intermediate_frequency` as kwarg. |
| `Measurement`             | One band's bundled buffer + rates. Single-band TS only. |
| `NamedTuple{<:Any,<:Tuple{Vararg{Measurement}}}` | Multi-band: one `Measurement` per band key. |

**Boundary validation** on entry:

1. The set of keys in the measurements NamedTuple equals the set of distinct
   bands used across `track_state.groups`. (Wrong keys → clear error.)
2. Each measurement's `samples` has the antenna shape that band's groups
   declare (`size(samples, 2) == num_ants`).
3. All measurements cover the exact same observation duration:
   `num_samples / sampling_frequency` is identical across bands. No
   tolerance — mismatched durations are an error.

These are O(num_bands) checks at the entry point — irrelevant next to the
inner loop.

### How "band key" maps from `Band` type to NamedTuple key

The measurements NamedTuple is keyed by `Symbol`s (`:l1`, `:l5`), but the
`SignalGroup`'s `band` field is a GNSSSignals `Band` instance (`L1()`,
`L5()`). The mapping is by name: `Symbol(lowercase(nameof(typeof(band))))`
yields `:l1` for `L1()`, `:l5` for `L5()`. Built once at TrackState
construction and cached; no per-call symbol munging in the hot loop.

Alternative considered: key the measurements NamedTuple by the group key
instead of by band key. Rejected because two groups on the same band would
then need duplicate buffers; keying by band makes the shared-measurement
case the natural default.

## Inner loop

The loop must keep the **two-phase structure** of today: first advance
every group's downconvert+correlate, *then* run estimation. The phase
boundary is what makes cross-band joint tracking possible — the estimator
needs all bands' correlator outputs in hand at once so a fused L1+L5
discriminator (for example) can see both sides.

```julia
while !all_groups_reached_end(track_state, measurements)
    # Phase 1 — downconvert + correlate every group.
    # Each group advances against its own band's measurement, walking to
    # its next integration boundary (or no-op if already at the chunk end).
    for group_key in keys(track_state.groups)        # unrolled at compile time
        g = track_state.groups[group_key]
        m = measurements[band_key_for(g.band)]
        downconvert_and_correlate!(dc, m, g, preferred_num_code_blocks_to_integrate)
    end

    # Phase 2 — estimation sees the whole TrackState.
    # Every group that completed an integration in Phase 1 has fresh
    # correlator output ready; the estimator decides whether to walk
    # groups independently (today's behaviour) or fuse across them
    # (cross-band joint tracking).
    estimate_dopplers_and_filter_prompt!(track_state, measurements,
                                         preferred_num_code_blocks_to_integrate)
end
```

Iterating `keys(...)` over a NamedTuple unrolls at compile time (same
machinery v2's `_foreach_system!` uses), so Phase 1 stays fully type-stable
and allocation-free.

### Integration cadence across groups

Per-group cadence (each group's "shortest signal" boundary, v2's MIN-cadence
rule) is unchanged. Across groups the loop uses the **same MIN rule**: the
outer `while` ticks at the rate of whichever group finishes its current
integration first. Groups that didn't hit a boundary this tick contribute
no fresh correlator output, and the estimator treats them as no-ops for
this tick.

This means:

- In the common case (L1 C/A and L5I both integrating 1 ms): the two
  groups hit a boundary in the same wall-clock moment, both contribute
  fresh correlator outputs, and the estimator fuses them every tick.
- In the mixed case (e.g. L1 C/A at 1 ms vs L1C-P at 20 ms): the L1 C/A
  group's correlator output is fresh on every tick; the L1C-P group
  contributes only on every 20th tick. The estimator sees this exactly
  as it sees today's multi-signal-per-sat case within one group.

The outer `while` exits once every group has consumed its measurement to
the end.

### Why estimation sees the whole TrackState

Today's `estimate_dopplers_and_filter_prompt!` takes one group's worth of
state because there was only ever one group's worth of correlator output
in flight at the phase boundary. With multiple groups, the estimator
needs visibility across them to do cross-band fusion — the same physical
satellite's L1 and L5 discriminators must be combinable in one
estimator step.

The signature change is mechanical: the function widens from "per-group"
to "per-TrackState," and the default implementation walks groups
independently (preserving today's behaviour). Cross-band fusion is added
later as an estimator method that joins across groups by matching PRNs
(see "Cross-band satellite linkage" below).

### Downconvert/correlate dispatch

The existing `_foreach_system!` walk
([downconvert_and_correlate_cpu.jl:616](src/downconvert_and_correlate_cpu.jl#L616))
iterates capabilities within one signal. The change is at the layer above:
the band-keyed measurement lookup happens in Phase 1 of the outer
`track` / `track!` loop, and each group is dispatched against *its*
measurement. **No changes to the fused/tile-share kernels** — they take a
single measurement plus a single group of sats, and that contract is
preserved.

## Cross-band satellite linkage (open item)

The Doppler estimator stays TrackState-level specifically so the same
physical satellite tracked on multiple bands (e.g. PRN 5 on L1 C/A and L5I)
can share a Doppler estimate, exploiting the known carrier-frequency ratio.
Realising this requires the estimator to recognise that two `TrackedSat`
entries — one in group `:legacy_gps_l1`, one in group `:gps_l5`, both with
`prn == 5` — are the same physical satellite.

**Not solved by this plan.** The plan must not close the door on it. Two
candidate mechanisms (decide later):

1. **Implicit by `(prn, band)`.** The estimator scans groups; sats sharing
   `prn` across different bands are linked. Simple but couples estimator
   logic to the group structure.
2. **Explicit handoff field.** `TrackedSat` gains an optional
   `physical_sat_id` (defaults to `prn`); the estimator groups by that ID.
   More flexible (lets two PRNs alias if needed) at the cost of one field.

The existing [`update_estimator_on_handoff`](src/sat_state.jl) hook is where
either mechanism plugs in. Today's tests don't exercise multi-band fusion,
so we land the structural plan first and revisit when a concrete fusion use
case appears.

## Backwards compatibility

| Caller pattern                                                | Stays working? | How                                  |
|---------------------------------------------------------------|----------------|--------------------------------------|
| `TrackState(; signal = GPSL1CA())`                            | yes            | wraps in single `:default` group     |
| `TrackState(; signals = (modern_gps = ...,))`                 | yes            | each entry wrapped into a SignalGroup with band inferred from `get_band(signals[1])` |
| `track(buf, state, fs)`                                       | yes            | scalar-signal shortcut, dispatched when state has one distinct band |
| `track!(buf, state, fs)`                                      | yes            | same shortcut                        |
| `add_satellite!(state; prn = …, capability = …)`              | **breaks**     | rename `capability` → `group` (unreleased on this branch); `band` is implicit |
| Legacy `TrackState(system, sats)` ([tracking_state.jl:115](src/tracking_state.jl#L115)) | yes | single-group shim, band inferred |
| Direct construction `TrackState(satellites::SatelliteDicts; …)` ([tracking_state.jl:168](src/tracking_state.jl#L168)) | **breaks** | rare, used in tests; migrate to `TrackState(; groups = ...)` or to a new `TrackState(groups::NamedTuple{<:Any,<:Tuple{Vararg{SignalGroup}}}; …)` ctor |

The one breaking surface is direct construction from a raw NamedTuple of
dicts. We'll migrate the small number of internal call sites.

## Implementation order

This is a substantial refactor; staging matters for keeping CI green.

1. **Introduce `SignalGroup` as a wrapper** that holds today's
   `(satellites_dict_for_one_capability, signals_tuple, num_ants, band)`.
   Make `TrackState` hold `groups :: NamedTuple` of `SignalGroup`s. The
   existing `satellites` and `signal_groups` fields become *views* /
   property forwards over `groups` so today's accessors keep working.
   Single-band construction auto-fills `band` via `get_band(signals[1])`.
   **No new functionality yet** — tests must stay green.

2. **Rename `capability` → `group`** on all public-facing kwargs. Clean
   break — the `capability` kwarg is unreleased on this branch. Update
   `add_satellite*`, `remove_satellite*`, `get_sat_state*`, `get_sat_states*`,
   `merge_sats`, `estimate_cn0`, `get_system`, and call sites in the test
   suite. Still no new functionality.

3. **Lift the inner-loop walk** to route by band. With one band the routing
   is trivial (every group reads the single measurement), so behaviour is
   unchanged. The downconvert/correlate kernels see exactly the same
   per-group calling convention.

4. **Open the API to multiple bands.** Allow groups with differing `band`
   fields. Allow `signals` and `sampling_frequencies` as NamedTuples keyed
   by band in `track` / `track!`. Add the boundary validation. Add a
   multi-band test (e.g. GPS L1 C/A + GPS L5I on a synthetic two-band
   signal).

5. **Per-band antenna counts.** Already wired by step 1 — verify by adding
   a test with different `NumAnts` per band, plus a validation test that
   two groups on the same band must agree on `num_ants`.

Each step is a self-contained commit. Steps 1–3 are non-functional
refactors and can land independently of step 4 if we want to defer the
public API surface change.

## Risks and open questions

- **Different sampling rates → different chunk lengths.** The constraint
  is exactly-equal *time* across bands, not equal *samples*. An L1 chunk
  at 4 MHz and an L5 chunk at 25 MHz over 1 ms is `4000` vs `25000`
  samples. Each group's `num_samples_signal` is read from its own band's
  measurement — trivially handled. The risk is users handing in
  mismatched-duration chunks; the boundary check catches that and errors
  on any inequality (no tolerance).
- **Two groups on the same band must agree on `num_ants`.** Same physical
  front-end → same antenna count. Validated at construction.
- **Acquisition handoff.** `TrackingLoopFilters` / acquisition code that
  builds `TrackedSat`s ahead of `add_satellite!` doesn't change — it still
  builds one sat per (prn, band) pair. The caller of `add_satellite!`
  decides which group the sat lands in.
- **Renamed accessor for direct walkers.** No new axis — the satellites
  layout is still a one-dimensional NamedTuple of per-group dicts, just
  longer (e.g. `:gps_legacy_l1, :gps_modern_l1, :galileo, :gps_l5`
  instead of `:gps_legacy, :gps_modern, :galileo`). Anything that walks
  `track_state.satellites[k]` directly switches to
  `track_state.groups[k].satellites`. Pure rename — grep covers it.
- **Cross-band single-sat fusion** (see "Cross-band satellite linkage"
  above) — deferred; this plan keeps the door open.

## Test plan (sketch)

- **Structural (steps 1–3):** existing test suite stays green with no
  changes to test code — all call sites still single-band.
- **Multi-band correctness (step 4):** generate a synthetic two-band
  signal (e.g. one GPS PRN on L1 C/A at 4 MHz and the same PRN on L5I at
  25 MHz, equal duration), track both, verify both bands converge to the
  injected Doppler / code phase.
- **Per-band antenna (step 5):** `NumAnts(2)` on L1, `NumAnts(1)` on L5;
  verify correlator widths match per band. Also a *negative* test: two
  groups on the same band with disagreeing `num_ants` must error at
  construction.
- **Boundary validation:** wrong keys, mismatched duration, mismatched
  antenna count → clear `ArgumentError`s.
- **No regression in single-band hot path:** rerun the v2 benchmarks; the
  group wrapper must add zero allocations and ≤ negligible runtime.
