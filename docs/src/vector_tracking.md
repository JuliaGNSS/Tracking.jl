# Vector Tracking

[`VectorPLLAndDLL`](@ref) is a Doppler estimator for **vector tracking**
(a vector delay/frequency-lock loop, VDFLL). Where the
[`ConventionalPLLAndDLL`](@ref) closes each satellite's carrier and code
loops independently with per-satellite loop filters, vector tracking closes
them *centrally*: an external navigation filter — living outside Tracking.jl,
e.g. in GNSSReceiver.jl — combines every satellite's measurements with the
receiver dynamics and feeds back per-satellite NCO corrections. A satellite
in a deep fade keeps being steered by the solution the healthy satellites
constrain, which is what gives vector tracking its robustness under weak
signal and high dynamics.

Tracking.jl owns only the signal-path half of that loop: it correlates,
forms discriminators, and applies the corrections the navigation filter
hands back. `VectorPLLAndDLL` is the interface between the two.

## The per-integration contract

`VectorPLLAndDLL` is configuration-only; the per-satellite state lives on
each [`TrackedSat`](@ref) as a `SatVectorPLLAndDLL` (seeded through
[`init_estimator_state`](@ref), like every estimator — see
[Custom Doppler Estimator](custom_doppler_estimator.md)). Each time a
satellite's estimator-driver signal (`signals[1]`) completes an integration,
`VectorPLLAndDLL` does one of two things depending on that satellite's
`vt_on` flag:

  - **`vt_on = false` — scalar fallback.** The satellite runs an ordinary
    FLL-assisted PLL and DLL, identical to
    [`ConventionalAssistedPLLAndDLL`](@ref) (same default loop-filter types
    and the same auto-sized loop bandwidths). This is the pull-in mode a
    freshly acquired satellite tracks in until [`enable_vt!`](@ref) puts it
    into the vector loop.

  - **`vt_on = true` — vector closure.** The navigation filter drives the
    NCOs instead of the local loop filters:

      * the **code Doppler follows `code_freq_update` directly** (the DLL
        loop filter is bypassed and holds its state);
      * the **carrier's FLL branch is driven by `carrier_freq_update`**
        while the PLL branch still runs on the satellite's own phase
        discriminator — so the carrier is a vector *frequency* lock loop
        with a retained local *phase* lock. (This needs the FLL-assisted
        `ThirdOrderAssistedBilinearLF`, the default carrier filter; with a
        plain PLL filter the `carrier_freq_update` has no input path and the
        vector carrier closure degrades to PLL-only.)

    In this mode the DLL and FLL discriminator outputs and the prompt
    magnitude are **accumulated** on the per-sat state for the navigation
    filter to read and reset — one per *signal*, see
    [Per-signal discriminator accumulators](#Per-signal-discriminator-accumulators).

A multi-signal satellite can fold its signals' discriminators into one
minimum-variance loop update here too. The setting is the same one the
conventional estimator reads — `discriminator_combining = true` on the
satellite's [`SignalGroup`](@ref), not a keyword of this estimator — and so are
the weighting and the driver-ordering precondition
[Multi-signal discriminator combining](tracking_state.md#Multi-signal-discriminator-combining)
describes. Which loops it reaches follows `vt_on`, because that is what decides
which loops this package still closes:

  - **`vt_on = false`** — every loop is local, so this *is* scalar tracking and
    all three discriminators combine. A satellite pulls in with the full
    combining gain.
  - **`vt_on = true`** — the carrier phase loop alone. The code and carrier
    frequency loops are the navigation filter's, and their discriminators reach
    it as [one raw measurement per signal](#Per-signal-discriminator-accumulators)
    for it to fuse.

Collecting those measurements is a separate matter from combining, and follows
`vt_on` **alone**: a satellite in the vector loop hands the filter every one of
its signals' measurements whether or not `discriminator_combining` is set, because the
flag says only what the loops this package still closes are fed. The two
settings are independent in both directions — combining without vector closure
is the scalar fallback above, and vector closure without combining collects
every signal's measurement while the loops stay the driver's own.

The satellite-shared carrier/code Doppler is always updated through the same
carrier-aiding (`aid_dopplers`) used by the conventional estimator, and the
same effective-bandwidth handling applies when a signal integrates `N` primary
code blocks coherently: `1/N` on the carrier loop, a stability cap against the
actual integration time on the code loop.

## The receiver-side loop

A navigation filter drives `VectorPLLAndDLL` through the exported state
managers. A typical iteration, after [`track!`](@ref) has produced new
correlations:

```julia
# 1. Put the satellites that have pulled in into the vector loop.
#    `prns_in_lock` is whatever the receiver decides is usable (e.g. a
#    C/N0 threshold); re-issuing the same set every iteration is a no-op.
enable_vt!(track_state, prns_in_lock)

# 2. Read the accumulated discriminator outputs the estimator collected
#    for each vector-loop satellite, then reset the accumulators so the
#    next block accumulates afresh.
for (prn, sat) in pairs(get_sat_states(track_state))
    get_doppler_estimator_state(sat).vt_on || continue
    for i in eachindex(get_signals(sat))
        code = mean_code_discr(sat, i)        # mean DLL discriminator, in chips
        carrier = mean_carrier_discr(sat, i)  # mean FLL discriminator, in Hz
        isnothing(code) && isnothing(carrier) && continue
        # Each mean covers the records since the last reset, all of one
        # coherent integration time — the one this receiver configured.
        # … feed this signal's measurements into the navigation filter …
    end
end
reset_code_discr_acc!(track_state)
reset_carrier_discr_acc!(track_state)

# 3. Run the navigation filter, then feed its per-satellite NCO
#    corrections back. Both setters take anything indexable by PRN
#    (e.g. a Dictionaries.Dictionary) with an entry for every vector-loop
#    satellite; satellites still in the scalar fallback are skipped.
set_code_freq_updates!(track_state, code_freq_updates)
set_carrier_freq_updates!(track_state, carrier_freq_updates)
```

Each accumulator is a `(count, sum)` pair. [`mean_code_discr`](@ref) /
[`mean_carrier_discr`](@ref) apply the averaging convention (`sum / count`,
returning `nothing` when nothing has accumulated) in one place, so consumers
don't each re-implement the divide and the `count == 0` guard. Reading and
resetting are deliberately separate calls so the filter can read at its own
(typically slower) rate than `track!`.

!!! note "What the pair assumes, and when it holds"

    A count and a sum carry no per-record timing, so the consumer places the
    mean itself: it covers `count` records of the coherent integration time the
    consumer configured, ending at the last record before it read. That is
    exact while every record the pair covers has the same length — which is what
    vector tracking gives it, on two counts.

    Every record is counted. `fll_disc` answers a placeholder 0 Hz when there is
    no previous prompt to difference against, and nothing rejects a record on
    that ground or any other, so the carrier `count` never falls behind the
    code one. The two cases that could produce a gap both sit before the
    satellite joins the vector loop: only a signal's very first record lacks a
    prompt pair (and the first after [`reset_loop_filters!`](@ref), which clears
    the prompt and the accumulators together), and a record is excluded from a
    passenger's combination only when it was correlated with a replica the
    bit/secondary-code sync had not yet corrected — the fold that *detects*
    sync, and no later one.

    Every record has the same length. A signal's coherent integration length
    changes exactly once, when its bit or secondary code is found, and nothing
    in this package changes it again. A navigation filter admits a satellite to
    the vector loop only once it is decoded, which implies that sync long
    since happened, so `vt_on` covers a stretch of constant-length records.

    The one caller-side way out of both is
    [`set_preferred_num_code_blocks_to_integrate!`](@ref): changing a signal's
    cadence while `vt_on` is set makes the mean cover records of two lengths
    with nothing in the pair to say so, and `count` stops mapping to an
    interval. Reset the accumulators across such a change, or make it while the
    satellite is out of the vector loop.

!!! note "Counting records is not weighting them"

    `count` is what a filter needs to weigh one signal against another, but the
    mean's variance is not the single-record variance over `count`, and for the
    FLL it is *smaller*. Successive FLL measurements are formed from overlapping
    prompt pairs — each record's prompt is the next one's reference — and the
    shared prompt enters the two differences with opposite signs, so their
    errors are negatively correlated and a contiguous run telescopes into one
    frequency estimate over the whole span: averaging `N` equal-length records
    buys `N²`, not `N` (independent prompt-phase noise assumed). Under vector
    tracking the run is contiguous by construction, per the note above. The
    DLL's records carry no shared prompt and no such term, so `count` is exactly
    right there.

### Per-signal discriminator accumulators

A satellite tracked on several signals accumulates **one `(count, sum)` pair per
signal and per loop**, in `sat.signals` order, and each is that signal's own discriminator —
never a mean across signals. [`mean_code_discr`](@ref) and
[`mean_carrier_discr`](@ref) take the same trailing signal selector as every
other per-signal accessor — an index or a signal type — which may be omitted
only when the satellite tracks one signal:

```julia
mean_code_discr(track_state, :gps_l1, 7, GPSL1C_D)   # by signal type
mean_code_discr(sat, 2)                              # by index
mean_code_discr(get_doppler_estimator_state(sat))    # single-signal satellite
```

The `SatVectorPLLAndDLL` form takes an index rather than a selector, since the
estimator state holds the accumulators but not the signals and so cannot resolve
a type. Asking a multi-signal satellite for "the" mean is an error, not a
driver-only read — the same refusal [`estimate_cn0`](@ref) and the other
per-signal accessors make.

Fusing them is deliberately left to the navigation filter, because that is where
the information to do it well already is. A filter that models its own
measurement noise — as GNSSReceiver.jl's does, from C/N₀, early-late spacing and
coherent integration length — can weigh each signal by what it actually
measured, which beats any weighting this package could apply from a nominal ICD
power split. It also keeps the variance honest: a measurement quietly turned
into a two-signal mean would enter such a filter carrying one signal's
variance.

The code measurements arrive **raw**, each on its own signal's code phase: the
[`get_group_delay`](@ref) difference against `signals[1]` is **not** subtracted.
That value belongs to the *combined code loop* — it exists so that a loop fed by
several signals still drives the satellite's shared `code_phase` to one datum —
and under `vt_on` there is no such loop here: the code loop is yours. Applying it
anyway would put Tracking in the business of correcting a measurement whose
consumer corrects it again downstream, since a receiver forming a pseudorange
already applies a per-signal group delay.

So a signal's accumulator is exactly what that signal measured, from its first
integration, whether or not a group delay has ever been supplied for it. The
scalar fallback is the other half of the same rule: while `vt_on` is unset this
package *does* close the combined code loop, so there it applies the difference
itself — see [Group delay](tracking_state.md#Group-delay).

What a consumer owes these measurements:

  - **The inter-signal bias, if it fuses them.** Several signals of one satellite
    do not share a code phase: their payload group delays differ, and two raw
    code discriminators therefore sit on two datums. Fusing them without
    referring both to one datum folds that difference into the result — the same
    moving bias the combined code loop's referral exists to prevent, arrived at
    from the filter's side. [`get_group_delay`](@ref) returns the per-signal
    value if you have given Tracking one; a filter holding its own ISC tables can
    equally use those. A filter that keeps one row per signal and lets the
    estimator absorb each signal's bias needs nothing.
  - **Room in the measurement model.** Several measurements from one satellite are
    several measurements along **one** line of sight. Whether they enter as
    separate rows sharing a geometry row, or are pre-fused into one row, is the
    filter's design choice; what they must not be is treated as independent
    satellites. Note too that their thermal noise is near-independent (the
    components' codes barely cross-correlate) while their orbit, clock and
    atmospheric errors are common-mode — so fusing sharpens the thermal part and
    nothing else.

Every signal's accumulators are cleared together by
[`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref); there is no
per-signal reset, because the filter reads a satellite's signals in one cycle.

### Multi-constellation addressing

Each state manager has an all-groups form and a **group-scoped** form that
takes a group selector (`Symbol`, `Integer`, or `Val`) as its first argument
after `track_state`:

```julia
enable_vt!(track_state, :gps, gps_prns_in_lock)
set_code_freq_updates!(track_state, :galileo, galileo_code_updates)
```

In a multi-constellation receiver a PRN alone is ambiguous — GPS PRN 5 and
Galileo PRN 5 are different satellites in different groups — so address each
constellation's group explicitly. The all-groups form matches a PRN in every
group and is only unambiguous for a single-group `TrackState`.

### Loop membership

[`enable_vt!`](@ref) and [`disable_vt!`](@ref) are the only way loop membership
changes; the estimator never joins or drops a satellite on its own. Membership
is not a lock indicator: a satellite in an outage stays in the loop and keeps
being steered by the navigation filter from the shared solution. Deciding
whether its (now uninformative) discriminator outputs should feed the
measurement update is the navigation filter's responsibility, tracked on the
receiver side rather than on the estimator state — [`disable_vt!`](@ref) is for
satellites the filter gives up on entirely, and wants handed back to their own
scalar loop (follow it with [`reset_loop_filters!`](@ref) for a
transient-free handoff).

## Resetting

[`reset_loop_filters!`](@ref) zeroes the loop-filter integrators, the
discriminator accumulators, **and** the NCO corrections, re-seeding the loop
from the satellite's current (converged) Doppler while preserving the
`vt_on` flag and any per-satellite bandwidth override. The
NCO corrections are zeroed because the re-seeded Doppler already contains the
last correction — keeping it would apply it twice. After a reset the
navigation filter must re-issue its corrections via
[`set_code_freq_updates!`](@ref) / [`set_carrier_freq_updates!`](@ref) before
the next vector-closed integration.

## API reference

```@docs
VectorPLLAndDLL
SatVectorPLLAndDLL
enable_vt!
disable_vt!
set_code_freq_updates!
set_carrier_freq_updates!
reset_code_discr_acc!
reset_carrier_discr_acc!
mean_code_discr
mean_carrier_discr
```
