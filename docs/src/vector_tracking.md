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
    filter to read and reset — one dump per *signal*, see
    [Per-signal discriminator dumps](#Per-signal-discriminator-dumps).

A multi-signal satellite can fold its signals' discriminators into one
minimum-variance loop update here too, with
`VectorPLLAndDLL(signal_combining = true)` — the same weighting and the same
driver-ordering precondition
[Multi-signal discriminator combining](tracking_state.md#Multi-signal-discriminator-combining)
describes. Which loops it reaches follows `vt_on`, because that is what decides
which loops this package still closes:

  - **`vt_on = false`** — every loop is local, so this *is* scalar tracking and
    all three discriminators combine, differential group delay included. A
    satellite pulls in with the full combining gain.
  - **`vt_on = true`** — the carrier phase loop alone. The code and carrier
    frequency loops are the navigation filter's, and their discriminators reach
    it as [one raw dump per signal](#Per-signal-discriminator-dumps), which is
    where they should be fused: the filter has each signal's measured C/N₀, tap
    spacing and dump length, so it can weigh them better than a nominal power
    split can, and fusing here as well would fuse the same measurements twice.
    For the same reason the differential group delay goes unread in this mode —
    referring a passenger's code measurement to one ranging datum belongs with
    the consumer that ranges on it.

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
    state = get_doppler_estimator_state(sat)
    state.vt_on || continue
    for i in eachindex(get_signals(sat))
        code_err = mean_code_discr(state, i)      # chips, or `nothing` if no data
        carrier_err = mean_carrier_discr(state, i) # Hz, or `nothing` if no data
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

The accumulators are stored as `(count, sum)` tuples;
[`mean_code_discr`](@ref) / [`mean_carrier_discr`](@ref) apply
the averaging convention (`sum / count`, returning `nothing` when nothing has
accumulated) in one place, so consumers don't each re-implement the divide and
the `count == 0` guard. Reading and resetting are deliberately separate calls
so the filter can read at its own (typically slower) rate than `track!`.

### Per-signal discriminator dumps

A satellite tracked on several signals accumulates **one `(count, sum)` pair per
signal**, in `sat.signals` order, and each is that signal's own discriminator —
never a mean across signals. `mean_code_discr(state, i)` and
`mean_carrier_discr(state, i)` read slot `i`; both default to `1`, the
estimator-driver signal, which is the ranging signal and the only slot a
single-signal satellite has.

Fusing them is deliberately left to the navigation filter, because that is where
the information to do it well already is. A filter that models its own
measurement noise — as GNSSReceiver.jl's does, from C/N₀, early-late spacing and
coherent dump length — can weigh each signal by what it actually measured, which
beats any weighting this package could apply from a nominal ICD power split. It
also keeps the variance honest: a dump quietly turned into a two-signal mean
would enter such a filter carrying one signal's variance.

Two things a consumer owes these measurements:

  - **A ranging datum for the code dumps.** Every signal shares the satellite's
    one `code_phase`, and their true code phases differ by the satellite's
    differential payload group delay, so a non-driver signal's code dump
    measures the driver's error *plus* that offset. Subtract it before fusing —
    the same quantity [`set_differential_group_delay!`](@ref) takes, which the
    scalar fallback applies for you and this mode deliberately does not. The
    carrier dumps need no counterpart: the signals of one satellite share a
    carrier, so their frequency dumps measure one Doppler.
  - **Room in the measurement model.** Several dumps from one satellite are
    several measurements along **one** line of sight. Whether they enter as
    separate rows sharing a geometry row, or are pre-fused into one row, is the
    filter's design choice; what they must not be is treated as independent
    satellites. Note too that their thermal noise is near-independent (the
    components' codes barely cross-correlate) while their orbit, clock and
    atmospheric errors are common-mode — so fusing sharpens the thermal part and
    nothing else.

Both dumps of every signal are cleared together by
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
