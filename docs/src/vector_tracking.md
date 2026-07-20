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
    freshly acquired satellite tracks in until it is promoted into the
    vector loop.

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
    filter to read and reset.

The satellite-shared carrier/code Doppler is always updated through the same
carrier-aiding (`aid_dopplers`) used by the conventional estimator, and
the same `1/N` loop-bandwidth scaling applies when a signal integrates `N`
primary code blocks coherently.

## The receiver-side loop

A navigation filter drives `VectorPLLAndDLL` through the exported state
managers. A typical iteration, after [`track!`](@ref) has produced new
correlations:

```julia
# 1. Promote satellites that have pulled in, and flag the ones currently
#    in an outage. `prns_in_lock` is whatever the receiver decides is
#    usable (e.g. a C/N0 threshold).
update_vt_states!(track_state, prns_in_lock)

# 2. Read the accumulated discriminator outputs the estimator collected
#    for each vector-loop satellite, then reset the accumulators so the
#    next block accumulates afresh.
for (prn, sat) in pairs(get_sat_states(track_state))
    state = get_doppler_estimator_state(sat)
    state.vt_on || continue
    code_err = mean_code_discriminator(state)      # chips, or `nothing` if no data
    carrier_err = mean_carrier_discriminator(state) # Hz, or `nothing` if no data
    # … feed the mean measurements into the navigation filter …
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
[`mean_code_discriminator`](@ref) / [`mean_carrier_discriminator`](@ref) apply
the averaging convention (`sum / count`, returning `nothing` when nothing has
accumulated) in one place, so consumers don't each re-implement the divide and
the `count == 0` guard. Reading and resetting are deliberately separate calls
so the filter can read at its own (typically slower) rate than `track!`.

### Multi-constellation addressing

Each state manager has an all-groups form and a **group-scoped** form that
takes a group selector (`Symbol`, `Integer`, or `Val`) as its first argument
after `track_state`:

```julia
update_vt_states!(track_state, :gps, gps_prns_in_lock)
set_code_freq_updates!(track_state, :galileo, galileo_code_updates)
```

In a multi-constellation receiver a PRN alone is ambiguous — GPS PRN 5 and
Galileo PRN 5 are different satellites in different groups — so address each
constellation's group explicitly. The all-groups form matches a PRN in every
group and is only unambiguous for a single-group `TrackState`.

Promotion is one-directional: [`update_vt_states!`](@ref) only ever *joins*
satellites to the vector loop. A satellite in an outage keeps being steered by
the navigation filter from the shared solution; deciding whether its (now
uninformative) discriminator outputs should feed the measurement update is the
navigation filter's responsibility, tracked on the receiver side rather than
on the estimator state. Use [`set_vt_on!`](@ref) to force loop membership on
or off manually.

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
Tracking.SatVectorPLLAndDLL
update_vt_states!
set_vt_on!
set_code_freq_updates!
set_carrier_freq_updates!
reset_code_discr_acc!
reset_carrier_discr_acc!
mean_code_discriminator
mean_carrier_discriminator
```
