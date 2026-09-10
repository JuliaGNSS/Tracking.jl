# Custom Doppler Estimator

Tracking.jl ships [`ConventionalPLLAndDLL`](@ref) (and the FLL-assisted
variant [`ConventionalAssistedPLLAndDLL`](@ref)), but you can plug in a
different Doppler-estimation algorithm — e.g. a Kalman filter or a
joint-channel estimator — by implementing a small set of methods.

## Where state lives

A [`TrackedSat`](@ref) carries everything per-satellite: the shared
carrier/code Doppler and phase, a tuple of [`TrackedSignal`](@ref)s
(one per signal being tracked), and the per-satellite estimator state
in `doppler_estimator_state`. The estimator object itself
(`<: AbstractDopplerEstimator`) is configuration plus any *shared* state
that spans satellites or systems. Cleanly separating per-sat state from
shared state makes the signal-path code (downconvert, correlate,
sample-bookkeeping) estimator-agnostic — it touches the per-signal and
per-sat fields directly and rewraps `doppler_estimator_state` unchanged.

## What you implement

1. **An estimator type** subtyping [`AbstractDopplerEstimator`](@ref). It
   carries configuration and any cross-satellite or cross-system shared
   state (filter parameters, joint-state vectors, …). Per-satellite state
   does *not* live here.

2. **A per-satellite state struct** (any name and shape you like). For
   the conventional estimator it is `SatConventionalPLLAndDLL` holding
   the loop filters and seed Dopplers; for a Kalman filter it might be
   the Kalman state vector.

3. **An [`init_estimator_state`](@ref) method** for your estimator.
   It produces the initial per-sat state for a sat entering the track
   set (acquisition → tracking handoff):

   ```julia
   Tracking.init_estimator_state(::MyEstimator, sat::TrackedSat) = MyPerSatState(...)
   ```

   It must be **pure** (no observable side effects): besides seeding
   real sats, it is also called on the throwaway PRN-0 template sat that
   fixes a group's slot type at `TrackState` construction, and as a type
   probe when validating pre-built sats. Register sats into shared state
   in `update_estimator_on_handoff` (below), never here.

   [`reset_loop_filters!`](@ref) re-seeds through
   `Tracking._reset_estimator_state(estimator, sat)`, which falls back to
   `init_estimator_state`. Specialize it if your per-sat state carries
   *configuration* a reset must keep — both shipped estimators do, to preserve
   a per-satellite loop bandwidth or `signal_combining` override that going
   through `init_estimator_state` would silently revert to the estimator's
   defaults.

4. **(Optional) An [`update_estimator_on_handoff`](@ref) method**, if
   your estimator carries cross-satellite or cross-system shared state
   that has to change when sats join the track set (joint-channel
   covariance, per-batch normalization terms, …):

   ```julia
   Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats) = ...
   ```

   It is called once per [`add_satellite!`](@ref) / [`merge_sats`](@ref)
   call with the dictionary of incoming sats and runs *after* per-sat
   seeding, so `init_estimator_state` sees the pre-update shared state.
   The default returns `est` unchanged, so estimators with no shared
   state need not implement it.

   **Type constraint:** the returned estimator must have the same
   concrete type as the input. [`TrackState`](@ref) is parameterized on
   the estimator type, and changing it would break inference. Keep the
   estimator an immutable `struct` and put any growing shared state in
   resizable containers (e.g. `Vector`, `Matrix`) you `push!`/`resize!`
   in place; if scalar fields need replacing, rebuild the estimator
   with `Setfield.@set` or a copying constructor.

   **Use the return value of the handoff functions.** Every entry point
   carries the estimator you return in the `TrackState` it gives back.
   Since `TrackState` is immutable, even [`add_satellite!`](@ref) can
   only honor a *rebuilt* estimator through its return value — write
   `track_state = add_satellite!(track_state; ...)`, not a bare
   `add_satellite!(track_state; ...)`, if your estimator rebuilds
   itself on handoff. (For in-place estimators, including the
   conventional ones, the returned `TrackState` is the input itself.)

5. **An `estimate_dopplers_and_filter_prompt` method** dispatched on
   `TrackState{<:Any, <:MyEstimator}`. This is where the actual update
   logic runs, once per **chunk** (see [Chunked Doppler updates](track.md#Chunked-Doppler-updates)).
   It walks each group in `track_state.groups`, reads the band's
   `BandMeasurement` from the `measurements::BandMeasurements` NamedTuple via
   `get_band_id(group.band)`, and produces new `TrackedSat`s with updated
   `carrier_doppler`/`code_doppler` and updated per-sat estimator state.

   Each signal's correlator outputs completed during the chunk are in its
   `correlator_outputs::Vector{`[`CorrelatorOutput`](@ref)`}` (a chunk may hold
   zero, one, or several per signal), each carrying the raw correlator, its
   integrated-sample count, and the end sample index.
   Fold over them in order — threading whatever filter state you carry — and
   write the NCO Doppler once (typically from the last output). Empty each
   signal's `correlator_outputs` when done, and remember the NCO Doppler is held
   fixed across the whole chunk, so a discriminator that needs the replica's
   Doppler should read the satellite's `code_doppler`/`carrier_doppler` (the
   value that generated the chunk), not an intermediate per-output estimate.

   **You own more than the Dopplers.** Nothing else in the pipeline touches a
   completed record, so per record and *per signal* your fold is also what
   normalizes the correlator by its sample count, runs the post-correlation
   filter, records the filtered prompt, advances the [bit
   buffer](bit_sync.md) and the [C/N₀ estimator](cn0_estimator.md), and moves
   the record to `last_fully_integrated_*`. Skip that and bit sync, C/N₀ and
   `get_filtered_prompts` are silently dead for every signal, the driver
   included. Two satellite-level duties come with it: the C/N₀ estimators read
   their signal's noise density out of `track_state.noise_estimators` (see
   [Noise Estimator](noise_estimator.md)), and the iteration on which a signal
   first reports bit or secondary-code sync is when the shared `code_phase` has
   to be anchored to that signal's secondary-chip window. The shipped
   `_apply_correlator_output` and `_update_tracked_sat_doppler` do all of this;
   reusing them (below) is much less work than restating them.

   The matching mutating method
   `estimate_dopplers_and_filter_prompt!(track_state, measurements)`
   is what [`track!`](@ref) calls. To support real-time loops, define
   both — the mutating version walks each group's
   `satellites.values::Vector{TrackedSat}` and reassigns slots in place.

   The estimator only needs the per-band **sampling frequency** out of
   `measurements` (to turn each output's `integrated_samples` into an
   integration time and to normalize the DLL discriminator). Both shipped
   estimators therefore also accept a bare per-band sampling-frequency source in
   place of `measurements` — a `NamedTuple`/`Dict` keyed by `get_band_id` — so
   an external correlator producer can run the estimator with no sample buffer
   (see
   [External correlator producers](track.md#External-correlator-producers)).
   A custom estimator that likewise reads only the rate is encouraged to offer
   the same overload.

## Skeleton

The skeleton below is the smallest possible working estimator: every
method returns a constant. Real estimators replace these bodies with
the actual algorithm, but the *structure* — five methods, two structs
— is what the rest of `Tracking.jl` dispatches on.

```jldoctest myestimator
julia> using Tracking, GNSSSignals

julia> using Tracking: AbstractDopplerEstimator, TrackedSat, TrackState,
                       BandMeasurements

julia> # 1. Estimator type — config + any shared state
       struct MyEstimator <: AbstractDopplerEstimator end

julia> # 2. Per-sat state struct
       struct SatMyEstimator end

julia> # 3. Seed each sat
       Tracking.init_estimator_state(::MyEstimator, ::TrackedSat) = SatMyEstimator();

julia> # 4. (Optional) shared-state update on handoff — default returns
       # `est` unchanged, so estimators with no shared state may skip this.
       Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats) = est;

julia> # 5a. Immutable form — walks each group, looks up its band's
       # `BandMeasurement`, returns a fresh `TrackState`. Real implementations
       # read `sat.signals[*].correlator` and compute new dopplers; here
       # we just return the state unchanged.
       function Tracking.estimate_dopplers_and_filter_prompt(
           track_state::TrackState{<:Any, <:MyEstimator},
           measurements::BandMeasurements,
       )
           return track_state
       end;

julia> # 5b. In-place form — what `track!` actually calls. Real
       # implementations write back into `g.satellites.values[i]`.
       function Tracking.estimate_dopplers_and_filter_prompt!(
           track_state::TrackState{<:Any, <:MyEstimator},
           measurements::BandMeasurements,
       )
           return track_state
       end;
```

Plug it into a `TrackState` like any other estimator:

```jldoctest myestimator
julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA(), doppler_estimator = MyEstimator());

julia> track_state = add_satellite!(track_state; prn = 1, code_phase = 0.0, carrier_doppler = 1000.0Hz);

julia> track_state = track(zeros(ComplexF64, 4000), track_state, 4e6Hz);

julia> get_doppler_estimator_state(get_sat_state(track_state, 1))
SatMyEstimator()
```

## Reusing the shipped fold

`src/conventional_pll_and_dll.jl` is worth reading as more than an example: most
of it is estimator-agnostic and is already shared with
[`VectorPLLAndDLL`](@ref), which is written as a per-satellite state type plus
one method. The seams, in the order a chunk passes through them:

  - `_update_tracked_sat_doppler(sat, sampling_frequency, noise, estimator)` —
    one satellite, pure, and the reason the immutable and in-place forms cannot
    drift. It does the per-signal walk, the sync-time code-phase snap and the
    rewrap.
  - `_process_signals(driver, passengers, sat, state, …)` — how the signals feed
    the loops. The generic method folds the passengers for their prompts, bits
    and C/N₀ only and hands the driver to
    `_process_estimator_driver_signal(driver, sat, state, …)`, dispatched on
    your per-satellite state type. Defining that one method is the cheap way in.
  - `_apply_correlator_output` / `_advance_driver_record` — everything one
    record does before a loop sees it, including the effective-bandwidth
    handling for an `N`-code-block integration.

None of these are public API and none carry a compatibility promise, but the
alternative is reimplementing rules — which record's blocks count toward a bit,
when a pre-sync-correlated prompt must be dropped — whose subtleties are the
reason they exist.

## What `signals[1]` means

Both shipped estimators give the [estimator-driver
signal](tracking_state.md#Estimator-driver-signal) the loop *cadence*, the loop
*bandwidths* and the *carrier-phase reference*, and the satellite ranges on its
`code_phase`. That is narrower than "the other signals are passengers": with
`signal_combining = true` every signal's discriminator is folded into the loop
update, weighted as
[Multi-signal discriminator combining](tracking_state.md#Multi-signal-discriminator-combining)
describes, and under [vector tracking](vector_tracking.md) every signal also
hands the navigation filter its own discriminator accumulator. All of it is convention
the shipped estimators chose; yours can use every signal's state any way it
likes.

## What stays generic

You do *not* override the downconvert/correlate path. The per-sat
`update` after correlation reads only estimator-agnostic fields (code/carrier
phase, integrated samples, correlator) and rewraps the existing
`doppler_estimator_state` unchanged. [`add_satellite!`](@ref) and
[`merge_sats`](@ref) call [`init_estimator_state`](@ref) with whatever
estimator the `TrackState` was built with.

## API reference

```@docs
AbstractDopplerEstimator
init_estimator_state
update_estimator_on_handoff
Tracking.estimate_dopplers_and_filter_prompt
Tracking.estimate_dopplers_and_filter_prompt!
```
