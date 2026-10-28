# Custom Doppler Estimator

The Doppler estimators are TrackingLoops'. [`track!`](@ref) drives any estimator
through TrackingLoops' estimator interface, so
[`ConventionalPLLAndDLL`](@extref TrackingLoops.ConventionalPLLAndDLL),
[`ConventionalAssistedPLLAndDLL`](@extref TrackingLoops.ConventionalAssistedPLLAndDLL-Union{Tuple{}, Tuple{Type{CO}}, Tuple{CO}} where CO<:TrackingLoopFilters.AbstractLoopFilter),
[`NCOReferencedPLLAndDLL`](@extref TrackingLoops.NCOReferencedPLLAndDLL), the vector
loop [`VectorPLLAndDLL`](@extref TrackingLoops.VectorPLLAndDLL) and your own (a Kalman
filter, a joint-channel estimator, …) are interchangeable as `doppler_estimator`.

## Where state lives

A [`TrackedSat`](@ref) carries everything per satellite: the shared carrier/code
Doppler and phase, one [`TrackedSignal`](@ref) per signal, and the estimator's
per-satellite state in `doppler_estimator_state`. The estimator object
(`<: AbstractDopplerEstimator`) holds configuration plus any state *shared* across
satellites. The signal path (downconvert, correlate) is estimator-agnostic and passes
`doppler_estimator_state` through unchanged.

## What you implement

1. **An estimator type** subtyping
   [`AbstractDopplerEstimator`](@extref TrackingLoops.AbstractDopplerEstimator):
   configuration and shared state, no per-satellite state.

2. **A per-satellite state struct**, e.g. `SatConventionalPLLAndDLL` (loop filters and
   seed Dopplers) or a Kalman state vector.

3. **An [`init_estimator_state`](@ref) method** seeding a satellite entering the
   track set:

   ```julia
   TrackingLoops.init_estimator_state(::MyEstimator, sat::TrackedSat) = MyPerSatState(...)
   ```

   It must be **pure**: it is also called on the PRN-0 template sat that fixes a
   group's slot type, as a type probe when validating pre-built sats, and by
   [`reset_loop_filters!`](@ref). Register sats into shared state in
   `update_estimator_on_handoff`, never here.

4. **(Optional) An [`update_estimator_on_handoff`](@ref) method**, if shared state
   has to change when sats join (joint-channel covariance, normalization terms, …):

   ```julia
   Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats) = ...
   ```

   It is called once per [`add_satellite!`](@ref) / [`merge_sats`](@ref) call with
   the dictionary of incoming sats, *after* per-sat seeding. The default returns
   `est` unchanged.

   The returned estimator must keep the input's **concrete type**
   ([`TrackState`](@ref) is parameterized on it): grow shared state in place in
   resizable containers, or rebuild scalar fields with `Setfield.@set` or a copying
   constructor. Since `TrackState` is immutable, a rebuilt estimator only reaches you
   through the return value: write `track_state = add_satellite!(track_state; ...)`.
   For in-place estimators the returned `TrackState` is the input itself.

5. **The per-record step**, in one of two ways.

   **A `step_loop` method** (TrackingLoops' per-record interface) is all `track!`
   needs. Every completed record of every signal of the satellite goes through
   `TrackingLoops.step_loop(estimator, state, record, words, landing_sample)` with
   the satellite's one state, the estimator-driver signal's (`signals[1]`) first. It
   returns the new state and carrier and code Dopplers; the satellite takes the
   driver's, as another signal's record returns the command in force. See
   [`LoopRecord`](@extref TrackingLoops.LoopRecord) for what a record carries.
   Tracking does the rest: the bit buffer, the C/N₀ estimators, the phase snap at
   sync. A `TrackingLoops.reset_estimator_state(estimator, state, carrier_doppler,
   code_doppler)` method chooses what [`reset_loop_filters!`](@ref) keeps; without
   one, the state is re-initialized. Such an estimator also runs unchanged in a
   hardware correlator's loop process.

   **An `estimate_dopplers_and_filter_prompt` method** on
   `TrackState{<:Any, <:MyEstimator}` instead replaces the whole estimate phase, for
   an estimator that does not fit a per-record step (e.g. one stepping all satellites
   jointly). It runs once per **chunk** (see
   [Chunked Doppler updates](track.md#Chunked-Doppler-updates)), walks
   `track_state.groups`, reads each band's `BandMeasurement` via
   `get_band_id(group.band)`, and produces `TrackedSat`s with new
   `carrier_doppler`/`code_doppler` and per-sat state.

   Each signal's records of the chunk are in its
   `correlator_outputs::Vector{`[`CorrelatorOutput`](@extref TrackingLoops.CorrelatorOutput)`}`
   (zero, one or several), each with the raw correlator, its sample count, its end
   sample index (from the current measurement's first sample) and the replica's code
   phase there. Fold over them in order, write the NCO Doppler once (typically from
   the last), and empty `correlator_outputs`. The NCO Doppler is fixed across the
   chunk, so a discriminator needing the replica's Doppler reads the satellite's
   `code_doppler`/`carrier_doppler`.

   [`track!`](@ref) calls the mutating
   `estimate_dopplers_and_filter_prompt!(track_state, measurements)`, which
   reassigns each group's `satellites.values` slots in place; define both.

   The generic estimate phase also accepts a `NamedTuple`/`Dict` of per-band
   sampling frequencies (keyed by `get_band_id`) in place of `measurements`, so an
   external correlator producer needs no sample buffer (see
   [External correlator producers](track.md#External-correlator-producers)). A
   custom estimator that reads only the rate should offer the same overload.

## Skeleton

The smallest working estimator that replaces the estimate phase, with constant
bodies; the structure (two structs, four methods) is what Tracking dispatches on.

```jldoctest myestimator
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking:
           TrackedSat,
           TrackState,
           SignalGroup,
           BandMeasurements,
           get_band_id

julia> using TrackingLoops: AbstractDopplerEstimator

julia> # 1. Estimator type — config + any shared state
       struct MyEstimator <: AbstractDopplerEstimator end

julia> # 2. Per-sat state struct
       struct SatMyEstimator end

julia> # 3. Seed each sat
       TrackingLoops.init_estimator_state(::MyEstimator, ::TrackedSat) = SatMyEstimator();

julia> # 4. (Optional) shared-state update on handoff
       Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats) = est;

julia> # 5a. Immutable form: a real one folds each signal's `correlator_outputs`
       # into new Dopplers and returns a fresh `TrackState`.
       function Tracking.estimate_dopplers_and_filter_prompt(
           track_state::TrackState{<:Any, <:MyEstimator},
           measurements::BandMeasurements,
       )
           return track_state
       end;

julia> # 5b. In-place form, what `track!` calls: writes `g.satellites.values[i]`.
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

The generic estimate phase in `src/conventional_pll_and_dll.jl` shows the full
pattern: the immutable form delegates to the in-place one, and the walk steps the
[estimator-driver signal](tracking_state.md#Estimator-driver-signal) (`signals[1]`,
whose records close the loops) before the passengers. That split is a convention of
the per-record path only.

## What stays generic

You do *not* override the downconvert/correlate path; it passes
`doppler_estimator_state` through unchanged. [`add_satellite!`](@ref) and
[`merge_sats`](@ref) call [`init_estimator_state`](@ref) with the `TrackState`'s
estimator.

## API reference

```@docs
init_estimator_state
update_estimator_on_handoff
Tracking.estimate_dopplers_and_filter_prompt
Tracking.estimate_dopplers_and_filter_prompt!
```

In the TrackingLoops manual:

- [`AbstractDopplerEstimator`](@extref TrackingLoops.AbstractDopplerEstimator)
- [`reset_estimator_state`](@extref TrackingLoops.reset_estimator_state-Tuple{ConventionalPLLAndDLL, TrackingLoops.SatConventionalPLLAndDLL, Any, Any})

### One record through a loop

The device-independent step lives in TrackingLoops, so the software receiver and a
hardware correlator's loop process run the same arithmetic. `track!` passes the
chunk's own replica word and no landing sample.

In the TrackingLoops manual:

- [`step_loop`](@extref TrackingLoops.step_loop-Tuple{ConventionalPLLAndDLL, TrackingLoops.SatConventionalPLLAndDLL, LoopRecord, Any, Int64})
- [`LoopRecord`](@extref TrackingLoops.LoopRecord)
- [`SatConventionalPLLAndDLL`](@extref TrackingLoops.SatConventionalPLLAndDLL)

### The delay-aware loop

`NCOReferencedPLLAndDLL` is for a correlator whose NCO word takes effect records
after the one that motivated it, as in hardware. Given a fixed word and no landing
sample, as `track!` supplies, it is the conventional loop to the bit at the same
bandwidths; by default it pulls in at the narrow 18 Hz, as a wider loop tolerates
less command delay.

In the TrackingLoops manual:

- [`NCOReferencedPLLAndDLL`](@extref TrackingLoops.NCOReferencedPLLAndDLL)
- [`SatNCOReferencedPLLAndDLL`](@extref TrackingLoops.SatNCOReferencedPLLAndDLL)
- [`NCOTimeline`](@extref TrackingLoops.NCOTimeline)
- [`scheduled_words`](@extref TrackingLoops.scheduled_words-Tuple{NCOTimeline})
- [`mean_nco_word`](@extref TrackingLoops.mean_nco_word-Tuple{NCOTimeline, Real, Real})
- [`FixedNCOWord`](@extref TrackingLoops.FixedNCOWord)

### The vector loop

`track!` drives `VectorPLLAndDLL(signals...)` through `step_loop` like any other
estimator. Inside, its navigation engine decodes every satellite's bits (of a plain
signal, or of the data component of a pilot + data pair, e.g.
`VectorPLLAndDLL(GalileoE1C() => GalileoE1B())` with
`signals = (GalileoE1C(), GalileoE1B())`), solves the PVT, and from its first fix
closes all loops with a navigation filter. Read `navigation_solution(estimator)` and
`navigation_status(estimator)`.

```julia
estimator = VectorPLLAndDLL(GPSL1CA())
track_state = TrackState(; signal = GPSL1CA(), doppler_estimator = estimator)
for chunk in chunks
    track!(chunk, track_state, sampling_frequency)
end
navigation_solution(estimator)
```

Each record gives it the satellite (`prn`), the replica's code phase at the
record's end, and a time grid shared by a band's satellites: the estimate phase adds
the band's running offset (`TrackState`'s `sample_offsets`, advanced by every
`track!`). Step every satellite at least every half `vt_update_interval` (50 ms by
default), which any `track!` on a shorter chunk does.

The engine and the offsets belong together: a `TrackState` built afresh around an
estimator that has run starts every band at `0`, so pass it the old state's
`sample_offsets` (see [`TrackState`](@ref)). To branch a state, `deepcopy` it; [`track`](@ref) shares the
estimator with its input.

In the vector loop the navigation filter alone sets a satellite's carrier NCO; a
residual PLL follows the carrier on the prompts without steering it, and
`carrier_loop_stage` reports `VECTOR_LOOP`. [`reset_loop_filters!`](@ref)
re-seeds the inner loop and keeps the satellite in the vector loop: it drops the
code correction (the converged Dopplers contain it) until the filter's next cycle,
and keeps the carrier command (they do not).

In the TrackingLoops manual:

- [`VectorPLLAndDLL`](@extref TrackingLoops.VectorPLLAndDLL)
- [`SatVectorPLLAndDLL`](@extref TrackingLoops.SatVectorPLLAndDLL)
- [`LoopRecord`](@extref TrackingLoops.LoopRecord)
