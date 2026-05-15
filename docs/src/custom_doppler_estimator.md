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
   Called once per satellite when a sat enters the track set
   (acquisition → tracking handoff), it produces the initial per-sat
   state for that sat:

   ```julia
   Tracking.init_estimator_state(::MyEstimator, sat::TrackedSat) = MyPerSatState(...)
   ```

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

5. **An `estimate_dopplers_and_filter_prompt` method** dispatched on
   `TrackState{<:..., <:..., <:MyEstimator}`. This is where the actual
   update logic runs, once per integration completion. It walks each
   capability's `Dictionary{Int, TrackedSat}` and produces new
   `TrackedSat`s with updated `carrier_doppler`/`code_doppler` and
   updated per-sat estimator state.

   The matching mutating method
   `estimate_dopplers_and_filter_prompt!(track_state, ...)` is what
   [`track!`](@ref) calls. To support real-time loops, define both — the
   mutating version walks `sats_dict.values` and reassigns slots in
   place.

## Skeleton

```julia
using Tracking
using Tracking: AbstractDopplerEstimator, TrackedSat, TrackedSignal, TrackState
using Tracking: init_estimator_state, update_estimator_on_handoff,
    estimate_dopplers_and_filter_prompt

# 1. Estimator type — config + any shared state
struct MyEstimator <: AbstractDopplerEstimator
    bandwidth_hz::Float64
    # ... shared state held in resizable containers, e.g. Vector{Float64}
    # for a growing joint-state vector — mutated in place on handoff ...
end

# 2. Per-sat state struct
struct SatMyEstimator
    # ... whatever your algorithm tracks per sat ...
    integrator::Float64
end

# 3. Seed each sat with its initial estimator state
function Tracking.init_estimator_state(est::MyEstimator, sat::TrackedSat)
    SatMyEstimator(0.0)
end

# 4. (Optional) update shared state when sats join the track set
function Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats)
    # `new_sats` is a `Dictionary{I, TrackedSat}` of just the incoming sats.
    # Grow the shared joint-state container in place (immutable estimator,
    # resizable field) and return the same `est`.
    for _ in 1:length(new_sats)
        push!(est.joint_state, 0.0)
    end
    return est
end

# 5. The update step (immutable form). Walks each capability's dict and
# produces new TrackedSats. The signals tuple is rebuilt per sat with
# updated per-signal state (cleared correlator, advanced bit buffer, etc.).
function Tracking.estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:Any, <:Any, <:MyEstimator},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    new_sats = map(track_state.satellites) do sats_dict
        map(sats_dict) do sat
            de_state = sat.doppler_estimator_state
            # ... compute new carrier_doppler, code_doppler, de_state ...
            new_signals = my_update_signals(sat.signals, ...)
            TrackedSat(
                sat.prn,
                sat.code_phase, sat.code_doppler,
                sat.carrier_phase, sat.carrier_doppler,
                sat.signal_start_sample, new_signals,
                new_de_state,
            )
        end
    end
    return TrackState(track_state; satellites = new_sats)
end
```

For full real-time support (zero-allocation steady state), additionally
provide an in-place version that walks each `sats_dict.values::Vector{TrackedSat}`
and reassigns slots:

```julia
function Tracking.estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:Any, <:Any, <:MyEstimator},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    for sats_dict in track_state.satellites
        vals = sats_dict.values
        @inbounds for i in eachindex(vals)
            vals[i] = my_per_sat_update(vals[i], ...)
        end
    end
    return track_state
end
```

The existing [`ConventionalPLLAndDLL`](@ref) implementation in
`src/conventional_pll_and_dll.jl` shows the full pattern, including how
the immutable and in-place forms share a `_update_tracked_sat_doppler`
helper so they cannot drift, and how the per-signal walk distinguishes
the Doppler-source signal (`signals[1]`, which drives the PLL/DLL) from
the other signals (which only have their prompts filtered).

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
```
