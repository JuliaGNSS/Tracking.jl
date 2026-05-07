# Custom Doppler Estimator

Tracking.jl ships [`ConventionalPLLAndDLL`](@ref) (and the FLL-assisted
variant [`ConventionalAssistedPLLAndDLL`](@ref)), but you can plug in a
different Doppler-estimation algorithm — e.g. a Kalman filter or a
joint-channel estimator — by implementing a small set of methods.

## Where state lives

The wrapper [`TrackedSat`](@ref) carries each `SatState` together with a
per-satellite estimator state. The estimator object itself
(`<: AbstractDopplerEstimator`) is configuration plus any *shared* state
that spans satellites or systems. Cleanly separating these makes the
signal-path code (downconvert, correlate, sample-bookkeeping) estimator-
agnostic — it touches `tracked_sat.sat_state` and rewraps the estimator
state unchanged.

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
   Tracking.init_estimator_state(::MyEstimator, sat_state::SatState) = MyPerSatState(...)
   ```

4. **(Optional) An [`update_estimator_on_handoff`](@ref) method**, if
   your estimator carries cross-satellite or cross-system shared state
   that has to change when sats join the track set (joint-channel
   covariance, per-batch normalization terms, …):

   ```julia
   Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats) = ...
   ```

   It is called once per [`merge_sats`](@ref) call with the dictionary
   of incoming `SatState`s and runs *after* per-sat seeding, so
   `init_estimator_state` sees the pre-update shared state. The default
   returns `est` unchanged, so estimators with no shared state need not
   implement it.

   **Type constraint:** the returned estimator must have the same
   concrete type as the input. [`TrackState`](@ref) is parameterized on
   the estimator type, and changing it would break inference. Keep the
   estimator an immutable `struct` and put any growing shared state in
   resizable containers (e.g. `Vector`, `Matrix`) you `push!`/`resize!`
   in place; if scalar fields need replacing, rebuild the estimator
   with `Setfield.@set` or a copying constructor.

5. **An `estimate_dopplers_and_filter_prompt` method** dispatched on
   `TrackState{<:..., <:MyEstimator}`. This is where the actual update
   logic runs, once per integration completion. It walks each system's
   `Dictionary{Int, TrackedSat}` and produces new `TrackedSat`s with
   updated `sat_state.carrier_doppler`/`code_doppler` and updated per-sat
   estimator state.

   The matching mutating method
   `estimate_dopplers_and_filter_prompt!(track_state, ...)` is what
   [`track!`](@ref) calls. To support real-time loops, define both — the
   mutating version walks `system_sats_state.states.values` and
   reassigns slots in place.

## Skeleton

```julia
using Tracking
using Tracking: AbstractDopplerEstimator, TrackedSat, SatState, TrackState
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
function Tracking.init_estimator_state(est::MyEstimator, sat::SatState)
    SatMyEstimator(0.0)
end

# 4. (Optional) update shared state when sats join the track set
function Tracking.update_estimator_on_handoff(est::MyEstimator, new_sats)
    # `new_sats` is a `Dictionary{I, SatState}` of just the incoming sats.
    # Grow the shared joint-state container in place (immutable estimator,
    # resizable field) and return the same `est`.
    for _ in 1:length(new_sats)
        push!(est.joint_state, 0.0)
    end
    return est
end

# 5. The update step (immutable form)
function Tracking.estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:Any, <:MyEstimator},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    new_mss = map(track_state.multiple_system_sats_state) do sss
        new_sats = map(sss.states) do tracked_sat
            sat = tracked_sat.sat_state
            est = tracked_sat.estimator_state
            # ... compute new carrier_doppler, code_doppler, est ...
            new_sat = SatState(sat; carrier_doppler = ..., code_doppler = ...)
            new_est = SatMyEstimator(est.integrator + ...)
            TrackedSat(new_sat, new_est)
        end
        Tracking.SystemSatsState(sss, new_sats)
    end
    return TrackState(track_state; multiple_system_sats_state = new_mss)
end
```

For full real-time support (zero-allocation steady state), additionally
provide an in-place version that walks
`sss.states.values::Vector{TrackedSat}` and reassigns slots:

```julia
function Tracking.estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:Any, <:MyEstimator},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    for sss in track_state.multiple_system_sats_state
        vals = sss.states.values
        @inbounds for i in eachindex(vals)
            vals[i] = my_per_sat_update(vals[i], sss.system, ...)
        end
    end
    return track_state
end
```

The existing [`ConventionalPLLAndDLL`](@ref) implementation in
`src/conventional_pll_and_dll.jl` shows the full pattern, including how
the immutable and in-place forms share a `_update_tracked_sat_doppler`
helper so they cannot drift.

## What stays generic

You do *not* override the downconvert/correlate path. The per-sat
`update` after correlation reads only estimator-agnostic fields
(code/carrier phase, integrated samples, correlator) and rewraps the
existing `tracked_sat.estimator_state` unchanged. Acquisition and
[`merge_sats`](@ref) already call [`init_estimator_state`](@ref) with
whatever estimator the `TrackState` was built with.

## API reference

```@docs
AbstractDopplerEstimator
init_estimator_state
update_estimator_on_handoff
```
