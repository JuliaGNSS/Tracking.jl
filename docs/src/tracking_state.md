# Tracking State

Tracking.jl uses a small hierarchy of state types to manage tracking across multiple satellites, multiple signals per satellite, and multiple signal groups.

## TrackState

The main container for all tracking state. It holds:

- `satellites` — a `NamedTuple` of `Dictionary{Int, TrackedSat}`, one per signal group.
- `signal_groups` — the per-group signal-instance tuple, used by [`add_satellite!`](@ref) to build new sats.
- `doppler_estimator` — the Doppler estimator configuration (e.g. PLL/DLL bandwidths).

```@docs
TrackState
```

## TrackedSat

Per-satellite tracking state. Carries the shared carrier and code Doppler and phase (one set of values per satellite, since all signals on a satellite share the same carrier), a tuple of per-signal [`TrackedSignal`](@ref) states, and the per-satellite Doppler-estimator state.

```@docs
TrackedSat
```

### TrackedSat accessor functions

These return per-satellite state. When invoked on a `TrackState` with `(track_state, group, prn)` arguments, they look up the sat first.

- `get_prn(sat)` — the PRN number.
- `get_code_phase(sat)` — the shared code phase (wraps at [`max_code_length`](@ref)).
- `get_code_doppler(sat)` — the shared code Doppler.
- `get_carrier_phase(sat)` — the shared carrier phase in radians.
- `get_carrier_doppler(sat)` — the shared carrier Doppler.
- `get_signal_start_sample(sat)` — index of the next sample to integrate.
- `get_signals(sat)` — the tuple of [`TrackedSignal`](@ref)s.
- `get_doppler_estimator_state(sat)` — the per-satellite Doppler-estimator state (e.g. the loop-filter integrators for the conventional PLL/DLL).

## TrackedSignal

Per-signal tracking state. Carries the correlator, post-correlation filter, CN0 estimator, bit buffer, and integration-progress flags for a single signal on a satellite. A multi-signal `TrackedSat` (e.g. one tracking GPS L1 C/A + L1C-D + L1C-P) carries one `TrackedSignal` per signal in `signals::Tuple{Vararg{TrackedSignal}}`.

The first signal in the tuple is the Doppler source — its correlator is what the PLL/DLL discriminator runs on.

```@docs
TrackedSignal
```

### TrackedSignal accessor functions

When called on a single-signal `TrackedSat`, these accessors forward to the sole `TrackedSignal`. For multi-signal sats, address a specific signal explicitly via `sat.signals[i]` (or `get_signals(sat)[i]`).

- `get_signal(tsig)` — the GNSS signal instance (e.g. `GPSL1CA()`).
- `get_correlator(tsig)` — the working correlator (the in-flight accumulator for the current integration).
- `get_last_fully_integrated_correlator(tsig)` — the correlator value from the last completed integration period.
- `get_last_fully_integrated_filtered_prompt(tsig)` — the filtered prompt value from the last completed integration period.
- `get_filtered_prompts(tsig)` — every filtered prompt produced during the most recent `track` call. The vector is reset at the start of each call and appended for every completed integration.
- `get_bit_buffer(tsig)` / `get_bits(tsig)` / `get_num_bits(tsig)` — the bit buffer, decoded bits, and bit count.
- `get_integrated_samples(tsig)` — number of samples accumulated into the current integration so far.
- `has_bit_or_secondary_code_been_found(tsig)` — `true` once bit or secondary-code synchronization has been achieved.

### Code-phase wrap period

The shared `TrackedSat.code_phase` wraps at the longest code period across all signals, including secondary code. For a sat tracking L1C-P (10230 primary chips × 1800 secondary chips ≈ 18.4 M chips ≈ 18 s at 1.023 Mcps), the wrap is 18.4 M; for a sat tracking only L1 C/A (1023 chips × 1) it's 1023. The helper [`max_code_length`](@ref) returns the wrap period for any signal tuple at compile time.

```@docs
max_code_length
```

## Capability accessors

```@docs
SatelliteDicts
get_sat_states
get_sat_state
get_system
```

## Adding and removing satellites

```@docs
add_satellite!
add_satellite
remove_satellite!
remove_satellite
merge_sats
```
