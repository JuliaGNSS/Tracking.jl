# Tracking State

Tracking.jl uses a hierarchy of state types to manage tracking across multiple satellites and GNSS systems.

## TrackState

The main container for all tracking state. It holds satellite states for one or more GNSS systems and the Doppler estimator (PLL/DLL).

```@docs
TrackState
```

## SatState

Per-satellite tracking state containing phase, Doppler, correlator, and bit buffer information.

```@docs
SatState
```

### SatState Accessor Functions

```@docs
get_prn
get_code_phase
get_code_doppler
get_carrier_phase
get_carrier_doppler
get_integrated_samples
get_signal_start_sample
get_correlator
get_last_fully_integrated_correlator
get_last_fully_integrated_filtered_prompt
get_bit_buffer
get_bits
get_num_bits
has_bit_or_secondary_code_been_found
```

## TrackedSat

Per-satellite wrapper that holds a [`SatState`](@ref) together with the
estimator-specific per-satellite state. The `Vector{TrackedSat}` backing
storage of [`TrackedSystem`](@ref) is what [`track!`](@ref) mutates
in place — the `TrackedSat` values themselves remain immutable.

```@docs
TrackedSat
get_estimator_state
```

## TrackedSystem

Container for multiple satellites of the same GNSS system.

```@docs
TrackedSystem
```

## Multi-System Support

```@docs
TrackedSystems
get_tracked_system
get_sat_states
get_sat_state
get_system
```

## Modifying Tracking State

```@docs
merge_sats
filter_out_sats
```
