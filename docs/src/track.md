# Track

The `track` function is the main entry point for processing GNSS signals. It takes an incoming signal and a `TrackState`, performs downconversion and correlation for all tracked satellites, estimates Doppler frequencies, and returns an updated `TrackState`.

For hard real-time SDR loops where GC pauses must be avoided, an in-place
counterpart [`track!`](@ref) is provided that mutates `TrackState` instead
of rebuilding new immutable wrappers. After one warmup call (to seat the
`filtered_prompts` buffer's capacity), `track!` runs with only a small
per-call residual that does not trigger GC pauses.

> **Allocation-free hot loop:** the
> [`CPUDownconvertAndCorrelator`](@ref) and [`CPUThreadedDownconvertAndCorrelator`](@ref)
> structs hold long-lived per-thread scratch buffers that grow on
> first use. **Construct the correlator once outside the `track!`
> loop** and pass it via `downconvert_and_correlator =` — the default
> kwarg value rebuilds the buffer matrix on every call.

## Function Reference

```@docs
track
track!
```

## Optional Parameters

- `downconvert_and_correlator`: The downconversion and correlation implementation to use. Defaults to `CPUDownconvertAndCorrelator`.
- `intermediate_frequency`: The intermediate frequency of the signal. Defaults to `0.0Hz`.
- `preferred_num_code_blocks_to_integrate`: The preferred number of code blocks to integrate. Defaults to `1`. Will only integrate more than one block once bit synchronization has been achieved.

## Real-time use

`track!` writes back into the existing `Vector{TrackedSat}` slots of each
system, so the tracking loop can run without producing GC pressure once the
sat set is steady. A typical setup looks like:

```julia
track_state = TrackState(system, initial_sats)
dc = CPUThreadedDownconvertAndCorrelator()  # hoist outside the loop

while got_signal_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_freq;
        downconvert_and_correlator = dc)
    # ... use track_state.sat_states ...
end
```

Hoisting `dc` is what makes the loop allocation-free: each correlator holds
long-lived per-thread scratch buffers that grow on first use and are reused
thereafter. Calling `track!` without the kwarg works, but rebuilds the
buffer matrix on every call.

The first `track!` call may also grow each satellite's `filtered_prompts`
buffer via `push!` reallocations; from the second call onwards the capacity
is settled and the loop allocates only the small `bit_buffer` find-bit
closure residual.

## Downconversion and Correlation

```@docs
CPUDownconvertAndCorrelator
CPUThreadedDownconvertAndCorrelator
AbstractDownconvertAndCorrelator
```

## Correlator Sample Shifts

```@docs
get_correlator_sample_shifts
get_early_late_sample_spacing
```
