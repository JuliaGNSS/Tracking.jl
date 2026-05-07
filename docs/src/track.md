# Track

The `track` function is the main entry point for processing GNSS signals. It takes an incoming signal and a `TrackState`, performs downconversion and correlation for all tracked satellites, estimates Doppler frequencies, and returns an updated `TrackState`.

For hard real-time SDR loops where GC pauses must be avoided, an in-place
counterpart [`track!`](@ref) is provided that mutates `TrackState` instead
of rebuilding new immutable wrappers. After a one-time call to
[`prewarm!`](@ref), each individual stage of `track!` is allocation-free in
steady state.

## Function Reference

```@docs
track
track!
prewarm!
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
prewarm!(track_state, max_prompts_per_track_call)  # one-time

while got_signal_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_freq;
        downconvert_and_correlator = dc)
    # ... use track_state.sat_states ...
end
```

`prewarm!` calls `sizehint!` on each satellite's `filtered_prompts` buffer
so the geometric `push!` reallocations don't fire during the hot loop.
`max_prompts_per_track_call` is the upper bound on how many integrations
can complete in one `track!` call (a few for short chunks).

## Downconversion and Correlation

```@docs
CPUDownconvertAndCorrelator
AbstractDownconvertAndCorrelator
```

## Correlator Sample Shifts

```@docs
get_correlator_sample_shifts
get_early_late_sample_spacing
```
