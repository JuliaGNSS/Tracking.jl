# Track

The `track` function is the main entry point for processing GNSS signals. It takes an incoming signal and a `TrackState`, performs downconversion and correlation for all tracked satellites, estimates Doppler frequencies, and returns an updated `TrackState`.

For hard real-time SDR loops where GC pauses must be avoided, an in-place
counterpart [`track!`](@ref) is provided that mutates `TrackState` instead
of rebuilding new immutable wrappers. After one warmup call (to seat the
`filtered_prompts` buffer's capacity), the single-threaded `track!` path
is **fully allocation-free**; the threaded path keeps a small irreducible
per-call residual from Polyester's `@batch` closure capture (160 B per
GNSS system).

!!! tip "Hoist the correlator outside the hot loop"
    Both [`CPUDownconvertAndCorrelator`](@ref) and
    [`CPUThreadedDownconvertAndCorrelator`](@ref) hold long-lived per-thread
    scratch buffers that grow on first use and are reused thereafter.
    **Construct the correlator once outside the `track`/`track!` loop**
    and pass it via the `downconvert_and_correlator =` keyword argument —
    relying on the default kwarg value rebuilds the scratch buffers on
    every call, which defeats the allocation-free design. This applies
    to both `track` (immutable) and `track!` (in-place); both default to
    `CPUThreadedDownconvertAndCorrelator()`.

## Function Reference

```@docs
track
track!
```

## Optional Parameters

- `downconvert_and_correlator`: The downconversion and correlation implementation to use. Defaults to `CPUThreadedDownconvertAndCorrelator()`. **For real-time loops, hoist this outside the loop** (see the tip above).
- `intermediate_frequency`: The intermediate frequency of the signal. Defaults to `0.0Hz`. Only accepted on the bare-buffer form `track!(buf, state, fs; intermediate_frequency = ...)`; for the [`Measurement`](@ref) and multi-band forms the IF lives on each `Measurement`.
- `preferred_num_code_blocks_to_integrate`: The preferred number of code blocks to integrate. Defaults to `1`. Will only integrate more than one block once bit synchronization has been achieved.

## Measurement

One band's incoming sample buffer plus the front-end metadata needed to process it. Bundles `samples` with `sampling_frequency` and `intermediate_frequency` — these are inseparable in practice, and the bundle removes the chance of mismatched parallel NamedTuples in a multi-band `track` call.

For multi-band tracking, build one `Measurement` per band and pass them as a NamedTuple keyed by [`band_key`](@ref):

```julia
track!((l1 = Measurement(buf_l1, 4e6Hz),
        l5 = Measurement(buf_l5, 25e6Hz)), track_state)
```

For single-band tracking, the bare-buffer form `track!(buf, state, fs)` and the single-`Measurement` form `track!(Measurement(buf, fs), state)` both auto-wrap into a one-entry NamedTuple internally.

```@docs
Measurement
Measurements
```

## Real-time use

`track!` writes back into the existing `Vector{TrackedSat}` slots of each
per-group dictionary, so the tracking loop can run without producing
GC pressure once the sat set is steady. A typical setup looks like:

```julia
track_state = TrackState(; signal = GPSL1CA())
add_satellite!(track_state; prn = 1, code_phase = 0.0, carrier_doppler = 1000.0Hz)
# ... more sats ...

dc = CPUThreadedDownconvertAndCorrelator()  # hoist outside the loop

while got_signal_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_freq;
        downconvert_and_correlator = dc)
    # ... read per-sat state e.g. via get_sat_state(track_state, group, prn) ...
end
```

Hoisting `dc` is what makes the loop allocation-free: each correlator holds
long-lived per-thread scratch buffers that grow on first use and are reused
thereafter. Calling `track!` without the kwarg works, but rebuilds the
buffers on every call.

The first `track!` call may also grow each signal's `filtered_prompts`
buffer via `push!` reallocations; from the second call onwards the capacity
is settled.

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
