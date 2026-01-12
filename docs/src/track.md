# Track

The `track` function is the main entry point for processing GNSS signals. It takes an incoming signal and a `TrackState`, performs downconversion and correlation for all tracked satellites, estimates Doppler frequencies, and returns an updated `TrackState`.

## Function Reference

```@docs
track
```

## Optional Parameters

- `downconvert_and_correlator`: The downconversion and correlation implementation to use. Defaults to `CPUDownconvertAndCorrelator`.
- `intermediate_frequency`: The intermediate frequency of the signal. Defaults to `0.0Hz`.
- `preferred_num_code_blocks_to_integrate`: The preferred number of code blocks to integrate. Defaults to `1`. Will only integrate more than one block once bit synchronization has been achieved.

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
