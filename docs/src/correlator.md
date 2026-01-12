# Correlator

The correlator computes the correlation between the incoming signal and the locally generated replica code at multiple code phase offsets (early, prompt, late, etc.).

## Default Correlators

The default correlator depends on the GNSS system and is returned by `get_default_correlator`:
- `EarlyPromptLateCorrelator` for GPS L1, GPS L5 (BPSK modulation)
- `VeryEarlyPromptLateCorrelator` for Galileo E1B (BOC modulation)

```@docs
get_default_correlator
```

## Correlator Types

```@docs
EarlyPromptLateCorrelator
VeryEarlyPromptLateCorrelator
```

## Accessing Correlator Values

```@docs
get_early
get_prompt
get_late
get_accumulators
get_num_accumulators
```

## Custom Correlators

You can implement your own correlator by creating a subtype of `AbstractCorrelator`
and implementing the required functions. See `src/correlators/correlator.jl` for the
interface that needs to be implemented:

- `get_accumulators(correlator)`: Return the accumulator values
- `get_num_accumulators(correlator)`: Return the number of accumulators
- `update_accumulator(correlator, accumulators)`: Create a new correlator with updated accumulators
- `get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)`: Return the sample shifts for each accumulator
