# Tracking.jl

*Modular tracking algorithm for various Global Navigation Satellite Systems (GNSS)*

Tracking.jl provides a modular and performant tracking algorithm.

```@contents
Pages = ["track.md", "tracking_state.md", "tracking_results.md", "correlator.md", "cn0_estimator.md"]
Depth = 1
```

## Installation

```julia-repl
julia> ]
pkg> add Tracking
```

## Usage

Tracking.jl processes each satellite individually. That means for each satellite a tracking state `TrackingState` must be initialized. From there, the `track` function needs to be called for every state and for every incoming signal to update the tracking parameters.
The required parameters to initialize the tracking state are the GNSS system, carrier doppler and the code phase, e.g.
```julia
state = TrackingState(GPSL1, carrier_doppler, code_phase)
```
These parameters are usually provided by the acquisition process of the satellite. Refer to [`TrackingState`](@ref) to find about other optional parameters.

The signal is tracked by
```julia
results = track(signal, state, prn, sample_frequency)
```
where `prn` is the PRN and `sample_frequency` the sample frequency. Refer to [`track`](@ref) to find about other optional parameters. The result contains the current state as well as some additional information such as the last valid correlator output, found data bits, etc. For each of those parameters a helper function exists to get the parameter (e.g. `get_prompt(results)`) - see [Tracking Results](@ref). The next track function needs the updated state:
```julia
next_results = track(next_signal, get_state(results), prn, sample_frequency)
```

Here is an example for a single PRN:
```julia
using Tracking
using Tracking: Hz, GPSL1
carrier_doppler = 1000Hz
code_phase = 50
sample_frequency = 2.5e6Hz
prn = 1
state = TrackingState(GPSL1, carrier_doppler, code_phase)
results = track(signal, state, prn, sample_frequency)
next_results = track(next_signal, get_state(results), prn, sample_frequency)
```

## Q/A

- Why are the correlator values zero?

The correlator output given in the tracking results is the correlation result after the
code phase has reached the full code length or multiple of the code length depending on the
maximal integration time `max_integration_time` (default: 1ms) you have set. If the current
tracked signal does not include the end of the PRN sequency (or multiples of that), the
correlator output will be zero. Moreover, a correlator output will only be valid, if the
integration time has been at least the miminimum configured integration time
`min_integration_time` (default: 0.75ms). Otherwise the correlator output will be zero as
well. This case occurs at the beginning of the tracking.
