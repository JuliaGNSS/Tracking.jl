# Tracking.jl

This package implements the tracking functionality of GNSS satellites that's part of the larger GNSS receiver.
Tracking.jl primarily consists of two main blocks:

1. Signal down-conversion and correlation
2. Code and carrier estimation to generate replicas and close the loop

Tracking.jl provides defaults for both blocks, but it provides mechanisms to hook in your own implementation (by using multiple dispatch). 
For signal down-conversion and correlation Tracking.jl provides a highly optimized CPU implementation. There is also a GPU implementation that is mainly there for reference, but is not yet as fast as the CPU implementation.
With respect to the second block Tracking.jl provides conventional DLLs and PLLs.
Down-conversion and correlation is done in full code blocks meaning from code start to code end or multiples of that (e.g. in GPS L1 from 0 to `N`*1023). The factor `N` can be specified, but will be `1` as long as the bit start is unknown in order to find the bit start. Once that is done for every tracked satellite the result will be handed over to the code and carrier estimation block.
Moreover, Tracking.jl allows tracking of signals from phased antenna arrays meaning that they are down-converted and correlated by the very same replica to conserve phase relationships.

```@contents
Pages = [
  "track.md",
  "tracking_state.md",
  "tracking_results.md",
  "loop_filter.md",
  "correlator.md",
  "cn0_estimator.md"
]
Depth = 1
```

## Installation

```julia-repl
julia> ]
pkg> add Tracking
```

## Usage

Tracking.jl processes each satellite individually. That means for each satellite a tracking
state `TrackingState` must be initialized. From there, the `track` function needs to be
called for every state and for every incoming signal to update the tracking parameters.
The required parameters to initialize the tracking state are the GNSS system, carrier
Doppler and the code phase, e.g.
```julia
state = TrackingState(GPSL1, carrier_doppler, code_phase)
```
These parameters are usually provided by the acquisition process of the satellite. Refer to [`TrackingState`](@ref) to find about other optional parameters.

The signal is tracked by
```julia
results = track(signal, state, prn, sampling_frequency)
```
where `prn` is the PRN and `sampling_frequency` the sampling frequency. Refer to [`track`](@ref)
to find about other optional parameters. The result contains the current state as well as
some additional information such as the last valid correlator output, found data bits, etc.
For each of those parameters a helper function exists to get the parameter
(e.g. `get_prompt(results)`) - see [Tracking Results](@ref). The next track function needs
the updated state:
```julia
next_results = track(next_signal, get_state(results), prn, sampling_frequency)
```

Here is an example for a single PRN:
```julia
using Tracking
using Tracking: Hz, GPSL1
carrier_doppler = 1000Hz
code_phase = 50
sampling_frequency = 2.5e6Hz
prn = 1
state = TrackingState(GPSL1, carrier_doppler, code_phase)
results = track(signal, state, prn, sampling_frequency)
next_results = track(next_signal, get_state(results), prn, sampling_frequency)
```

## Track multiple signals coherently

Tracking.jl provides a way to track multiple signals coherently, e.g. to track signals from
a phased array. In that case the input `signal` should be a Matrix instead of a Vector,
where the number of rows is equal to the number of antenna elements and the number of
columns is equal to the number of samples. Furthermore, you need to specify the number
of antenna elements to the tracking state:
```julia
state = TrackingState(GPSL1, carrier_doppler, code_phase, num_ants = NumAnts(4))
```
By default the track function will use the first antenna channel as the reference signal to
drive the discriminators etc. However, an appropiate beamforming algorithm will probably
suit better. For that, you'll have to pass a function `post_corr_filter` to the track
function like the following:
```julia
results = track(signal, state, prn, sampling_frequency, post_corr_filter = x -> x[end])
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
