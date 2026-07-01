# Tracking.jl

This package implements the tracking functionality of GNSS satellites that's part of the larger GNSS receiver.
Tracking.jl primarily consists of two main blocks:

1. Signal down-conversion and correlation
2. Code and carrier estimation to generate replicas and close the loop

Tracking.jl provides defaults for both blocks, but it provides mechanisms to hook in your own implementation (by using multiple dispatch).
For signal down-conversion and correlation Tracking.jl provides a highly optimized CPU implementation using hand-tuned SIMD intrinsics.
With respect to the second block Tracking.jl provides conventional DLLs and PLLs with FLL-assisted carrier tracking as the default.

Down-conversion and correlation is done in full code blocks meaning from code start to code end or multiples of that (e.g. in GPS L1 C/A from 0 to `N`*1023).
The factor `N` can be specified, but will be `1` as long as the bit start is unknown — Tracking.jl uses single-code-period integrations to locate the bit edge. Once the bit start is known for every tracked satellite, longer coherent integrations become available and the result is handed over to the code-and-carrier estimation block.

Moreover, Tracking.jl allows tracking of signals from phased antenna arrays meaning that they are down-converted and correlated by the very same replica to conserve phase relationships.
Multi-signal tracking is supported: a single satellite can be tracked on several signals at once (e.g. GPS L1 C/A together with L1C-D and L1C-P) sharing one carrier downconvert per outer iteration, with a per-signal correlator each.

## Supported signals

- GPS L1 C/A
- GPS L1C-D (data)
- GPS L1C-P (pilot, with 1800-chip overlay-code sync)
- GPS L2CM (data)
- GPS L2CL (pilot, dataless — no secondary code)
- GPS L5I
- GPS L5Q (pilot, with NH20 secondary-code sync)
- Galileo E1B (and its BOC(1,1) approximation)
- Galileo E1C (pilot, with CS25 secondary-code sync; and its BOC(1,1) approximation)
- Galileo E5a-I (data, with CS20 secondary-code sync)
- Galileo E5a-Q (pilot, with per-PRN CS100 secondary-code sync)

```@contents
Pages = [
  "track.md",
  "tracking_state.md",
  "bit_sync.md",
  "loop_filter.md",
  "custom_doppler_estimator.md",
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

## Quick start

The minimum to track a single GPS L1 C/A satellite: build a [`TrackState`](@ref), seed it from an acquisition handoff via [`add_satellite!`](@ref), then call [`track`](@ref) on each measurement.

```jldoctest quickstart; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_code_frequency, get_code_center_frequency_ratio

julia> track_state = TrackState(; signal = GPSL1CA());

julia> track_state = add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);

julia> sampling_frequency = 4e6Hz;

julia> num_samples = 4000;

julia> signal = GPSL1CA();

julia> code_frequency = 1000.0Hz * get_code_center_frequency_ratio(signal) + get_code_frequency(signal);

julia> measurement = cis.(2π .* 1000.0Hz .* (0:num_samples-1) ./ sampling_frequency) .*
           gen_code(num_samples, signal, 1, sampling_frequency, code_frequency, 50.0);

julia> track_state = track(measurement, track_state, sampling_frequency);

julia> get_carrier_doppler(track_state, 1)
999.9999883655299 Hz

julia> get_code_phase(track_state, 1)
50.00064935064897
```

`estimate_cn0(track_state, prn)` returns the CN0 estimate in dB-Hz. With a noise-free test signal it is `Inf dB-Hz`; real signals typically land in 30–50 dB-Hz.

Everything beyond this minimal case — multi-satellite, multi-system, multi-signal, multi-band, phased arrays, acquisition handoff, removing satellites, the accessor ladder — is covered in [Tracking State](tracking_state.md). Real-time loop patterns (hoisting the correlator, `track!`, allocation behavior) are covered in [Track](track.md).

!!! tip "Real-time loops: hoist the correlator"
    For loops processing many chunks in sequence, construct
    [`CPUThreadedDownconvertAndCorrelator`](@ref) (or
    [`CPUDownconvertAndCorrelator`](@ref)) **once outside** the loop and
    pass it via `downconvert_and_correlator =`. See [Track](track.md#Real-time-use).

## Q/A

- Why are the correlator values zero?

The correlator output given by `get_last_fully_integrated_correlator` is the correlation result after the
code phase has reached the full code length or multiples of the code length. If the current
tracked signal does not include the end of the PRN sequence (or multiples of that), the
correlator from the last complete integration will be returned. At the very start of tracking,
before any complete integration has occurred, the correlator values will be zero.
