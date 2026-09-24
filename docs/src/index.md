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

The per-record loop arithmetic — correlators, discriminators, the bit buffer,
the C/N₀ estimators and the Doppler estimators — lives in
[TrackingLoops.jl](https://github.com/JuliaGNSS/TrackingLoops.jl), so that a
hardware correlator's loop process can run the same code without this
package's sample-domain half. Reading results needs only Tracking.jl: it
re-exports the TrackingLoops functions that read them — [`estimate_cn0`](@ref),
[`get_prompt`](@ref), `get_soft_bits`, `has_bit_or_secondary_code_been_found`
and a few more. Load
TrackingLoops next to it to configure a correlator, a C/N₀ or noise estimator,
a Doppler estimator or a post-correlation filter. Its API is listed in the
[Loop core reference](trackingloops.md).

```julia
using Tracking, TrackingLoops, GNSSSignals
```

## Supported signals

Tracking.jl tracks **every** concrete signal type GNSSSignals.jl defines — GPS
L1 C/A, L1C, L2C and L5, Galileo E1, E5a (including the E5a-QP acquisition aid),
E5b and E6, and BeiDou B1I, B3I, B1C, B2a and B2b.

The [signal support matrix](signals.md) lists each one with its coherent
integration length, its sync feature, whether it carries navigation data, and
the PRNs its code table defines — plus the known per-signal limitations.
`test/signal_coverage.jl` checks that matrix against the live GNSSSignals type
tree, so a release that adds a signal fails the test suite rather than a user's
`get_default_correlator` call.

Signals sharing a band can be tracked in one [`BandMeasurement`](@ref)
regardless of constellation: Galileo E5b and BeiDou B2b share the 1207.14 MHz
carrier (`E5b` for both), BeiDou B2a shares `L5` with GPS L5 and Galileo E5a,
and BeiDou B1C shares `L1` with GPS L1 C/A, GPS L1C and Galileo E1.

```@contents
Pages = [
  "signals.md",
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
