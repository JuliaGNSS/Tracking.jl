# Tracking.jl

Tracking.jl tracks GNSS satellites as part of a receiver in two blocks:

1. downconversion and correlation, with SIMD-optimized CPU kernels;
2. code and carrier estimation, closing the loop (by default FLL-assisted PLL
   and DLL).

Both can be replaced through multiple dispatch.

Correlation runs over whole code blocks: `N` primary code periods (e.g. `N`·1023
chips for GPS L1 C/A), with `N = 1` until the bit edge or secondary code is
found. Antenna-array samples are correlated with one replica, preserving their
phase relationships. One satellite can be tracked on several signals at once
(e.g. GPS L1 C/A, L1C-D and L1C-P), sharing the carrier downconvert, with a
correlator per signal.

The per-record loop arithmetic (correlators, discriminators, bit buffer, C/N₀
and Doppler estimators) lives in
[TrackingLoops.jl](https://github.com/JuliaGNSS/TrackingLoops.jl), so a hardware
correlator's loop can run the same code. Tracking.jl re-exports the functions
that read results:
[`estimate_cn0`](@extref TrackingLoops.estimate_cn0-Tuple{MomentsCN0Estimator, Any}),
[`get_prompt`](@extref TrackingLoops.get_prompt), `get_soft_bits`,
`has_bit_or_secondary_code_been_found` and a few more. Load TrackingLoops to
configure a correlator, a C/N₀ or noise estimator, a Doppler estimator or a
post-correlation filter, or for vector tracking; see the
[TrackingLoops.jl manual](https://juliagnss.github.io/TrackingLoops.jl/stable/).

```julia
using Tracking, TrackingLoops, GNSSSignals
```

## Supported signals

Tracking.jl tracks **every** concrete signal type GNSSSignals.jl defines: GPS
L1 C/A, L1C, L2C and L5, Galileo E1, E5a (including the E5a-QP acquisition aid),
E5b and E6, and BeiDou B1I, B3I, B1C, B2a and B2b. The
[signal support matrix](signals.md) lists each one's integration length, sync
feature, navigation data, PRNs and known limitations.

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

`estimate_cn0(track_state, prn)` returns the C/N₀ in dB-Hz: `Inf dB-Hz` for a noise-free test signal, typically 30–50 dB-Hz for real ones.

Multiple satellites, systems, signals and bands, antenna arrays, acquisition handoff and the accessors are covered in [Tracking State](tracking_state.md); real-time loops in [Track](track.md).

!!! tip "Real-time loops: hoist the correlator"
    Construct [`CPUThreadedDownconvertAndCorrelator`](@ref) (or
    [`CPUDownconvertAndCorrelator`](@ref)) **once outside** the loop and pass it
    via `downconvert_and_correlator =`. See [Track](track.md#Real-time-loops).

## Q/A

- Why are the correlator values zero?

`get_last_fully_integrated_correlator` returns the last *completed* integration,
one that reached the end of its code blocks. Before the first one completes it
is zero.
