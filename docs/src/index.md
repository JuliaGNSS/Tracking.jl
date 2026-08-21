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

### GPS

- GPS L1 C/A
- GPS L1C-D (data)
- GPS L1C-P (pilot, with 1800-chip overlay-code sync)
- GPS L2CM (data)
- GPS L2CL (pilot, dataless — no secondary code)
- GPS L5I (data, with NH10 secondary-code sync)
- GPS L5Q (pilot, with NH20 secondary-code sync)

### Galileo

- Galileo E1B (and its BOC(1,1) approximation)
- Galileo E1C (pilot, with CS25 secondary-code sync; and its BOC(1,1) approximation)
- Galileo E5a-I (data, with CS20 secondary-code sync)
- Galileo E5a-Q (pilot, with per-PRN CS100 secondary-code sync)
- Galileo E5b-I (data, with CS4 secondary-code sync)
- Galileo E5b-Q (pilot, with per-PRN CS100 secondary-code sync)
- Galileo E6-B (data, C/NAV at 1000 sym/s — one symbol per primary period)
- Galileo E6-C (pilot, with per-PRN CS100 secondary-code sync)

### BeiDou

- BeiDou B1I (data, with per-PRN NH20 secondary-code sync)
- BeiDou B3I (data, with per-PRN NH20 secondary-code sync)
- BeiDou B2b-I (data, B-CNAV3 at 1000 sym/s — one symbol per primary period)
- BeiDou B2a data (with 5-chip secondary-code sync)
- BeiDou B2a pilot (with per-PRN 100-chip secondary-code sync)
- BeiDou B1C data (BOC(1,1), one symbol per primary period)
- BeiDou B1C pilot (BOC(1,1), with 1800-chip overlay-code sync)

Galileo E5b and BeiDou B2b share the 1207.14 MHz carrier, which GNSSSignals.jl
names the `E5b` band for both; BeiDou B2a shares `L5` with GPS L5 and Galileo
E5a, and BeiDou B1C shares `L1` with GPS L1 C/A, GPS L1C and Galileo E1. Signals
on one band can be tracked in one [`BandMeasurement`](@ref) regardless of
constellation.

Galileo E5a-QP is **not** tracked, by design: it is an acquisition aid rather
than a tracking signal — a dataless BPSK(5) quasi-pilot whose 330-chip code
repeats every 64.5 µs, there to let a receiver find E5a cheaply and hand the
satellite over to E5a-I / E5a-Q (OS SIS ICD v2.2 §2.3.1.4). It carries neither
data nor a secondary code, so there is no sync feature to lock, and one primary
code period — the unit Tracking integrates in and derives its loop-bandwidth
defaults from — is 64.5 µs, far too short for either to be meaningful.
`test/signal_coverage.jl` records that decision, and fails if GNSSSignals.jl
adds a further signal Tracking.jl has no support for.

On the BeiDou GEO satellites (PRN 1-5 and 59-63) the B1I/B3I ranging codes carry
no NH20 overlay and the navigation message is D2 at 500 sym/s rather than D1 at
50 sym/s. Those two together leave the secondary-code detector nothing to lock
onto, so a GEO satellite tracks and can be *ranged* on but never syncs, and its
data is not decoded: a 2-block D2 bit-edge search would need a per-PRN data rate,
while `get_data_frequency` is a per-signal-type accessor reporting the D1 rate
for every PRN — see `src/beidou/b1i.jl`.

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
