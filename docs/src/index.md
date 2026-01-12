# Tracking.jl

This package implements the tracking functionality of GNSS satellites that's part of the larger GNSS receiver.
Tracking.jl primarily consists of two main blocks:

1. Signal down-conversion and correlation
2. Code and carrier estimation to generate replicas and close the loop

Tracking.jl provides defaults for both blocks, but it provides mechanisms to hook in your own implementation (by using multiple dispatch).
For signal down-conversion and correlation Tracking.jl provides a highly optimized CPU implementation using SIMD vectorization via LoopVectorization.jl.
With respect to the second block Tracking.jl provides conventional DLLs and PLLs with FLL-assisted carrier tracking as the default.
Down-conversion and correlation is done in full code blocks meaning from code start to code end or multiples of that (e.g. in GPS L1 from 0 to `N`*1023). The factor `N` can be specified, but will be `1` as long as the bit start is unknown in order to find the bit start. Once that is done for every tracked satellite the result will be handed over to the code and carrier estimation block.
Moreover, Tracking.jl allows tracking of signals from phased antenna arrays meaning that they are down-converted and correlated by the very same replica to conserve phase relationships.

```@contents
Pages = [
  "track.md",
  "tracking_state.md",
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

Tracking.jl can track one or multiple satellites simultaneously. The tracking state is held in a `TrackState` which contains `SatState` objects for each tracked satellite.

### Single Satellite Tracking

To track a single satellite, create a `SatState` with the GNSS system, PRN, code phase, and carrier Doppler (typically from acquisition):

```jldoctest single_sat
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> system = GPSL1();

julia> prn = 1;

julia> code_phase = 50.0;

julia> carrier_doppler = 1000.0Hz;

julia> sat_state = SatState(system, prn, code_phase, carrier_doppler);

julia> track_state = TrackState(system, sat_state);

julia> get_prn(track_state)
1

julia> get_code_phase(track_state)
50.0

julia> get_carrier_doppler(track_state)
1000.0 Hz
```

To process a signal, call the `track` function. The signal should be a vector of complex samples:

```jldoctest single_sat
julia> using GNSSSignals: gen_code, get_code_frequency, get_code_center_frequency_ratio

julia> sampling_frequency = 4e6Hz;

julia> num_samples = 4000;

julia> code_frequency = carrier_doppler * get_code_center_frequency_ratio(system) + get_code_frequency(system);

julia> signal = cis.(2π .* carrier_doppler .* (0:num_samples-1) ./ sampling_frequency) .*
           gen_code(num_samples, system, prn, sampling_frequency, code_frequency, code_phase);

julia> track_state = track(signal, track_state, sampling_frequency);
```

After tracking, you can retrieve the updated tracking parameters and results:

```jldoctest single_sat; filter = r"[0-9]+\.[0-9]+" => "***"
julia> get_carrier_doppler(track_state)
999.9999883655299 Hz

julia> get_code_phase(track_state)
50.00064935064897

julia> get_prompt(get_last_fully_integrated_correlator(track_state))
3805.0 - 0.000841952976770699im
```

Note: `estimate_cn0(track_state)` returns the CN0 estimate in dB-Hz. With a noise-free
test signal, this will be `Inf dB-Hz`. With real signals containing noise, you'll get
finite values typically in the range of 30-50 dB-Hz.

### Multiple Satellite Tracking

To track multiple satellites of the same system:

```jldoctest multi_sat
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> system = GPSL1();

julia> sat_state1 = SatState(system, 1, 50.0, 1000.0Hz);

julia> sat_state2 = SatState(system, 5, 120.0, -500.0Hz);

julia> sat_state3 = SatState(system, 17, 890.0, 2000.0Hz);

julia> track_state = TrackState(system, [sat_state1, sat_state2, sat_state3]);

julia> get_carrier_doppler(track_state, 5)
-500.0 Hz

julia> get_code_phase(track_state, 17)
890.0
```

### Multi-System Tracking

Tracking.jl supports tracking satellites from different GNSS systems simultaneously:

```jldoctest multi_system
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> gps_sat = SatState(GPSL1(), 1, 50.0, 1000.0Hz);

julia> galileo_sat = SatState(GalileoE1B(), 11, 200.0, -300.0Hz);

julia> gps_system_state = SystemSatsState(GPSL1(), gps_sat);

julia> galileo_system_state = SystemSatsState(GalileoE1B(), galileo_sat);

julia> track_state = TrackState((gps_system_state, galileo_system_state));

julia> get_carrier_doppler(track_state, 1, 1)  # system 1, PRN 1
1000.0 Hz

julia> get_carrier_doppler(track_state, 2, 11)  # system 2, PRN 11
-300.0 Hz
```

You can also use a named tuple for easier access:
```jldoctest multi_system_named
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> gps_sat = SatState(GPSL1(), 1, 50.0, 1000.0Hz);

julia> galileo_sat = SatState(GalileoE1B(), 11, 200.0, -300.0Hz);

julia> gps_system_state = SystemSatsState(GPSL1(), gps_sat);

julia> galileo_system_state = SystemSatsState(GalileoE1B(), galileo_sat);

julia> track_state = TrackState((gps = gps_system_state, galileo = galileo_system_state));

julia> get_carrier_doppler(track_state, :gps, 1)
1000.0 Hz

julia> get_carrier_doppler(track_state, :galileo, 11)
-300.0 Hz
```

### Adding and Removing Satellites

Satellites can be dynamically added or removed during tracking:

```jldoctest merge_filter
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> system = GPSL1();

julia> sat_state = SatState(system, 1, 50.0, 1000.0Hz);

julia> track_state = TrackState(system, sat_state);

julia> new_sat = SatState(system, 23, 500.0, 1500.0Hz);

julia> track_state = merge_sats(track_state, new_sat);

julia> get_carrier_doppler(track_state, 23)
1500.0 Hz

julia> track_state = filter_out_sats(track_state, [1]);

julia> haskey(get_sat_states(track_state), 23)
true

julia> haskey(get_sat_states(track_state), 1)
false
```

## Track Multiple Signals Coherently (Phased Arrays)

Tracking.jl provides a way to track multiple signals coherently, e.g. to track signals from
a phased array. In that case the input `signal` should be a Matrix instead of a Vector,
where the number of rows is equal to the number of samples and the number of
columns is equal to the number of antenna elements. Furthermore, you need to specify the number
of antenna elements to the satellite state:

```jldoctest phased_array
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> system = GPSL1();

julia> sat_state = SatState(system, 1, 50.0, 1000.0Hz, num_ants = NumAnts(4));

julia> track_state = TrackState(system, sat_state);

julia> get_num_ants(track_state)
4
```

By default the track function will use the first antenna channel as the reference signal to
drive the discriminators. However, an appropriate beamforming algorithm will probably
suit better. For that, you can pass a custom `post_corr_filter` to the `SatState`.

## GPU Support

To use GPU acceleration features, you need to explicitly load CUDA:

```julia
using Tracking
using CUDA  # Activates GPU functionality

# Access GPU types via the extension
ext = Base.get_extension(Tracking, :TrackingCUDAExt)
gpu_correlator = ext.GPUDownconvertAndCorrelator(...)
```

Note: The GPU implementation is available for reference but is not yet as optimized as the
CPU implementation which uses SIMD vectorization via LoopVectorization.jl.

## Q/A

- Why are the correlator values zero?

The correlator output given by `get_last_fully_integrated_correlator` is the correlation result after the
code phase has reached the full code length or multiples of the code length. If the current
tracked signal does not include the end of the PRN sequence (or multiples of that), the
correlator from the last complete integration will be returned. At the very start of tracking,
before any complete integration has occurred, the correlator values will be zero.
