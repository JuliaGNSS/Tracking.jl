# Tracking.jl

This package implements the tracking functionality of GNSS satellites that's part of the larger GNSS receiver.
Tracking.jl primarily consists of two main blocks:

1. Signal down-conversion and correlation
2. Code and carrier estimation to generate replicas and close the loop

Tracking.jl provides defaults for both blocks, but it provides mechanisms to hook in your own implementation (by using multiple dispatch).
For signal down-conversion and correlation Tracking.jl provides a highly optimized CPU implementation using hand-tuned SIMD intrinsics.
With respect to the second block Tracking.jl provides conventional DLLs and PLLs with FLL-assisted carrier tracking as the default.
Down-conversion and correlation is done in full code blocks meaning from code start to code end or multiples of that (e.g. in GPS L1 C/A from 0 to `N`*1023). The factor `N` can be specified, but will be `1` as long as the bit start is unknown in order to find the bit start. Once that is done for every tracked satellite the result will be handed over to the code and carrier estimation block.
Moreover, Tracking.jl allows tracking of signals from phased antenna arrays meaning that they are down-converted and correlated by the very same replica to conserve phase relationships.
Multi-signal tracking is supported: a single satellite can be tracked on several signals at once (e.g. GPS L1 C/A together with L1C-D and L1C-P) sharing one carrier downconvert per outer iteration, with a per-signal correlator each.

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

Tracking.jl can track one or multiple satellites simultaneously, each potentially on multiple signals. The tracking state is held in a [`TrackState`](@ref). Satellites are organized into **capabilities** — named groups that share the same signal-tuple shape.

### Constructing a `TrackState`

Declare which signals each capability tracks. Each entry is a tuple of `AbstractGNSSSignal` instances; the first signal in the tuple is the Doppler source (its correlator drives the PLL/DLL).

```jldoctest single_signal
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());
```

The singular `signal = GPSL1CA()` keyword is the shortcut for the common one-capability, one-signal case. It desugars internally to `signals = (default = (GPSL1CA(),),)` so the rest of the API can stay uniform. With one capability, [`add_satellite!`](@ref) may omit the `capability=` keyword. Use the plural `signals = (...)` keyword for multi-signal or multi-capability tracking — see below.

### Adding satellites

Satellites are added to a `TrackState` via [`add_satellite!`](@ref). The acquisition handoff values (`prn`, `code_phase`, `carrier_doppler`, optionally `code_doppler` and `carrier_phase`) get wired into a fresh [`TrackedSat`](@ref) with the library's default correlator and post-correlation filter.

```jldoctest single_signal
julia> add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);

julia> get_prn(track_state, :default, 1)
1

julia> get_code_phase(track_state, :default, 1)
50.0

julia> get_carrier_doppler(track_state, :default, 1)
1000.0 Hz
```

Adding a satellite with the same PRN again overwrites the existing entry (matching [`merge_sats`](@ref) semantics — no error).

### Tracking a signal

To process a signal, call the [`track`](@ref) function. The signal should be a vector of complex samples:

```jldoctest single_signal
julia> using GNSSSignals: gen_code, get_code_frequency, get_code_center_frequency_ratio

julia> sampling_frequency = 4e6Hz;

julia> num_samples = 4000;

julia> system = GPSL1CA();

julia> code_frequency = 1000.0Hz * get_code_center_frequency_ratio(system) + get_code_frequency(system);

julia> signal = cis.(2π .* 1000.0Hz .* (0:num_samples-1) ./ sampling_frequency) .*
           gen_code(num_samples, system, 1, sampling_frequency, code_frequency, 50.0);

julia> track_state = track(signal, track_state, sampling_frequency);
```

After tracking, you can retrieve the updated tracking parameters and results:

```jldoctest single_signal; filter = r"[0-9]+\.[0-9]+" => "***"
julia> get_carrier_doppler(track_state, :default, 1)
999.9999883655299 Hz

julia> get_code_phase(track_state, :default, 1)
50.00064935064897
```

`estimate_cn0(track_state, capability, prn)` returns the CN0 estimate in dB-Hz. With a noise-free test signal, this will be `Inf dB-Hz`. With real signals containing noise, you'll get finite values typically in the range of 30-50 dB-Hz.

!!! tip "Real-time loops: hoist the correlator"
    For loops that process many signal chunks in sequence (e.g. an SDR
    capture stream), construct
    [`CPUThreadedDownconvertAndCorrelator`](@ref) (or
    [`CPUDownconvertAndCorrelator`](@ref)) **once outside** the loop and
    pass it via `downconvert_and_correlator =`:

    ```julia
    dc = CPUThreadedDownconvertAndCorrelator()
    while got_chunk(rx)
        chunk = read_chunk!(rx)
        track_state = track(chunk, track_state, sampling_frequency;
                            downconvert_and_correlator = dc)
    end
    ```

    The correlator holds long-lived per-thread scratch buffers; relying on
    the default kwarg value rebuilds them on every call. See also
    [`track!`](@ref) for the in-place variant.

### Multi-satellite tracking

To track several satellites, simply call `add_satellite!` repeatedly:

```jldoctest multi_sat
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> add_satellite!(track_state; prn = 1,  code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> add_satellite!(track_state; prn = 5,  code_phase = 120.0, carrier_doppler = -500.0Hz);

julia> add_satellite!(track_state; prn = 17, code_phase = 890.0, carrier_doppler = 2000.0Hz);

julia> get_carrier_doppler(track_state, :default, 5)
-500.0 Hz

julia> get_code_phase(track_state, :default, 17)
890.0
```

### Multi-system tracking (different signals on different sats)

When different satellites carry different signal types, use multiple named capabilities. Each capability has its own concrete `TrackedSat` value type, so type inference stays sharp across the heterogeneous mix.

```jldoctest multi_system
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (
               gps     = (GPSL1CA(),),
               galileo = (GalileoE1B(),),
           ),
       );

julia> add_satellite!(track_state; prn = 1,  capability = :gps,     code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> add_satellite!(track_state; prn = 11, capability = :galileo, code_phase = 200.0, carrier_doppler = -300.0Hz);

julia> get_carrier_doppler(track_state, :gps, 1)
1000.0 Hz

julia> get_carrier_doppler(track_state, :galileo, 11)
-300.0 Hz
```

### Multi-signal tracking (one satellite, several signals)

A modern GPS satellite transmits L1 C/A, L1C-D, and L1C-P simultaneously on the same carrier. Tracking.jl can track all three together on one satellite, sharing a single carrier downconvert per outer iteration:

```julia
track_state = TrackState(;
    signals = (
        modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
    ),
)
add_satellite!(track_state;
    prn = 11, capability = :modern_gps,
    code_phase = 0.0, carrier_doppler = 1234.0Hz,
)
```

The first signal in each capability's tuple is the **Doppler source** — its correlator is what the PLL/DLL discriminator runs on. Putting a pilot signal first (e.g. `GPSL1C_P()`) is encouraged when one is available: pilot signals carry no data-bit modulation, which lets the PLL run longer coherent integrations and reach lower phase-noise floors. The data-bearing signals (L1C-D, L1 C/A) still recover their navigation bits independently — each [`TrackedSignal`](@ref) carries its own `bit_buffer` regardless of which signal drives the loop filter.

When a satellite tracks signals with different primary-code lengths (e.g. L1 C/A at 1 ms vs L1C-P at 10 ms), each outer iteration integrates to the **shortest** signal's next primary-code boundary. The shorter signal's correlator completes every iteration; the longer signal's correlator accumulates across multiple iterations and only marks `is_integration_completed = true` on its own boundary. Doppler updates therefore happen at the shortest signal's cadence (1 ms in this example), and longer signals see their integration windows spanned by piecewise Doppler updates — the natural per-iteration-Doppler-correction behaviour of a real receiver.

### Removing satellites

```jldoctest remove_sats
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> add_satellite!(track_state; prn = 1,  code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> add_satellite!(track_state; prn = 23, code_phase = 500.0, carrier_doppler = 1500.0Hz);

julia> remove_satellite!(track_state; prn = 1);

julia> haskey(get_sat_states(track_state, :default), 23)
true

julia> haskey(get_sat_states(track_state, :default), 1)
false
```

## Track Multiple Signals Coherently (Phased Arrays)

Tracking.jl provides a way to track multiple signals coherently, e.g. to track signals from
a phased array. In that case the input `signal` should be a Matrix instead of a Vector,
where the number of rows is equal to the number of samples and the number of
columns is equal to the number of antenna elements. Furthermore, you need to specify the number
of antenna elements at `TrackState` construction:

```jldoctest phased_array
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signal = GPSL1CA(),
           num_ants = NumAnts(4),
       );

julia> add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);

julia> get_num_ants(track_state, :default, 1)
4
```

By default the track function will use the first antenna channel as the reference signal to
drive the discriminators. However, an appropriate beamforming algorithm will probably
suit better. For that, construct a [`TrackedSat`](@ref) with a custom `post_corr_filter` and
hand it to `add_satellite!`'s escape-hatch overload:

```julia
sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                 num_ants = NumAnts(4),
                 post_corr_filter = MyBeamformer())
add_satellite!(track_state, :default, sat)
```

## Q/A

- Why are the correlator values zero?

The correlator output given by `get_last_fully_integrated_correlator` is the correlation result after the
code phase has reached the full code length or multiples of the code length. If the current
tracked signal does not include the end of the PRN sequence (or multiples of that), the
correlator from the last complete integration will be returned. At the very start of tracking,
before any complete integration has occurred, the correlator values will be zero.
