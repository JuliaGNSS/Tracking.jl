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

Tracking.jl can track one or multiple satellites simultaneously, each potentially on multiple signals. The tracking state is held in a [`TrackState`](@ref). Satellites are organized into **signal groups** — named groups of sats that share the same signal-tuple shape (and therefore the same concrete `TrackedSat` type, which is what gives the hot loop type stability).

### Constructing a `TrackState`

Declare which signals each group tracks. Each entry is a tuple of `AbstractGNSSSignal` instances; the first signal in the tuple is the Doppler source (its correlator drives the PLL/DLL).

```jldoctest single_signal
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());
```

The singular `signal = GPSL1CA()` keyword is the shortcut for the common one-group, one-signal case. It desugars internally to `signals = (default = (GPSL1CA(),),)` so the rest of the API can stay uniform. With one group, [`add_satellite!`](@ref) may omit the `group=` keyword. Use the plural `signals = (...)` keyword for multi-signal or multi-group tracking — see below.

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

`estimate_cn0(track_state, group, prn)` returns the CN0 estimate in dB-Hz. With a noise-free test signal, this will be `Inf dB-Hz`. With real signals containing noise, you'll get finite values typically in the range of 30-50 dB-Hz.

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

When different satellites carry different signal types, use multiple named groups. Each group has its own concrete `TrackedSat` value type, so type inference stays sharp across the heterogeneous mix.

```jldoctest multi_system
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (
               gps     = (GPSL1CA(),),
               galileo = (GalileoE1B(),),
           ),
       );

julia> add_satellite!(track_state; prn = 1,  group = :gps,     code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> add_satellite!(track_state; prn = 11, group = :galileo, code_phase = 200.0, carrier_doppler = -300.0Hz);

julia> get_carrier_doppler(track_state, :gps, 1)
1000.0 Hz

julia> get_carrier_doppler(track_state, :galileo, 11)
-300.0 Hz
```

### Multi-signal tracking (one satellite, several signals)

A modern GPS satellite transmits L1 C/A, L1C-D, and L1C-P simultaneously on the same carrier. Tracking.jl can track all three together on one satellite, sharing a single carrier downconvert per outer iteration:

```jldoctest modern_gps
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (
               modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
           ),
       );

julia> add_satellite!(track_state;
           prn = 11, group = :modern_gps,
           code_phase = 0.0, carrier_doppler = 1234.0Hz,
       );

julia> get_carrier_doppler(track_state, :modern_gps, 11)
1234.0 Hz
```

The first signal in each group's tuple is the **Doppler source** — its correlator is what the PLL/DLL discriminator runs on. Putting a pilot signal first (e.g. `GPSL1C_P()`) is encouraged when one is available: pilot signals carry no data-bit modulation, which lets the PLL run longer coherent integrations and reach lower phase-noise floors. The data-bearing signals (L1C-D, L1 C/A) still recover their navigation bits independently — each [`TrackedSignal`](@ref) carries its own `bit_buffer` regardless of which signal drives the loop filter.

When a satellite tracks signals with different primary-code lengths (e.g. L1 C/A at 1 ms vs L1C-P at 10 ms), each outer iteration integrates to the **shortest** signal's next primary-code boundary. The shorter signal's correlator completes every iteration; the longer signal's correlator accumulates across multiple iterations and only marks `is_integration_completed = true` on its own boundary. Doppler updates therefore happen at the shortest signal's cadence (1 ms in this example), and longer signals see their integration windows spanned by piecewise Doppler updates — the natural per-iteration-Doppler-correction behaviour of a real receiver.

### Multi-band tracking (different RF carriers)

A satellite often broadcasts on more than one RF band — GPS broadcasts on L1 (1575.42 MHz) and L5 (1176.45 MHz); Galileo broadcasts on E1 (L1) and E5a (L5). In a multi-band receiver these arrive from separate front-ends, generally at different sample rates, and need to be downconverted and correlated against their own carrier replicas. Tracking.jl exposes this as a multi-band `TrackState` where each group declares which RF band it sits on.

**Why this matters.** A single physical satellite tracked on two bands gives the receiver two near-independent observations of the same path. The classic uses:

- **Ionospheric correction** via dual-frequency (iono-free) pseudorange combinations.
- **Wider effective bandwidth** for code-phase observations (L5 carries far more chip-rate bandwidth than L1 C/A).
- **Cross-band-aided tracking**: the carrier Doppler ratio between L1 and L5 is exactly the ratio of their RF carrier frequencies. A joint estimator can fuse the two bands' discriminators and produce a more accurate Doppler estimate than either band alone — particularly valuable at low CN0 where the wider data-aided integration on L5 helps the noisier L1 C/A.

This release ships the **structural enablers** for multi-band: per-band groups, per-band measurement routing, an estimation barrier that sees every band's correlator outputs at once. The cross-band joint-tracking *algorithm* (e.g. linking PRN-X-on-L1 with PRN-X-on-L5 in one estimator step) is a follow-up — see [docs/plans/2026-05-15-multi-band-tracking-design.md](https://github.com/JuliaGNSS/Tracking.jl/blob/master/docs/plans/2026-05-15-multi-band-tracking-design.md) for the design and the open mechanism question.

#### Declaring bands

The `band` field of each [`SignalGroup`](@ref) is inferred from `get_band(signals[1])`, so you don't normally type it. Mix signals from different bands in the same `signals = (...)` keyword and the bands fall out:

```jldoctest multi_band_groups
julia> using Tracking, GNSSSignals

julia> track_state = TrackState(;
           signals = (
               legacy_gps_l1 = (GPSL1CA(),),
               modern_gps_l1 = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
               galileo       = (GalileoE1B(),),
               gps_l5        = (GPSL5I(),),
           ),
       );

julia> keys(track_state.groups)
(:legacy_gps_l1, :modern_gps_l1, :galileo, :gps_l5)
```

Four groups, two distinct bands — the first three groups all sit on L1 (GPS L1 and Galileo E1 share the 1575.42 MHz carrier), the fourth sits on L5. Two groups sharing a band is fine; the grouping partitions satellites by *signal-tuple shape* (the type-stability axis), not by band.

#### Tracking against multiple measurements

For multi-band tracking, build one [`Measurement`](@ref) per band — bundling sample buffer and front-end metadata — and pass them as a NamedTuple keyed by [`band_key`](@ref):

```jldoctest multi_band_track; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_code_frequency, get_code_center_frequency_ratio

julia> function make_signal(sys, prn, carrier_doppler, num_samples, fs)
           code_freq = carrier_doppler * get_code_center_frequency_ratio(sys) + get_code_frequency(sys)
           range = 0:num_samples-1
           cis.(2π .* carrier_doppler .* range ./ fs) .*
               gen_code(num_samples, sys, prn, fs, code_freq, 0.0)
       end;

julia> track_state = TrackState(;
           signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)),
       );

julia> add_satellite!(track_state; prn = 1, group = :legacy_gps_l1, code_phase = 0.0, carrier_doppler = 200Hz);

julia> add_satellite!(track_state; prn = 1, group = :gps_l5, code_phase = 0.0, carrier_doppler = -150Hz);

julia> buf_l1 = make_signal(GPSL1CA(),  1, 200Hz,  4000,  4e6Hz);  # 1 ms at 4 MHz

julia> buf_l5 = make_signal(GPSL5I(),   1, -150Hz, 25000, 25e6Hz); # 1 ms at 25 MHz

julia> track!((l1 = Measurement(buf_l1, 4e6Hz),
               l5 = Measurement(buf_l5, 25e6Hz)), track_state);

julia> get_carrier_doppler(track_state, :legacy_gps_l1, 1)
200.00000359633913 Hz
```

The keys (`:l1`, `:l5`) come from `band_key(L1())` and `band_key(L5())`. All measurements must cover the **exact same observation duration** — `num_samples / sampling_frequency` must compare equal across bands. An L1 chunk of 4000 samples at 4 MHz and an L5 chunk of 25000 samples at 25 MHz both cover 1 ms, so they're compatible; an L5 chunk of 25001 samples is rejected.

#### Per-band antenna counts

Different bands often come from different front-ends with different antenna arrangements. To declare per-band antenna counts, pass [`SignalGroup`](@ref) instances directly as the entries — the bare-tuple shortcut uses the constructor's single `num_ants` kwarg for all groups, but the `SignalGroup` form lets each group set its own:

```jldoctest per_band_ants
julia> using Tracking, GNSSSignals

julia> using Tracking: NumAnts

julia> track_state = TrackState(;
           signals = (
               legacy_gps_l1 = SignalGroup((GPSL1CA(),); num_ants = NumAnts(2)),
               gps_l5        = SignalGroup((GPSL5I(),);  num_ants = NumAnts(1)),
           ),
       );

julia> track_state.groups[:legacy_gps_l1].num_ants
NumAnts{2}()

julia> track_state.groups[:gps_l5].num_ants
NumAnts{1}()
```

Two groups on the same band must declare the same `num_ants` — they share a physical front-end. The constructor errors at TrackState construction if they disagree.

#### Bare-buffer compatibility

A single-band receiver doesn't need to type any of this. The bare-buffer call `track!(buf, state, fs)` keeps working for any `TrackState` that spans exactly one band — internally it wraps the buffer into a one-entry NamedTuple keyed by the lone band. Pass `intermediate_frequency` via the same kwarg as before, or move it onto a [`Measurement`](@ref) when you migrate to multi-band.

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
