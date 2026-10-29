# Tracking State

The tracking state nests as **TrackState → SignalGroup → TrackedSat → TrackedSignal**:

- [`TrackState`](@ref): a `NamedTuple` of [`SignalGroup`](@ref)s and the Doppler estimator.
- [`SignalGroup`](@ref): satellites with the same signal tuple, and so one concrete `TrackedSat` type (the hot loop's type stability), on one RF band and antenna array.
- [`TrackedSat`](@ref): the carrier/code Doppler and phase its signals share, a tuple of [`TrackedSignal`](@ref)s, and the Doppler-estimator state.
- [`TrackedSignal`](@ref): correlator, post-correlation filter, C/N₀ estimator, bit buffer and integration progress.

### Estimator-driver signal

The first signal of each group's tuple is the **estimator-driver signal**, the one the Doppler estimator updates the shared carrier and code Doppler from. The conventional estimators ([`ConventionalPLLAndDLL`](@extref TrackingLoops.ConventionalPLLAndDLL), [`ConventionalAssistedPLLAndDLL`](@extref TrackingLoops.ConventionalAssistedPLLAndDLL-Union{Tuple{}, Tuple{Type{CO}}, Tuple{CO}} where CO<:TrackingLoopFilters.AbstractLoopFilter)) run their discriminators on it and take its default loop bandwidths (a carrier loop pulling in at 50 Hz and narrowing to 18 Hz, 1 Hz code). A custom [`AbstractDopplerEstimator`](@extref TrackingLoops.AbstractDopplerEstimator) may use any signal.

Bit sync, the post-correlation filter and the **C/N₀ estimator** run per signal, so a multi-signal satellite has one C/N₀ per signal; see [CN0 Estimator](cn0_estimator.md) and [`NoCN0Estimator`](@extref TrackingLoops.NoCN0Estimator), the per-signal opt-out.

## Choosing a `TrackState` constructor

The right constructor depends on **when you know which satellites you'll track**.

### One-shot scripts: acquire once, then track

The [Acquisition](https://github.com/JuliaGNSS/Acquisition.jl) extension's constructors derive signals, loop bandwidths and satellites from `AcquisitionResults`:

```julia
using Acquisition  # loads the extension; required for TrackState(acq...)

# Single satellite
acq = acquire(GPSL1CA(), data, sampling_frequency, 7)
track_state = TrackState(acq)

# Many satellites, one signal
acqs = acquire(GPSL1CA(), data, sampling_frequency, 1:32)
track_state = TrackState(filter(is_detected, acqs))

# Many satellites, multiple signals
acqs = vcat(
    acquire(GPSL1CA(),    data, sampling_frequency, 1:32),
    acquire(GalileoE1B(), data, sampling_frequency, 1:36),
)
track_state = TrackState(filter(is_detected, acqs);
    signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(),)),
)
```

[`add_satellite!`](@ref) can later add satellites only for the signals those acqs declared, as the groups are frozen at construction. To track more signals, build empty and populate instead.

### Real-time / repeating loops: build empty, then populate

A receiver that re-acquires periodically builds an **empty** `TrackState` once with `TrackState(; signal = ...)` or `TrackState(; signals = (...))` and adds and removes satellites with [`add_satellite!`](@ref) and [`remove_satellite!`](@ref):

```julia
# Build once
track_state = TrackState(;
    signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(),)),
)

# In your acquisition loop
while running
    acqs = acquire(GPSL1CA(), latest_chunk, sampling_frequency, candidate_prns)
    track_state = add_satellite!(track_state, filter(is_detected, acqs))   # routes each acq to the matching group
    track_state = track(latest_chunk, track_state, sampling_frequency)
end
```

The `TrackState`'s concrete type stays fixed as satellites come and go. `signal = GPSL1CA()` is short for `signals = (default = (GPSL1CA(),),)`; with one group, [`add_satellite!`](@ref) may omit `group_key =`.

### Power-user: pre-built `TrackedSat`s

To change the correlator or post-correlation filter *type*, build the `TrackedSat`s yourself, pass them to `TrackState(signal, sats)`, and add more with `add_satellite!(track_state, group, sat)`. The single-signal constructor:

```julia
TrackedSat(
    signal,                # e.g. GPSL1CA()
    prn::Int,
    code_phase,
    carrier_doppler;
    # all kwargs below are optional and have signal-derived defaults
    doppler_estimator     = ConventionalAssistedPLLAndDLL(...),
    num_ants              = NumAnts(1),
    correlator            = get_default_correlator(signal, num_ants),
    carrier_phase         = 0.0,
    code_doppler          = carrier_doppler * get_chip_rate_to_center_frequency_ratio(signal),
    num_prompts_for_cn0_estimation = 100,
    cn0_estimator         = TrackingLoops.default_cn0_estimator(signal, num_prompts_for_cn0_estimation),
    post_corr_filter      = DefaultPostCorrFilter(),
)
```

`cn0_estimator` takes any [`AbstractCN0Estimator`](@extref TrackingLoops.AbstractCN0Estimator);
as a type parameter of [`TrackedSignal`](@ref), unlike the correlator and filter,
it may change the satellite's type (see [CN0 Estimator](cn0_estimator.md)).

An example with a narrower correlator, a trivial mean-of-antennas beamformer and
a larger C/N₀ buffer. A filter declares its combining weights rather than
combining the taps itself, so the C/N₀ path can reduce the measured noise
covariance through the same weights:

```jldoctest power_user
julia> using Tracking, TrackingLoops, GNSSSignals, StaticArrays

julia> using Tracking: Hz

julia> using TrackingLoops: NumAnts, AbstractPostCorrFilter

julia> # Trivial beamformer — averages across antenna elements
       struct MyBeamformer <: AbstractPostCorrFilter end

julia> TrackingLoops.update(f::MyBeamformer, prompt) = f;

julia> TrackingLoops.get_weights(::MyBeamformer, ::NumAnts{1}) = 1.0 + 0.0im;

julia> TrackingLoops.get_weights(::MyBeamformer, ::NumAnts{M}) where {M} =
           SVector{M,ComplexF64}(ntuple(_ -> 1 / M + 0.0im, M));

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
           num_ants                       = NumAnts(4),
           correlator                     = EarlyPromptLateCorrelator(
               num_ants = NumAnts(4),
               preferred_early_late_to_prompt_code_shift = 0.1,
           ),
           post_corr_filter               = MyBeamformer(),
           num_prompts_for_cn0_estimation = 200,
       );

julia> track_state = TrackState(GPSL1CA(), sat);

julia> get_num_ants(track_state, 1)
4

julia> get_correlator(track_state, 1).preferred_early_late_to_prompt_code_shift
0.1
```

For a multi-signal satellite, customize individual `TrackedSignal`s with the
kwarg-update constructor `TrackedSat(sat; signals = (...))`.

```@docs
TrackState
```

## Adding satellites

[`add_satellite!`](@ref) builds a [`TrackedSat`](@ref) with the default correlator and post-correlation filter from the acquisition handoff (`prn`, `code_phase`, `carrier_doppler`, optionally `code_doppler` and `carrier_phase`). Adding a PRN again overwrites it, as [`merge_sats`](@ref) does.

### Multi-satellite tracking

```jldoctest multi_sat
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> track_state = add_satellite!(track_state; prn = 1,  code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> track_state = add_satellite!(track_state; prn = 5,  code_phase = 120.0, carrier_doppler = -500.0Hz);

julia> track_state = add_satellite!(track_state; prn = 17, code_phase = 890.0, carrier_doppler = 2000.0Hz);

julia> get_carrier_doppler(track_state, 5)
-500.0 Hz

julia> get_code_phase(track_state, 17)
890.0
```

### Multi-system tracking (different signals on different sats)

Satellites with different signals go into different named groups:

```jldoctest multi_system
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (
               gps     = (GPSL1CA(),),
               galileo = (GalileoE1B(),),
           ),
       );

julia> track_state = add_satellite!(track_state; prn = 1,  group_key = :gps,     code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> track_state = add_satellite!(track_state; prn = 11, group_key = :galileo, code_phase = 200.0, carrier_doppler = -300.0Hz);

julia> get_carrier_doppler(track_state, :gps, 1)
1000.0 Hz

julia> get_carrier_doppler(track_state, :galileo, 11)
-300.0 Hz
```

### Multi-signal tracking (one satellite, several signals)

A modern GPS satellite transmits L1 C/A, L1C-D and L1C-P on one carrier; all three can be tracked together, sharing the carrier downconvert:

```jldoctest modern_gps
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (
               modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),
           ),
       );

julia> track_state = add_satellite!(track_state;
           prn = 11, group_key = :modern_gps,
           code_phase = 0.0, carrier_doppler = 1234.0Hz,
       );

julia> get_carrier_doppler(track_state, :modern_gps, 11)
1234.0 Hz
```

Put a pilot first where there is one: without data bits the PLL can integrate longer. The data signals still decode their bits, each with its own `bit_buffer`.

Signals with different code periods (L1 C/A 1 ms, L1C-P 10 ms) each complete records on their own boundaries; the conventional estimators update the Doppler from the driver's records (see [Chunked Doppler updates](track.md#Chunked-Doppler-updates)).

### Phased-array tracking

Pass a `Matrix` measurement (rows = samples, columns = antennas) and declare the antenna count:

```jldoctest phased_array
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signal = GPSL1CA(),
           num_ants = NumAnts(4),
       );

julia> track_state = add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);

julia> get_num_ants(track_state, 1)
4
```

By default the last antenna drives the discriminators. For a beamformer, build the [`TrackedSat`](@ref) with a custom `post_corr_filter`, whose [`get_weights`](@extref TrackingLoops.get_weights) are applied to the correlator and to the measured noise covariance alike, so the C/N₀ stays correct:

```jldoctest beamformer_array
julia> using Tracking, TrackingLoops, GNSSSignals, StaticArrays

julia> using Tracking: Hz

julia> using TrackingLoops: NumAnts, AbstractPostCorrFilter

julia> # Same trivial mean-of-antennas filter as the power-user example above
       struct MyBeamformer <: AbstractPostCorrFilter end

julia> TrackingLoops.update(f::MyBeamformer, prompt) = f;

julia> TrackingLoops.get_weights(::MyBeamformer, ::NumAnts{1}) = 1.0 + 0.0im;

julia> TrackingLoops.get_weights(::MyBeamformer, ::NumAnts{M}) where {M} =
           SVector{M,ComplexF64}(ntuple(_ -> 1 / M + 0.0im, M));

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        num_ants = NumAnts(4),
                        post_corr_filter = MyBeamformer());

julia> track_state = TrackState(GPSL1CA(), sat);

julia> get_num_ants(track_state, 1)
4
```

### Acquisition handoff

With `using Acquisition`, [`add_satellite!`](@ref) / [`add_satellite`](@ref) accept `AcquisitionResults` (or a vector of them, each routed independently). Without `group_key =`, an acq goes to the group whose longest-primary-code signal matches `acq.signal`.

```julia
using Acquisition  # loads the extension

# Single acq
ts = add_satellite!(ts, acq)                       # auto-route
ts = add_satellite!(ts, acq; group_key = :legacy_gps)  # explicit group, asserts match

# Vector of acqs (mixed constellations OK)
ts = add_satellite!(ts, filter(is_detected, acqs))
```

Only the longest code's phase is unambiguous for the whole group, so hand over an L1C-P acq (not L1 C/A) for `(GPSL1C_P(), GPSL1C_D(), GPSL1CA())`.

### Removing satellites

```jldoctest remove_sats
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> track_state = add_satellite!(track_state; prn = 1,  code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> track_state = add_satellite!(track_state; prn = 23, code_phase = 500.0, carrier_doppler = 1500.0Hz);

julia> track_state = remove_satellite!(track_state; prn = 1);

julia> haskey(get_sat_states(track_state, :default), 23)
true

julia> haskey(get_sat_states(track_state, :default), 1)
false
```

```@docs
add_satellite!
add_satellite
remove_satellite!
remove_satellite
merge_sats
```

## Multi-band tracking

Satellites broadcast on several RF bands (GPS L1 and L5, Galileo E1 and E5a), which arrive from separate front ends, generally at different sample rates. Each group of a multi-band `TrackState` sits on one band. Two bands give near-independent observations of one path: for ionosphere-free combinations, wider code bandwidth, and cross-band aiding, as the L1/L5 Doppler ratio is that of the carrier frequencies.

Tracking.jl provides per-band groups, measurement routing and an estimate phase that sees every band's records at once; joint cross-band estimation is up to the Doppler estimator, which receives them (see [Custom Doppler Estimator](custom_doppler_estimator.md)).

### Declaring bands

A group's band is `get_band(signals[1])`:

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

Three groups on L1 (GPS L1 and Galileo E1 share the carrier), one on L5.

### Tracking against multiple measurements

Pass one [`BandMeasurement`](@ref) per band, in a NamedTuple keyed by `GNSSSignals.get_band_id`:

```jldoctest multi_band_track; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_chip_rate, get_chip_rate_to_center_frequency_ratio

julia> function make_signal(sys, prn, carrier_doppler, num_samples, fs)
           code_freq = carrier_doppler * get_chip_rate_to_center_frequency_ratio(sys) + get_chip_rate(sys)
           range = 0:num_samples-1
           cis.(2π .* carrier_doppler .* range ./ fs) .*
               gen_code(num_samples, sys, prn, fs, code_freq, 0.0)
       end;

julia> track_state = TrackState(;
           signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)),
       );

julia> track_state = add_satellite!(track_state; prn = 1, group_key = :legacy_gps_l1, code_phase = 0.0, carrier_doppler = 200Hz);

julia> track_state = add_satellite!(track_state; prn = 1, group_key = :gps_l5, code_phase = 0.0, carrier_doppler = -150Hz);

julia> buf_l1 = make_signal(GPSL1CA(),  1, 200Hz,  4000,  4e6Hz);  # 1 ms at 4 MHz

julia> buf_l5 = make_signal(GPSL5I(),   1, -150Hz, 25000, 25e6Hz); # 1 ms at 25 MHz

julia> track!((L1 = BandMeasurement(buf_l1, 4e6Hz),
               L5 = BandMeasurement(buf_l5, 25e6Hz)), track_state);

julia> get_carrier_doppler(track_state, :legacy_gps_l1, 1)
200.00000359633913 Hz
```

`get_band_id` is `nameof` the band type, so user-defined bands have keys too. All measurements must cover **exactly the same duration**: 4000 samples at 4 MHz and 25000 at 25 MHz are both 1 ms; 25001 at 25 MHz is rejected.

### Per-band antenna counts

Pass [`SignalGroup`](@ref)s as the entries to give each group its own `num_ants`:

```jldoctest per_band_ants
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using TrackingLoops: NumAnts

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

Groups on one band share a front end and must declare the same `num_ants`; the
constructor throws otherwise.

A single-band `TrackState` can still be tracked with a bare buffer, `track!(buf, state, fs; intermediate_frequency)`.

## SignalGroup

```@docs
SignalGroup
Tracking.SignalGroups
```

### Band routing

[`get_band_ids`](@ref) lists the band ids a `TrackState` expects as measurement keys.

```@docs
get_band_ids
```

## Addressing satellites and signals

`track_state.groups[:legacy_gps].satellites` reaches a group directly. The accessors take the group key as an argument; their argument count says where the lookup stops:

| Form | Meaning |
|------|---------|
| `f(track_state)` | Single-group, single-sat — folds via `only(...)` at each level. |
| `f(track_state, prn)` | Single-group, multi-sat — picks the sat by PRN. |
| `f(track_state, group_key, prn)` | Multi-group — picks the sat in the named group. |
| `f(track_state, group_key, prn, sig)` | Per-signal — picks one [`TrackedSignal`](@ref) within a multi-signal sat. |

The trailing `sig` selector is either:

- an **`Integer`** index into `signals` (`1` is the [estimator-driver signal](#Estimator-driver-signal)), or
- a **signal type** like `GPSL1CA` (not `GPSL1CA()`), which must match exactly one signal.

### What you can read

**Sat-level — shared across all signals on the sat (no `sig` selector):**

| Accessor | Returns |
|----------|---------|
| `get_prn` | PRN number. |
| `get_num_ants` | Number of antenna elements. |
| `get_code_phase` | Shared code phase (wraps at [`max_code_wrap`](@ref Tracking.max_code_wrap)). |
| `get_code_doppler` | Shared code Doppler. |
| `get_carrier_phase` | Shared carrier phase in radians. |
| [`get_carrier_phase_correction`](@ref) | What the estimator adds to `get_carrier_phase` for the carrier-phase measurement (vector loop only; takes the band's sampling frequency). |
| `get_carrier_doppler` | Shared carrier Doppler. |
| `get_signal_start_sample` | Index of the next sample to integrate. |

**Per-signal — pass `sig` (index or type) on multi-signal sats:**

| Accessor | Returns |
|----------|---------|
| `get_correlator` | The working correlator (in-flight accumulator). |
| `get_last_fully_integrated_correlator` | Correlator value from the last completed integration. |
| `get_last_fully_integrated_filtered_prompt` | Filtered prompt value from the last completed integration. |
| `get_last_fully_integrated_num_code_blocks` | Primary-code block count of the last completed integration. |
| `get_last_fully_integrated_integration_time` | Integration time `T` of the last completed integration; the post-integration SNR is C/N₀ × `T`. |
| `get_post_corr_filter` | The post-correlation filter. |
| `get_cn0_estimator` | The CN0 estimator. |
| `estimate_cn0` | CN0 estimate in dB-Hz. |
| `get_bit_buffer` / `get_soft_bits` / `get_num_bits` | Bit buffer, decoded soft bits (sign = the hard bit), bit count. |
| `get_integrated_samples` | Number of samples accumulated into the current integration so far. |
| `has_bit_or_secondary_code_been_found` | `true` once bit/secondary-code synchronization has been achieved. |

The per-signal form always names the group, `:default` on a single-group `TrackState`: `estimate_cn0(track_state, :default, 11, GPSL1C_P)`.

```jldoctest addressing
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),),
       );

julia> track_state = add_satellite!(track_state; prn = 11, group_key = :modern_gps,
                                   code_phase = 0.0, carrier_doppler = 1234.0Hz);

julia> get_carrier_doppler(track_state, :modern_gps, 11)  # sat-level: same for all signals
1234.0 Hz

julia> get_num_bits(track_state, :modern_gps, 11, 1)  # per-signal by index (estimator-driver signal)
0

julia> get_num_bits(track_state, :modern_gps, 11, GPSL1CA)  # per-signal by type
0
```

### Group accessors

```@docs
Tracking.SatelliteDicts
get_sat_states
get_sat_state
get_signal
get_carrier_phase_correction
```

## TrackedSat

```@docs
TrackedSat
```

The sat-level accessors also take a `TrackedSat`:

```@docs
get_prn(::TrackedSat)
get_code_phase(::TrackedSat)
get_code_doppler(::TrackedSat)
get_carrier_phase(::TrackedSat)
get_carrier_doppler(::TrackedSat)
get_signal_start_sample(::TrackedSat)
get_signals(::TrackedSat)
get_doppler_estimator_state(::TrackedSat)
```

## TrackedSignal

```@docs
TrackedSignal
get_correlator(::TrackedSignal)
get_last_fully_integrated_correlator(::TrackedSignal)
get_last_fully_integrated_filtered_prompt(::TrackedSignal)
get_last_fully_integrated_num_code_blocks(::TrackedSignal)
get_last_fully_integrated_integration_time(::TrackedSignal)
get_filtered_prompts(::TrackedSignal)
get_bit_buffer(::TrackedSignal)
get_num_bits(::TrackedSignal)
get_integrated_samples(::TrackedSignal)
get_preferred_num_code_blocks_to_integrate(::TrackedSignal)
```

The other per-signal accessors also take a `TrackedSignal`, as does `get_signal(tsig)`, the signal instance.

Bit and secondary-code sync, and how it widens the code-phase wrap, have their own page: [Bit and Secondary-Code Sync](bit_sync.md). Most users only need `has_bit_or_secondary_code_been_found` to gate `get_soft_bits`.

## Handing over to PositionVelocityTime

With the navigation data GNSSDecoder.jl decodes from their soft bits, tracked
satellites become the measurements `PositionVelocityTime.calc_pvt` solves. Pass a
`NamedTuple` of per-group PRN → `GNSSDecoderState` mappings:

```julia
using PositionVelocityTime: PositionVelocityTime, calc_pvt
pvt = calc_pvt(PositionVelocityTime.signal_groups(track_state, decoders))
```

```@docs
Tracking.PositionVelocityTime.signal_groups(::TrackState, ::NamedTuple)
Tracking.PositionVelocityTime.SignalGroup(::Tracking.SignalGroup, ::Any)
Tracking.PositionVelocityTime.SatelliteState(::Tracking.GNSSDecoderState, ::Tracking.AbstractGNSSSignal, ::TrackedSat)
```
