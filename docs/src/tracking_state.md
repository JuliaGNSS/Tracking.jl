# Tracking State

Tracking.jl uses a small hierarchy of state types to manage tracking across multiple satellites, multiple signals per satellite, multiple signal groups, and multiple RF bands.

The tracking state nests as **TrackState → SignalGroup → TrackedSat → TrackedSignal**:

- [`TrackState`](@ref) — top-level container; holds a `NamedTuple` of [`SignalGroup`](@ref)s plus the Doppler-estimator configuration.
- [`SignalGroup`](@ref) — named group of sats that share the same signal-tuple shape (and therefore the same concrete `TrackedSat` value type, which is what gives the hot loop type stability). Each group also carries its RF band and antenna count.
- [`TrackedSat`](@ref) — per-satellite state: shared carrier/code Doppler and phase (one set of values per satellite, since all signals on a satellite share the same carrier), a tuple of [`TrackedSignal`](@ref)s, and the per-satellite Doppler-estimator state.
- [`TrackedSignal`](@ref) — per-signal state: correlator, post-correlation filter, CN0 estimator, bit buffer, and integration-progress flags.

### Estimator-driver signal

The first signal in each group's tuple is the **estimator-driver signal** — the one the Doppler estimator uses to update the satellite-shared carrier and code Doppler. With the default [`ConventionalPLLAndDLL`](@ref) / [`ConventionalAssistedPLLAndDLL`](@ref), `signals[1]`'s correlator is the input to the PLL/DLL discriminator, and the per-signal default loop bandwidths are sized off this signal's primary-code period. A user-supplied [`AbstractDopplerEstimator`](@ref) is free to use the other signals' state too — `signals[1]`'s privileged role is a convention of the conventional estimators, not a structural constraint of `TrackedSat`.

## Choosing a `TrackState` constructor

`TrackState` has several constructors. The right choice depends on **when you know which satellites you'll track**.

### One-shot scripts: acquire once, then track

If you acquire once at the start of a script and hand the results straight to tracking, use the `AcquisitionResults`-aware constructors from the [Acquisition](https://github.com/JuliaGNSS/Acquisition.jl) extension. They derive everything (signal, default loop bandwidths, satellite parameters) from the acquisition results, so the whole acquire→track handoff is one line.

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

You can still call [`add_satellite!`](@ref) on a `TrackState` built this way later — but only for signals already declared by the acqs you handed to the constructor (the groups, and therefore the slot types, are frozen at construction). If you anticipate tracking a wider set of signals than the initial acquisition produced, use the empty-construct-then-populate pattern below instead.

### Real-time / repeating loops: build empty, then populate

If you re-acquire periodically (a typical receiver: re-search PRNs every few seconds, hand new detections off to tracking without rebuilding the whole `TrackState`), build an **empty** `TrackState` once with `TrackState(; signal = ...)` or `TrackState(; signals = (...))`, then add satellites later with [`add_satellite!`](@ref) (or remove them with [`remove_satellite!`](@ref)):

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

This pattern keeps the `TrackState`'s concrete type fixed across the loop — the satellite-dict's slot type is frozen at construction, so the tracking hot path stays type-stable as sats come and go.

The singular `signal = GPSL1CA()` keyword is the shortcut for the common one-group, one-signal case. It desugars internally to `signals = (default = (GPSL1CA(),),)`, so the rest of the API can stay uniform. With one group, [`add_satellite!`](@ref) may omit the `group =` keyword.

### Power-user: pre-built `TrackedSat`s

If you need to customize the correlator or post-correlation filter type (the slot type itself), build the `TrackedSat`s yourself and hand them to the positional constructor `TrackState(signal, sats)` or the `add_satellite!(track_state, group, sat)` escape hatch. The kwarg-based constructors only let you customize the satellite's *values*, not its concrete type.

The single-signal `TrackedSat` constructor surface:

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
    code_doppler          = carrier_doppler * get_code_center_frequency_ratio(signal),
    num_prompts_for_cn0_estimation = 100,
    post_corr_filter      = DefaultPostCorrFilter(),
)
```

A worked example combining a narrower-than-default correlator, a custom
post-correlation filter (beamformer), and a larger CN0 buffer. The
beamformer here is a trivial mean-of-antennas — a real receiver would
plug in an actual beamforming algorithm:

```jldoctest power_user
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz, NumAnts, AbstractPostCorrFilter

julia> # Trivial beamformer — averages across antenna elements
       struct MyBeamformer <: AbstractPostCorrFilter end

julia> Tracking.update(f::MyBeamformer, prompt) = f;

julia> (::MyBeamformer)(x::AbstractVector) = sum(x) / length(x);

julia> (::MyBeamformer)(x) = x;

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

For a multi-signal satellite, the empty `TrackState(; signals = (group =
(sig1, sig2, …),))` path will build a default template sat the first time
you `add_satellite!` to that group; customize individual `TrackedSignal`s
afterwards via the [`TrackedSat`](@ref) kwarg-update constructor
(`TrackedSat(sat; signals = (...))`).

```@docs
TrackState
```

## Adding satellites

Satellites are added to a `TrackState` via [`add_satellite!`](@ref). The acquisition handoff values (`prn`, `code_phase`, `carrier_doppler`, optionally `code_doppler` and `carrier_phase`) get wired into a fresh [`TrackedSat`](@ref) with the library's default correlator and post-correlation filter. Adding a satellite with the same PRN again overwrites the existing entry (matching [`merge_sats`](@ref) semantics — no error).

### Multi-satellite tracking

To track several satellites on the same signal, simply call `add_satellite!` repeatedly:

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

julia> track_state = add_satellite!(track_state; prn = 1,  group = :gps,     code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> track_state = add_satellite!(track_state; prn = 11, group = :galileo, code_phase = 200.0, carrier_doppler = -300.0Hz);

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

julia> track_state = add_satellite!(track_state;
           prn = 11, group = :modern_gps,
           code_phase = 0.0, carrier_doppler = 1234.0Hz,
       );

julia> get_carrier_doppler(track_state, :modern_gps, 11)
1234.0 Hz
```

Putting a pilot signal first (e.g. `GPSL1C_P()`) is encouraged with the conventional estimators when one is available: pilot signals carry no data-bit modulation, which lets the PLL run longer coherent integrations and reach lower phase-noise floors. The data-bearing signals (L1C-D, L1 C/A) still recover their navigation bits independently — each [`TrackedSignal`](@ref) carries its own `bit_buffer` regardless of which signal drives the estimator.

When a satellite tracks signals with different primary-code lengths (e.g. L1 C/A at 1 ms vs L1C-P at 10 ms), each outer iteration integrates to the **shortest** signal's next primary-code boundary. The shorter signal's correlator completes every iteration; the longer signal's correlator accumulates across multiple iterations and only marks `is_integration_completed = true` on its own boundary. Doppler updates therefore happen at the shortest signal's cadence (1 ms in this example), and longer signals see their integration windows spanned by piecewise Doppler updates — the natural per-iteration-Doppler-correction behaviour of a real receiver.

### Phased-array tracking

To track signals coherently across an antenna array, pass a `Matrix` measurement (rows = samples, columns = antenna elements) and declare the number of antennas at `TrackState` construction:

```jldoctest phased_array
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signal = GPSL1CA(),
           num_ants = NumAnts(4),
       );

julia> track_state = add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);

julia> get_num_ants(track_state, 1)
4
```

By default the track function uses the first antenna channel as the reference signal to drive the discriminators. An appropriate beamforming algorithm will probably suit better — construct a [`TrackedSat`](@ref) with a custom `post_corr_filter` and build the `TrackState` from it (so the slot type takes the custom filter type rather than the default):

```jldoctest beamformer_array
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz, NumAnts, AbstractPostCorrFilter

julia> # Same trivial mean-of-antennas filter as the power-user example above
       struct MyBeamformer <: AbstractPostCorrFilter end

julia> Tracking.update(f::MyBeamformer, prompt) = f;

julia> (::MyBeamformer)(x::AbstractVector) = sum(x) / length(x);

julia> (::MyBeamformer)(x) = x;

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        num_ants = NumAnts(4),
                        post_corr_filter = MyBeamformer());

julia> track_state = TrackState(GPSL1CA(), sat);

julia> get_num_ants(track_state, 1)
4
```

### Acquisition handoff

When the [Acquisition.jl](https://github.com/JuliaGNSS/Acquisition.jl) extension is loaded (via `using Acquisition`), [`add_satellite!`](@ref) / [`add_satellite`](@ref) gain `AcquisitionResults` overloads that read `prn` / `code_phase` / `carrier_doppler` straight off the acq result. With `group = nothing` (the default) the routing is inferred by matching `acq.system` against each group's longest-primary-code signal; pass an explicit `group =` to bypass the inference. The batch form takes an `AbstractVector{<:AcquisitionResults}` and routes each entry independently — convenient for the `filter(is_detected, acquire(...))` pipeline.

```julia
using Acquisition  # loads the extension

# Single acq
ts = add_satellite!(ts, acq)                       # auto-route
ts = add_satellite!(ts, acq; group = :legacy_gps)  # explicit group, asserts match

# Vector of acqs (mixed constellations OK)
ts = add_satellite!(ts, filter(is_detected, acqs))
```

`acq.system` must match the **longest-primary-code** signal in the target group's tuple — its code phase is the only one that's unambiguous when the group tracks multiple signals on shared chips. Hand over an L1C-P acq (not L1 C/A) for a group tracking `(GPSL1C_P(), GPSL1C_D(), GPSL1CA())`.

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

A satellite often broadcasts on more than one RF band — GPS broadcasts on L1 (1575.42 MHz) and L5 (1176.45 MHz); Galileo broadcasts on E1 (L1) and E5a (L5). In a multi-band receiver these arrive from separate front-ends, generally at different sample rates, and need to be downconverted and correlated against their own carrier replicas. Tracking.jl exposes this as a multi-band `TrackState` where each group declares which RF band it sits on.

**Why this matters.** A single physical satellite tracked on two bands gives the receiver two near-independent observations of the same path. The classic uses:

- **Ionospheric correction** via dual-frequency (iono-free) pseudorange combinations.
- **Wider effective bandwidth** for code-phase observations (L5 carries far more chip-rate bandwidth than L1 C/A).
- **Cross-band-aided tracking**: the carrier Doppler ratio between L1 and L5 is exactly the ratio of their RF carrier frequencies. A joint estimator can fuse the two bands' discriminators and produce a more accurate Doppler estimate than either band alone — particularly valuable at low CN0 where the wider data-aided integration on L5 helps the noisier L1 C/A.

This release ships the **structural enablers** for multi-band: per-band groups, per-band measurement routing, an estimation barrier that sees every band's correlator outputs at once. The cross-band joint-tracking *algorithm* (e.g. linking PRN-X-on-L1 with PRN-X-on-L5 in one estimator step) is a follow-up — see [docs/plans/2026-05-15-multi-band-tracking-design.md](https://github.com/JuliaGNSS/Tracking.jl/blob/master/docs/plans/2026-05-15-multi-band-tracking-design.md) for the design and the open mechanism question.

### Declaring bands

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

### Tracking against multiple measurements

For multi-band tracking, build one [`BandMeasurement`](@ref) per band — bundling sample buffer and front-end metadata — and pass them as a NamedTuple keyed by the band's `GNSSSignals.get_band_id` (e.g. `:L1`, `:L5`):

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

julia> track_state = add_satellite!(track_state; prn = 1, group = :legacy_gps_l1, code_phase = 0.0, carrier_doppler = 200Hz);

julia> track_state = add_satellite!(track_state; prn = 1, group = :gps_l5, code_phase = 0.0, carrier_doppler = -150Hz);

julia> buf_l1 = make_signal(GPSL1CA(),  1, 200Hz,  4000,  4e6Hz);  # 1 ms at 4 MHz

julia> buf_l5 = make_signal(GPSL5I(),   1, -150Hz, 25000, 25e6Hz); # 1 ms at 25 MHz

julia> track!((L1 = BandMeasurement(buf_l1, 4e6Hz),
               L5 = BandMeasurement(buf_l5, 25e6Hz)), track_state);

julia> get_carrier_doppler(track_state, :legacy_gps_l1, 1)
200.00000359633913 Hz
```

The keys (`:L1`, `:L5`) come from `GNSSSignals.get_band_id(L1())` and `GNSSSignals.get_band_id(L5())` — `nameof` of the band type, so any band (including user-defined ones) has a key without Tracking-side registration. All measurements must cover the **exact same observation duration** — `num_samples / sampling_frequency` must compare equal across bands. An L1 chunk of 4000 samples at 4 MHz and an L5 chunk of 25000 samples at 25 MHz both cover 1 ms, so they're compatible; an L5 chunk of 25001 samples is rejected.

### Per-band antenna counts

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

### Bare-buffer compatibility

A single-band receiver doesn't need to type any of this. The bare-buffer call `track!(buf, state, fs)` keeps working for any `TrackState` that spans exactly one band — internally it wraps the buffer into a one-entry NamedTuple keyed by the lone band. Pass `intermediate_frequency` via the same kwarg as before, or move it onto a [`BandMeasurement`](@ref) when you migrate to multi-band.

## SignalGroup

A group of satellites that all track the same tuple of GNSS signal types, on the same RF band, observed by the same antenna array. Groups are the unit of type stability — every `TrackedSat` inside a `SignalGroup` shares the same concrete signal-tuple shape, so the satellites dictionary has a concrete value type and the hot loop sees no dynamic dispatch.

Two groups may share a band: e.g. a `:legacy_gps` group tracking `(GPSL1CA(),)` and a `:galileo` group tracking `(GalileoE1B(),)` both report `band = L1()`. The grouping is by signal-tuple shape, not by band — `band` is metadata each group carries so `track` can route the right measurement to it.

```@docs
SignalGroup
SignalGroups
```

### Band routing

The `Symbol` keys used in multi-band measurement collections (see [`BandMeasurement`](@ref)) are the bands' ids as reported by `GNSSSignals.get_band_id`; [`band_keys`](@ref) lists the ids a given `TrackState` expects.

```@docs
band_keys
```

## Addressing satellites and signals

To reach per-group state, index `track_state.groups` by the group's key (e.g. `track_state.groups[:legacy_gps].satellites`). The high-level accessors below ([`get_sat_states`](@ref), [`get_sat_state`](@ref), …) take the group key as an argument and fold to compile-time constants when the groups type is known.

The accessor's argument count tells the lookup where to stop:

| Form | Meaning |
|------|---------|
| `f(track_state)` | Single-group, single-sat — folds via `only(...)` at each level. |
| `f(track_state, prn)` | Single-group, multi-sat — picks the sat by PRN. |
| `f(track_state, group, prn)` | Multi-group — picks the sat in the named group. |
| `f(track_state, group, prn, sig)` | Per-signal — picks one [`TrackedSignal`](@ref) within a multi-signal sat. |

The trailing `sig` selector is either:

- an **`Integer`** index into the sat's `signals` tuple (`1` = first signal = [estimator-driver signal](#Estimator-driver-signal)) — the canonical form, unambiguous even when the same signal type appears twice in the tuple, or
- a **signal type** like `GPSL1CA` (the bare type, not `GPSL1CA()`) — readable sugar that errors if the type appears zero or more than once in the tuple.

### What you can read

**Sat-level — shared across all signals on the sat (no `sig` selector):**

| Accessor | Returns |
|----------|---------|
| `get_prn` | PRN number. |
| `get_num_ants` | Number of antenna elements. |
| `get_code_phase` | Shared code phase (wraps at [`max_code_length`](@ref)). |
| `get_code_doppler` | Shared code Doppler. |
| `get_carrier_phase` | Shared carrier phase in radians. |
| `get_carrier_doppler` | Shared carrier Doppler. |
| `get_signal_start_sample` | Index of the next sample to integrate. |

**Per-signal — pass `sig` (index or type) on multi-signal sats:**

| Accessor | Returns |
|----------|---------|
| `get_correlator` | The working correlator (in-flight accumulator). |
| `get_last_fully_integrated_correlator` | Correlator value from the last completed integration. |
| `get_last_fully_integrated_filtered_prompt` | Filtered prompt value from the last completed integration. |
| `get_post_corr_filter` | The post-correlation filter. |
| `get_cn0_estimator` | The CN0 estimator. |
| `estimate_cn0` | CN0 estimate in dB-Hz. |
| `get_bit_buffer` / `get_soft_bits` / `get_num_bits` | Bit buffer, decoded complex soft bits, bit count. |
| `get_integrated_samples` | Number of samples accumulated into the current integration so far. |
| `has_bit_or_secondary_code_been_found` | `true` once bit/secondary-code synchronization has been achieved. |

The per-signal form always names the group explicitly, even on a single-group `TrackState`. Use `:default` as the group key in that case — `estimate_cn0(track_state, :default, 11, GPSL1C_P)`.

```jldoctest addressing
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),),
       );

julia> track_state = add_satellite!(track_state; prn = 11, group = :modern_gps,
                                   code_phase = 0.0, carrier_doppler = 1234.0Hz);

julia> get_carrier_doppler(track_state, :modern_gps, 11)  # sat-level: same for all signals
1234.0 Hz

julia> get_soft_bits(track_state, :modern_gps, 11, 1)  # per-signal by index (estimator-driver signal)
ComplexF32[]

julia> get_soft_bits(track_state, :modern_gps, 11, GPSL1CA)  # per-signal by type
ComplexF32[]
```

### Group accessors

```@docs
SatelliteDicts
get_sat_states
get_sat_state
get_signal
```

## TrackedSat

```@docs
TrackedSat
```

The sat-level accessors listed under [Addressing satellites and signals](#Addressing-satellites-and-signals) all have a single-argument form that dispatches on a `TrackedSat` directly:

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
```

The per-signal accessors in the table under [Addressing satellites and signals](#Addressing-satellites-and-signals) all dispatch directly on a `TrackedSignal` too. Additionally:

- `get_signal(tsig)` — the GNSS signal instance (e.g. `GPSL1CA()`).
- `get_filtered_prompts(tsig)` — every filtered prompt produced during the most recent `track` call. The vector is reset at the start of each call and appended for every completed integration.

### Bit sync, secondary-code sync, and the code-phase wrap

The mechanics of bit / secondary-code synchronization — including the
runtime widening of `TrackedSat.code_phase` from one primary code period
to one full symbol period, the per-signal `BitBuffer` lifecycle, and the
code-phase seeding from a recovered secondary-code phase — have a
dedicated page: [Bit and Secondary-Code Sync](bit_sync.md).

For day-to-day use, the only thing most users need is
`has_bit_or_secondary_code_been_found` (per signal) to gate calls to
`get_soft_bits` — both are listed in the [per-signal accessor table](#What-you-can-read).
