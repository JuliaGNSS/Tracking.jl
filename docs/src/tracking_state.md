# Tracking State

Tracking.jl uses a small hierarchy of state types to manage tracking across multiple satellites, multiple signals per satellite, multiple signal groups, and multiple RF bands.

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
    add_satellite!(track_state, filter(is_detected, acqs))   # routes each acq to the matching group
    track_state = track(latest_chunk, track_state, sampling_frequency)
end
```

This pattern keeps the `TrackState`'s concrete type fixed across the loop — the satellite-dict's slot type is frozen at construction, so the tracking hot path stays type-stable as sats come and go.

### Power-user: pre-built `TrackedSat`s

If you need to customize the correlator or post-correlation filter type (the slot type itself), build the `TrackedSat`s yourself and hand them to the legacy constructor `TrackState(signal, sats)` or the `add_satellite!(track_state, group, sat)` escape hatch. The kwarg-based constructors only let you customize the satellite's *values*, not its concrete type.

## TrackState

The main container for all tracking state. It holds:

- `groups` — a `NamedTuple` of [`SignalGroup`](@ref)s. Each group bundles its satellites dictionary, signal-instance tuple, band, and antenna count.
- `doppler_estimator` — the Doppler estimator configuration (e.g. PLL/DLL bandwidths).

To reach per-group state, index `track_state.groups` by the group's key (e.g. `track_state.groups[:legacy_gps].satellites`). The high-level accessors below ([`get_sat_states`](@ref), [`get_sat_state`](@ref), …) take the group key as an argument and fold to compile-time constants when the groups type is known.

```@docs
TrackState
```

## SignalGroup

A group of satellites that all track the same tuple of GNSS signal types, on the same RF band, observed by the same antenna array. Groups are the unit of type stability — every `TrackedSat` inside a `SignalGroup` shares the same concrete signal-tuple shape, so the satellites dictionary has a concrete value type and the hot loop sees no dynamic dispatch.

Two groups may share a band: e.g. a `:legacy_gps` group tracking `(GPSL1CA(),)` and a `:galileo` group tracking `(GalileoE1B(),)` both report `band = L1()`. The grouping is by signal-tuple shape, not by band — `band` is metadata each group carries so `track` can route the right measurement to it.

```@docs
SignalGroup
SignalGroups
```

## Band routing

The mapping between GNSSSignals `Band` instances and the `Symbol` keys used in multi-band measurement collections (see [`Measurement`](@ref)).

```@docs
band_key
band_keys
```

## Addressing satellites and signals

Tracking state nests as **TrackState → SignalGroup → TrackedSat → TrackedSignal**, and the accessor functions accept a matching ladder of arguments. The accessor's argument count tells the lookup where to stop:

| Form | Meaning |
|------|---------|
| `f(track_state)` | Single-group, single-sat — folds via `only(...)` at each level. |
| `f(track_state, prn)` | Single-group, multi-sat — picks the sat by PRN. |
| `f(track_state, group, prn)` | Multi-group — picks the sat in the named group. |
| `f(track_state, group, prn, sig)` | Per-signal — picks one [`TrackedSignal`](@ref) within a multi-signal sat. |

The trailing `sig` selector is either:

- an **`Integer`** index into the sat's `signals` tuple (`1` = first signal = Doppler source) — the canonical form, unambiguous even when the same signal type appears twice in the tuple, or
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
| `get_bit_buffer` / `get_bits` / `get_num_bits` | Bit buffer, decoded bits, bit count. |
| `get_integrated_samples` | Number of samples accumulated into the current integration so far. |
| `has_bit_or_secondary_code_been_found` | `true` once bit/secondary-code synchronization has been achieved. |

The per-signal form always names the group explicitly, even on a single-group `TrackState`. Use `:default` as the group key in that case — `estimate_cn0(track_state, :default, 11, GPSL1C_P)`.

```jldoctest addressing
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (modern_gps = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),),
       );

julia> add_satellite!(track_state; prn = 11, group = :modern_gps,
                                   code_phase = 0.0, carrier_doppler = 1234.0Hz);

julia> get_carrier_doppler(track_state, :modern_gps, 11)  # sat-level: same for all signals
1234.0 Hz

julia> get_bits(track_state, :modern_gps, 11, 1)  # per-signal by index (Doppler source)
0x00000000000000000000000000000000

julia> get_bits(track_state, :modern_gps, 11, GPSL1CA)  # per-signal by type
0x00000000000000000000000000000000
```

## TrackedSat

Per-satellite tracking state. Carries the shared carrier and code Doppler and phase (one set of values per satellite, since all signals on a satellite share the same carrier), a tuple of per-signal [`TrackedSignal`](@ref) states, and the per-satellite Doppler-estimator state.

```@docs
TrackedSat
```

In addition to the accessors listed under [Addressing satellites and signals](#Addressing-satellites-and-signals), a `TrackedSat` exposes:

- `get_signals(sat)` — the tuple of [`TrackedSignal`](@ref)s.
- `get_doppler_estimator_state(sat)` — the per-satellite Doppler-estimator state (e.g. the loop-filter integrators for the conventional PLL/DLL).

## TrackedSignal

Per-signal tracking state. Carries the correlator, post-correlation filter, CN0 estimator, bit buffer, and integration-progress flags for a single signal on a satellite. A multi-signal `TrackedSat` (e.g. one tracking GPS L1 C/A + L1C-D + L1C-P) carries one `TrackedSignal` per signal in `signals::Tuple{Vararg{TrackedSignal}}`.

The first signal in the tuple is the Doppler source — its correlator is what the PLL/DLL discriminator runs on.

```@docs
TrackedSignal
```

The per-signal accessors in the table under [Addressing satellites and signals](#Addressing-satellites-and-signals) all dispatch directly on a `TrackedSignal` too. Additionally:

- `get_signal(tsig)` — the GNSS signal instance (e.g. `GPSL1CA()`).
- `get_filtered_prompts(tsig)` — every filtered prompt produced during the most recent `track` call. The vector is reset at the start of each call and appended for every completed integration.

### Code-phase wrap period

The shared `TrackedSat.code_phase` wraps at two distinct values depending on the per-signal sync state:

- **Before any signal has synced** — wrap is the largest primary code length across the sat's signals. For a sat tracking only L1 C/A this is 1023 chips. The wrap is intentionally narrow at this stage because we don't yet know which data bit or secondary-code chip the current primary period belongs to.

- **After a signal syncs** — its contribution to the wrap widens to one full *symbol period*: `primary × secondary_code_length` for a pilot, or `primary × blocks_per_data_bit` for a data-bearing signal. The shared wrap is then the `max` across signals, so the longest synced signal pins it.

Concrete values:

| Sat tracks         | Before sync | After full sync                          |
|--------------------|-------------|------------------------------------------|
| GPS L1 C/A only    | 1023        | 1023 × 20 = **20460** (one 50 Hz data bit) |
| GPS L5I only       | 10230       | 10230 × 10 = **102300** (one NH10 cycle / one data bit) |
| GPS L1C-P only     | 10230       | 10230 × 1800 ≈ **18.4 M** (one overlay-code cycle, ≈ 18 s) |
| L1C-P + L1C-D + L1CA | 10230 (longest primary) | 18.4 M (L1C-P dominates) |

The post-sync widening is what lets downstream consumers (e.g. PositionVelocityTime.jl) distinguish which primary-code period within the symbol they're currently in — `mod(code_phase, primary)` gives the per-signal replica phase, while `div(code_phase, primary)` gives the symbol-internal position.

[`max_code_length`](@ref) returns the *upper bound* (the post-full-sync value) at compile time. [`current_code_wrap`](@ref) returns the *runtime* value honoring the current per-signal sync state — this is what the inner loop's `mod` actually uses.

```@docs
max_code_length
current_code_wrap
```

### Bit-sync and secondary-code-sync detection

Each per-signal `BitBuffer` runs an `is_upcoming_integration_new_bit` detector against the running buffer of primary-code-block signs. The detector returns a `SyncResult` containing whether sync was found, the secondary-code phase (chip offset within the secondary code, used for code-phase seeding — see below), and the locked polarity (±1).

Per-signal contract. *Min-to-fire* is the smallest `num_code_blocks` the detector accepts before it tries a match (smaller windows return `SyncResult(false, …)` without searching). *Buffer width* is the sliding window's container type, picked by `get_code_block_buffer_type(signal)`; it's at least `2 × min-to-fire` for the symmetric-template signals so a full bit/symbol period can slide past either polarity.

| Signal | Min-to-fire | Buffer width | Template | Tolerance | Phase | Blocks per symbol |
|--------|-------------|--------------|----------|-----------|-------|---------------------|
| GPS L1 C/A | 40 blocks | `UInt64` (40 bits used) | bit-edge `0xfffff` (20+20) | 3 errors | 0 | 20 |
| Galileo E1B | 8 blocks | `UInt8` (8 bits) | bit-edge `0x0f` (4+4) | 1 error | 0 | 1 |
| GPS L5I | 10 blocks | `UInt32` (20 bits used, 2 × NH10) | NH10 `0x035` shared | 2 errors | 0 | 10 (secondary code) |
| GPS L1C-D | n/a | `UInt8` (unused) | trivial (1 block per symbol) | n/a | 0 | 1 |
| GPS L1C-P | 1800 blocks | `UInt1800` (exact width) | per-PRN overlay | 36 errors (2 %) | `0..1799` | n/a (pilot) |

The buffer-width type threads through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF}` as a type parameter. The L1C-P case uses an exact-width `UInt1800` defined via `BitIntegers.@define_integers 1800`; the other signals use built-in `UInt8` / `UInt32` / `UInt64`.

#### Lifecycle of a `BitBuffer`

Two distinct phases, separated by the `found::Bool` flag:

1. **Pre-sync search** (`found = false`). Each completed integration shifts one bit (the sign of the prompt's real part) into `code_block_buffer::B`. `is_upcoming_integration_new_bit` is called on every shift; while it returns `SyncResult(false, ...)` the loop keeps integrating one primary code period at a time. By construction, the detector only matches at a true bit/symbol boundary — so the very call that transitions `found` to `true` lands exactly on the edge between two data symbols.

2. **Post-sync accumulation** (`found = true`). The post-sync branch in `Tracking.buffer` ignores `code_block_buffer` and instead accumulates the complex prompt into `prompt_accumulator`. Each integration also bumps `prompt_accumulator_integrated_code_blocks`. Once that counter reaches the per-signal "blocks per symbol" value above — `_calc_num_code_blocks_that_form_a_bit(signal) = get_code_frequency(signal) / (get_code_length(signal) * get_data_frequency(signal))` — one decoded bit is committed to `buffer::UInt128` and the accumulator resets to zero. For Galileo E1B and GPS L1C-D the counter is `1`, so one symbol commits per integration; for GPS L1 C/A it's `20`, so the loop counts 20 primary-code periods (≈ 20 × 1023 chips = 20460 chips of `code_phase` advance, modulo wrap) per data bit; for GPS L5I it's `10`. The polarity flag flips the accumulator's sign at commit time when the detector locked at negative polarity, so downstream consumers always see `1 = data symbol 0`.

Pilot signals (`get_data_frequency = 0 Hz`, e.g. GPS L1C-P) never enter the post-sync accumulation branch — their `bit_buffer` carries the recovered secondary-code phase but no decoded bits, and the post-sync work is purely the [code-phase seeding](#Code-phase-wrap-period) described next.

#### Code-phase seeding from the secondary-code phase

When a signal whose detector exposes a non-zero `secondary_phase` syncs — currently only GPS L1C-P — that phase is used to seed `TrackedSat.code_phase` so subsequent wrap-mod-[`current_code_wrap`](@ref) arithmetic gives the absolute position in the longest secondary-code cycle. The seeding follows a fallback chain: the synced signal with the largest `(primary × secondary)` code length wins.

Signals with `secondary_code_length == 1` (bit-edge only, e.g. GPS L1 C/A) do **not** carry an explicit `secondary_phase` to snap — there's no per-PRN overlay to recover a chip offset from. Instead, their post-sync bit-edge alignment is captured by **two complementary mechanisms**:

1. The `BitBuffer.prompt_accumulator_integrated_code_blocks` counter tracks "how many of the next N primary periods have I integrated since the last bit commit." `reset(bit_buffer)` preserves it, so the bit cadence survives intra-call resets without re-syncing.

2. The wrap returned by [`current_code_wrap`](@ref) widens from `primary` to `primary × blocks_per_data_bit` (e.g. 1023 → 20460 for L1 C/A) the moment `bit_buffer.found` flips to `true`. From that point on `mod(code_phase, primary)` continues to give the replica phase, while `div(code_phase, primary)` reads off which primary period within the data bit we're in. The transition is one-shot: on the call that flips `found`, `code_phase` is implicitly the start of a fresh data bit because the detector only matches at a true bit boundary (the `[0, primary)` range of `code_phase` then represents primary period 0 of the new bit), so no explicit snap is needed.

## Group accessors

```@docs
SatelliteDicts
get_sat_states
get_sat_state
get_signal
```

## Adding and removing satellites

```@docs
add_satellite!
add_satellite
remove_satellite!
remove_satellite
merge_sats
```

### Acquisition handoff

When the [Acquisition.jl](https://github.com/JuliaGNSS/Acquisition.jl) extension is loaded (via `using Acquisition`), [`add_satellite!`](@ref) / [`add_satellite`](@ref) gain `AcquisitionResults` overloads that read `prn` / `code_phase` / `carrier_doppler` straight off the acq result. With `group = nothing` (the default) the routing is inferred by matching `acq.system` against each group's longest-primary-code signal; pass an explicit `group =` to bypass the inference. The batch form takes an `AbstractVector{<:AcquisitionResults}` and routes each entry independently — convenient for the `filter(is_detected, acquire(...))` pipeline.

```julia
using Acquisition  # loads the extension

# Single acq
add_satellite!(ts, acq)                       # auto-route
add_satellite!(ts, acq; group = :legacy_gps)  # explicit group, asserts match

# Vector of acqs (mixed constellations OK)
add_satellite!(ts, filter(is_detected, acqs))
```

`acq.system` must match the **longest-primary-code** signal in the target group's tuple — its code phase is the only one that's unambiguous when the group tracks multiple signals on shared chips. Hand over an L1C-P acq (not L1 C/A) for a group tracking `(GPSL1C_P(), GPSL1C_D(), GPSL1CA())`.
