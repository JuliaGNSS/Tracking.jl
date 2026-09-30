# Tracking State

Tracking.jl uses a small hierarchy of state types to manage tracking across multiple satellites, multiple signals per satellite, multiple signal groups, and multiple RF bands.

The tracking state nests as **TrackState → SignalGroup → TrackedSat → TrackedSignal**:

- [`TrackState`](@ref) — top-level container; holds a `NamedTuple` of [`SignalGroup`](@ref)s plus the Doppler-estimator configuration.
- [`SignalGroup`](@ref) — named group of sats that share the same signal-tuple shape (and therefore the same concrete `TrackedSat` value type, which is what gives the hot loop type stability). Each group also carries its RF band and antenna count.
- [`TrackedSat`](@ref) — per-satellite state: shared carrier/code Doppler and phase (one set of values per satellite, since all signals on a satellite share the same carrier), a tuple of [`TrackedSignal`](@ref)s, and the per-satellite Doppler-estimator state.
- [`TrackedSignal`](@ref) — per-signal state: correlator, post-correlation filter, CN0 estimator, bit buffer, and integration-progress flags.

### Estimator-driver signal

The first signal in each group's tuple is the **estimator-driver signal**. With the default [`ConventionalPLLAndDLL`](@ref) / [`ConventionalAssistedPLLAndDLL`](@ref) it sets the loop *cadence* (Dopplers are updated when it completes an integration), the loop *bandwidths* (sized off its primary-code period), the *carrier-phase reference* the loops lock onto, and the datum of the satellite-shared `code_phase`.

What the driver does **not** have to do is monopolise the measurements: with `discriminator_combining = true` the other signals' discriminator outputs are folded into the driver's loop update wherever their records coincide with the driver's — see [Multi-signal discriminator combining](#Multi-signal-discriminator-combining). A user-supplied [`AbstractDopplerEstimator`](@ref) is free to use the signals' state any way it likes; `signals[1]`'s privileged role is a convention of the shipped estimators, not a structural constraint of `TrackedSat`.

The driver signal is privileged for the Doppler estimator only. Bit synchronisation, the post-correlation filter and the **CN0 estimator** all run per signal, so a multi-signal satellite produces one C/N₀ per signal rather than one for the driver — see [CN0 Estimator](cn0_estimator.md) for what that costs and for [`NoCN0Estimator`](@ref), the per-signal opt-out.

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
    cn0_estimator         = Tracking.default_cn0_estimator(signal, num_prompts_for_cn0_estimation),
    post_corr_filter      = DefaultPostCorrFilter(),
)
```

`cn0_estimator` takes any [`AbstractCN0Estimator`](@ref) — unlike the correlator
and the post-corr filter it *is* free to change the sat's concrete type, since it
is a type parameter of [`TrackedSignal`](@ref). See
[CN0 Estimator](cn0_estimator.md) for the estimators that ship with Tracking and
for writing your own.

A worked example combining a narrower-than-default correlator, a custom
post-correlation filter (beamformer), and a larger CN0 buffer. The
beamformer here is a trivial mean-of-antennas — a real receiver would
plug in an actual beamforming algorithm. A filter declares its combining
weights rather than combining the taps itself, which is what lets the C/N₀
path reduce the measured noise covariance through those same weights (see
[CN0 Estimator](cn0_estimator.md)):

```jldoctest power_user
julia> using Tracking, GNSSSignals, StaticArrays

julia> using Tracking: Hz, NumAnts, AbstractPostCorrFilter

julia> # Trivial beamformer — averages across antenna elements
       struct MyBeamformer <: AbstractPostCorrFilter end

julia> Tracking.update(f::MyBeamformer, prompt) = f;

julia> Tracking.get_weights(::MyBeamformer, ::NumAnts{1}) = 1.0 + 0.0im;

julia> Tracking.get_weights(::MyBeamformer, ::NumAnts{M}) where {M} =
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

### Multi-signal discriminator combining

Tracking several signals of one satellite and closing the loops on only one of them throws away most of the information. A signal group can therefore combine them: with `discriminator_combining = true`, a passenger signal's (`signals[2:end]`) discriminator outputs are folded into the driver's before the loop filters see them, as a minimum-variance weighted mean.

It is declared on the group, next to its band and antenna count, because the assumption it rests on — `signals[1]` integrates at least as long as every passenger — is a property of the group's signal tuple, which every satellite of the group shares:

```julia
# Every coincident record drives the loops, for every satellite of the group.
TrackState(; signals = (galileo_e1 = (GalileoE1C(), GalileoE1B()),), discriminator_combining = true)

# Or on a hand-built group, which then keeps its own setting.
SignalGroup((GalileoE1C(), GalileoE1B()); discriminator_combining = true)

# The default: driver-only.
TrackState(; signals = (galileo_e1 = (GalileoE1C(), GalileoE1B()),))
```

Combining is opt-in because it changes a multi-signal satellite's carrier/code Doppler, code phase and decoded-bit timing. With it off — and on every single-signal satellite either way — the loops are exactly the driver's own, bit for bit.

**Which records combine: coincident ones only.** For each driver record, a passenger record contributes to that loop update if and only if it ended on the **same sample** — its `sample_index` equals the driver record's (see [`CorrelatorOutput`](@ref)). Every other passenger record is still applied to its own signal in full — prompt, post-correlation filter, C/N₀ estimator, bit buffer, `last_fully_integrated_*` — but its discriminators reach no loop. Nothing is carried from one driver record, or one processing chunk, to the next: the combination for a driver record is formed from that record and the passenger records that ended with it, and is spent on that one update.

That makes the integration lengths the thing to get right. For full aiding, driver and passengers integrate equally long, so every passenger record ends with a driver record. That is the default: every signal starts at one primary code block (except Galileo E5a-QP, which is tracked alone), and every intra-band pilot/data pair GNSSSignals defines — GPS L5, GPS L1C, Galileo E1, E5a, E5b and E6, BeiDou B1C and B2a — shares one primary code period between its components. What is lost where the lengths differ:

  - A passenger integrating **shorter** than the driver — `k` of its records per driver record, like a 1 ms L1 C/A passenger beside a 10 ms L1C-P driver — contributes one record in `k`: the one that ends with the driver's. The other `k − 1` are applied and dropped.
  - A caller lengthening the **driver alone** with [`set_preferred_num_code_blocks_to_integrate!`](@ref) puts its passengers in exactly that position, until they are lengthened to match.
  - A passenger integrating **longer** than the driver breaks the assumption; it aids only the driver records it happens to end with.

None of this is enforced: a mis-ordered or mismatched group is built and tracked like any other, and the setter accepts any length its own signal allows. What a mismatch costs is gain, never a wrong measurement, and the receiver is the one holding the reason a group is configured as it is — so it is a documented assumption, which a receiver such as GNSSReceiver.jl checks before it turns combining on.

**How records are weighted.** Each record contributes with weight `P · N / G`: the component's **nominal** power share `P` (`GNSSSignals.get_relative_power` — the ICD power split: L1C-P `0.75` against L1C-D `0.25`, E1B/E1C `0.5` each) times the record's sample count `N`, over the discriminator's noise gain `G`. For the code loop `G` is [`dll_disc_noise_gain`](@ref Tracking.dll_disc_noise_gain); for the carrier phase loop it is the same for every Costas discriminator and cancels, leaving `P · N`; for the carrier frequency loop the weight is `P · N³`, because `fll_disc` divides an inter-prompt phase difference by `2π·T`, which scales its noise by `1/T`. `P · N` is the post-integration SNR up to factors common to the satellite's signals, and every such factor — front-end noise density, sampling frequency, code amplitude — cancels and is never formed.

The power is nominal rather than measured because within one constellation and band the ICD ratio is the *stable* quantity: elevation, satellite block, free-space loss, antenna gain and front-end gain are common to the components and cancel out of it, whereas a measured `|P|²` carries the estimator's own noise (biased upward by `σ²/N`) into the loop gain. It also costs nothing — `get_relative_power` folds to a compile-time constant — and needs no C/N₀ estimate, so it works for a signal tracked with [`NoCN0Estimator`](@ref).

The noise gain is what lets a BOC pilot's very-early-minus-late discriminator outvote a legacy BPSK early-late one at equal C/N₀, as it should. Mis-weighting costs combining efficiency but never introduces bias, since each discriminator is individually calibrated in chips. There is no generic gain, though: a custom correlator type with its own `dll_disc` method must define `dll_disc_noise_gain` as well before a satellite using it can combine.

The gains are evaluated at the origin of the S-curve, so "never introduces bias" holds while every contributor is inside its own linear range. A sharper discriminator has a *narrower* one — roughly ±0.15 chips for the default VEML taps on a BOC(1,1) peak against ±0.5 for a 1-chip early-late layout — and it is also the one carrying ~8× the weight, so a large starting error is the one condition where the weighting works against you. It does not arise for the pilot/data pairs GNSSSignals defines, which are either all-BOC (Galileo E1, GPS L1C, BeiDou B1C) or all-BPSK (GPS L5, Galileo E5a, BeiDou B2a). For a mixed group, put the BOC signal first or hand off from acquisition inside a fraction of a chip.

**A mean, not a sum.** The combined discriminator per loop is `Σ wᵢ dᵢ / Σ wᵢ` over the driver's record and the coincident passenger records, so the loop gain is the same however many signals reported and only the measurement noise changes. Where no passenger contributed to a loop — nothing coincided, or every coincident record carried zero weight for that loop — the loop reads the driver's own discriminator unchanged, bit for bit. Zero total weight reads as zero.

Three records are held out, and each only from the loops the reason applies to:

  - **A pre-sync passenger record on a secondary-coded signal.** On a signal with a secondary code — every pilot component, plus the Neumann-Hoffman-coded data components (GPS L5I, Galileo E5a-I and E5b-I, BeiDou B2a-I, B1I and B3I) — the record still correlated against the pre-sync replica in the chunk where that signal's secondary-code sync was detected has a coherent sum that partially cancels and measures nothing, so it is excluded from combining altogether: one record per signal per sync event. The driver's own record is never excluded — it is the record the update is formed on.
  - **A record with no previous prompt** — a signal's first record, or its first after [`reset_loop_filters!`](@ref) — has no frequency measurement, and contributes zero weight to the FLL.
  - **A passenger without a known group delay**, or any passenger while `signals[1]`'s is unknown, aids the two carrier loops only — see [Group delay](#Group-delay). The driver's own code discriminator is never withheld.

**Declare only signals the satellite transmits.** Nominal weights say what a component's power share *should* be; they cannot tell you whether it is being received. That is by design — a signal is in a satellite's tuple only because you put it there — but it does mean a component the satellite does not broadcast contributes noise at its full nominal weight. A pre-Block-III GPS satellite carries no L1C at all; route such satellites to a C/A-only group instead.

**Why discriminators and not prompts.** `pll_disc` and `fll_disc` are blind to the unknown ±1 navigation-data sign, so a data component's discriminator can be averaged with a pilot's; summing their *prompts* would let a bit flip cancel the pilot. A passenger transmitted in carrier quadrature with the driver (GPS L5I against L5Q, Galileo E5a-I against E5a-Q) has its correlator de-rotated onto the driver's phase frame before `pll_disc` reads it; `dll_disc` and `fll_disc` are rotation-invariant and read it as it is.

That is the reason for the two *carrier* loops. The code loop's is separate, and it holds even where the data sign is known: `dll_disc` is noncoherent — a function of the tap magnitudes, not of a prompt — and every method calibrates itself against its own signal's correlation shape, a triangular BPSK autocorrelation for the early-late form and a sine-BOC(1,1) envelope for the VEML one. Taps summed across a BOC signal and a BPSK one produce an envelope with no single slope to divide by, so the result would not be calibrated in chips at all — and being calibrated in chips is exactly what makes two signals' code errors averageable.

The setting is independent of the estimator: [`VectorPLLAndDLL`](@ref) reads the same group flag, with one restriction — a satellite already in the vector loop (`vt_on = true`) combines the **carrier phase** discriminator only, because its other two loops are the navigation filter's and each signal's discriminator reaches the filter as its own raw measurement for it to weigh (see [Vector tracking](vector_tracking.md)). In its scalar fallback, which is where a satellite pulls in, all three combine exactly as described here.

The design, and what the coincidence rule gives up against accumulating every passenger record inside a driver record's window, is recorded in `docs/plans/2026-09-30-simple-discriminator-combining.md`.

### Group delay

The carrier loops combine from the first integration: every component rides the same carrier, so once de-rotated each signal's Costas discriminator measures the same phase error the driver's does. The **code** loop needs one number from you first.

A satellite's signals leave it through different payload paths, so their code phases differ by the differences between their **group delays** — nanoseconds, i.e. decimetres to metres of range. A combined DLL left uncorrected drives the satellite-shared `code_phase` towards a weighted average of the signals' code phases, at an offset that *moves as the weights move*. The offset is not negligible where combining helps most: an L1C-P passenger can outweigh an L1 C/A driver by roughly 80:1, so a 1–3 ns difference lands as 0.3–0.9 m of bias — more than the jitter the combining just bought. So each [`TrackedSignal`](@ref) carries its `group_delay`, the signal's own payload delay, larger for the signal that leaves the satellite through the longer path and so arrives later and sits at the *smaller* code phase, and each passenger's code discriminator is referred to the driver's code phase with it at the loop update — `(delay(signals[1]) − delay(i)) · f_code` chips are subtracted, at the code frequency in effect. The satellite's `code_phase` then keeps meaning the driver's. That happens wherever this package closes the combined code loop itself — under [`VectorPLLAndDLL`](@ref) that is its scalar fallback only, since a vector-closed satellite's code loop belongs to the navigation filter and the per-signal accumulators that filter reads ([`mean_code_discr`](@ref)) are raw, for it to refer as it sees fit (see [Vector tracking](vector_tracking.md)).

The value is a **time** and carries its unit, like every other dimensioned quantity here; a bare number is refused rather than assumed to be seconds, and a value in metres gets a sentence naming the conversion, since a realistic value is sub-nanosecond and an assumed unit costs metres. Set it with [`set_group_delay!`](@ref), or pass `group_delay` to the `TrackedSignal` constructor:

```julia
# Every slot, `signals[1]` included — the datum is yours to choose.
set_group_delay!(track_state, :modern_gps, 11, 1, 0.0s)
set_group_delay!(track_state, :modern_gps, 11, GPSL1CA, -1.0e-9s)   # or -1.0u"ns"
```

**Only differences between a satellite's signals are ever used.** What is applied to passenger `i` is `delay(signals[1]) − delay(i)` — positive where the passenger is *less* delayed than the driver, and so sits at the larger code phase — which means the datum your values are stated against cancels and is yours to choose. Two usages follow, and they are the same rule:

  - **Driver as datum.** Put `0.0s` on `signals[1]` and state each passenger's delay relative to it; a passenger's stored value is then directly its bias.
  - **Your own datum.** Give every signal its own payload delay on whatever reference you hold, `signals[1]` included — derived as **Where the value comes from** below describes, which for GPS means `−ISC_x` rather than the ISC as it comes. Any term common to the satellite's signals (GPS's `T_GD`) cancels in the difference, so nothing has to be referred on the way in.

[`get_group_delay`](@ref) reads back what you stored — this signal's own delay, not the difference; the difference is formed inside the fold and never has to be built by a caller. A **single-signal** satellite has no difference to form, so a value stored there is kept and never used.

**Every slot starts at `nothing`, `signals[1]` included** — the datum is a statement about the satellite that only you can make, and no slot is treated differently from any other. So the code loop combines once you have supplied a value for the driver *and* for each passenger you want in it: `nothing` on a passenger withholds that one signal's code contribution, `nothing` on `signals[1]` withholds every passenger's, since nothing can be referred to an unknown datum. The carrier loops combine throughout. `nothing` and `0.0s` are deliberately different states. `0.0s` on every slot asserts that the components share a code phase, which is true by construction for a pilot/data pair out of one payload chain and is yours to state. `nothing` says the difference is unknown, and nothing is assumed on your behalf, because assuming zero is the unsafe direction: fused code measurements would then drive the shared `code_phase` to a weighted average of the signals' phases that wanders as the weights move, and a downstream consumer such as `PositionVelocityTime.jl` — which applies the group-delay correction of whichever ranging signal you name — would apply it to a phase that is no longer that signal's. An unknown difference is also indistinguishable from a value you meant to supply and forgot.

**Where the value comes from is deliberately not this package's business.** Tracking neither knows which constellation broadcasts what nor parses navigation messages. Since only differences within a satellite are read, any term the signals share — a per-SV `T_GD`, say — cancels and need not be resolved. A value may be:

  - **The same number on every signal, usually `0.0s`**, justified by the ICD. Every Galileo pair — E1B/E1C, E5aI/E5aQ, E5bI/E5bQ, E6B/E6C — leaves the satellite as one composite modulation out of one payload chain, and the broadcast group delays are cross-band only (E1-to-E5a/E5b), so there is no intra-band difference to state.
  - **A broadcast inter-signal correction, per signal.** GPS broadcasts a *per-component* ISC (`ISC_L1CA`, `ISC_L1CD`, `ISC_L1CP`, `ISC_L2C`, `ISC_L5I5`, `ISC_L5Q5`) — precisely a statement that the components are not assumed to share a group delay — and BeiDou does the same for its B1C and B2a pairs (`ISC_B1Cd`, `ISC_B2ad`). Mind the sign, which differs between the two families because they reference their corrections differently. IS-GPS-705/800 has the receiver correct a range on signal `x` by `−T_GD + ISC_x`, so that signal's *delay* is `T_GD − ISC_x`: write **`−ISC_x`** on each signal and let `T_GD` cancel. BeiDou states the pilot's delay as `T_GD` and the data component's as `T_GD + ISC`, so write **`0.0s` on the pilot and `+ISC` on the data component**.
  - **A ground calibration**, for a signal pair whose ICD broadcasts nothing useful.

Which ISCs a decoder can give you depends on the *message*, not the signal it came from: CNAV-2 (from an L1C-D decoder) carries all six GPS terms; CNAV (L5I / L2C, message type 30) carries everything except the L1C pair; LNAV carries none, so a legacy-only receiver has nothing to derive an L1 C/A + L1C difference from.

```@docs
set_group_delay!
get_group_delay
Tracking.dll_disc_noise_gain
```

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

By default the track function uses the last antenna channel as the reference signal to drive the discriminators. An appropriate beamforming algorithm will probably suit better — construct a [`TrackedSat`](@ref) with a custom `post_corr_filter` and build the `TrackState` from it (so the slot type takes the custom filter type rather than the default). A filter supplies its combining weights through [`get_weights`](@ref); Tracking applies them to the correlator and reduces the measured noise covariance through the same weights, so the C/N₀ stays correct for whatever the beamformer does:

```jldoctest beamformer_array
julia> using Tracking, GNSSSignals, StaticArrays

julia> using Tracking: Hz, NumAnts, AbstractPostCorrFilter

julia> # Same trivial mean-of-antennas filter as the power-user example above
       struct MyBeamformer <: AbstractPostCorrFilter end

julia> Tracking.update(f::MyBeamformer, prompt) = f;

julia> Tracking.get_weights(::MyBeamformer, ::NumAnts{1}) = 1.0 + 0.0im;

julia> Tracking.get_weights(::MyBeamformer, ::NumAnts{M}) where {M} =
           SVector{M,ComplexF64}(ntuple(_ -> 1 / M + 0.0im, M));

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
| `get_last_fully_integrated_num_code_blocks` | Primary-code block count of the last completed integration. |
| `get_last_fully_integrated_integration_time` | Integration time `T` of the last completed integration. Post-integration SNR is `estimate_cn0` × this; C/N₀ alone is processing-independent. |
| `get_post_corr_filter` | The post-correlation filter. |
| `get_cn0_estimator` | The CN0 estimator. |
| `estimate_cn0` | CN0 estimate in dB-Hz. |
| `get_bit_buffer` / `get_soft_bits` / `get_num_bits` | Bit buffer, decoded soft bits (sign = the hard bit), bit count. |
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

julia> get_num_bits(track_state, :modern_gps, 11, 1)  # per-signal by index (estimator-driver signal)
0

julia> get_num_bits(track_state, :modern_gps, 11, GPSL1CA)  # per-signal by type
0
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
get_last_fully_integrated_integration_time(::TrackedSignal)
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
