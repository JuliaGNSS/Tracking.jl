# Loop Filter

The loop filters are provided by
[TrackingLoopFilters.jl](https://github.com/JuliaGNSS/TrackingLoopFilters.jl). This includes:
- first order loop filter `FirstOrderLF`
- second order bilinear loop filter `SecondOrderBilinearLF`
- second order boxcar loop filter `SecondOrderBoxcarLF`
- third order bilinear loop filter `ThirdOrderBilinearLF`
- third order boxcar loop filter `ThirdOrderBoxcarLF`
- third order assisted bilinear loop filter `ThirdOrderAssistedBilinearLF` (combines PLL and FLL)

## Default Configuration

The default Doppler estimator is `ConventionalAssistedPLLAndDLL` which uses:
- `ThirdOrderAssistedBilinearLF` for the carrier loop (FLL-assisted PLL for improved dynamics)
- `SecondOrderBilinearLF` for the code loop

When [`TrackState`](@ref) builds the default estimator implicitly from a
signal-tuple declaration, the **carrier** bandwidth is a flat 18 Hz and the
**code** bandwidth a flat 1 Hz for every signal.

Both are one-sided noise bandwidths `BL` in the sense of Kaplan & Hegarty, so
they plug into the usual PLL jitter and dynamic-stress formulas. (Before
[#244](https://github.com/JuliaGNSS/Tracking.jl/issues/244) the carrier filter
was fed radians, so a configured 18 Hz loop behaved like one of ~100 Hz.)

18 Hz is the third-order PLL bandwidth of the literature (Kaplan & Hegarty
Table 5.6, Pany Table 3.3; Borre et al. quote about 20 Hz). None of them scales
it with the primary code period: thermal jitter and dynamic stress, which set
the bandwidth, do not depend on it.

### Stability cap

Stability does depend on the loop update interval `Δt`, so at filter time each
bandwidth is capped against the record's actual integration time
(`TrackingLoops.effective_carrier_loop_filter_bandwidth`,
[`effective_code_loop_filter_bandwidth`](@extref TrackingLoops.effective_code_loop_filter_bandwidth)). An explicit bandwidth is
capped the same way; the cap only ever narrows a loop.

```julia
BL_carrier = min(BL_carrier_configured, 0.09  / Δt)
BL_code    = min(BL_code_configured,    0.018 / Δt)
```

The default FLL-assisted third-order carrier filter diverges at `BL · Δt ≈ 0.4`
(the plain third-order one at ≈ 0.43). Below that, the loop's actual noise
bandwidth runs wider than configured:

| `BL · Δt`                    | 0.018 | 0.036 | 0.072 | 0.09  | 0.18  | 0.36 |
|------------------------------|------:|------:|------:|------:|------:|-----:|
| `ThirdOrderAssistedBilinearLF` | 1.05× | 1.09× | 1.19× | 1.25× | 1.65× | 5.2× |
| `ThirdOrderBilinearLF`         | 0.97× | 1.00× | 1.07× | 1.12× | 1.38× | 2.8× |

The carrier cap at 0.09 keeps the loop within 25 % of its configured bandwidth
with about 4× stability margin. It is the product of Kaplan & Hegarty's
third-order design example (18 Hz at 5 ms) and matches GNSS-SDR's narrow
post-sync bandwidths (5 Hz at 20 ms). Kaplan & Hegarty (Fig. 5.24) and Pany
also run 18 Hz at 20 ms, but warn that analog-derived loops reach their design
bandwidth only for `BL · Δt` well below unity — this filter is 5× wider there.
The resulting defaults:

| Integration | Signals                                                   | Carrier BL |  Code BL |
|-------------|-----------------------------------------------------------|-----------:|---------:|
| 1 ms        | GPS L1 C/A, GPS L5, Galileo E5a, …                        |      18 Hz |     1 Hz |
| 2 ms        | Galileo E5a-QP                                            |      18 Hz |     1 Hz |
| 4 ms        | Galileo E1B / E1C                                         |      18 Hz |     1 Hz |
| 10 ms       | GPS L1C-D / L1C-P, BeiDou B1C; GPS L5I synced at 10 ms    |       9 Hz |     1 Hz |
| 20 ms       | GPS L2 CM; GPS L1 C/A, L5Q, Galileo E5a-I synced at 20 ms |     4.5 Hz |   0.9 Hz |
| 1.5 s       | GPS L2 CL                                                 |    0.06 Hz | 0.012 Hz |

The code cap of 0.018 is conservative for the second-order code filter, which
destabilizes only around `BL · Δt ≈ 0.4` (S. A. Stephens and J. B. Thomas,
"Controlled-Root Formulation for Digital Phase-Locked Loops", IEEE Trans.
Aerospace and Electronic Systems 31(1), 1995). Carrier-aided, the DLL has
almost no dynamics to track, so its 1 Hz is a thermal-noise-versus-pull-in
choice, applied in full up to 18 ms.

Override per signal by defining methods of
[`default_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_carrier_loop_filter_bandwidth) /
[`default_code_loop_filter_bandwidth`](@extref TrackingLoops.default_code_loop_filter_bandwidth), or override at construction
time by passing your own `doppler_estimator =` to `TrackState`.

## Carrier loop staging

Every satellite's carrier loop starts as an FLL-assisted PLL and drops the FLL
once its frequency has converged, along Kaplan & Hegarty's closure sequence
(§5.3, §5.5): "apply the error inputs from both discriminators as an
FLL-assisted PLL until phase lock is achieved, then convert to pure PLL". There
is no down-staging: the pure PLL stays until [`reset_loop_filters!`](@ref)
restarts the staging. The staging is TrackingLoops' (its
`FrequencyLockIndicator`, held in the per-satellite estimator state); what
Tracking adds is the per-record information it needs.

  - **Frequency lock** is declared once the mean FLL discriminator, the
    residual frequency error, stays below `frequency_lock_threshold` (3 Hz, at
    most 1/(16T) on long records, a quarter of the two-quadrant FLL's range)
    over a `frequency_lock_window` (0.5 s, at least four records). The window
    mean is the phase advance across the window over its length, so its noise
    falls with the window length rather than with each record's SNR. Both are
    overridable per signal type by defining a method of the TrackingLoops
    function. Kaplan & Hegarty's phase lock indicator (§5.11.2), specified
    unnormalised for 20 ms updates, never declares lock at 1 ms and 30 dB-Hz,
    which would leave the noisy FLL branch in.
  - **The pure PLL** is the FLL-assisted filter with a zero FLL input, which is
    exactly the third-order PLL `ThirdOrderBilinearLF` (same state, same
    coefficients), so the switch costs nothing and leaves no transient.
  - **Four-quadrant discriminators** apply where the replica wipes every sign
    modulation off the prompt: a dataless signal
    (`get_data_frequency(signal) == 0`), synced to its secondary code where it
    has one. Tracking hands each record this state of the driver's bit buffer
    from before the fold (TrackingLoops' `is_wiped_off` and `sync_polarity`).
      * The FLL is four-quadrant (twice the pull-in range) from the first
        record correlated after that sync, and from the start for the pilots
        without a secondary code (GPS L2 CL, Galileo E5a-QP). It needs no sign,
        only that consecutive prompts share it.
      * The PLL is four-quadrant (linear over ±180°, worth up to 6 dB of
        threshold) from the same record on, reading the prompt with the sign
        the secondary-code sync found. A short overlay (GPS L5Q's 20 ms) can
        sync while the loop is still pulling in, and a Costas slip after it
        makes the switch a half-cycle jump of the carrier phase: the start of
        the resolved phase ([`get_carrier_phase_polarity`](@ref)), not part of
        a continuous one.
      * The pilots without a secondary code (GPS L2 CL, Galileo E5a-QP) keep
        the Costas PLL. That is a choice, not a necessity: their prompt keeps
        its sign too, so the PLL could turn four-quadrant with the sign the
        Costas loop holds. But without a sync that sign is arbitrary, so the
        switch would leave the carrier phase unresolved, and it would have to
        be taken off a single noisy prompt; the wider range alone was not
        considered worth that.

    Data signals stay on the two-quadrant (Costas) discriminators.

The FLL compares a record's prompt with the previous record's, so Tracking
leaves a record out of it where that comparison does not hold: the first record
after [`add_satellite!`](@ref) or [`reset_loop_filters!`](@ref); a record whose
length differs from the previous one's, as at a data signal's bit sync or after
[`set_preferred_num_code_blocks_to_integrate!`](@ref), since the reading divides
the rotation by the record's own integration time; and a pilot's first record
after the sync that wipes its prompt off, which may differ in sign from the one
before.

With a carrier filter other than the FLL-assisted one the loop is a PLL from the
start and runs no frequency lock indicator.

## Signal combining

A satellite tracked on several signals of one band (e.g. Galileo E1C and E1B)
closes its loops on the driver, `signals[1]`. With `combine_signals = true` on
the estimator, e.g. `ConventionalAssistedPLLAndDLL(; combine_signals = true)`,
the discriminators of the other signals, the passengers, are combined with the
driver's into a weighted mean before the loop filters read it. A mean rather
than a sum, so the loop gain does not change with the number of signals. The
combining is TrackingLoops'. Tracking steps every passenger record with the
estimator's `step_loop`, in sample order alongside the driver's, as it does for
every estimator; an estimator that does not combine leaves the state as it is.

  - **Weights** are each signal's ICD power share
    (`GNSSSignals.get_relative_power`) times the record's integration time, and
    for the FLL times the integration time squared on top, as its noise variance
    falls with its cube. All signals of a satellite share one antenna and one
    path, so the power split fixes their C/N₀ ratio; no C/N₀ estimate is read.
    Discriminators are calibrated (PLL in cycles, FLL in Hz, DLL in chips), so
    the mean is unbiased whatever the weights; they only decide how much noise
    is removed.
  - **Time alignment:** every passenger record ending within a driver
    integration is combined into that driver record. Records after the driver's
    last one in a chunk stay pending, carried in the satellite's state into the
    next chunk; a sync phase snap drops them with the driver's in-flight
    integration. Where no passenger record ended, the driver's loops close on its
    own discriminators. Each passenger is assumed to integrate no longer than
    the driver, as a longer record would dominate the one driver update it is
    combined into: make the longest-integrating signal (typically the pilot) the
    driver.
  - **PLL and FLL:** passengers read the two-quadrant (Costas) discriminators,
    which are blind to a sign flip of a whole record, so data passengers and
    records correlated before a passenger's own sync are combined like any
    other. They are combined into a carrier loop while it is formed, the FLL
    until frequency lock. Into a four-quadrant driver discriminator, after the
    sync, they are combined only while its reading lies within their own
    two-quadrant range, ±1/4 cycle for the PLL and ±1/(4T) for the FLL, where
    both read the error alike. In simulation (Galileo E1C driving, E1B combined)
    this lowers the carrier-phase jitter after the sync by about 30 % from 20 to
    40 dB-Hz, and less at 18 dB-Hz, without adding cycle slips. A passenger's
    prompt is rotated onto the driver's phase frame by the signals' nominal
    carrier phase offsets (`get_carrier_phase_offset`); a residual phase bias
    between the components is not modelled and shifts the combined lock point.
  - **DLL:** a passenger is combined into the code loop only where both its own
    and the driver's group delay are known ([`set_group_delay!`](@ref)); zero is
    not assumed. Its discriminator is referred to the driver's code phase by the
    difference of the two, in chips at the driver's code frequency.
  - **Vector tracking:** out of the vector loop, `VectorPLLAndDLL` combines as
    its inner loop does. In the vector loop passengers are combined into the PLL
    only: the code loop and the FLL branch belong to the navigation filter.
  - `NCOReferencedPLLAndDLL` does not combine: it steps the driver's phase error
    predicted to the landing sample, which the passengers' records are not.

## Doppler Estimators

In the TrackingLoops manual:

- [`ConventionalPLLAndDLL`](@extref TrackingLoops.ConventionalPLLAndDLL)
- [`ConventionalAssistedPLLAndDLL`](@extref TrackingLoops.ConventionalAssistedPLLAndDLL)
- [`default_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_carrier_loop_filter_bandwidth)
- [`default_code_loop_filter_bandwidth`](@extref TrackingLoops.default_code_loop_filter_bandwidth)
- `effective_carrier_loop_filter_bandwidth` (TrackingLoops 4)
- [`effective_code_loop_filter_bandwidth`](@extref TrackingLoops.effective_code_loop_filter_bandwidth)

## Resetting loop filters

When you change a signal's coherent-integration length mid-track with
[`set_preferred_num_code_blocks_to_integrate!`](@ref), reset the affected
loop filters for a clean handoff so the previous integration length's filter
state does not leak into the new one.

```@docs
reset_loop_filters!
```

## Custom Configuration

You can customize the loop filters and bandwidths when creating the Doppler estimator:

```jldoctest custom_loop
julia> using Tracking, TrackingLoops, GNSSSignals, TrackingLoopFilters

julia> using Tracking: Hz

julia> # Use non-assisted PLL with custom loop filter types
       doppler_estimator = ConventionalPLLAndDLL(
           ThirdOrderBilinearLF,      # carrier loop filter type
           SecondOrderBilinearLF;     # code loop filter type
           carrier_loop_filter_bandwidth = 15.0Hz,
           code_loop_filter_bandwidth = 0.5Hz
       );

julia> track_state = TrackState(; signal = GPSL1CA(), doppler_estimator);

julia> track_state = add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);
```

## Custom Loop Filters

You can implement a custom loop filter `MyLoopFilter <: AbstractLoopFilter`. In this
case, a specialized `filter_loop` function is needed. For more information
refer to [TrackingLoopFilters.jl](https://github.com/JuliaGNSS/TrackingLoopFilters.jl).

## Custom Doppler Estimator

To replace the loop-filter-based estimator with a different algorithm
(Kalman filter, joint-channel estimator, …), see the dedicated guide
in [`custom_doppler_estimator.md`](custom_doppler_estimator.md).
