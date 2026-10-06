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
