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
**code** bandwidth a flat 1 Hz for every signal. Both are one-sided noise
bandwidths `BL` in the sense of Kaplan & Hegarty, so they plug into the usual
PLL jitter and dynamic-stress formulas: the carrier filter is fed the phase
error in cycles and the FLL error in Hz. 18 Hz is the third-order PLL bandwidth
of the literature (Kaplan & Hegarty, *Understanding GPS*, Table 5.6).

Neither default is scaled with the signal's code period or the number of
integrated blocks. Each is instead capped at filter time by its stability
product `BL · T_int` against the record's actual integration time: the carrier
by `0.09` (`TrackingLoops.MAX_CARRIER_LOOP_BANDWIDTH_TIME_PRODUCT`, Kaplan's
18 Hz at a 5 ms update), the code by `0.018`
(`TrackingLoops.MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT`). The bandwidths that end
up in the loop are:

| Integration time            | Carrier BL | Code BL  |
|-----------------------------|-----------:|---------:|
| ≤ 5 ms (L1 C/A, L5I, E1B)   |    18 Hz   |    1 Hz  |
| 10 ms (L1C-D/P, B1C)        |     9 Hz   |    1 Hz  |
| 20 ms (L1 C/A or L2 CM after sync) | 4.5 Hz |  0.9 Hz |
| 1.5 s (L2 CL)               |   0.06 Hz  | 0.012 Hz |

The two caps differ because the loops do. The carrier loop carries the
dynamic stress, so it should stay as wide as stability allows. The DLL is
carrier-aided and has almost no dynamic stress to track, so its bandwidth is a
thermal-noise-versus-pull-in choice, and its tighter cap only binds past 18 ms.

Both caps keep well clear of instability: transform-designed digital loops of
this kind only destabilize around `BL · Δt ≈ 0.4` (S. A. Stephens and
J. B. Thomas, "Controlled-Root Formulation for Digital Phase-Locked Loops",
IEEE Trans. Aerospace and Electronic Systems 31(1), 1995). As `BL · Δt` grows,
the true noise bandwidth runs wider than configured: about 5 % at `0.018` and
25 % at `0.09` for the FLL-assisted carrier filter.

Before TrackingLoops 4, the carrier filter was fed the phase error in radians
while its output was read as Hz. The loop gain was 2π too high and the loop
about 5.6× wider than configured (≈ 100 Hz for the 18 Hz L1 C/A default). To
get approximately that loop back, configure
`carrier_loop_filter_bandwidth = 85.0Hz`.

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
- [`effective_carrier_loop_filter_bandwidth`](@extref TrackingLoops.effective_carrier_loop_filter_bandwidth)
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
