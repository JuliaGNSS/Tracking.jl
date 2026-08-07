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
signal-tuple declaration, the **carrier** bandwidth is sized **per signal**
from the signal's primary code period `T`, at `BL · T ≈ 0.018` — ~10× margin
from the `BL · T < 0.18` stability edge of the bilinear third-order filter.
The **code** bandwidth is a flat 1 Hz for every signal. The values fall out to:

| Signal      | Primary period | Carrier BL | Code BL |
|-------------|----------------|-----------:|--------:|
| GPS L1 C/A  | 1 ms           |    18 Hz   |    1 Hz |
| GPS L5I     | 1 ms           |    18 Hz   |    1 Hz |
| Galileo E1B | 4 ms           |   4.5 Hz   |    1 Hz |
| GPS L1C-D   | 10 ms          |   1.8 Hz   |    1 Hz |
| GPS L1C-P   | 10 ms          |   1.8 Hz   |    1 Hz |
| GPS L2 CM   | 20 ms          |   0.9 Hz   |    1 Hz |
| GPS L2 CL   | 1.5 s          | 0.012 Hz   |    1 Hz |

The 1-ms-primary-period signals (L1 C/A, L5I) keep the historical 18 Hz /
1 Hz default; longer-period signals get appropriately tighter carrier loops
so the PLL stays stable. The DLL does **not** follow the carrier loop down:
being carrier-aided it has almost no dynamic stress to track, so its
bandwidth is a thermal-noise-versus-pull-in choice that does not scale with
the symbol rate.

The two are also treated differently at filter time, since only the carrier
bandwidth is a per-code-period reference. Integrating `N` primary blocks
coherently scales the carrier bandwidth to `BL/N`, holding its `BL · Δt`
product at the single-period value, while the code bandwidth is left alone and
merely capped by the same product against the record's actual integration time
(`Tracking.effective_code_loop_filter_bandwidth`). That cap is what pulls the
L2 C primaries down in practice — 0.9 Hz for a 20 ms L2 CM integration,
0.012 Hz for a 1.5 s L2 CL one — while every integration shorter than 18 ms,
whatever its signal or block count, runs the DLL at the full 1 Hz.

Reusing the third-order carrier filter's `0.018` product to cap the
*second*-order code filter is conservative: transform-designed digital loops
of this kind only destabilize around `BL · Δt ≈ 0.4` (S. A. Stephens and
J. B. Thomas, "Controlled-Root Formulation for Digital Phase-Locked Loops",
IEEE Trans. Aerospace and Electronic Systems 31(1), 1995 — the standard
treatment of digital-loop stability at large `BL · Δt`), so the code loop's
cap carries even more stability margin than the carrier loop's.

Override per signal by defining methods of
[`default_carrier_loop_filter_bandwidth`](@ref) /
[`default_code_loop_filter_bandwidth`](@ref), or override at construction
time by passing your own `doppler_estimator =` to `TrackState`.

## Doppler Estimators

```@docs
ConventionalPLLAndDLL
ConventionalAssistedPLLAndDLL
default_carrier_loop_filter_bandwidth
default_code_loop_filter_bandwidth
Tracking.effective_code_loop_filter_bandwidth
```

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
julia> using Tracking, GNSSSignals, TrackingLoopFilters

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
