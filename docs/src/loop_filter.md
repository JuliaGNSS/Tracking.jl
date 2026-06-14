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
signal-tuple declaration, the loop bandwidths are sized **per signal** at
`BL · T ≈ 0.018`, where `T` is the signal's primary code period — that's
~10× margin from the `BL · T < 0.18` stability edge of the bilinear
third-order filter. The values fall out to:

| Signal      | Primary period | Carrier BL | Code BL  |
|-------------|----------------|-----------:|---------:|
| GPS L1 C/A  | 1 ms           |    18 Hz   |   1 Hz   |
| GPS L5I     | 1 ms           |    18 Hz   |   1 Hz   |
| Galileo E1B | 4 ms           |   4.5 Hz   |  0.25 Hz |
| GPS L1C-D   | 10 ms          |   1.8 Hz   |  0.1 Hz  |
| GPS L1C-P   | 10 ms          |   1.8 Hz   |  0.1 Hz  |

The 1-ms-primary-period signals (L1 C/A, L5I) keep the historical 18 Hz /
1 Hz default; longer-period signals get appropriately tighter loops so the
PLL stays stable. Override per signal by defining methods of
[`default_carrier_loop_filter_bandwidth`](@ref) /
[`default_code_loop_filter_bandwidth`](@ref), or override at construction
time by passing your own `doppler_estimator =` to `TrackState`.

## Doppler Estimators

```@docs
ConventionalPLLAndDLL
ConventionalAssistedPLLAndDLL
default_carrier_loop_filter_bandwidth
default_code_loop_filter_bandwidth
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
