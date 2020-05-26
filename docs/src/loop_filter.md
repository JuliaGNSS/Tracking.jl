# Loop Filter

The default loop filter are provided by
[TrackingLoopFilters.jl](https://github.com/JuliaGNSS/TrackingLoopFilters.jl). This includes
- first order loop filter `FirstOrderLF`
- second order biliner loop filter `SecondOrderBilinearLF`
- second order boxcar loop filter `SecondOrderBoxcarLF`
- third order bilinear loop filter `ThirdOrderBilinarLF`
- third order boxcar loop filter `ThirdOrderBoxcarLF`

The default loop filter for the carrier loop is `ThirdOrderBilinarLF` and the default loop
filter for the code loop is `SecondOrderBilinearLF`. This is set by the initialization of
the `TrackingState` and can be changed by:
```julia
state = TrackingState(
  GPSL1,
  carrier_doppler,
  code_phase,
  carrier_loop_filter = ThirdOrderBilinearLF(),
  code_loop_filter = SecondOrderBilinearLF()
)
```

The bandwidth of the loop filter is set by the `track` function:
```julia
next_results = track(
  next_signal,
  get_state(results),
  prn,
  sampling_frequency,
  carrier_loop_filter_bandwidth = 18Hz,
  code_loop_filter_bandwidth = 1Hz
)
```

You can easily implement a custom loop filter `MyLoopFilter <: AbstractLoopFilter`. In this
case, a specialized `propagate` and `get_filter_output` is needed. For more information
refer to [TrackingLoopFilters.jl](https://github.com/JuliaGNSS/TrackingLoopFilters.jl).
