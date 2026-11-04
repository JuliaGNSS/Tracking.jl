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

The default tracking loop is `ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)`, a
scalar PLL/DLL with
- `ThirdOrderAssistedBilinearLF` for the carrier loop (FLL-assisted PLL)
- `SecondOrderBilinearLF` for the code loop

When [`TrackState`](@ref) builds the default tracking loop from a signal-tuple
declaration, the carrier loop is **staged** (see
[`CarrierLoopStage`](@extref TrackingLoops.CarrierLoopStage)): it pulls in as an
FLL-assisted PLL at a **wide** 50 Hz with a 5 Hz **FLL-assist** bandwidth, drops the
FLL once a phase-lock indicator reads lock, and later narrows to a **narrow** 18 Hz.
The **code** bandwidth is a flat 1 Hz. All are one-sided noise bandwidths `BL`
(Kaplan & Hegarty), for every signal. The stages only advance, and
[`reset_loop_filters!`](@ref) keeps the stage reached; a satellite added anew starts
over. Read a satellite's stage with
`get_carrier_loop_stage(get_loop_state(sat))`.

Where the replica wipes every sign modulation off the prompt (a dataless signal
synced to its secondary code) the carrier discriminators turn four-quadrant.

The defaults are not scaled with the integration time; each is capped by its
stability product `BL · T_int` at filter time: the PLL by `0.09`, the FLL path by `0.02`, the
code loop by `0.018`
([`effective_carrier_loop_filter_bandwidth`](@extref TrackingLoops.effective_carrier_loop_filter_bandwidth-Tuple{Any, Any})
and siblings). The bandwidths that end up in the loop are:

| Integration time                   |  Wide BL | Narrow BL | FLL-assist BL |  Code BL |
|------------------------------------|---------:|----------:|--------------:|---------:|
| 1 ms (L1 C/A, L5I, …)              |    50 Hz |     18 Hz |          5 Hz |     1 Hz |
| 4 ms (E1B / E1C)                   |  22.5 Hz |     18 Hz |          5 Hz |     1 Hz |
| 10 ms (L1C-D/P, B1C)               |     9 Hz |      9 Hz |          2 Hz |     1 Hz |
| 20 ms (L1 C/A or L2 CM after sync) |   4.5 Hz |    4.5 Hz |          1 Hz |   0.9 Hz |
| 1.5 s (L2 CL)                      |  0.06 Hz |   0.06 Hz |      0.013 Hz | 0.012 Hz |

The [TrackingLoops manual](https://juliagnss.github.io/TrackingLoops.jl/stable/)
("Loop-filter bandwidths") gives the rationale for the bandwidths and caps.

Before TrackingLoops 4 the carrier filter was fed radians but its output read as Hz,
making the loop about 5.6× wider than configured (≈ 100 Hz for 18 Hz); to get
approximately that loop back, set both `wide_carrier_loop_filter_bandwidth` and
`narrow_carrier_loop_filter_bandwidth` to `85.0Hz`. Before TrackingLoops 5 the loop
ran FLL-assisted at one carrier bandwidth throughout; pass equal wide and narrow
bandwidths for one PLL bandwidth.

Override per signal with methods of
[`default_wide_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_wide_carrier_loop_filter_bandwidth-Tuple{AbstractGNSSSignal}),
[`default_narrow_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_narrow_carrier_loop_filter_bandwidth-Tuple{AbstractGNSSSignal}),
[`default_fll_assist_loop_filter_bandwidth`](@extref TrackingLoops.default_fll_assist_loop_filter_bandwidth-Tuple{AbstractGNSSSignal}) /
[`default_code_loop_filter_bandwidth`](@extref TrackingLoops.default_code_loop_filter_bandwidth-Tuple{AbstractGNSSSignal}), or per `TrackState` by passing your own `tracking_loop =`.

## Tracking Loops

In the TrackingLoops manual:

- [`ScalarTrackingLoop`](@extref TrackingLoops.ScalarTrackingLoop)
- [`DelayCompensatedTrackingLoop`](@extref TrackingLoops.DelayCompensatedTrackingLoop)
- [`default_wide_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_wide_carrier_loop_filter_bandwidth-Tuple{AbstractGNSSSignal})
- [`default_narrow_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_narrow_carrier_loop_filter_bandwidth-Tuple{AbstractGNSSSignal})
- [`default_fll_assist_loop_filter_bandwidth`](@extref TrackingLoops.default_fll_assist_loop_filter_bandwidth-Tuple{AbstractGNSSSignal})
- [`default_code_loop_filter_bandwidth`](@extref TrackingLoops.default_code_loop_filter_bandwidth-Tuple{AbstractGNSSSignal})
- [`effective_carrier_loop_filter_bandwidth`](@extref TrackingLoops.effective_carrier_loop_filter_bandwidth-Tuple{Any, Any})
- [`effective_fll_assist_loop_filter_bandwidth`](@extref TrackingLoops.effective_fll_assist_loop_filter_bandwidth-Tuple{Any, Any})
- [`effective_code_loop_filter_bandwidth`](@extref TrackingLoops.effective_code_loop_filter_bandwidth-Tuple{Any, Any})
- [`LoopBandwidths`](@extref TrackingLoops.LoopBandwidths)
- [`CarrierLoopStage`](@extref TrackingLoops.CarrierLoopStage)
- [`get_carrier_loop_stage`](@extref TrackingLoops.get_carrier_loop_stage-Tuple{TrackingLoops.ScalarLoopState})
- [`get_phase_lock_indicator`](@extref TrackingLoops.get_phase_lock_indicator-Tuple{TrackingLoops.ScalarLoopState})

## Resetting loop filters

After changing a signal's integration length with
[`set_preferred_num_code_periods_to_integrate!`](@ref), reset the affected loop
filters so the old integration length's filter state does not leak into the new one.

```@docs
reset_loop_filters!
```

## Custom Configuration

You can customize the loop filters and bandwidths when creating the tracking loop:

```jldoctest custom_loop
julia> using Tracking, TrackingLoops, GNSSSignals, TrackingLoopFilters

julia> using Tracking: Hz

julia> # Use non-assisted PLL with custom loop filter types
       tracking_loop = ScalarTrackingLoop(
           ThirdOrderBilinearLF,      # carrier loop filter type
           SecondOrderBilinearLF;     # code loop filter type
           wide_carrier_loop_filter_bandwidth = 30.0Hz,
           narrow_carrier_loop_filter_bandwidth = 15.0Hz,
           code_loop_filter_bandwidth = 0.5Hz
       );

julia> track_state = TrackState(; signal = GPSL1CA(), tracking_loop);

julia> track_state = add_satellite!(track_state; prn = 1, code_phase = 50.0, carrier_doppler = 1000.0Hz);
```

## Custom Loop Filters

You can implement a custom loop filter `MyLoopFilter <: AbstractLoopFilter`. In this
case, a specialized `filter_loop` function is needed. For more information
refer to [TrackingLoopFilters.jl](https://github.com/JuliaGNSS/TrackingLoopFilters.jl).

## Custom Tracking Loop

To replace the loop-filter-based tracking loop with a different algorithm
(Kalman filter, joint-channel loop, …), see the dedicated guide
in [`custom_tracking_loop.md`](custom_tracking_loop.md).
