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
signal-tuple declaration, the carrier loop is **staged** on phase lock: it pulls
in as an FLL-assisted PLL at a **wide** 50 Hz, its FLL path at a 5 Hz
**FLL-assist** bandwidth, drops the FLL once a phase-lock indicator has read
lock, and later narrows the PLL to a **narrow** 18 Hz. The **code** bandwidth is
a flat 1 Hz. All of them hold for every signal and are one-sided noise
bandwidths `BL` in the sense of Kaplan & Hegarty, so they plug into the usual
jitter and dynamic-stress formulas: the carrier filter is fed the phase error in
cycles and the FLL error in Hz. 18 Hz is the third-order PLL bandwidth of the
literature (Kaplan & Hegarty, *Understanding GPS*, Table 5.6). The stages only
advance, and [`reset_loop_filters!`](@ref) keeps the stage a loop has reached:
it hands over to a new integration length, say, rather than re-acquiring. A
satellite added anew starts over. Read a satellite's stage
with `carrier_loop_stage(get_doppler_estimator_state(sat))` (see "Carrier-loop
staging" in the
[TrackingLoops manual](https://juliagnss.github.io/TrackingLoops.jl/stable/)).

Where the replica wipes every sign modulation off the prompt — a dataless
signal synced to its secondary code — the carrier discriminators turn
four-quadrant, over twice the range of the two-quadrant (Costas) ones.

None of the defaults is scaled with the signal's code period or the number of
integrated blocks. Each is instead capped at filter time by its stability
product `BL · T_int` against the record's actual integration time: the PLL, wide
and narrow alike, by `0.09` (`TrackingLoops.MAX_CARRIER_LOOP_BANDWIDTH_TIME_PRODUCT`,
Kaplan's 18 Hz at a 5 ms update), the FLL path by `0.02`
(`TrackingLoops.MAX_FLL_ASSIST_LOOP_BANDWIDTH_TIME_PRODUCT`), the code by `0.018`
(`TrackingLoops.MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT`). The bandwidths that end
up in the loop are:

| Integration time                   |  Wide BL | Narrow BL | FLL-assist BL |  Code BL |
|------------------------------------|---------:|----------:|--------------:|---------:|
| 1 ms (L1 C/A, L5I, …)              |    50 Hz |     18 Hz |          5 Hz |     1 Hz |
| 4 ms (E1B / E1C)                   |  22.5 Hz |     18 Hz |          5 Hz |     1 Hz |
| 10 ms (L1C-D/P, B1C)               |     9 Hz |      9 Hz |          2 Hz |     1 Hz |
| 20 ms (L1 C/A or L2 CM after sync) |   4.5 Hz |    4.5 Hz |          1 Hz |   0.9 Hz |
| 1.5 s (L2 CL)                      |  0.06 Hz |   0.06 Hz |      0.013 Hz | 0.012 Hz |

The caps differ because the loops do. The carrier loop carries the
dynamic stress, so it should stay as wide as stability allows. The DLL is
carrier-aided and has almost no dynamic stress to track, so its bandwidth is a
thermal-noise-versus-pull-in choice, and its tighter cap only binds past 18 ms.

The caps keep well clear of instability: transform-designed digital loops of
this kind only destabilize around `BL · Δt ≈ 0.4` (S. A. Stephens and
J. B. Thomas, "Controlled-Root Formulation for Digital Phase-Locked Loops",
IEEE Trans. Aerospace and Electronic Systems 31(1), 1995). As `BL · Δt` grows,
the true noise bandwidth runs wider than configured: about 5 % at `0.018` and
25 % at `0.09` for the FLL-assisted carrier filter.

Before TrackingLoops 4, the carrier filter was fed the phase error in radians
while its output was read as Hz. The loop gain was 2π too high and the loop
about 5.6× wider than configured (≈ 100 Hz for the 18 Hz L1 C/A default). To
get approximately that loop back, configure
`wide_carrier_loop_filter_bandwidth = 85.0Hz` and
`narrow_carrier_loop_filter_bandwidth = 85.0Hz`. Before TrackingLoops 5, the
loop ran FLL-assisted at one carrier bandwidth throughout; pass the same wide
and narrow bandwidth for one PLL bandwidth.

Override per signal by defining methods of
[`default_wide_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_wide_carrier_loop_filter_bandwidth),
[`default_narrow_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_narrow_carrier_loop_filter_bandwidth),
[`default_fll_assist_loop_filter_bandwidth`](@extref TrackingLoops.default_fll_assist_loop_filter_bandwidth) /
[`default_code_loop_filter_bandwidth`](@extref TrackingLoops.default_code_loop_filter_bandwidth), or override at construction
time by passing your own `doppler_estimator =` to `TrackState`.

## Doppler Estimators

In the TrackingLoops manual:

- [`ConventionalPLLAndDLL`](@extref TrackingLoops.ConventionalPLLAndDLL)
- [`ConventionalAssistedPLLAndDLL`](@extref TrackingLoops.ConventionalAssistedPLLAndDLL)
- [`default_wide_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_wide_carrier_loop_filter_bandwidth)
- [`default_narrow_carrier_loop_filter_bandwidth`](@extref TrackingLoops.default_narrow_carrier_loop_filter_bandwidth)
- [`default_fll_assist_loop_filter_bandwidth`](@extref TrackingLoops.default_fll_assist_loop_filter_bandwidth)
- [`default_code_loop_filter_bandwidth`](@extref TrackingLoops.default_code_loop_filter_bandwidth)
- [`effective_carrier_loop_filter_bandwidth`](@extref TrackingLoops.effective_carrier_loop_filter_bandwidth)
- [`effective_fll_assist_loop_filter_bandwidth`](@extref TrackingLoops.effective_fll_assist_loop_filter_bandwidth)
- [`effective_code_loop_filter_bandwidth`](@extref TrackingLoops.effective_code_loop_filter_bandwidth)
- [`LoopBandwidths`](@extref TrackingLoops.LoopBandwidths)
- [`CarrierLoopStage`](@extref TrackingLoops.CarrierLoopStage)
- [`carrier_loop_stage`](@extref TrackingLoops.carrier_loop_stage)
- [`phase_lock_indicator`](@extref TrackingLoops.phase_lock_indicator)

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
           wide_carrier_loop_filter_bandwidth = 30.0Hz,
           narrow_carrier_loop_filter_bandwidth = 15.0Hz,
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
