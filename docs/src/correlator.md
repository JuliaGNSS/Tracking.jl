# Correlator

A correlator correlates the incoming signal with the local code replica at one
code-phase offset; a correlator bank holds several (early, prompt, late, …). Each
[`TrackedSignal`](@ref) holds its own correlator bank, so a multi-signal satellite
carries one per signal.

## Why early/late spacing matters

The PLL discriminator and C/N₀ estimator consume the prompt; the early and late
accumulators give the DLL the slope of the code-correlation triangle. Early minus
late crosses zero when the replica is aligned, and its sign tells the DLL which way
to push.

Notation: `d` is the early-to-late spacing in chips (`d = 0.5`: E and L ½ chip
apart), `ε` the code-phase error in chips, `B_fe` the front-end (analog + ADC)
bandwidth, `T_chip = 1 / f_code` the chip period, `f_code` the chip rate (1.023 MHz
for GPS L1 C/A) and `Fs` the sampling frequency.

- **Thermal-noise jitter** decreases as `d` shrinks, because the noise on E and L is
  correlated and partially cancels in their difference (Van Dierendonck, Fenton &
  Ford, "Theory and Performance of Narrow Correlator Spacing in a GPS Receiver",
  NAVIGATION 39(3), 1992). A wider spacing only smooths a time-multiplexed
  (dither) DLL; Tracking.jl computes E, P and L simultaneously from the same samples.
- **Multipath error** also decreases as `d` shrinks, as the discriminator samples
  the correlation function nearer the peak (same paper).
- **Robustness to dynamics decreases as `d` shrinks.** The discriminator is linear
  only over `|ε| < d/2` (the pull-in range `T_DLL`, Springer Handbook of GNSS,
  §14.4.4); beyond it the slope drops to zero with the sign still correct. The lock
  rule of thumb `3σ_DLL + ε_D ≤ T_DLL ≈ d/2` (eq. 14.57) leaves a narrow `d` less
  budget for dynamic stress error `ε_D`, which is why NovAtel's narrow-correlator
  design widens `d` during acquisition and narrows it once locked.
- **Front-end bandwidth** caps the gains: below roughly `d ≈ 1 / (B_fe · T_chip)`
  the band-limited peak is rounded, the slope flattens and jitter grows again. Van
  Dierendonck reports `d ≈ 0.1` as optimal at 8 MHz on L1 C/A.

### The sample rate sets the minimum `d`

`preferred_early_late_to_prompt_spacing` is the early-to-prompt (and
prompt-to-late) shift, so `d` is twice it. It is rounded to the nearest whole sample
shift, at least one, so the realized spacing is
`d = 2 · max(1, round(shift · Fs / f_code)) · f_code / Fs` chips. A `d ≈ 0.1` on
L1 C/A needs `Fs ≥ 20 · f_code ≈ 20.46 MHz`; at 4 MHz `d` saturates at
`2 · 1.023e6 / 4e6 ≈ 0.51` chips whatever the preferred value.

The `VeryEarlyPromptLateCorrelatorBank` adds a very-early / very-late pair for the
multi-peaked correlation function of BOC signals like Galileo E1B: the inner pair
locks the central peak, the outer pair disambiguates it from the side peaks.

## Default Correlators

[`get_default_correlator_bank`](@extref TrackingLoops.get_default_correlator_bank) returns, per
signal type:

- [`VeryEarlyPromptLateCorrelatorBank`](@extref TrackingLoops.VeryEarlyPromptLateCorrelatorBank)
  for the BOC-family signals GPS L1C-D, GPS L1C-P, Galileo E1B, Galileo E1C and
  BeiDou B1C;
- [`EarlyPromptLateCorrelatorBank`](@extref TrackingLoops.EarlyPromptLateCorrelatorBank) for
  every other signal.

## Correlator Types

In the TrackingLoops manual:

- [`EarlyPromptLateCorrelatorBank`](@extref TrackingLoops.EarlyPromptLateCorrelatorBank)
- [`VeryEarlyPromptLateCorrelatorBank`](@extref TrackingLoops.VeryEarlyPromptLateCorrelatorBank)

## Accessing Correlator Values

In the TrackingLoops manual:

- [`get_early`](@extref TrackingLoops.get_early-Tuple{TrackingLoops.AbstractEarlyPromptLateCorrelatorBank})
- [`get_prompt`](@extref TrackingLoops.get_prompt-Tuple{TrackingLoops.AbstractCorrelatorBank})
- [`get_late`](@extref TrackingLoops.get_late-Tuple{TrackingLoops.AbstractEarlyPromptLateCorrelatorBank})
- [`get_correlator_outputs`](@extref TrackingLoops.get_correlator_outputs-Tuple{TrackingLoops.AbstractCorrelatorBank})
- [`get_num_correlators`](@extref TrackingLoops.get_num_correlators-Tuple{TrackingLoops.AbstractCorrelatorBank})
- [`get_num_antennas(::AbstractCorrelatorBank)`](@extref TrackingLoops.get_num_antennas-Union{Tuple{TrackingLoops.AbstractCorrelatorBank{M}}, Tuple{M}} where M)

## Correlator outputs

Each completed integration within a processing chunk is recorded as a
`CorrelatorBankOutput`, collected per signal and consumed by the tracking loop
after the chunk (see [Chunked Doppler updates](track.md#Chunked-Doppler-updates)).
An external correlator producer (e.g. an FPGA) can build these itself and feed
them straight to the tracking loop — see
[External correlator producers](track.md#External-correlator-producers).

```@docs
get_correlator_bank_outputs
append_correlator_bank_output!
```

In the TrackingLoops manual:

- [`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput)

## Sample shifts

`calc_correlator_sample_offsets` places each correlator relative to the prompt, in
samples; `calc_early_late_sample_spacing`, the spacing of the outermost pair, scales
the DLL discriminator back to chips.

In the TrackingLoops manual:

- [`calc_correlator_sample_offsets`](@extref TrackingLoops.calc_correlator_sample_offsets-Tuple{EarlyPromptLateCorrelatorBank, Any, Any})
- [`calc_early_late_sample_spacing`](@extref TrackingLoops.calc_early_late_sample_spacing-Tuple{TrackingLoops.AbstractEarlyPromptLateCorrelatorBank, Any, Any})
- [`with_correlator_outputs(::EarlyPromptLateCorrelatorBank, ::Any)`](@extref TrackingLoops.with_correlator_outputs-Tuple{EarlyPromptLateCorrelatorBank, Any})
- [`with_correlator_outputs(::VeryEarlyPromptLateCorrelatorBank, ::Any)`](@extref TrackingLoops.with_correlator_outputs-Tuple{VeryEarlyPromptLateCorrelatorBank, Any})

## Antenna and correlator counts

In the TrackingLoops manual:

- [`NumAntennas`](@extref TrackingLoops.NumAntennas)
- [`NumCorrelators`](@extref TrackingLoops.NumCorrelators)

## Beamformer

A beamformer reduces the per-antenna correlator outputs to the single channel the
Doppler and C/N₀ estimators see. It declares its combining weights through
[`get_beamforming_weights`](@extref TrackingLoops.get_beamforming_weights) and Tracking applies them, so the
measured noise covariance is reduced through the *same* weights and the C/N₀'s
numerator and denominator describe one channel (see [CN0 Estimator](cn0_estimator.md)).

[`DefaultBeamformer`](@extref TrackingLoops.DefaultBeamformer) is the identity
for a single antenna and selects the last element for an array. For an array a
beamformer is the natural override; see [Tracking State](tracking_state.md).

In the TrackingLoops manual:

- [`AbstractBeamformer`](@extref TrackingLoops.AbstractBeamformer)
- [`DefaultBeamformer`](@extref TrackingLoops.DefaultBeamformer)
- [`update!(::DefaultBeamformer, ::Any)`](@extref TrackingLoops.update!-Tuple{DefaultBeamformer, Any})
- [`get_beamforming_weights`](@extref TrackingLoops.get_beamforming_weights)

## Adjusting the spacing of the stock correlator

Changing `d` needs no custom correlator: pass `d / 2` in chips as
`preferred_early_late_to_prompt_spacing` (here `d = 0.2`):

```jldoctest narrow_corr
julia> using Tracking, TrackingLoops, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        correlator_bank = EarlyPromptLateCorrelatorBank(
                            preferred_early_late_to_prompt_spacing = 0.1,
                        ));

julia> track_state = add_satellite!(track_state, :default, sat);

julia> get_correlator_bank(track_state, :default, 1, 1).preferred_early_late_to_prompt_spacing
0.1
```

`VeryEarlyPromptLateCorrelatorBank` exposes the inner
(`preferred_early_late_to_prompt_spacing`, default `0.15`) and the outer
(`preferred_very_early_late_to_prompt_spacing`, default `0.6`) shift. Either is
subject to the sample-rate floor above.

### Changing the spacing of a satellite that's already tracking

The sat, signal and correlator are immutable, so build a replacement sat with
`Accessors.@set` and write it into the slot with `Dictionaries.set!`:

```jldoctest narrow_corr
julia> using Accessors, Dictionaries

julia> sat     = track_state.groups[:default].satellites[1];

julia> new_sat = @set sat.signals[1].correlator_bank.preferred_early_late_to_prompt_spacing = 0.7;

julia> Dictionaries.set!(track_state.groups[:default].satellites, 1, new_sat);

julia> get_correlator_bank(track_state, :default, 1, 1).preferred_early_late_to_prompt_spacing
0.7
```

`@set` keeps everything else, including the in-flight correlator outputs, and the
`TrackState`'s concrete type does not change.

## Custom Correlator Banks

Write a correlator bank for a *layout* the stock types lack: extra correlators, asymmetric
shifts, or a different prompt index. Subtype `AbstractCorrelatorBank` (or
`AbstractEarlyPromptLateCorrelatorBank` to inherit the early/late accessors) and
implement:

| Method | Purpose |
|--------|---------|
| `get_correlator_outputs(c)` | Return the backing storage of the correlator outputs (e.g. an `SVector` of `ComplexF64`). |
| `get_num_correlators(c)` | Return how many correlators the bank holds. Defaults to `size(get_correlator_outputs(c), 1)` — only override if your layout differs. |
| `with_correlator_outputs(c, correlator_outputs)` | Construct a fresh correlator bank with the given correlator outputs (immutable update). |
| `calc_correlator_sample_offsets(c, sampling_frequency, code_frequency)` | Return an `SVector{N,Int}` of sample offsets relative to the prompt, ordered from latest replica to earliest. |

`get_prompt_index` defaults to the middle correlator; override it for an
asymmetric layout.

### Skeleton

A bank of seven correlators, e.g. for a multipath-mitigation discriminator that fits a
slope over several early/late correlators:

```jldoctest multitap
julia> using Tracking, TrackingLoops, GNSSSignals, StaticArrays

julia> using TrackingLoops:
           AbstractEarlyPromptLateCorrelatorBank,
           zero_correlator_outputs,
           calc_correlator_sample_offset,
           NumAntennas,
           NumCorrelators

julia> struct SevenCorrelatorBank{M,T} <: AbstractEarlyPromptLateCorrelatorBank{M}
           correlator_outputs::SVector{7,T}
           correlator_spacing_chips::Float64
           # Pin M from the correlator output element type, as TrackingLoops' EPL does.
           function SevenCorrelatorBank(
               correlator_outputs::AbstractVector{SVector{M,T}}, correlator_spacing_chips,
           ) where {M,T<:Complex}
               new{M,SVector{M,T}}(correlator_outputs, correlator_spacing_chips)
           end
           function SevenCorrelatorBank(
               correlator_outputs::AbstractVector{T}, correlator_spacing_chips,
           ) where {T<:Complex}
               new{1,T}(correlator_outputs, correlator_spacing_chips)
           end
       end

julia> # Keyword constructor, like the EPL one
       function SevenCorrelatorBank(; num_antennas::NumAntennas = NumAntennas(1),
                                     correlator_spacing_chips = 0.1)
           SevenCorrelatorBank(
               zero_correlator_outputs(num_antennas, NumCorrelators(7)),
               correlator_spacing_chips,
           )
       end;

julia> TrackingLoops.with_correlator_outputs(c::SevenCorrelatorBank, correlator_outputs) =
           SevenCorrelatorBank(correlator_outputs, c.correlator_spacing_chips);

julia> function TrackingLoops.calc_correlator_sample_offsets(
           c::SevenCorrelatorBank, sampling_frequency, code_frequency,
       )
           s = calc_correlator_sample_offset(
               c.correlator_spacing_chips, sampling_frequency, code_frequency,
           )
           SVector(-3s, -2s, -s, 0, s, 2s, 3s)
       end;
```

Build the `TrackState` from a sat carrying it, so the slot type is the custom
correlator bank's:

```jldoctest multitap
julia> using Tracking: Hz

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        correlator_bank = SevenCorrelatorBank(correlator_spacing_chips = 0.1));

julia> track_state = TrackState(GPSL1CA(), sat);

julia> TrackingLoops.get_num_correlators(get_correlator_bank(track_state, 1))
7
```

The default DLL discriminator reads only the outermost early/late pair, so pair a
multi-correlator bank with a [Custom Tracking Loop](custom_tracking_loop.md).
TrackingLoops' `src/correlators/early_prompt_late.jl` and
`very_early_prompt_late.jl` show the full pattern, including multi-antenna
correlator outputs.
