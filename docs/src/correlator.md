# Correlator

The correlator computes the correlation between the incoming signal and the
locally generated replica code at multiple code-phase offsets (early,
prompt, late, …). Each [`TrackedSignal`](@ref) on a [`TrackedSat`](@ref)
holds its own correlator, so a multi-signal satellite carries one
correlator per signal.

## Why early/late spacing matters

The prompt accumulator is what the PLL discriminator and CN0 estimator
consume; the surrounding early/late accumulators exist purely to give the
DLL a slope estimate of the code-correlation triangle. The early-minus-late
amplitude crosses zero exactly when the replica is aligned with the
incoming code, and its sign on either side tells the DLL which direction
to push.

Notation used below:

- `d` — early-to-late spacing, expressed as a fraction of one PRN chip
  (so `d = 0.5` means E and L sit ½ chip apart).
- `ε` — code-phase error: the offset between the locally generated
  replica and the incoming code, in chips. `ε = 0` is perfect lock.
- `E`, `P`, `L` — early, prompt, and late correlator outputs (complex
  accumulators).
- `B_fe` — pre-correlation (front-end) bandwidth in Hz: the analog +
  ADC bandwidth that limits how sharp the correlation peak can be.
- `T_chip` — one PRN chip period in seconds (`1 / f_code`); equivalently
  `B_fe · T_chip` is the front-end bandwidth measured in chip-rate
  units.
- `Fs` — sampling frequency in Hz.
- `f_code` — chip rate in Hz (1.023 MHz for GPS L1 C/A).

With those names settled, the spacing `d` determines several DLL
properties together:

- **Thermal-noise jitter** decreases as `d` shrinks, because the noise on
  the early and late accumulators is correlated and partially cancels in
  the early-minus-late difference. This is the central result of Van
  Dierendonck, Fenton & Ford ("Theory and Performance of Narrow
  Correlator Spacing in a GPS Receiver", NAVIGATION 39(3), 1992) and is
  contrary to the intuition that wider should give a smoother
  discriminator — it only does so when E and L are time-multiplexed
  (dither-DLL), which Tracking.jl does **not** do; here E, P, L are
  computed simultaneously from the same sample stream.
- **Multipath error** also decreases as `d` shrinks: a narrower
  discriminator samples the correlation function nearer the peak, where
  multipath distortion is smallest. Same source, same paper.
- **Robustness to dynamics decreases as `d` shrinks.** The
  discriminator is only linear in `ε` over `|ε| < d/2` (the Springer
  Handbook of GNSS, §14.4.4, calls this the pull-in range `T_DLL`);
  outside it the slope drops to zero with the sign still correct, so
  the loop converges but no longer behaves as designed. The lock
  rule of thumb (Springer eq. 14.57) is `3σ_DLL + ε_D ≤ T_DLL ≈ d/2`,
  so a narrow `d` shrinks the budget for dynamic stress error `ε_D`
  and the receiver is more prone to losing lock under fast LOS
  acceleration or scintillation. This is why NovAtel's original
  narrow-correlator design widens `d` during acquisition and only
  narrows once in steady-state lock.
- **Front-end bandwidth interaction.** The benefits above saturate (and
  eventually reverse) once `d` drops below roughly `1 / (B_fe · T_chip)`:
  a band-limited correlation peak is rounded, so the early-minus-late
  slope flattens and noise jitter starts to grow again. Van Dierendonck
  reports `d ≈ 0.1` is optimal at 8 MHz front-end bandwidth on L1 C/A.

So "narrower is better, up to a point that the front-end bandwidth
allows" — there is no fundamental noise penalty until the band-limit
kicks in.

### Tracking.jl-specific constraint: sample rate sets the minimum `d`

Tracking.jl expresses the spacing as a
`preferred_early_late_to_prompt_code_shift` fraction-of-chip; the actual
sample shift is the nearest integer at the current sampling rate, with a
**floor of one sample**. The realized `d` is therefore
`max(1, round(d_pref · Fs / f_code)) · f_code / Fs` chips. To realize the
recommended narrow `d ≈ 0.1` on L1 C/A (1.023 MHz chip rate) the
sampling rate has to be at least `10 · f_code ≈ 10.23 MHz`; at 4 MHz the
nearest integer shift saturates the spacing at `1 / (4e6 / 1.023e6) ≈
0.26` chips regardless of what `preferred_early_late_to_prompt_code_shift`
is set to. So narrow-correlator operation in Tracking.jl is bottlenecked
by both the front-end bandwidth (the physical argument above) and the
ADC/sampling rate (this implementation detail).

The `VeryEarlyPromptLateCorrelator` adds a second outer pair (very-early /
very-late) to handle the multi-peaked code-correlation function of
BOC-modulated signals like Galileo E1B — the inner pair locks the central
peak and the outer pair disambiguates from the side-peaks.

## Default Correlators

The default correlator depends on the GNSS signal type and is returned by
[`get_default_correlator`](@ref):

- [`EarlyPromptLateCorrelator`](@ref) for GPS L1 C/A, GPS L1C-D, GPS L1C-P,
  GPS L2CM, GPS L2CL, GPS L5I, GPS L5Q, Galileo E5a-I, Galileo E5a-Q
  (BPSK / TMBOC modulation)
- [`VeryEarlyPromptLateCorrelator`](@ref) for Galileo E1B and Galileo E1C
  (CBOC / BOC modulation)

```@docs
get_default_correlator
```

## Correlator Types

```@docs
EarlyPromptLateCorrelator
VeryEarlyPromptLateCorrelator
```

## Accessing Correlator Values

```@docs
get_early
get_prompt
get_late
get_accumulators
get_num_accumulators
get_num_ants(::AbstractCorrelator)
```

## Sample shifts

The downconvert/correlate inner loop walks the incoming sample buffer once
and asks the correlator where to lay down each accumulator relative to the
prompt — that mapping is `get_correlator_sample_shifts`. The spacing in
*samples* between the outermost shifts is `get_early_late_sample_spacing`,
which is what the DLL uses to scale its discriminator output back to chips.

```@docs
get_correlator_sample_shifts
get_early_late_sample_spacing
update_accumulator(::EarlyPromptLateCorrelator, ::Any)
update_accumulator(::VeryEarlyPromptLateCorrelator, ::Any)
```

## Antenna and accumulator counts

```@docs
NumAnts
NumAccumulators
```

## Post-correlation filter

A post-correlation filter is applied to the correlator output before the
Doppler estimator sees it. The default is a passthrough; for multi-antenna
tracking a beamformer is the natural override.

```@docs
AbstractPostCorrFilter
DefaultPostCorrFilter
```

## Integration sizing

```@docs
calc_signal_samples_to_integrate
```

## Adjusting the spacing of the stock correlator

The stock `EarlyPromptLateCorrelator` already exposes its spacing as a
constructor knob — you do **not** need a custom correlator type just to
change `d`. Pass the desired chip fraction via
`preferred_early_late_to_prompt_code_shift` and hand it to a `TrackedSat`:

```jldoctest narrow_corr
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        correlator = EarlyPromptLateCorrelator(
                            preferred_early_late_to_prompt_code_shift = 0.1,
                        ));

julia> track_state = add_satellite!(track_state, :default, sat);

julia> get_correlator(track_state, :default, 1, 1).preferred_early_late_to_prompt_code_shift
0.1
```

Likewise `VeryEarlyPromptLateCorrelator` exposes both the inner
(`preferred_early_late_to_prompt_code_shift`, default `0.15`) and the
outer (`preferred_very_early_late_to_prompt_code_shift`, default `0.6`)
spacing.

Remember the sample-rate floor: the actual realized `d` is the nearest
integer sample shift at the current sampling rate (minimum 1), so on
low-`Fs` configurations the *preferred* value may be rounded to a coarser
realized value (see the Tracking.jl-specific constraint above).

### Changing the spacing of a satellite that's already tracking

`TrackedSat`, `TrackedSignal`, and `EarlyPromptLateCorrelator` are all
immutable structs, so changing the spacing means producing a replacement
sat with the new value and writing it into the satellites-dictionary slot.
The cleanest way is `Accessors.@set` for the deeply-nested field update,
followed by `Dictionaries.set!` for the dictionary write:

```jldoctest narrow_corr
julia> using Accessors, Dictionaries

julia> sat     = track_state.groups[:default].satellites[1];

julia> new_sat = @set sat.signals[1].correlator.preferred_early_late_to_prompt_code_shift = 0.7;

julia> Dictionaries.set!(track_state.groups[:default].satellites, 1, new_sat);

julia> get_correlator(track_state, :default, 1, 1).preferred_early_late_to_prompt_code_shift
0.7
```

`@set` preserves everything else on the correlator (including the
in-flight accumulator), only swapping the one field. The dictionary slot
is reassigned in place, so the `TrackState`'s concrete type doesn't
change and the hot path stays type-stable across the swap.

## Custom Correlators

Roll your own correlator when you need a *layout* the stock types don't
provide — extra taps, asymmetric shifts, a non-standard
`get_correlator_sample_shifts` mapping, or a different `get_prompt_index`
arrangement. Subtype `AbstractCorrelator` (or
`AbstractEarlyPromptLateCorrelator` to inherit the early/late accessors)
and implement a small interface.

### What to implement

| Method | Purpose |
|--------|---------|
| `get_accumulators(c)` | Return the backing storage of accumulators (e.g. an `SVector` of `ComplexF64`). |
| `get_num_accumulators(c)` | Return how many accumulators the correlator holds. Defaults to `size(get_accumulators(c), 1)` — only override if your layout differs. |
| `update_accumulator(c, accumulators)` | Construct a fresh correlator instance with the given accumulators (immutable update). |
| `get_correlator_sample_shifts(c, sampling_frequency, code_frequency)` | Return an `SVector{N,Int}` of sample offsets relative to the prompt, ordered from latest replica to earliest. |

The `get_prompt_index` default places the prompt at the middle of the
accumulator vector — if you want a different layout (e.g. asymmetric
shifts) override that too.

### Skeleton

A 7-tap correlator with three early and three late accumulators (e.g. as
a building block for a custom multipath-mitigation discriminator that
fits a slope over multiple early/late samples rather than using a single
E−L difference):

```jldoctest multitap
julia> using Tracking, GNSSSignals, StaticArrays

julia> using Tracking: AbstractEarlyPromptLateCorrelator,
                       get_initial_accumulator,
                       calc_preferred_code_shift_to_sample_shift,
                       NumAnts, NumAccumulators

julia> struct MultiTapCorrelator{M,T} <: AbstractEarlyPromptLateCorrelator{M}
           accumulators::SVector{7,T}
           tap_spacing_chips::Float64
           # Dispatch on accumulator-vector element type to pin M, mirroring
           # the pattern in `src/correlators/early_prompt_late.jl`.
           function MultiTapCorrelator(
               accumulators::AbstractVector{SVector{M,T}}, tap_spacing_chips,
           ) where {M,T<:Complex}
               new{M,SVector{M,T}}(accumulators, tap_spacing_chips)
           end
           function MultiTapCorrelator(
               accumulators::AbstractVector{T}, tap_spacing_chips,
           ) where {T<:Complex}
               new{1,T}(accumulators, tap_spacing_chips)
           end
       end

julia> # Convenience kwarg-only constructor — matches the EPL surface
       function MultiTapCorrelator(; num_ants::NumAnts = NumAnts(1),
                                     tap_spacing_chips = 0.1)
           MultiTapCorrelator(
               get_initial_accumulator(num_ants, NumAccumulators(7)),
               tap_spacing_chips,
           )
       end;

julia> Tracking.update_accumulator(c::MultiTapCorrelator, accumulators) =
           MultiTapCorrelator(accumulators, c.tap_spacing_chips);

julia> function Tracking.get_correlator_sample_shifts(
           c::MultiTapCorrelator, sampling_frequency, code_frequency,
       )
           s = calc_preferred_code_shift_to_sample_shift(
               c.tap_spacing_chips, sampling_frequency, code_frequency,
           )
           SVector(-3s, -2s, -s, 0, s, 2s, 3s)
       end;
```

To use it, plug it into a `TrackedSat` and build the `TrackState` from
that sat (so the slot type is taken from the custom correlator rather
than the default `EarlyPromptLateCorrelator`):

```jldoctest multitap
julia> using Tracking: Hz

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        correlator = MultiTapCorrelator(tap_spacing_chips = 0.1));

julia> track_state = TrackState(GPSL1CA(), sat);

julia> get_num_accumulators(get_correlator(track_state, 1))
7
```

The default DLL discriminator only reads the outermost early/late pair,
so a multi-tap correlator like this is most useful in combination with a
custom Doppler estimator (see [Custom Doppler Estimator](custom_doppler_estimator.md))
that walks every accumulator. The reference implementations in
`src/correlators/early_prompt_late.jl` and
`src/correlators/very_early_prompt_late.jl` show the full pattern
including the multi-antenna `SVector` element type.
