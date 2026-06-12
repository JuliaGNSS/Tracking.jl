# CN0 Estimator

The CN0 (Carrier-to-Noise density ratio) estimator provides a measure of
signal quality in dB-Hz. Each [`TrackedSignal`](@ref) on a
[`TrackedSat`](@ref) holds its own CN0 estimator, so a multi-signal
satellite produces one CN0 value per signal.

## Default Estimator

The default CN0 estimator is the
[Moments Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1),
implemented as `MomentsCN0Estimator`. It buffers prompt correlation values
and estimates CN0 from the first and second moments of their magnitudes.
The buffer holds `num_prompts_for_cn0_estimation` prompts (default `100`,
~100 ms for L1 C/A) — the estimator returns `0.0 dB-Hz` until at least one
prompt has been pushed in, and the estimate sharpens as the buffer fills.

```@docs
MomentsCN0Estimator
Tracking.update(::MomentsCN0Estimator, ::Any)
```

You can shrink or grow the buffer by building the `TrackedSat` yourself
and handing it to `add_satellite!`'s escape-hatch overload:

```jldoctest cn0_buffer
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        num_prompts_for_cn0_estimation = 200);

julia> track_state = add_satellite!(track_state, :default, sat);

julia> get_prn(track_state, 1)
1
```

## What the estimator returns

The estimator's output depends on what's been pushed into its buffer:

- **Before any `track` call has completed an integration** — the buffer is
  empty and `estimate_cn0` returns `0.0 dB-Hz`. This is the value seen
  immediately after `add_satellite!`.
- **Driving a noise-free signal through `track`** — the noise-moment
  estimate goes to zero, so the C/N₀ ratio diverges and the result is
  `Inf dB-Hz`. The [Quick start](index.md#Quick-start) shows this
  behavior.
- **Driving a noisy signal** — the estimate converges to the underlying
  CN0 as more prompts fill the buffer. Real-world signals typically land
  in 30–50 dB-Hz; values below ~25 dB-Hz are unusable for conventional
  tracking.

The doctest below builds a noisy L1 C/A signal at a known 45 dB-Hz CN0
and drives 25 1-ms tracking cycles through it. The estimator's output
converges to the input within a few dB:

```jldoctest cn0_example; filter = r"[0-9]+\.[0-9]+" => "***"
julia> using Tracking, GNSSSignals, Random

julia> using Tracking: Hz

julia> using GNSSSignals: gen_code, get_code_frequency

julia> function run_cn0_demo()
           sys = GPSL1CA()
           fs = 4e6Hz
           num_samples = 4000  # 1 ms at 4 MHz
           prn = 1
           cn0_db_hz = 45.0
           sigma = sqrt((fs/Hz) / 10^(cn0_db_hz/10) / 2)
           code_freq = get_code_frequency(sys)
           rng = MersenneTwister(0)
           track_state = TrackState(; signal = GPSL1CA())
           track_state = add_satellite!(track_state; prn, code_phase = 0.0, carrier_doppler = 0.0Hz)
           for _ in 1:25
               clean = gen_code(num_samples, sys, prn, fs, code_freq, 0.0)
               noise = sigma .* randn(rng, ComplexF64, num_samples)
               track_state = track(clean .+ noise, track_state, fs)
           end
           estimate_cn0(track_state, 1)
       end;

julia> run_cn0_demo()  # converges toward the 45 dB-Hz input as the buffer fills
49.1 dB-Hz
```

The 1-σ noise amplitude is derived from the target CN0 by inverting
`C/N₀ = signal_power / (σ² / Fs)`. With unit-amplitude code samples this
collapses to `σ² = Fs / 10^(CN0_dB/10)`, split evenly across the complex
sample's real and imaginary parts. The seeded `MersenneTwister(0)` keeps
the doctest deterministic; the bounds above (40–50 dB-Hz) leave room for
the natural spread of the moment estimator with a 25-prompt fill of the
default 100-sample buffer.

```@docs
estimate_cn0
```

## Custom CN0 Estimators

```@docs
AbstractCN0Estimator
```

You can implement your own estimator by creating a subtype of
[`AbstractCN0Estimator`](@ref) and implementing:

- `Tracking.update(cn0_estimator::MyCN0Estimator, prompt)` — return a new
  estimator with the latest prompt added (immutable update).
- `estimate_cn0(cn0_estimator::MyCN0Estimator, integration_time)` —
  return the CN0 estimate as a `dB-Hz` quantity.

To plug your custom estimator into a satellite, build a
[`TrackedSignal`](@ref) with your `cn0_estimator`, splice it into a fresh
[`TrackedSat`](@ref), and add it via the escape-hatch overload — see the
[`TrackedSat`](@ref) docstring for the full constructor surface.
