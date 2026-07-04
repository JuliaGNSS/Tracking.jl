# Track

[`track`](@ref) is the main entry point: it takes an incoming measurement and a [`TrackState`](@ref), runs downconversion + correlation for every tracked satellite, drives the Doppler estimator, and returns the updated `TrackState`. [`track!`](@ref) is the in-place counterpart for hard real-time loops where GC pauses must be avoided — after one warmup call (to seat the `filtered_prompts` buffer's capacity), the single-threaded `track!` path is **fully allocation-free**; the threaded path keeps a small irreducible per-call residual from Polyester's `@batch` closure capture (160 B per GNSS system).

```@docs
track
track!
```

## Optional parameters

  - `downconvert_and_correlator` — the downconversion and correlation implementation. Defaults to `CPUThreadedDownconvertAndCorrelator()`. **For real-time loops, hoist this outside the loop** (see below).
  - `intermediate_frequency` — the IF of the signal. Defaults to `0.0Hz`. Only accepted on the bare-buffer form `track!(buf, state, fs; intermediate_frequency = ...)`; on the [`BandMeasurement`](@ref) and multi-band forms the IF lives on each `BandMeasurement`.
    The **coherent-integration length** is not a `track!` argument — it is a per-signal setting on each [`TrackedSignal`](@ref) (its `preferred_num_code_blocks_to_integrate` field), changed with [`set_preferred_num_code_blocks_to_integrate!`](@ref). It defaults to `1`, is capped per integration by the signal's bit/secondary-code period, and only takes effect once bit/secondary-code synchronization has been achieved. For data-bearing signals the length must evenly divide the number of code blocks that form one bit (e.g. a divisor of 20 for GPS L1 C/A, of 10 for GPS L5I) so integrations stay aligned to bit boundaries; other values throw an `ArgumentError`. With the conventional estimator the loop bandwidth auto-scales by `1/N` so longer integration stays stable without re-tuning.

```@docs
set_preferred_num_code_blocks_to_integrate!
```

## Real-time loops

A typical receiver loop builds the `TrackState` once, hoists the correlator outside the loop, and calls `track!` per chunk:

```julia
track_state = TrackState(; signal = GPSL1CA())
track_state =
    add_satellite!(track_state; prn = 1, code_phase = 0.0, carrier_doppler = 1000.0Hz)
# ... more sats ...

dc = CPUThreadedDownconvertAndCorrelator()  # hoist outside the loop

while got_signal_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
    # ... read per-sat state e.g. via get_sat_state(track_state, group, prn) ...
end
```

Both [`CPUDownconvertAndCorrelator`](@ref) and [`CPUThreadedDownconvertAndCorrelator`](@ref) hold long-lived per-thread scratch buffers that grow on first use and are reused thereafter. Relying on the default kwarg value rebuilds them on every call, which defeats the allocation-free design. This applies to both `track` (immutable) and `track!` (in-place).

`track!` writes back into the existing `Vector{TrackedSat}` slots of each per-group dictionary, so the tracking loop runs without GC pressure once the sat set is steady. The first `track!` call may grow each signal's `filtered_prompts` buffer via `push!`; from the second call onwards the capacity is settled.

## BandMeasurement

One band's incoming sample buffer plus the front-end metadata needed to process it. Bundles `samples` with `sampling_frequency` and `intermediate_frequency` — these are inseparable in practice, and the bundle removes the chance of mismatched parallel NamedTuples in a multi-band `track` call.

For single-band tracking, the bare-buffer form `track!(buf, state, fs)` and the single-`BandMeasurement` form `track!(BandMeasurement(buf, fs), state)` both auto-wrap into a one-entry NamedTuple internally — see the [Quick start](index.md#Quick-start) for the bare-buffer form.

For multi-band tracking, build one `BandMeasurement` per band and pass them as a NamedTuple keyed by [`band_key`](@ref):

```jldoctest multi_band_call
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)));

julia> track_state = add_satellite!(
           track_state;
           prn = 1,
           group = :legacy_gps_l1,
           code_phase = 0.0,
           carrier_doppler = 0.0Hz,
       );

julia> track_state = add_satellite!(
           track_state;
           prn = 1,
           group = :gps_l5,
           code_phase = 0.0,
           carrier_doppler = 0.0Hz,
       );

julia> buf_l1 = zeros(ComplexF64, 4000);   # 1 ms at  4 MHz

julia> buf_l5 = zeros(ComplexF64, 25000);  # 1 ms at 25 MHz

julia> track!(
           (l1 = BandMeasurement(buf_l1, 4e6Hz), l5 = BandMeasurement(buf_l5, 25e6Hz)),
           track_state,
       );

```

See [Multi-band tracking](tracking_state.md#Multi-band-tracking) for the full setup (group declaration, per-band antenna counts, duration matching).

```@docs
BandMeasurement
BandMeasurements
```

## Downconversion and correlation

The default [`CPUThreadedDownconvertAndCorrelator`](@ref) runs a Float32 pipeline
and accepts any complex sample type. For `Complex{Int16}` (integer ADC) sample
buffers there is an opt-in integer backend,
[`Int16ThreadedDownconvertAndCorrelator`](@ref) (and its single-threaded sibling
[`Int16DownconvertAndCorrelator`](@ref)), which is typically ~1.3–2.9× faster.
Select it explicitly via the `downconvert_and_correlator` keyword; it errors on a
non-`Complex{Int16}` measurement.

For an even faster **bit-wise** option on `Complex{Int16}` captures of **BPSK**
signals, [`OneBitThreadedDownconvertAndCorrelator`](@ref) (and its single-threaded
sibling [`OneBitDownconvertAndCorrelator`](@ref)) hard-limits the measurement,
carrier and code to a single sign bit, so downconversion becomes XOR and the tap
accumulate becomes popcount — measured ~1.5–5.6× faster than the Float32 backend
(the gap grows with sampling rate). It trades ≈2–3 dB of SNR for that speed; since
the discriminators, C/N0 and bit buffer are ratio-normalised, the coarse amplitude
is immaterial. Bit-wise correlation is awkward for non-binary modulations, so this
backend is BPSK-only and errors on CBOC/BOC code types.

Between the one-bit and Int16 backends sits a **two-bit** sign+magnitude option,
[`TwoBitThreadedDownconvertAndCorrelator`](@ref) (and its single-threaded sibling
[`TwoBitDownconvertAndCorrelator`](@ref)). It keeps a magnitude bit as well as a
sign bit for the measurement — the near-optimal 4-level `{±1, ±3}` quantiser — so
it still correlates with XOR + popcount but recovers ≈1.4 dB of the SNR the
one-bit backend gives up, while staying in the same performance class as the
Int16 backend (≈2–3× faster than Float32).
The measurement is always two-bit; the carrier precision is selectable via the
`carrier_bits` keyword (`1`, the default sign-only carrier, or `2` for a
sign+magnitude carrier that recovers a further ≈0.6 dB). The measurement magnitude
split point is the `threshold` keyword (ADC counts, set near 1σ of your front
end's input). Like the one-bit backend it is BPSK-only and requires
`Complex{Int16}` samples.

```@docs
CPUDownconvertAndCorrelator
CPUThreadedDownconvertAndCorrelator
Int16DownconvertAndCorrelator
Int16ThreadedDownconvertAndCorrelator
OneBitDownconvertAndCorrelator
OneBitThreadedDownconvertAndCorrelator
TwoBitDownconvertAndCorrelator
TwoBitThreadedDownconvertAndCorrelator
AbstractDownconvertAndCorrelator
```

Correlator sample shifts and the early/late spacing are documented in [Correlator](correlator.md#Sample-shifts).
