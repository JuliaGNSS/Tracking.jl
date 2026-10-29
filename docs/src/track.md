# Track

```@meta
CurrentModule = Tracking
```

[`track`](@ref) takes a measurement and a [`TrackState`](@ref), downconverts and correlates every tracked satellite, drives the Doppler estimator, and returns the updated `TrackState`. [`track!`](@ref) is the in-place counterpart for real-time loops; after one warmup call its single-threaded path is **allocation-free** (see [`track!`](@ref) for the threaded path's residual).

```@docs
track
track!
```

## Optional parameters

- `downconvert_and_correlator`: the backend, by default `CPUThreadedDownconvertAndCorrelator()`. **Hoist it outside real-time loops** (see below).
- `intermediate_frequency`: the IF, by default `0.0Hz`. Bare-buffer form only; otherwise it lives on each [`BandMeasurement`](@ref).
- `doppler_update_interval`: the NCO-update interval, a time (e.g. `1u"ms"`). `nothing` means the smallest primary-code period of the tracked signals (1 ms for GPS L1 C/A); see [Chunked Doppler updates](#Chunked-Doppler-updates). Longer intervals cost less estimation at a lower update rate.

The **coherent-integration length** is set per signal with [`set_preferred_num_code_blocks_to_integrate!`](@ref). It starts at [`default_num_code_blocks_to_integrate`](@extref TrackingLoops.default_num_code_blocks_to_integrate-Tuple{AbstractGNSSSignal}) (`1`, except a 31-block, 2 ms code cycle for Galileo E5a-QP, whose 64.5 µs code period is too short to run a loop on), is capped by [`max_num_code_blocks_to_integrate`](@extref TrackingLoops.max_num_code_blocks_to_integrate-Tuple{AbstractGNSSSignal}) (one bit or secondary-code period), and takes effect only after bit or secondary-code sync.

```@docs
set_preferred_num_code_blocks_to_integrate!
```

In the TrackingLoops manual:

- [`default_num_code_blocks_to_integrate`](@extref TrackingLoops.default_num_code_blocks_to_integrate-Tuple{AbstractGNSSSignal})
- [`max_num_code_blocks_to_integrate`](@extref TrackingLoops.max_num_code_blocks_to_integrate-Tuple{AbstractGNSSSignal})

## Chunked Doppler updates

`track` / `track!` walk each measurement in chunks of `doppler_update_interval`, each with one correlate pass and one estimate:

1. **Correlate** each satellite up to its last completed integration boundary in the chunk, collecting every completed integration into the signal's `correlator_outputs` as a [`CorrelatorOutput`](@extref TrackingLoops.CorrelatorOutput) tagged with its end sample. A 1 ms code in a 1 ms chunk yields 0, 1 or 2 outputs.
2. **Estimate**: the Doppler estimator folds the outputs in order and writes every satellite's NCO **once per chunk**, at a common epoch.

The next chunk's pass starts at that boundary, so every integration runs at one NCO Doppler and each correction takes effect at the boundary where its integration completed, as with a per-code-period update. A final pass drains each satellite's trailing partial into its accumulator, carried into the next `track!` call. [`get_correlator_outputs`](@ref) reads a chunk's outputs before the estimate clears them.

## Real-time loops

A typical receiver loop builds the `TrackState` once, hoists the correlator outside the loop, and calls `track!` per chunk:

```julia
track_state = TrackState(; signal = GPSL1CA())
track_state = add_satellite!(track_state; prn = 1, code_phase = 0.0, carrier_doppler = 1000.0Hz)
# ... more sats ...

dc = CPUThreadedDownconvertAndCorrelator()  # hoist outside the loop

while got_signal_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_frequency;
           downconvert_and_correlator = dc)
    # ... read per-sat state e.g. via get_sat_state(track_state, group_key, prn) ...
end
```

The backends hold per-thread scratch buffers that grow on first use; the keyword's default rebuilds them on every call, for `track` and `track!` alike. `track!` writes into the existing satellite slots, and from its second call on the buffers' capacity is settled.

## External correlator producers

The Doppler estimator reads only each signal's `correlator_outputs`, never the samples, so an **external producer** (e.g. an FPGA streaming correlator dumps) can supply them and the host runs only the loops. (Transport, wire format and NCO marshaling live in [GNSSReceiver.jl](https://github.com/JuliaGNSS/GNSSReceiver.jl).)

Per chunk:

1. **Ingest.** Build a [`CorrelatorOutput`](@extref TrackingLoops.CorrelatorOutput) from the producer's raw accumulator and append it with [`append_correlator_output!`](@ref), **per signal in `sample_index` order**:

   ```julia
   append_correlator_output!(track_state, output, group_key, prn, sig)
   ```

2. **Estimate.** Call [`estimate_dopplers_and_filter_prompt!`](@ref) with the **per-band sampling frequencies** (a `NamedTuple` or `Dict` keyed by `get_band_id`) in place of the samples:

   ```julia
   estimate_dopplers_and_filter_prompt!(track_state, (L1 = 25e6Hz, L5 = 25e6Hz))
   # or: estimate_dopplers_and_filter_prompt!(track_state, Dict(:L1 => 25e6Hz))
   ```

   [`BandMeasurements`](@ref Tracking.BandMeasurements) work too. The call **consumes and clears** each signal's `correlator_outputs`.

3. **Marshal** the updated [`get_code_doppler`](@ref) / [`get_carrier_doppler`](@ref) back to the hardware NCOs.

### Caller contract

- `CorrelatorOutput.correlator`: the **raw** accumulator (sum of products over `integrated_samples`, as `normalize` expects), e.g. built with [`EarlyPromptLateCorrelator`](@extref TrackingLoops.EarlyPromptLateCorrelator) / `update_accumulator`.
- `CorrelatorOutput.integrated_samples`: the producer's true sample count.
- `CorrelatorOutput.sample_index`: the end sample relative to the start of the current measurement, as the software path writes it (`signal_start_sample` returns to 1 each `track!`). A producer with a free-running **global** counter subtracts the measurement's origin; the estimate phase adds the band's [`TrackState`](@ref) `sample_offsets` before handing the record to the Doppler estimator.

### Transport delay

The DLL discriminator takes the **pre-update** `code_doppler` as the one in effect during the integration, right for a one-chunk feedback pipeline. With a deeper delay, tag outputs with their epoch and schedule each NCO update for a known future epoch.

## BandMeasurement

One band's samples with their `sampling_frequency` and `intermediate_frequency`. A single band may pass a bare buffer, `track!(buf, state, fs)`, or one `BandMeasurement`; several bands pass a NamedTuple keyed by `GNSSSignals.get_band_id`:

```jldoctest multi_band_call
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> using Random: Xoshiro

julia> track_state = TrackState(;
           signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)),
       );

julia> track_state = add_satellite!(track_state; prn = 1, group_key = :legacy_gps_l1, code_phase = 0.0, carrier_doppler = 0.0Hz);

julia> track_state = add_satellite!(track_state; prn = 1, group_key = :gps_l5,        code_phase = 0.0, carrier_doppler = 0.0Hz);

julia> buf_l1 = randn(Xoshiro(1), ComplexF64, 4000);   # 1 ms at  4 MHz

julia> buf_l5 = randn(Xoshiro(2), ComplexF64, 25000);  # 1 ms at 25 MHz

julia> track!((L1 = BandMeasurement(buf_l1, 4e6Hz),
               L5 = BandMeasurement(buf_l5, 25e6Hz)), track_state);
```

Noise, not `zeros`: the default C/N₀ estimator measures the noise floor from the samples, and `track!` warns, once per signal, about an all-zero buffer (see [`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)). See [Multi-band tracking](tracking_state.md#Multi-band-tracking) for the full setup.

```@docs
BandMeasurement
Tracking.BandMeasurements
```

## Downconversion and correlation

The default [`CPUThreadedDownconvertAndCorrelator`](@ref) runs in Float32 and
accepts any complex sample type. For `Complex{Int16}` samples there are three
opt-in backends, each with a threaded and a single-threaded variant, selected via
the `downconvert_and_correlator` keyword:

- [`Int16ThreadedDownconvertAndCorrelator`](@ref) /
  [`Int16DownconvertAndCorrelator`](@ref): integer arithmetic, typically
  ~1.3–2.9× faster. Its **required** `max_meas` is the front end's full scale
  (the largest `|real|`/`|imag|`, e.g. `2^11` for a 12-bit ADC), which sizes the
  carrier replica against overflow; see [`Int16DownconvertAndCorrelator`](@ref)
  for why under-declaring it is catastrophic.
- [`OneBitThreadedDownconvertAndCorrelator`](@ref) /
  [`OneBitDownconvertAndCorrelator`](@ref): measurement, carrier and code
  hard-limited to a sign bit, so downconversion is XOR and accumulation popcount;
  ~1.5–5.6× faster than Float32 (more at higher rates) for ≈2–3 dB of SNR. The
  ratio-normalised discriminators, C/N₀ and bit buffer do not mind the coarse
  amplitude. BPSK only.
- [`TwoBitThreadedDownconvertAndCorrelator`](@ref) /
  [`TwoBitDownconvertAndCorrelator`](@ref): adds a **magnitude** bit to
  measurement and carrier (`{±1, ±3}`), recovering ≈2 dB of the one-bit loss
  (≈0.8 dB short of Float32) at about Int16 speed. `threshold` is the magnitude
  split in ADC counts (≈1σ of the input is near-optimal; the default 512 suits an
  AGC'd 12-bit capture). BPSK only.

```@docs
CPUDownconvertAndCorrelator
CPUThreadedDownconvertAndCorrelator
Int16DownconvertAndCorrelator
Int16ThreadedDownconvertAndCorrelator
OneBitDownconvertAndCorrelator
OneBitThreadedDownconvertAndCorrelator
TwoBitDownconvertAndCorrelator
TwoBitThreadedDownconvertAndCorrelator
Tracking.AbstractDownconvertAndCorrelator
Tracking.downconvert_and_correlate
Tracking.downconvert_and_correlate!
```

Correlator sample shifts and the early/late spacing are documented in [Correlator](correlator.md#Sample-shifts).
