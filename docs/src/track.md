# Track

```@meta
CurrentModule = Tracking
```

[`track`](@ref) is the main entry point: it takes an incoming measurement and a [`TrackState`](@ref), runs downconversion + correlation for every tracked satellite, drives the Doppler estimator, and returns the updated `TrackState`. [`track!`](@ref) is the in-place counterpart for hard real-time loops where GC pauses must be avoided — after one warmup call (to seat the `filtered_prompts` buffer's capacity), the single-threaded `track!` path is **fully allocation-free**; the threaded path keeps a small irreducible per-call residual from Polyester's `@batch` closure capture (160 B per GNSS system).

```@docs
track
track!
```

## Optional parameters

- `downconvert_and_correlator` — the downconversion and correlation implementation. Defaults to `CPUThreadedDownconvertAndCorrelator()`. **For real-time loops, hoist this outside the loop** (see below).
- `intermediate_frequency` — the IF of the signal. Defaults to `0.0Hz`. Only accepted on the bare-buffer form `track!(buf, state, fs; intermediate_frequency = ...)`; on the [`BandMeasurement`](@ref) and multi-band forms the IF lives on each `BandMeasurement`.
- `doppler_update_interval` — the Doppler-estimation / NCO-update interval, a time (e.g. `1u"ms"`). Defaults to `nothing` ⇒ auto = the smallest primary-code period across all tracked signals (1 ms for GPS L1 C/A). Each measurement is processed in fixed-size chunks of this length: within a chunk the NCO Doppler is held fixed and every correlator output that completes is collected, then the estimator processes them in order and updates **every** satellite's NCO once, at a common epoch (see [Chunked Doppler updates](#Chunked-Doppler-updates)). Pick a longer interval to reduce Doppler-estimation cost at the expense of update rate.

The **coherent-integration length** is not a `track!` argument — it is a per-signal setting on each [`TrackedSignal`](@ref) (its `preferred_num_code_blocks_to_integrate` field), changed with [`set_preferred_num_code_blocks_to_integrate!`](@ref). It defaults to `1`, is capped per integration by the signal's bit/secondary-code period, and only takes effect once bit/secondary-code synchronization has been achieved. For data-bearing signals the length must evenly divide the number of code blocks that form one bit (e.g. a divisor of 20 for GPS L1 C/A, of 10 for GPS L5I) so integrations stay aligned to bit boundaries; other values throw an `ArgumentError`. With the conventional estimator the loop bandwidth auto-scales by `1/N` so longer integration stays stable without re-tuning.

```@docs
set_preferred_num_code_blocks_to_integrate!
```

## Chunked Doppler updates

`track` / `track!` walk each measurement in fixed-size time chunks of length `doppler_update_interval` (default: the smallest code period across all signals). Each chunk runs one correlate pass and one estimate:

1. **Correlate to the last completed boundary** — each satellite integrates from wherever it stands up to its last coherent-integration boundary inside the chunk; every completed integration is collected into that signal's `correlator_outputs` buffer, tagged with the sample index at which it ended (a [`CorrelatorOutput`](@ref); the sample index is important for vector tracking). A 1 ms-code signal in a 1 ms chunk yields 0, 1, or 2 outputs; a signal whose coherent integration is longer than the chunk yields outputs only on the chunks where it completes.
2. **Estimate** — the Doppler estimator processes the collected outputs **in order** (threading the loop-filter state across them) and writes the resulting Doppler to the NCO **once per chunk** — all satellites' NCOs update at the same point in the processing, a common epoch.

The chunk's trailing partial — from each satellite's last completed boundary to the chunk end — is *not* integrated separately: the **next** chunk's pass starts right at that boundary, so each integration runs boundary → boundary in one kernel window, entirely at the freshly updated Doppler. Every completed integration is therefore produced by a single NCO Doppler and each correction takes effect right at the boundary where its integration completed — the same loop timing as a classic per-code-period update. A final pass after the last chunk drains the buffer's trailing partial into each satellite's live accumulator so it carries into the next `track!` call.

Read the collected outputs for the most recent chunk with [`get_correlator_outputs`](@ref) (they are cleared after each chunk's estimate).

Choosing a larger `doppler_update_interval` batches more correlator outputs per NCO update, trading Doppler-tracking bandwidth for lower estimation cost. At the default interval the behavior is numerically very close to updating the NCO once per code period.

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
    # ... read per-sat state e.g. via get_sat_state(track_state, group, prn) ...
end
```

Both [`CPUDownconvertAndCorrelator`](@ref) and [`CPUThreadedDownconvertAndCorrelator`](@ref) hold long-lived per-thread scratch buffers that grow on first use and are reused thereafter. Relying on the default kwarg value rebuilds them on every call, which defeats the allocation-free design. This applies to both `track` (immutable) and `track!` (in-place).

`track!` writes back into the existing `Vector{TrackedSat}` slots of each per-group dictionary, so the tracking loop runs without GC pressure once the sat set is steady. The first `track!` call may grow each signal's `filtered_prompts` buffer via `push!`; from the second call onwards the capacity is settled.

## External correlator producers

The correlate phase and the Doppler estimator are decoupled: the estimator consumes each signal's `correlator_outputs` buffer and never touches the sample buffer directly. That lets an **external producer** — e.g. an FPGA/hardware correlator that streams completed correlator dumps — supply the outputs and run only the loop filters on the host. The FPGA does downconversion + correlation; the host folds the dumps through the estimator and streams the resulting NCO Dopplers back. (The DMA transport, wire format and NCO marshaling live in [GNSSReceiver.jl](https://github.com/JuliaGNSS/GNSSReceiver.jl); this section is only the Tracking.jl-side contract.)

The offload loop, per processing chunk (epoch):

1. **Ingest.** For each satellite/signal, build a [`CorrelatorOutput`](@ref) from the producer's raw accumulator and append it with [`append_correlator_output!`](@ref) — appended **per signal in `sample_index` order**:

   ```julia
   append_correlator_output!(track_state, output, group, prn, sig)
   ```

   Prefer this over mutating the vector [`get_correlator_outputs`](@ref) returns: it documents intent and type-checks the correlator against the signal's.

2. **Estimate.** Fold the batch and update every satellite's NCO once by calling [`estimate_dopplers_and_filter_prompt!`](@ref) with a **per-band sampling-frequency source** instead of a sample buffer — a `NamedTuple`/`Dict` keyed by the band's `get_band_id` (e.g. `:L1`, `:L5`):

   ```julia
   estimate_dopplers_and_filter_prompt!(track_state, (L1 = 25e6Hz, L5 = 25e6Hz))
   # or: estimate_dopplers_and_filter_prompt!(track_state, Dict(:L1 => 25e6Hz))
   ```

   This skips `downconvert_and_correlate!` entirely. The rate must be **per band** (the estimator walks groups that may be on different bands); passing the original [`BandMeasurements`](@ref) works too and reads the rate off each band's `BandMeasurement`, so the CPU `track!` path is unchanged. The estimator **consumes and clears** each signal's `correlator_outputs` as part of this call — it is empty again afterwards, ready for the next chunk. This is the public contract external producers rely on.

3. **Marshal** the updated `code_doppler`/`carrier_doppler` (read with [`get_code_doppler`](@ref)/[`get_carrier_doppler`](@ref)) back to the hardware NCOs.

### Caller contract

- `CorrelatorOutput.correlator` — the **raw** accumulator (sum-of-products over `integrated_samples`, matching what `normalize` expects). Reusing [`EarlyPromptLateCorrelator`](@ref) / `update_accumulator` on the producer side satisfies this by construction.
- `CorrelatorOutput.integrated_samples` — the producer's true sample count for that integration.
- `CorrelatorOutput.sample_index` — the chunk-relative end sample. The software path writes it buffer-relative (`signal_start_sample` returns to 1 each `track!`); a producer with a free-running **global** sample counter must subtract the current chunk/epoch origin so every satellite reads a consistent per-chunk time grid (the estimator itself does not read it — it is preserved for downstream vector/Kalman tracking).

### Transport delay

The DLL discriminator uses the satellite's **pre-update** `code_doppler` as the value in effect during the integration — correct for a one-epoch feedback pipeline (the NCO written from this chunk's estimate takes effect on the next chunk). If your hardware pipeline has deeper feedback delay, that is the caller's responsibility: tag outputs with their epoch and schedule each NCO update to a known future epoch so the loop stays consistent with when the Doppler actually reaches the correlator.

## BandMeasurement

One band's incoming sample buffer plus the front-end metadata needed to process it. Bundles `samples` with `sampling_frequency` and `intermediate_frequency` — these are inseparable in practice, and the bundle removes the chance of mismatched parallel NamedTuples in a multi-band `track` call.

For single-band tracking, the bare-buffer form `track!(buf, state, fs)` and the single-`BandMeasurement` form `track!(BandMeasurement(buf, fs), state)` both auto-wrap into a one-entry NamedTuple internally — see the [Quick start](index.md#Quick-start) for the bare-buffer form.

For multi-band tracking, build one `BandMeasurement` per band and pass them as a NamedTuple keyed by the band's `GNSSSignals.get_band_id` (e.g. `:L1`, `:L5`):

```jldoctest multi_band_call
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(;
           signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)),
       );

julia> track_state = add_satellite!(track_state; prn = 1, group = :legacy_gps_l1, code_phase = 0.0, carrier_doppler = 0.0Hz);

julia> track_state = add_satellite!(track_state; prn = 1, group = :gps_l5,        code_phase = 0.0, carrier_doppler = 0.0Hz);

julia> buf_l1 = zeros(ComplexF64, 4000);   # 1 ms at  4 MHz

julia> buf_l5 = zeros(ComplexF64, 25000);  # 1 ms at 25 MHz

julia> track!((L1 = BandMeasurement(buf_l1, 4e6Hz),
               L5 = BandMeasurement(buf_l5, 25e6Hz)), track_state);
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
non-`Complex{Int16}` measurement. Its constructor takes one **required** positional
argument, `max_meas` — your front end's full-scale (the largest `|real|`/`|imag|` a
sample can take, e.g. `2^11` for a 12-bit ADC) — from which it sizes the carrier
replica so the integer carrier wipe cannot overflow. There is no default: see
[`Int16DownconvertAndCorrelator`](@ref) for why under-declaring it is catastrophic.

For an even faster **bit-wise** option on `Complex{Int16}` captures of **BPSK**
signals, [`OneBitThreadedDownconvertAndCorrelator`](@ref) (and its single-threaded
sibling [`OneBitDownconvertAndCorrelator`](@ref)) hard-limits the measurement,
carrier and code to a single sign bit, so downconversion becomes XOR and the tap
accumulate becomes popcount — measured ~1.5–5.6× faster than the Float32 backend
(the gap grows with sampling rate). It trades ≈2–3 dB of SNR for that speed; since
the discriminators, C/N0 and bit buffer are ratio-normalised, the coarse amplitude
is immaterial. Bit-wise correlation is awkward for non-binary modulations, so this
backend is BPSK-only and errors on CBOC/BOC code types.

Between the two sits [`TwoBitThreadedDownconvertAndCorrelator`](@ref) (and its
single-threaded sibling [`TwoBitDownconvertAndCorrelator`](@ref)): the same
bit-plane XOR + popcount machinery, but with a second **magnitude** bit for the
measurement and the carrier, making both 4-level `{±1, ±3}` quantities (the
carrier's sign and magnitude bit planes come straight off SinCosLUT's 2-bit NCO).
That recovers ≈2 dB of the one-bit backend's SNR loss (leaving only ≈0.8 dB vs
Float32) at roughly Int16 speed — so pick **one-bit** when raw speed matters most,
**two-bit** when you want near-Int16 sensitivity with the bit-wise memory/layout
advantages, and **Int16** when quantisation loss must be negligible. The `threshold`
keyword sets the measurement's magnitude split point in ADC counts (≈1σ of the
front end's input is the classic near-optimal choice; the default 512 suits a
properly-AGC'd 12-bit capture). Same scope as one-bit: `Complex{Int16}` samples,
binary (BPSK) codes only.

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
