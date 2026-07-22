# Chunked Doppler-update tracking loop

**Date:** 2026-07-21
**Status:** Implemented

## Goal

Decouple the Doppler-estimation (NCO-update) rate from the code period. Process
each measurement in fixed-size **time chunks**; within a chunk the NCO Doppler
is held fixed and every correlator output that completes is collected — tagged
with the sample index at which it was taken — into a per-signal array. After the
chunk, the estimator folds over those outputs in order and updates every
satellite's NCO once, at a common epoch.

## Motivation

Previously `track!` advanced each satellite to its *next code-block boundary*
and immediately re-estimated and rewrote that satellite's NCO Doppler. Three
limitations:

- The update rate was welded to the code period — no way to trade
  Doppler-estimation compute for a slower update rate.
- Correlator outputs were consumed one at a time and never collected; the
  estimator could only ever see the single just-completed correlator. Vector
  tracking needs a *batch* of outputs across channels, each tagged with its
  sample epoch, on a common time grid.
- NCO updates were staggered — each satellite updated its NCO the instant *it*
  finished a code period, so different satellites updated at different times.

## Design

### Chunk grid

The chunk length is the Doppler-update interval, a **time** (`update_interval`
kwarg on [`track`](@ref) / [`track!`](@ref)). Default `nothing` ⇒ auto = the
smallest primary-code period across all tracked signals (`_smallest_code_period`;
e.g. 1 ms for GPS L1 C/A, or for a mixed L1 C/A + Galileo E1B state). For band
`b` at sampling frequency `fs_b`, chunk `k` (0-based) ends at sample
`min(round((k+1)·Δt·fs_b), num_samples_b)` (`_chunk_last_sample`). The boundary
is re-anchored to the absolute chunk index each call (not accumulated), so
rounding never drifts and different bands stay time-aligned to within a sample.
A one-sample floor is validated up front (`_validate_update_interval`) so the
loop always makes progress.

### Correlate phase (records outputs, two passes per chunk)

`_update_tracked_sat_correlator` runs an inner loop over the chunk. Each
sub-step integrates to the next code-block boundary of any signal (or the chunk
end), reusing `_calc_min_samples_and_completed` with the chunk end substituted
for the buffer end (the true buffer length is still used for replica sizing).
When a signal's integration completes, `update` / `_build_new_signals` snapshots
its raw accumulator into a [`CorrelatorOutput`](@ref) — `(correlator,
integrated_samples, sample_index)` — appended to the reused
per-signal `correlator_outputs` buffer, and resets the accumulator so the next
integration in the chunk starts fresh. A partial integration at the chunk (or
buffer) boundary carries in the accumulator. The NCO Doppler is untouched here.
A chunk yields 0, 1, or several outputs per signal depending on code phase and
Doppler.

`track!` runs this correlate step **twice per chunk**, around the estimate:

1. with `stop_before_partial = true` — integrate every completion inside the
   chunk, but stop at the last code-block boundary instead of integrating the
   trailing chunk-clamped partial (a sub-step that completes no signal is
   exactly that partial);
2. after the estimate has written the new NCO Doppler — integrate the residue
   (last boundary → chunk end), which is the head of the *next* integration.

This keeps every completed integration on a single NCO Doppler and applies each
Doppler correction right at the boundary where its integration completed — the
same loop timing as the classic per-completion update — while still estimating
once per chunk at a common epoch. (Without the split, the residue was
integrated at the pre-update Doppler, which measurably tightened the FLL
pull-in edge.) The updated code Doppler can, rarely, pull a boundary that pass
1 measured as beyond the chunk back inside it; such an output is folded by the
next chunk's estimate, and `track!` runs one trailing fold after the last chunk
so nothing is dropped at the buffer end.

### Estimate phase (folds, one NCO write)

`_apply_correlator_output` applies one record: normalize by its sample count,
post-corr filter, push the filtered prompt, advance CN0 and bit buffer, move the
record to `last_fully_integrated_*` (it no longer resets the accumulator — the
correlate phase already did). The driver processor
`_process_estimator_driver_signal` folds over `correlator_outputs` in order,
threading the loop-filter state and the FLL `previous_prompt` across records
(per-record `integration_time` and `1/N` bandwidth scaling), computing a Doppler
each record and keeping the **last**; the NCO is written once. Passenger signals
apply the per-record advance only. Each signal's `correlator_outputs` is
`empty!`-ed after folding. With no outputs in a chunk the Doppler holds.

### Phase-snap at sync

Sync is detected during the estimate pass, *after* the whole chunk was
correlated. The one-time secondary-code phase snap
(`_snap_code_phase_from_synced_signal`) still runs against the end-of-chunk
`code_phase`, preserving the within-primary-block phase (issue #117). Because any
in-flight partial for the sync chunk was accumulated at the *pre-snap* phase (and
pre-sync without a secondary overlay), the snap now also resets every signal's
in-flight accumulator (`_reset_inflight_integration`): the phase bookkeeping is
kept, but the next chunk re-integrates cleanly from the snapped phase. Otherwise
a flipped overlay chip can cancel the partial to zero and feed a 0/0 into the
discriminators. (With the two-pass chunk the estimate usually runs exactly on a
boundary, so the reset is a no-op there; it still matters for signals whose
integration spans multiple chunks and thus carry a partial at estimate time.)

Records that *follow* the syncing record within the same fold were correlated
with pre-sync replicas (no secondary-code wipe-off), so accumulating them into
the bit buffer as if wiped would feed sign-corrupted prompts into the first
post-sync bits. The folds detect the mid-fold `found` transition and apply
those records with `skip_bit_buffer = true` (`_apply_correlator_output`):
prompt filter, CN0, and loop-filter processing still run, only the bit-buffer
update is skipped. Irrelevant at the default interval (a chunk then holds at
most one record past the sync), it protects enlarged `update_interval`s.

### Data model

- [`CorrelatorOutput`](@ref)`{C}` in `src/correlators/correlator.jl`.
- [`TrackedSignal`](@ref) gains `correlator_outputs::Vector{CorrelatorOutput{C}}`,
  preallocated and reused like `filtered_prompts` (`sizehint!`-ed, `empty!`-ed at
  each `track` start and after each chunk's fold), so steady-state tracking stays
  allocation-free.

## Backends

All new logic lives in the shared `_update_tracked_sat_correlator` / `update` /
estimator layer, above the backend boundary `_correlate_signals`. The CPU-fused
and Int16 backends (single + threaded) need no changes; the per-`_correlate_signals`
call count is unchanged from before.

## Behavior compatibility

With the two-pass chunk, every completed integration is produced by a single
NCO Doppler and each correction takes effect at its completing boundary — the
same loop timing as the previous per-completion update. The suite's edge
probes match the pre-chunking behavior exactly: the FLL pull-in edge sits at
its historical 240 Hz (converges) / 250 Hz (fails), and the 1800-block L1C-P
overlay sync locks after exactly one overlay cycle. (An earlier single-pass
design integrated the chunk residue at the pre-update Doppler, which tightened
the FLL edge to ~210 Hz and could land the overlay sync one block late.)

## Deferred follow-ups

- A global (cross-`track!`-call) sample index for streamed real-time chunks
  (currently `sample_index` is relative to the current measurement buffer).
- The vector-tracking consumer of the batched, sample-indexed correlator outputs.
