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

### Correlate phase (records outputs)

`_update_tracked_sat_correlator` runs an inner loop over the chunk. Each
sub-step integrates to the next code-block boundary of any signal (or the chunk
end), reusing `_calc_min_samples_and_completed` with the chunk end substituted
for the buffer end (the true buffer length is still used for replica sizing).
When a signal's integration completes, `update` / `_build_new_signals` snapshots
its raw accumulator into a [`CorrelatorOutput`](@ref) — `(correlator,
integrated_samples, sample_index)` — appended to the reused
per-signal `correlator_outputs` buffer, and resets the accumulator so the next
integration in the chunk starts fresh. A partial integration at the chunk (or
buffer) boundary carries in the accumulator. The NCO Doppler is untouched here,
so it stays fixed for the whole chunk. So a chunk yields 0, 1, or several
outputs per signal depending on code phase and Doppler.

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
discriminators.

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

At the default (chunk = smallest code period) the result is numerically very
close to the previous per-code-period behavior, but **not** bit-identical: a
code period straddling a chunk boundary has its leading residue integrated with
the pre-update Doppler, and — unlike the old per-completion update — a period
occasionally skipped by a one-sample floating-point drift is not "caught up" by
a second completion in a later buffer. The observable effects are tiny and
confined to convergence edges: the FLL pull-in edge tightened from ~240 Hz to
~210 Hz for the flagship convergence test, and a 1800-block overlay sync can land
one block late. Both are recalibrated in the test suite.

## Deferred follow-ups

- A global (cross-`track!`-call) sample index for streamed real-time chunks
  (currently `sample_index` is relative to the current measurement buffer).
- The vector-tracking consumer of the batched, sample-indexed correlator outputs.
