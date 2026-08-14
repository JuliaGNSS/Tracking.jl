"""
$(SIGNATURES)

Abstract supertype for CN0 (carrier-to-noise-density ratio) estimators.
Each [`TrackedSignal`](@ref) holds one estimator instance, stored in a type
parameter — pass any subtype instance as the `cn0_estimator` keyword of
[`TrackedSignal`](@ref) / [`TrackedSat`](@ref) to replace the default
[`NoiseRefCN0Estimator`](@ref); see [`default_cn0_estimator`](@ref) for which to
pick when. Custom estimators subtype this and implement `Tracking.update` and
[`estimate_cn0`](@ref), plus [`requires_noise_density`](@ref) if they read a
measured noise floor.
"""
abstract type AbstractCN0Estimator end

"""
$(SIGNATURES)

Per-record side information handed to `Tracking.update` alongside the prompt:
everything the tracking loop knows about the record that an estimator cannot
recover from the prompt stream alone. Its reason for existing is the
**navigation-bit grid** — where the data-bit boundaries sit and whether they are
known yet — which is what lets an estimator sum prompts coherently over exactly
one symbol ([`NWPRCN0Estimator`](@ref) does) instead of straddling a bit flip.

Fields:

  - `signal` — the signal the record belongs to.
  - `num_code_blocks` — primary-code blocks this record spanned. One for the
    default configuration; more when the correlate step was lengthened by
    [`set_preferred_num_code_blocks_to_integrate!`](@ref) or an external
    producer handed over longer records.
  - `num_code_blocks_per_bit` — blocks that form one navigation bit (symbol) of
    `signal`; `20` for GPS L1 C/A, `1` for GPS L1C-D / Galileo E1B, and `0` for a
    pilot, which carries no data and therefore has **no** bit grid — post-sync
    its prompts stay coherent for as long as the loops hold them.
  - `bit_code_block_index` — blocks already accumulated into the currently open
    navigation bit *before* this record, i.e. this record's own offset inside
    the bit; `0` means the record starts a fresh bit. It is `-1` whenever this
    record's prompt cannot be summed coherently against the grid: before the
    signal's bit / secondary-code sync has been found, and — for a
    secondary-coded signal only — for a record of the fold that follows a sync
    detected in that same fold, which was correlated without the overlay
    wipe-off the sync just established (see `_apply_correlator_output`). Such a
    record's *blocks* still count towards the grid, so the offset stays right for
    the records after it.
  - `bit_buffer` — the signal's `BitBuffer` as of *before* this record:
    the soft bits decoded so far, the open coherent bit accumulator, the lock
    polarity and the secondary-code phase.
  - `noise_density` — **this signal's** measured noise density `N₀` (dimension
    `1/Hz`), for an estimator that divides the prompt's power by a *measured*
    floor instead of inferring one from the prompt stream
    ([`NoiseRefCN0Estimator`](@ref) does). Per signal and not per band because
    the floor is the post-correlation one — see [`AbstractNoiseEstimator`](@ref).
    It is a **type parameter, not a sentinel**: `Nothing` says no
    [`AbstractNoiseEstimator`](@ref) is configured for this signal at all, which
    is a static property of the setup, so every call site monomorphises and the
    estimators stay allocation-free. A configured
    source whose window is merely still empty never reaches here — the fold
    skips the update instead, so this field is unconditionally a plain scalar
    whenever it is not `Nothing`.
  - `integration_time` — this record's own `T`, so an estimator that works per
    record does not have to be told one `T` for a whole ring of records
    integrated at different lengths. `nothing` when the caller did not supply
    one (the bare-prompt-stream path).
"""
struct CN0UpdateContext{S<:AbstractGNSSSignal,B<:Unsigned,N,T}
    signal::S
    num_code_blocks::Int
    num_code_blocks_per_bit::Int
    bit_code_block_index::Int
    bit_buffer::BitBuffer{B}
    noise_density::N
    integration_time::T
end

# Positional convenience for the five load-bearing fields, defaulting the two
# that only a noise-referenced estimator reads. Keeps the pre-6.1 call shape
# working: an estimator that ignores them is unaffected by their arrival.
@inline CN0UpdateContext(
    signal::AbstractGNSSSignal,
    num_code_blocks::Integer,
    num_code_blocks_per_bit::Integer,
    bit_code_block_index::Integer,
    bit_buffer::BitBuffer,
) = CN0UpdateContext(
    signal,
    Int(num_code_blocks),
    Int(num_code_blocks_per_bit),
    Int(bit_code_block_index),
    bit_buffer,
    nothing,
    nothing,
)

# Build the context from the state the correlate/fold step has at hand. The
# bit buffer must be the one from *before* this record was folded in, so
# `bit_code_block_index` is the record's own offset inside the open bit rather
# than the next record's. `bit_sync_usable = false` forces the "no bit grid"
# marker for a record whose prompt cannot be summed coherently against the grid
# — a record correlated with a pre-sync replica whose secondary-code wipe-off the
# sync changed under it (see `drop_prompt` in `_apply_correlator_output`).
@inline function CN0UpdateContext(
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer,
    num_code_blocks::Integer,
    bit_sync_usable::Bool = true;
    noise_density = nothing,
    integration_time = nothing,
)
    usable = bit_sync_usable && has_bit_or_secondary_code_been_found(bit_buffer)
    CN0UpdateContext(
        signal,
        Int(num_code_blocks),
        _calc_num_code_blocks_that_form_a_bit(signal),
        usable ? bit_buffer.prompt_accumulator_integrated_code_blocks : -1,
        bit_buffer,
        noise_density,
        integration_time,
    )
end

"""
$(SIGNATURES)

Does `estimator` read its signal's noise density out of its
[`CN0UpdateContext`](@ref)? `false` for every estimator that infers its noise
floor from the prompt stream ([`MomentsCN0Estimator`](@ref),
[`NWPRCN0Estimator`](@ref), [`NoCN0Estimator`](@ref)); `true` only for
[`NoiseRefCN0Estimator`](@ref).

It is a trait rather than a hard-coded type check because
[`update(::AbstractCN0Estimator, ::Any, ::CN0UpdateContext)`](@ref) is a
documented extension point: a third-party estimator that wants a density must be
able to say so, and one that does not must not pay for it.

Two things key off it, and both are compile-time constants on the estimator's
type:

  - **Provisioning.** [`TrackState`](@ref) gives a signal a
    [`CorrelatorNoiseEstimator`](@ref) only where that signal's estimator returns
    `true`. A signal that does not ask gets no entry at all, so its despread
    never runs and costs exactly zero — which is the answer for anyone who
    deliberately stays on NWPR.
  - **The warm-up skip.** While a configured source's window is still empty, the
    fold skips the C/N₀ update for the requiring signal *only* — its co-residents
    on the same band and the same satellite are untouched, because each has its
    own source. An `NWPRCN0Estimator` beside it would otherwise be corrupted: a
    record missing from the bit grid makes `_update_nwpr` drop its open
    narrowband window, so NWPR would silently degrade to its fallback.

The trait's real home is the **type**, and the instance method forwards to it.
That is what lets provisioning be decided from a group's already-fixed slot type
rather than from a satellite value: the whole `noise_estimators` NamedTuple —
its keys included — then folds out of `TrackState`'s type parameters instead of
inferring as a union of "provisioned" and "not". A custom estimator may define
either form; defining the type form is the one that keeps that folding.
"""
requires_noise_density(::Type{<:AbstractCN0Estimator}) = false
requires_noise_density(estimator::AbstractCN0Estimator) =
    requires_noise_density(typeof(estimator))

"""
$(SIGNATURES)

Fold one completed record's `prompt` into `estimator` and return the updated
estimator (immutable update). This three-argument form is the one the tracking
loop calls and the **extension point** for custom
[`AbstractCN0Estimator`](@ref)s that profit from the navigation-bit information
in `context` (see [`CN0UpdateContext`](@ref)).

The default implementation drops `context` and calls the two-argument
`Tracking.update(estimator, prompt)`, so an estimator that only needs the
prompt stream — like [`MomentsCN0Estimator`](@ref) — implements just that one.
"""
update(estimator::AbstractCN0Estimator, prompt, ::CN0UpdateContext) =
    update(estimator, prompt)

"""
$(SIGNATURES)

The default CN0 estimator for `signal`: a [`NoiseRefCN0Estimator`](@ref)
averaging over `num_prompts_for_cn0_estimation` records, against that signal's
own measured noise density.

It reads a density, so [`TrackState`](@ref) provisions the signal a
[`CorrelatorNoiseEstimator`](@ref) automatically (see
[`requires_noise_density`](@ref)) and `track!` fills it from the samples — the
sample-driven path needs no configuration at all.

# Why this and not NWPR

[`NWPRCN0Estimator`](@ref) is accurate where it applies, but it needs a coherent
narrowband window and so does not apply uniformly. Measured through `track!` on
a data-modulated GPS L1 C/A signal, 1200 code blocks, median over 9 seeds, with
the fraction of runs reporting `-Inf dB-Hz` in brackets:

| true | NWPR, 1-block records | NoiseRef    | NWPR, 20-block records | NoiseRef    |
|:---- | ---------------------:| -----------:| ----------------------:| -----------:|
| 25   | 11.2 (56 %)           | 23.4 (0)    | —                      | —           |
| 30   | 29.3 ± 1.40           | 29.6 ± 0.62 | —                      | —           |
| 40   | 39.7 ± 0.35           | 39.9 ± 0.19 | 27.1 ± 5.79            | 39.8 ± 0.37 |
| 45   | 44.7 ± 0.30           | 44.9 ± 0.13 | 32.9 ± 2.00            | 44.7 ± 0.30 |

Three things in that table decided the default. At **25 dB-Hz** NWPR's window
lands outside `1 < μ̂ < M` in over half the runs and the surviving estimates read
14 dB low, which is exactly the regime a lock detector has to work in. At **long
coherent records** NWPR collapses — a record as long as its own window has
`NBP ≡ WBP` and no window exists at all, so it falls back — while the
non-coherent reference is **immune to the phase-noise wash** and barely moves
between 1- and 20-block records. And on GPS L1C-D, Galileo E1B and any
secondary-coded signal before sync, NWPR admits no window ever and defers
permanently to a fallback with a different bias (issue #217).

What NWPR is still better at is the top of the range, where the reference
carries a self-leakage bias it does not: ≈0.13 dB at 45 dB-Hz and ≈0.40 at 50 on
L1 C/A. See [`NoiseRefCN0Estimator`](@ref).

# When to pass something else

`NWPRCN0Estimator` remains exported and is the estimator to configure explicitly
for **externally supplied correlator outputs without a noise observation** — a
correlator-ingest path that cannot also report `Σ|B|²` for an untracked PRN. It
is the one place it is still necessary; everywhere else, prefer appending a
[`NoiseObservation`](@ref) per signal with [`append_noise_observation!`](@ref).
"""
function default_cn0_estimator(
    signal::AbstractGNSSSignal,
    num_prompts_for_cn0_estimation::Int,
)
    NoiseRefCN0Estimator(; num_records = num_prompts_for_cn0_estimation)
end
