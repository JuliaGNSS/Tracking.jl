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
what the tracking loop knows about the record that an estimator cannot recover
from the prompt stream alone — chiefly the **navigation-bit grid**, which lets an
estimator sum prompts coherently without straddling a bit flip
([`NWPRCN0Estimator`](@ref) does).

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
  - `bit_code_block_index` — this record's offset inside the open navigation
    bit (`0` = starts a fresh bit). `-1` when the prompt cannot be summed
    coherently against the grid: before bit / secondary-code sync, and, for a
    secondary-coded signal, for a record correlated without the overlay wipe-off
    a sync in the same fold just established (see `_apply_correlator_output`).
    Such a record's blocks still count towards the grid.
  - `bit_buffer` — the signal's `BitBuffer` as of *before* this record:
    the soft bits decoded so far, the open coherent bit accumulator, the lock
    polarity and the secondary-code phase.
  - `noise_density` — **this signal's** measured noise density `N₀` (dimension
    `1/Hz`), read by [`NoiseRefCN0Estimator`](@ref). `Nothing` (in the type, not
    as a sentinel) means no [`AbstractNoiseEstimator`](@ref) is configured for
    the signal, so call sites monomorphise. A configured source whose window is
    still empty never reaches here — the fold skips the update — so otherwise
    this is always a plain scalar.
  - `integration_time` — this record's own `T`, so records of different lengths
    in one ring are each handled correctly. `nothing` on the bare-prompt-stream
    path.
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

# Positional convenience for the first five fields, defaulting the two only a
# noise-referenced estimator reads.
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

# Build the context from the correlate/fold step's state. `bit_buffer` must be the
# one from *before* this record was folded in, so `bit_code_block_index` is this
# record's offset. `bit_sync_usable = false` forces the "no bit grid" marker (see
# `drop_prompt` in `_apply_correlator_output`).
#
# The keyword form is for callers building a context by hand; the per-record fold
# calls the positional core below, because a keyword call allocates per call on
# Julia 1.10 (still supported) and this is the hottest construction site.
@inline CN0UpdateContext(
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer,
    num_code_blocks::Integer,
    bit_sync_usable::Bool = true;
    noise_density = nothing,
    integration_time = nothing,
) = CN0UpdateContext(
    signal,
    bit_buffer,
    num_code_blocks,
    bit_sync_usable,
    noise_density,
    integration_time,
)

@inline function CN0UpdateContext(
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer,
    num_code_blocks::Integer,
    bit_sync_usable::Bool,
    noise_density,
    integration_time,
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

A trait, so a custom estimator can opt in. Two things key off it:

  - **Provisioning.** [`TrackState`](@ref) gives a signal a
    [`CorrelatorNoiseEstimator`](@ref) only where its estimator returns `true`;
    any other signal runs no despread at all.
  - **The warm-up skip.** While a configured source's window is still empty, the
    fold skips the C/N₀ update for the requiring signal *only*. Skipping its
    neighbours too would make an `NWPRCN0Estimator` beside it drop its open
    window on every skipped record (see `_update_nwpr`).

The trait's home is the **type** and the instance method forwards to it, so
provisioning folds out of `TrackState`'s type parameters. A custom estimator
should define the type form to keep that.
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

[`NWPRCN0Estimator`](@ref) needs a coherent narrowband window, so it does not
apply uniformly. Measured through `track!` on a data-modulated GPS L1 C/A signal,
1200 code blocks, median over 9 seeds, fraction of runs reporting `-Inf dB-Hz`
in brackets:

| true | NWPR, 1-block records | NoiseRef    | NWPR, 20-block records | NoiseRef    |
|:---- | ---------------------:| -----------:| ----------------------:| -----------:|
| 25   | 11.2 (56 %)           | 23.4 (0)    | —                      | —           |
| 30   | 29.3 ± 1.40           | 29.6 ± 0.62 | —                      | —           |
| 40   | 39.7 ± 0.35           | 39.9 ± 0.19 | 27.1 ± 5.79            | 39.8 ± 0.37 |
| 45   | 44.7 ± 0.30           | 44.9 ± 0.13 | 32.9 ± 2.00            | 44.7 ± 0.30 |

At **25 dB-Hz**, where lock decisions are made, NWPR's ratio leaves
`1 < μ̂ < M` in over half the runs. At **long coherent records** a record as long
as its window has `NBP ≡ WBP` and NWPR falls back, while the non-coherent
reference is immune to phase noise. On GPS L1C-D, Galileo E1B and secondary-coded
signals before sync NWPR never gets a window at all (issue #217). NWPR is better
only at the top of the range, where the reference carries a self-leakage bias
(see [`NoiseRefCN0Estimator`](@ref)).

# When to pass something else

Configure `NWPRCN0Estimator` explicitly for **externally supplied correlator
outputs without a noise observation**. Everywhere else, prefer appending a
[`NoiseObservation`](@ref) per signal with [`append_noise_observation!`](@ref).
"""
function default_cn0_estimator(
    signal::AbstractGNSSSignal,
    num_prompts_for_cn0_estimation::Int,
)
    NoiseRefCN0Estimator(; num_records = num_prompts_for_cn0_estimation)
end
