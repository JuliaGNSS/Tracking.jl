"""
$(SIGNATURES)

Abstract supertype for CN0 (carrier-to-noise-density ratio) estimators.
Each [`TrackedSignal`](@ref) holds one estimator instance, stored in a type
parameter — pass any subtype instance as the `cn0_estimator` keyword of
[`TrackedSignal`](@ref) / [`TrackedSat`](@ref) to replace the default
[`NWPRCN0Estimator`](@ref). Custom estimators subtype this and implement
`Tracking.update` and [`estimate_cn0`](@ref).
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
  - `noise_density` — the band's measured noise density `N₀` (dimension `1/Hz`),
    for an estimator that divides the prompt's power by a *measured* floor
    instead of inferring one from the prompt stream
    ([`NoiseRefCN0Estimator`](@ref) does). It is a **type parameter, not a
    sentinel**: `Nothing` says no [`AbstractNoiseEstimator`](@ref) is configured
    on the band at all, which is a static property of the setup, so every call
    site monomorphises and the estimators stay allocation-free. A configured
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

Does `estimator` read the band's noise density out of its
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

  - **Provisioning.** [`TrackState`](@ref) gives a band a
    [`CorrelatorNoiseEstimator`](@ref) only where some signal on it returns
    `true`. A band with no requiring signal gets no entry at all, so the
    per-band despread never runs and costs exactly zero — which is the answer
    for anyone who deliberately stays on NWPR.
  - **The warm-up skip.** While a configured source's window is still empty, the
    fold skips the C/N₀ update for the requiring signals *only*. Skipping the
    whole band would corrupt a co-resident `NWPRCN0Estimator`: a record missing
    from the bit grid makes `_update_nwpr` drop its open narrowband window, so
    NWPR would silently degrade to its fallback.
"""
requires_noise_density(::AbstractCN0Estimator) = false

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

The default CN0 estimator for `signal`: an [`NWPRCN0Estimator`](@ref) averaging
over `num_prompts_for_cn0_estimation` records, with a
`MomentsCN0Estimator(num_prompts_for_cn0_estimation)` as its fallback.

The narrowband window is the whole code blocks covering about 5 ms, at least two
of them: 5 blocks for a 1 ms code, 2 for GPS L1C-P's 10 ms code. It caps the
coherent sum for a data-bearing signal (whose window tiles the navigation bit)
and is the window outright for a pilot. About 5 ms is what a coherent sum
survives with the default 18 Hz carrier loop at the low C/N₀ this estimator
exists for; see [`NWPRCN0Estimator`](@ref) for the measurements and for when to
raise it.
"""
function default_cn0_estimator(
    signal::AbstractGNSSSignal,
    num_prompts_for_cn0_estimation::Int,
)
    code_period = get_code_length(signal) / get_code_frequency(signal)
    num_narrowband_code_blocks = max(2, round(Int, 5ms / code_period))
    NWPRCN0Estimator(;
        num_records = num_prompts_for_cn0_estimation,
        num_narrowband_code_blocks,
        fallback = MomentsCN0Estimator(num_prompts_for_cn0_estimation),
    )
end
