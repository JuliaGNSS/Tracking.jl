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
"""
struct CN0UpdateContext{S<:AbstractGNSSSignal,B<:Unsigned}
    signal::S
    num_code_blocks::Int
    num_code_blocks_per_bit::Int
    bit_code_block_index::Int
    bit_buffer::BitBuffer{B}
end

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
    bit_sync_usable::Bool = true,
)
    usable = bit_sync_usable && has_bit_or_secondary_code_been_found(bit_buffer)
    CN0UpdateContext(
        signal,
        Int(num_code_blocks),
        _calc_num_code_blocks_that_form_a_bit(signal),
        usable ? bit_buffer.prompt_accumulator_integrated_code_blocks : -1,
        bit_buffer,
    )
end

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

MomentsCN0Estimator to estimate the CN0
"""
struct MomentsCN0Estimator <: AbstractCN0Estimator
    prompt_buffer::Vector{ComplexF64}
    buffer_current_index::Int
    filled_buffer_length::Int
end

function MomentsCN0Estimator(N)
    MomentsCN0Estimator(zeros(ComplexF64, N), 0, 0)
end

length(estimator::MomentsCN0Estimator) = estimator.filled_buffer_length
get_prompt_buffer(estimator::MomentsCN0Estimator) = estimator.prompt_buffer
get_current_index(estimator::MomentsCN0Estimator) = estimator.buffer_current_index

"""
$(SIGNATURES)

Buffers the prompts such that they can be used to estimate the CN0.
Returns a new estimator with the latest prompt added (immutable update).
The moment ratio needs no side information, so this estimator implements only
the prompt-only form of `Tracking.update`; see
[`update(::AbstractCN0Estimator, ::Any, ::CN0UpdateContext)`](@ref) for the
extension point that also receives the navigation-bit context.
"""
function update(estimator::MomentsCN0Estimator, prompt)
    buffer_length = length(estimator.prompt_buffer)
    next_index = mod(get_current_index(estimator), buffer_length) + 1
    estimator.prompt_buffer[next_index] = prompt
    filled_buffer_length = min(length(estimator) + 1, buffer_length)
    MomentsCN0Estimator(estimator.prompt_buffer, next_index, filled_buffer_length)
end

"""
$(SIGNATURES)

Estimates the CN0 based on the struct `MomentsCN0Estimator`.
"""
function estimate_cn0(estimator::MomentsCN0Estimator, integration_time)
    length(estimator) == 0 && return 0.0dBHz
    prompt_buffer = get_prompt_buffer(estimator)
    M₂ = 1 / length(estimator) * sum(abs2, prompt_buffer)
    M₄ = 1 / length(estimator) * sum(x -> abs2(x)^2, prompt_buffer)
    Pd = sqrt(abs(2 * M₂^2 - M₄))
    noise_power = abs(M₂ - Pd)
    SNR = Pd / noise_power
    dBHz(SNR / integration_time)
end

"""
$(SIGNATURES)

Van Dierendonck's **narrowband/wideband power ratio** (NWPR) CN0 estimator, the
default estimator of every [`TrackedSignal`](@ref).

Over a narrowband window of `M` consecutive records it forms

```
NBP = |Σ_M prompt|²          (coherent — narrowband)
WBP =  Σ_M |prompt|²         (incoherent — wideband)
```

combines the windows that fit in `num_records` records into a mean power ratio
`μ̂`, and reports

```
μ̂    = Σ_K NBP_k / Σ_K WBP_k
Ĉ/N₀ = (1 / T) · (μ̂ − 1) / (M − μ̂)
```

with `T` the record's own integration time (the `integration_time` argument of
[`estimate_cn0`](@ref)). Reference: A. J. Van Dierendonck, "GPS Receivers",
ch. 8 in *Global Positioning System: Theory and Applications*, Vol. I, ed.
B. W. Parkinson & J. J. Spilker Jr.; the formulas above are reproduced on
[ESA Navipedia](https://gssc.esa.int/navipedia/index.php/Lock_Detectors).

# Why the ratio of the sums, and not the mean of the ratios

The reference spells the mean ratio as `μ̂ = (1/K) Σ NBP_k / WBP_k`, the mean of
the per-window ratios. The inversion `(μ̂ − 1) / (M − μ̂)` is derived from
`μ = E[NBP] / E[WBP]`, which the ratio of the sums estimates consistently and
the mean of the ratios does not: `E[NBP/WBP] < E[NBP]/E[WBP]` at finite `M`, so
the mean of the ratios reads low. Measured on synthetic prompts with `K → ∞`, so
that only the bias is left:

| true C/N₀ | mean of ratios, `M = 2` | `M = 5` | `M = 20` | ratio of sums, any `M` |
|:--------- | -----------------------:| -------:| --------:| ----------------------:|
| 20 dB-Hz  | 18.4                    | 19.3    | 19.8     | 20.0                   |
| 25 dB-Hz  | 23.6                    | 24.4    | 24.9     | 25.0                   |
| 30 dB-Hz  | 28.9                    | 29.6    | 29.9     | 30.0                   |

At the `M = 20` of the classic GPS L1 C/A configuration the two agree to
~0.15 dB, which is why the literature's form is good enough there; the
difference only matters at the short windows this estimator uses while the
navigation-bit grid is unknown or the loop's coherence time is short. The
spread is unchanged (the correction is a shift, not a variance trade), and at
equal false-alarm rate the ratio of the sums is the slightly better *detector*
too — at `M = 2` and a 0.1 % false-alarm rate it detects a true 25 dB-Hz signal
in 22.5 % of updates against the mean of the ratios' 15.0 %.

# Why it is the default

NWPR is usable as a *detection* statistic near threshold, which
[`MomentsCN0Estimator`](@ref) is not: the moment ratio's sample moments
manufacture signal power out of noise at finite window length, so at its
default 100-prompt window M2M4 reports a median of ~27.6 dB-Hz on **pure
noise** and cannot separate a true 20 dB-Hz signal from noise at all (a code
lock threshold below ~30 dB-Hz can never trip, and one at 30 dB-Hz has a ~19 %
per-update false-alarm rate). NWPR on the identical prompt streams tracks the
truth from 20 dB-Hz up and reports "no signal" on noise, with a visibly tighter
spread below 32 dB-Hz (Falletti, Pini & Lo Presti, *IEEE T-AES* 47(1):420–437,
2011). See JuliaGNSS/Tracking.jl#217 for the measurements.

# The coherence constraint, and why this belongs in Tracking

`NBP` is a **coherent** sum, so a window must not straddle a navigation-bit
flip (a mid-window flip costs ~7 dB) and must stay short against the residual
Doppler (`M·T ≪ 1/(2·Δf)`; 25 Hz over a 20 ms window is half a cycle and costs
~9 dB). Where the bit flips sit is something the tracking loop knows and a
consumer of `get_filtered_prompts` does not, so the window is taken
from the navigation-bit grid in [`CN0UpdateContext`](@ref) — that is the whole
reason this estimator lives here rather than on top of the prompt stream:

| signal state                                            | narrowband window                                       |
|:------------------------------------------------------- |:------------------------------------------------------- |
| bit / secondary sync found, data-bearing signal         | one navigation bit, started on a bit boundary           |
| bit / secondary sync found, pilot (no data)             | `num_narrowband_code_blocks` (no bit grid to respect)   |
| sync not found yet, data-bearing without secondary code | `num_presync_narrowband_code_blocks`, unaligned         |
| sync not found yet, signal with a secondary code        | none — the unknown overlay flips sign every code block  |
| one symbol per code block (GPS L1C-D, Galileo E1B)      | none — no coherent window longer than one record exists |

The pre-sync window matters more than it may look: the CFAR bit-edge detector
needs seconds to lock at 35 dB-Hz and does not lock at all below ~30 dB-Hz, so
a bit-aligned window alone would leave exactly the regime this estimator exists
for on the fallback estimator. An unaligned window of `M` records inside an
`L`-block bit straddles a flip with probability `(M−1)/L`, which at the default
two blocks of a 20-block GPS L1 C/A bit costs ~1 dB of bias and still removes
the moment ratio's noise floor entirely (0 % false alarm at a 25 dB-Hz
threshold on pure noise, against M2M4's ~80 %).

Where no window is admissible at all the `fallback` estimator's value is
reported instead, so such signals keep exactly the behaviour they had before
NWPR became the default.

# Fields / configuration

  - `num_records` — how many records the estimate averages over, i.e. the
    memory of the estimator (100 by default, ~100 ms at GPS L1 C/A). The
    ring buffer of ratios is sized from it: `num_records ÷ M` ratios are
    averaged, so the memory stays put when `M` changes.
  - `num_narrowband_code_blocks` — window length in primary-code blocks for a
    signal with **no** navigation-bit grid: a pilot, or a bare prompt stream
    fed through the two-argument `Tracking.update`.
  - `num_presync_narrowband_code_blocks` — window length in primary-code blocks
    used while the bit grid is still unknown (see the table above); `0`
    disables the pre-sync window and reports the `fallback` until sync.
  - `buffered_narrowband_powers`, `buffered_wideband_powers`,
    `ratio_current_index`, `filled_ratio_length`, `num_records_per_ratio` — the
    ring buffers of completed windows' `NBP` and `WBP` and the record count `M`
    they were formed with. The two powers are buffered separately because the
    estimate divides their sums; `M` enters the estimate, so a window completing
    with a different record count (records lengthened at bit sync, say) restarts
    the buffers.
  - `narrowband_sum`, `wideband_power`, `num_accumulated_records`,
    `num_accumulated_code_blocks` — the currently open window.
  - `fallback` — the estimator reported while no window has completed yet, and
    for the signals of the table above that never get one. Defaults to a
    [`MomentsCN0Estimator`](@ref) and is fed every prompt.
"""
struct NWPRCN0Estimator{F<:AbstractCN0Estimator} <: AbstractCN0Estimator
    num_records::Int
    num_narrowband_code_blocks::Int
    num_presync_narrowband_code_blocks::Int
    buffered_narrowband_powers::Vector{Float64}
    buffered_wideband_powers::Vector{Float64}
    ratio_current_index::Int
    filled_ratio_length::Int
    num_records_per_ratio::Int
    narrowband_sum::ComplexF64
    wideband_power::Float64
    num_accumulated_records::Int
    num_accumulated_code_blocks::Int
    fallback::F
end

"""
$(SIGNATURES)

Construct a fresh [`NWPRCN0Estimator`](@ref) averaging over the last
`num_records` records. `num_narrowband_code_blocks` /
`num_presync_narrowband_code_blocks` are the narrowband window lengths for the
signal states that have no navigation-bit period to take them from — see
[`NWPRCN0Estimator`](@ref) for the full table and for what the defaults cost.

`fallback` is the estimator reported until the first window completes; it
defaults to a [`MomentsCN0Estimator`](@ref) over the same `num_records`
prompts.
"""
function NWPRCN0Estimator(;
    num_records::Int = 100,
    num_narrowband_code_blocks::Int = 20,
    num_presync_narrowband_code_blocks::Int = 2,
    fallback::AbstractCN0Estimator = MomentsCN0Estimator(num_records),
)
    num_records >= 2 ||
        throw(ArgumentError("num_records must be at least 2, got $num_records"))
    num_narrowband_code_blocks >= 1 || throw(
        ArgumentError(
            "num_narrowband_code_blocks must be at least 1, got " *
            "$num_narrowband_code_blocks",
        ),
    )
    num_presync_narrowband_code_blocks >= 0 || throw(
        ArgumentError(
            "num_presync_narrowband_code_blocks must not be negative, got " *
            "$num_presync_narrowband_code_blocks",
        ),
    )
    # A window holds at least two records (one carries no information), so at
    # most `num_records ÷ 2` windows can ever be combined.
    NWPRCN0Estimator(
        num_records,
        num_narrowband_code_blocks,
        num_presync_narrowband_code_blocks,
        zeros(Float64, div(num_records, 2)),
        zeros(Float64, div(num_records, 2)),
        0,
        0,
        0,
        complex(0.0, 0.0),
        0.0,
        0,
        0,
        fallback,
    )
end

length(estimator::NWPRCN0Estimator) = estimator.filled_ratio_length
get_buffered_narrowband_powers(estimator::NWPRCN0Estimator) =
    estimator.buffered_narrowband_powers
get_buffered_wideband_powers(estimator::NWPRCN0Estimator) =
    estimator.buffered_wideband_powers
get_current_index(estimator::NWPRCN0Estimator) = estimator.ratio_current_index
get_fallback_cn0_estimator(estimator::NWPRCN0Estimator) = estimator.fallback

# How many of the ring buffer's slots are in use at the current window size:
# enough to cover `num_records` records, capped by the buffer. Recomputed from
# the fields rather than stored, since it only changes together with
# `num_records_per_ratio` (which restarts the ring anyway).
@inline function _num_ratios(estimator::NWPRCN0Estimator)
    capacity = Base.length(estimator.buffered_narrowband_powers)
    estimator.num_records_per_ratio < 1 && return capacity
    clamp(div(estimator.num_records, estimator.num_records_per_ratio), 1, capacity)
end

# Narrowband window for this record, as `(window_code_blocks, bit_code_block_index)`:
# the window length in primary-code blocks (`0` = no admissible window, report
# the fallback) and the record's offset inside the navigation bit when the window
# has to follow the bit grid (`-1` = free-running, any alignment). See the table
# in [`NWPRCN0Estimator`](@ref) for the reasoning behind each case.
@inline function _narrowband_window(estimator::NWPRCN0Estimator, context::CN0UpdateContext)
    num_code_blocks_per_bit = context.num_code_blocks_per_bit
    if context.bit_code_block_index < 0
        # Bit grid unknown. A signal carrying a secondary code has no coherence
        # to exploit at all before sync — the unknown overlay flips sign from
        # code block to code block — while a purely data-modulated signal stays
        # coherent inside a bit, so a window well short of the bit period only
        # occasionally straddles a flip.
        num_code_blocks_per_bit > 1 && get_secondary_code_length(context.signal) == 1 ||
            return (0, -1)
        return (
            min(estimator.num_presync_narrowband_code_blocks, num_code_blocks_per_bit),
            -1,
        )
    end
    # One symbol per code block: every record may flip sign, nothing to sum.
    num_code_blocks_per_bit == 1 && return (0, -1)
    # Pilot: no data modulation, and post-sync the secondary code is wiped off
    # in the replica, so the window is free-running at the configured length.
    num_code_blocks_per_bit == 0 && return (estimator.num_narrowband_code_blocks, -1)
    # Data-bearing and synced: one bit per window, aligned to the bit boundary.
    (num_code_blocks_per_bit, context.bit_code_block_index)
end

"""
$(SIGNATURES)

Accumulate one record's `prompt` into the open narrowband window, taking the
window length and its alignment from the navigation-bit grid in `context` (see
[`NWPRCN0Estimator`](@ref) for the per-signal-state table). A window that
follows the bit grid is only started on a bit boundary, and an open window is
dropped rather than closed whenever it is no longer admissible — no window at
all applies to this record (sync lost, or never established), or the grid moved
under it when sync was found.

A closed window contributes a ratio only if it held at least two records; with
one record `NBP == WBP` identically and the ratio carries no information.
"""
function update(estimator::NWPRCN0Estimator, prompt, context::CN0UpdateContext)
    window_code_blocks, alignment_offset = _narrowband_window(estimator, context)
    _update_nwpr(
        estimator,
        prompt,
        update(estimator.fallback, prompt, context),
        context.num_code_blocks,
        window_code_blocks,
        alignment_offset,
    )
end

"""
$(SIGNATURES)

Advance the estimator on a bare prompt stream, with no navigation-bit context:
each prompt counts as a one-block record and windows run back to back at
`num_narrowband_code_blocks`, aligned to the first prompt. This is the form to
use when folding a captured prompt stream by hand; inside `track` the
three-argument form above is called, which aligns the window to the
navigation-bit grid instead.
"""
update(estimator::NWPRCN0Estimator, prompt) = _update_nwpr(
    estimator,
    prompt,
    update(estimator.fallback, prompt),
    1,
    estimator.num_narrowband_code_blocks,
    -1,
)

# Shared NWPR accumulation core. `fallback` is the already-advanced fallback
# estimator (the caller picks the arity to advance it with), and the trailing
# arguments describe the record and its window: the record's block count, the
# window length in blocks (0 = no admissible window) and the record's offset
# inside the navigation bit when the window must follow the bit grid (-1 = free
# running).
@inline function _update_nwpr(
    estimator::NWPRCN0Estimator,
    prompt,
    fallback::AbstractCN0Estimator,
    num_code_blocks::Int,
    window_code_blocks::Int,
    bit_code_block_index::Int,
)
    # No admissible window: drop whatever was open (its remainder would either
    # straddle the sync transition or was never coherent to begin with).
    window_code_blocks < 1 && return _with_window_state(estimator, fallback)
    # A bit-aligned window that has accumulated `k` blocks must sit exactly `k`
    # blocks into the bit. Anything else means the grid moved under the open
    # window — sync was just found, so the window that was opened at the free
    # pre-sync length is now inside a bit at an unknown offset — and its partial
    # sum has to go. Only a record that starts a bit may then open a new window,
    # so the coherent sum can never straddle a data-bit transition.
    if bit_code_block_index >= 0 &&
       bit_code_block_index != estimator.num_accumulated_code_blocks
        bit_code_block_index == 0 || return _with_window_state(estimator, fallback)
        estimator = _with_window_state(estimator, fallback)
    end
    narrowband_sum = estimator.narrowband_sum + prompt
    wideband_power = estimator.wideband_power + abs2(prompt)
    num_records = estimator.num_accumulated_records + 1
    num_code_blocks_accumulated = estimator.num_accumulated_code_blocks + num_code_blocks
    num_code_blocks_accumulated < window_code_blocks && return _with_window_state(
        estimator,
        fallback;
        narrowband_sum,
        wideband_power,
        num_accumulated_records = num_records,
        num_accumulated_code_blocks = num_code_blocks_accumulated,
    )
    # Window complete. `NBP / WBP` needs at least two records to say anything
    # (with one, `NBP == WBP` by construction), and `M` enters the estimate, so
    # a changed record count invalidates the powers buffered at the old one.
    if num_records < 2 || iszero(wideband_power)
        return _with_window_state(estimator, fallback)
    end
    narrowband_powers = estimator.buffered_narrowband_powers
    wideband_powers = estimator.buffered_wideband_powers
    narrowband_power = abs2(narrowband_sum)
    if num_records == estimator.num_records_per_ratio
        num_ratios = _num_ratios(estimator)
        ratio_current_index = mod(estimator.ratio_current_index, num_ratios) + 1
        narrowband_powers[ratio_current_index] = narrowband_power
        wideband_powers[ratio_current_index] = wideband_power
        return _with_window_state(
            estimator,
            fallback;
            ratio_current_index,
            filled_ratio_length = min(estimator.filled_ratio_length + 1, num_ratios),
            num_records_per_ratio = num_records,
        )
    end
    # First window at this record count: the rings' slots have to be re-zeroed,
    # since `estimate_cn0` sums the whole buffers (and their in-use length
    # shrinks as `M` grows).
    fill!(narrowband_powers, 0.0)
    fill!(wideband_powers, 0.0)
    narrowband_powers[1] = narrowband_power
    wideband_powers[1] = wideband_power
    _with_window_state(
        estimator,
        fallback;
        ratio_current_index = 1,
        filled_ratio_length = 1,
        num_records_per_ratio = num_records,
    )
end

# Rebuild an NWPR estimator with a new ring-buffer / open-window state, reusing
# the configuration and the (in-place updated) power vectors. The defaults are
# "ring unchanged, window closed", which is what most branches of
# `_update_nwpr` want. Immutable rebuild of a bits-and-pointers struct, so this
# allocates nothing.
@inline _with_window_state(
    estimator::NWPRCN0Estimator,
    fallback::AbstractCN0Estimator;
    ratio_current_index::Int = estimator.ratio_current_index,
    filled_ratio_length::Int = estimator.filled_ratio_length,
    num_records_per_ratio::Int = estimator.num_records_per_ratio,
    narrowband_sum::ComplexF64 = complex(0.0, 0.0),
    wideband_power::Float64 = 0.0,
    num_accumulated_records::Int = 0,
    num_accumulated_code_blocks::Int = 0,
) = NWPRCN0Estimator(
    estimator.num_records,
    estimator.num_narrowband_code_blocks,
    estimator.num_presync_narrowband_code_blocks,
    estimator.buffered_narrowband_powers,
    estimator.buffered_wideband_powers,
    ratio_current_index,
    filled_ratio_length,
    num_records_per_ratio,
    narrowband_sum,
    wideband_power,
    num_accumulated_records,
    num_accumulated_code_blocks,
    fallback,
)

"""
$(SIGNATURES)

Estimate the CN0 from the buffered narrowband and wideband powers, dividing by
`integration_time` — the *record's* integration time, which is the predetection
integration time `T` of Van Dierendonck's formula (see
[`NWPRCN0Estimator`](@ref)).

The mean power ratio is `μ̂ = Σ NBP_k / Σ WBP_k` — the ratio of the sums, not the
mean of the per-window ratios, which is biased low at short windows; see
[`NWPRCN0Estimator`](@ref).

Until the first narrowband window has completed — and for good on a signal that
admits no window at all — the `fallback` estimator's value is returned instead.

The mean ratio `μ̂` lies in `[1, M]` by construction and both ends are reported
as the limits of the expression rather than clamped to a magic number: `μ̂ ≤ 1`
means the coherent sum holds no more power than the incoherent one, i.e. **no
detectable signal**, and yields `-Inf dB-Hz`; `μ̂ ≥ M` means no detectable
noise and yields `Inf dB-Hz` (the same value a noise-free signal gets out of
[`MomentsCN0Estimator`](@ref)). A consumer thresholding the estimate needs no
special case for either; one *averaging* it does, which is why the infinities
are documented rather than hidden behind an epsilon.
"""
function estimate_cn0(estimator::NWPRCN0Estimator, integration_time)
    length(estimator) == 0 && return estimate_cn0(estimator.fallback, integration_time)
    num_records = estimator.num_records_per_ratio
    total_wideband_power = sum(get_buffered_wideband_powers(estimator))
    iszero(total_wideband_power) &&
        return estimate_cn0(estimator.fallback, integration_time)
    mean_ratio = sum(get_buffered_narrowband_powers(estimator)) / total_wideband_power
    mean_ratio <= 1 && return dBHz(0.0 / integration_time)
    mean_ratio >= num_records && return dBHz(Inf / integration_time)
    SNR = (mean_ratio - 1) / (num_records - mean_ratio)
    dBHz(SNR / integration_time)
end

"""
$(SIGNATURES)

The default CN0 estimator for `signal`: an [`NWPRCN0Estimator`](@ref) averaging
over `num_prompts_for_cn0_estimation` records, with a
`MomentsCN0Estimator(num_prompts_for_cn0_estimation)` as its fallback.

The free-running narrowband window only has to be chosen here for signals that
have no navigation-bit period to take it from (pilots — a data-bearing signal's
window is one bit, see [`CN0UpdateContext`](@ref)). It is set to the whole code
blocks covering about 20 ms, the GPS L1 C/A bit period, which is a reasonable
proxy for how long a coherent sum survives typical residual Doppler: 20 blocks
for a 1 ms code, 2 for GPS L1C-P's 10 ms code.
"""
function default_cn0_estimator(
    signal::AbstractGNSSSignal,
    num_prompts_for_cn0_estimation::Int,
)
    code_period = get_code_length(signal) / get_code_frequency(signal)
    num_narrowband_code_blocks = max(1, round(Int, 20ms / code_period))
    NWPRCN0Estimator(;
        num_records = num_prompts_for_cn0_estimation,
        num_narrowband_code_blocks,
        fallback = MomentsCN0Estimator(num_prompts_for_cn0_estimation),
    )
end
