"""
$(SIGNATURES)

Van Dierendonck's **narrowband/wideband power ratio** (NWPR) CN0 estimator.

The estimator to configure on a **correlator-ingest path with no noise
observation**; see [`default_cn0_estimator`](@ref) for why it is no longer the
default, and [`NWPRCN0Estimator(::AbstractGNSSSignal)`](@ref) for the constructor
that sizes its window for the signal.

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

The reference spells `μ̂` as the mean of the per-window ratios, but the inversion
is derived from `μ = E[NBP] / E[WBP]`, which only the ratio of the sums estimates
consistently: `E[NBP/WBP] < E[NBP]/E[WBP]` at finite `M`. Measured with `K → ∞`,
so only the bias is left:

| true C/N₀ | mean of ratios, `M = 2` | `M = 5` | `M = 20` | ratio of sums, any `M` |
|:--------- | -----------------------:| -------:| --------:| ----------------------:|
| 20 dB-Hz  | 18.4                    | 19.3    | 19.8     | 20.0                   |
| 25 dB-Hz  | 23.6                    | 24.4    | 24.9     | 25.0                   |
| 30 dB-Hz  | 28.9                    | 29.6    | 29.9     | 30.0                   |

At the default five-block window the two forms agree to a few tenths of a dB in
the loop; at two blocks the ratio of the sums reads 1–2 dB higher. It is
preferred because it is consistent at *every* window length. It shifts the
estimate without changing its spread.

Unlike [`MomentsCN0Estimator`](@ref), NWPR reports "no signal" on pure noise
rather than a ~27.6 dB-Hz floor (issue #217; see "Why not the moment method" in
docs/src/cn0_estimator.md).

# The coherence constraint

`NBP` is a **coherent** sum, so a window must not straddle a navigation-bit
flip (~7 dB loss) and must stay short against the residual Doppler
(`M·T ≪ 1/(2·Δf)`). The tracking loop knows where the bit flips sit, so the
window is taken from the navigation-bit grid in [`CN0UpdateContext`](@ref):

| signal state                                            | narrowband window                                                      |
|:------------------------------------------------------- |:---------------------------------------------------------------------- |
| bit / secondary sync found, data-bearing signal         | `num_narrowband_code_blocks`, tiling the navigation bit from its start |
| bit / secondary sync found, pilot (no data)             | `num_narrowband_code_blocks` (no bit grid to respect)                  |
| sync not found yet, data-bearing without secondary code | `num_presync_narrowband_code_blocks`, unaligned                        |
| sync not found yet, signal with a secondary code        | none — the unknown overlay flips sign every code block                 |
| one symbol per code block (GPS L1C-D, Galileo E1B)      | none — no coherent window longer than one record exists                |
| record at least as long as its own window               | none — a one-record window has `NBP == WBP` by construction            |

The unaligned pre-sync window exists because bit sync takes seconds and fails
below ~30 dB-Hz; it straddles a flip with probability `(M−1)/L`, costing ~0.6 dB
at the defaults (see docs/src/cn0_estimator.md).

# The window is capped by the loop's coherence time, not by the bit period

Windows of `num_narrowband_code_blocks` tile the bit from its start. A longer
window buys little spread (the same `num_records` records are only partitioned
differently), while residual phase noise costs a *bias* that averaging cannot
remove. GPS L1 C/A, 1 ms records, locked at 45 dB-Hz and faded to a true
25 dB-Hz, median / 10th percentile over 96 runs:

| window                    | 2 records   | 5           | 10          | 20 (one bit) |
|:------------------------- | -----------:| -----------:| -----------:| ------------:|
| default FLL-assisted loop | 25.0 / 20.6 | 24.3 / 21.5 | 23.5 / 19.6 | 21.5 / 9.5   |
| plain PLL                 | 24.9 / 20.9 | 24.8 / 22.6 | 24.8 / 22.8 | 24.7 / 22.6  |

The plain 18 Hz PLL (`ConventionalPLLAndDLL`) holds phase over a whole bit; the
default loop's FLL assist does not, so the default cap is short (~5 ms). Raise
it for a plain PLL, a pilot or a signal that is never weak.

Where no window is admissible the `fallback` is reported instead — at the default
[`MomentsCN0Estimator`](@ref), with its ≈27.6 dB-Hz noise floor (issue #217).
That includes records at least as long as the window, e.g. with
[`set_preferred_num_code_blocks_to_integrate!`](@ref) at one navigation bit.

# Fields / configuration

  - `num_records` — how many records the estimate averages over (100 by
    default, ~100 ms at GPS L1 C/A); `num_records ÷ M` windows are averaged, so
    the memory is independent of `M`.
  - `num_narrowband_code_blocks` — window length in primary-code blocks: the cap
    on the coherent sum, and the window length itself for a signal with **no**
    navigation-bit grid (a pilot, or a bare prompt stream fed through the
    two-argument `Tracking.update`).
  - `num_presync_narrowband_code_blocks` — window length in primary-code blocks
    used while the bit grid is still unknown (see the table above); `0`
    disables the pre-sync window and reports the `fallback` until sync.
  - `buffered_narrowband_powers`, `buffered_wideband_powers`,
    `ratio_current_index`, `filled_ratio_length`, `num_records_per_ratio`,
    `ratios_are_bit_aligned` — the ring buffers of completed windows' `NBP` and
    `WBP` (separate, since the estimate divides their sums), the record count
    `M` they were formed with, and whether they followed the bit grid. A change
    of either restarts the ring; see `_update_nwpr`.
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
    ratios_are_bit_aligned::Bool
    narrowband_sum::ComplexF64
    wideband_power::Float64
    num_accumulated_records::Int
    num_accumulated_code_blocks::Int
    fallback::F
end

"""
$(SIGNATURES)

Construct a fresh [`NWPRCN0Estimator`](@ref) averaging over the last
`num_records` records; see there for the window parameters.

`fallback` is the estimator reported until the first window completes; it
defaults to a [`MomentsCN0Estimator`](@ref) over the same `num_records`
prompts.
"""
function NWPRCN0Estimator(;
    num_records::Int = 100,
    num_narrowband_code_blocks::Int = 5,
    num_presync_narrowband_code_blocks::Int = 5,
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
        false,
        complex(0.0, 0.0),
        0.0,
        0,
        0,
        fallback,
    )
end

"""
$(SIGNATURES)

Construct an [`NWPRCN0Estimator`](@ref) whose coherent window is sized for
`signal`: the whole code blocks covering about 5 ms, at least two of them — 5
blocks for a 1 ms code, 2 for GPS L1C-P's 10 ms one.

About 5 ms is what a coherent sum survives with the default FLL-assisted
carrier loop at low C/N₀ (see [`NWPRCN0Estimator`](@ref)); a window counted in
blocks without the code period would be off by the period's ratio. This is the
form to use on a correlator-ingest path with no noise observation:

```julia
TrackedSat(
    GPSL1C_P(),
    prn,
    code_phase,
    doppler;
    cn0_estimator = NWPRCN0Estimator(GPSL1C_P()),
)
```

Every other keyword is forwarded unchanged.
"""
function NWPRCN0Estimator(
    signal::AbstractGNSSSignal;
    num_records::Int = 100,
    num_narrowband_code_blocks::Int = _default_narrowband_code_blocks(signal),
    kwargs...,
)
    NWPRCN0Estimator(; num_records, num_narrowband_code_blocks, kwargs...)
end

# Whole code blocks covering ~5 ms, at least two.
@inline function _default_narrowband_code_blocks(signal::AbstractGNSSSignal)
    code_period = get_code_length(signal) / get_code_frequency(signal)
    max(2, round(Int, 5ms / code_period))
end

length(estimator::NWPRCN0Estimator) = estimator.filled_ratio_length
get_buffered_narrowband_powers(estimator::NWPRCN0Estimator) =
    estimator.buffered_narrowband_powers
get_buffered_wideband_powers(estimator::NWPRCN0Estimator) =
    estimator.buffered_wideband_powers
get_current_index(estimator::NWPRCN0Estimator) = estimator.ratio_current_index
get_fallback_cn0_estimator(estimator::NWPRCN0Estimator) = estimator.fallback

# Ring slots in use at the current window size: enough to cover `num_records`
# records, capped by the buffer. Derived, since it only changes together with
# `num_records_per_ratio` (which restarts the ring).
@inline function _num_ratios(estimator::NWPRCN0Estimator)
    capacity = Base.length(estimator.buffered_narrowband_powers)
    estimator.num_records_per_ratio < 1 && return capacity
    clamp(div(estimator.num_records, estimator.num_records_per_ratio), 1, capacity)
end

# Narrowband window for this record, as `(window_code_blocks, window_block_index)`:
# the window length in primary-code blocks (`0` = no admissible window, report
# the fallback) and the record's offset inside its window when the window has to
# follow the bit grid (`-1` = free-running, any alignment). See the table in
# [`NWPRCN0Estimator`](@ref) for the reasoning behind each case.
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
    # Data-bearing and synced: windows tile the bit from its start. The bit's
    # trailing `num_code_blocks_per_bit % window_code_blocks` blocks are left
    # out — a shorter window has a different `M` and would restart the ring
    # every bit.
    window_code_blocks = min(estimator.num_narrowband_code_blocks, num_code_blocks_per_bit)
    window_start =
        div(context.bit_code_block_index, window_code_blocks) * window_code_blocks
    window_start + window_code_blocks > num_code_blocks_per_bit && return (0, -1)
    (window_code_blocks, context.bit_code_block_index - window_start)
end

"""
$(SIGNATURES)

Accumulate one record's `prompt` into the open narrowband window, taking the
window length and alignment from the navigation-bit grid in `context` (see the
table in [`NWPRCN0Estimator`](@ref)). A grid-following window only starts at a
multiple of the window length inside the bit, and an open window is dropped
whenever it is no longer admissible (no window applies, or the grid moved under
it at sync). A window of one record carries no information (`NBP == WBP`) and is
not buffered.
"""
function update(estimator::NWPRCN0Estimator, prompt, context::CN0UpdateContext)
    window_code_blocks, window_block_index = _narrowband_window(estimator, context)
    _update_nwpr(
        estimator,
        prompt,
        update(estimator.fallback, prompt, context),
        context.num_code_blocks,
        window_code_blocks,
        window_block_index,
    )
end

"""
$(SIGNATURES)

Advance the estimator on a bare prompt stream, with no navigation-bit context:
each prompt counts as a one-block record and windows run back to back at
`num_narrowband_code_blocks`. For folding a captured prompt stream by hand;
`track` calls the three-argument form.
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
# inside its window when the window must follow the bit grid (-1 = free running).
@inline function _update_nwpr(
    estimator::NWPRCN0Estimator,
    prompt,
    fallback::AbstractCN0Estimator,
    num_code_blocks::Int,
    window_code_blocks::Int,
    window_block_index::Int,
)
    # No admissible window: drop whatever was open (its remainder would either
    # straddle the sync transition or was never coherent to begin with).
    window_code_blocks < 1 && return _with_window_state(estimator, fallback)
    # A record spanning no whole code block (the fractional one after a sync
    # phase-snap) has no position on the grid; folding it in would break the
    # inversion's assumption of `M` equal records. Skip it, window untouched.
    num_code_blocks < 1 && return _with_open_window(estimator, fallback)
    # A grid-following window holding `k` blocks must sit `k` blocks in.
    # Otherwise the grid moved under it (sync was just found) and its partial sum
    # goes; only a record that starts a window may open a new one.
    if window_block_index >= 0 &&
       window_block_index != estimator.num_accumulated_code_blocks
        window_block_index == 0 || return _with_window_state(estimator, fallback)
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
    # Window complete. With one record `NBP == WBP`: the records have grown to
    # the window length, a lasting configuration, so the ring is emptied too.
    # Otherwise the estimate would freeze on the old windows (`estimate_cn0`
    # only consults the `fallback` while the ring is empty), reported at the new
    # record's `integration_time`.
    num_records < 2 && return _with_window_state(
        estimator,
        fallback;
        ratio_current_index = 0,
        filled_ratio_length = 0,
        num_records_per_ratio = 0,
    )
    # All-zero window: nothing to buffer (division guard for degenerate input).
    iszero(wideband_power) && return _with_window_state(estimator, fallback)
    narrowband_powers = estimator.buffered_narrowband_powers
    wideband_powers = estimator.buffered_wideband_powers
    narrowband_power = abs2(narrowband_sum)
    bit_aligned = window_block_index >= 0
    if num_records == estimator.num_records_per_ratio &&
       bit_aligned == estimator.ratios_are_bit_aligned
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
            ratios_are_bit_aligned = bit_aligned,
        )
    end
    # First window at this `M`, or the first across the sync transition: re-zero
    # the rings, since `estimate_cn0` sums the whole buffers. Pre-sync windows go
    # even at unchanged `M` — some straddled a bit flip and would drag the
    # estimate down for `num_records` after sync.
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
        ratios_are_bit_aligned = bit_aligned,
    )
end

# Rebuild with a new ring / open-window state, reusing the configuration and the
# in-place power vectors. Defaults: ring unchanged, window closed. Allocation-free.
@inline _with_window_state(
    estimator::NWPRCN0Estimator,
    fallback::AbstractCN0Estimator;
    ratio_current_index::Int = estimator.ratio_current_index,
    filled_ratio_length::Int = estimator.filled_ratio_length,
    num_records_per_ratio::Int = estimator.num_records_per_ratio,
    ratios_are_bit_aligned::Bool = estimator.ratios_are_bit_aligned,
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
    ratios_are_bit_aligned,
    narrowband_sum,
    wideband_power,
    num_accumulated_records,
    num_accumulated_code_blocks,
    fallback,
)

# Advance nothing but the fallback: the open window is carried over untouched,
# where `_with_window_state`'s defaults would close it.
@inline _with_open_window(estimator::NWPRCN0Estimator, fallback::AbstractCN0Estimator) =
    _with_window_state(
        estimator,
        fallback;
        narrowband_sum = estimator.narrowband_sum,
        wideband_power = estimator.wideband_power,
        num_accumulated_records = estimator.num_accumulated_records,
        num_accumulated_code_blocks = estimator.num_accumulated_code_blocks,
    )

"""
$(SIGNATURES)

Estimate the CN0 from the buffered narrowband and wideband powers, dividing by
`integration_time` — the *record's* integration time, which is the predetection
integration time `T` of Van Dierendonck's formula (see
[`NWPRCN0Estimator`](@ref)).

`μ̂` is the ratio of the sums (see [`NWPRCN0Estimator`](@ref)). Until the first
window completes, and for good on a signal that admits none, the `fallback`'s
value is returned.

`μ̂ ≤ 1` (**no detectable signal**) yields `-Inf dB-Hz` and `μ̂ ≥ M` (no
detectable noise) yields `Inf dB-Hz`, the limits of the expression rather than a
clamped value. Thresholding needs no special case; averaging the estimate does.
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
