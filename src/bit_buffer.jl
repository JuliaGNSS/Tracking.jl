"""
    SyncResult

Outcome of a per-signal bit-sync / secondary-code-sync detector call.

Fields:

  - `found::Bool` — whether the detector locked on this update.
  - `phase::Int` — when `found = true`, the secondary-code chip the
    *upcoming* integration aligns to, in `0:secondary_code_length-1`
    (recovered by the hard [`_secondary_code_search`](@ref)). Always zero for the
    soft CFAR detectors and signals without a secondary code, which fire on a
    bit / period boundary.
  - `polarity::Int8` — `+1` or `-1`; which match orientation the detector
    locked. Applied to the post-sync prompt accumulator so a negative-polarity
    lock doesn't invert the decoded bits.
"""
struct SyncResult
    found::Bool
    phase::Int
    polarity::Int8
end

"""
$(SIGNATURES)

Standard-normal quantile (inverse CDF) `Φ⁻¹(probability)` for
`probability ∈ (0, 1)`, as `√2 · erfinv(2·probability − 1)` (`erfinv` from
SpecialFunctions.jl). Used by [`_t_quantile`](@ref) as the `dof → ∞` anchor
that its tail expansion corrects. Returns `±Inf` at `probability = 1` / `0`;
callers keep the argument in the open interval.
"""
@inline _norm_quantile(probability::Float64) = sqrt(2.0) * erfinv(2 * probability - 1)

"""
$(SIGNATURES)

`probability`-quantile of a Student-t distribution with `dof` degrees of
freedom, i.e. `t` such that `P(T ≤ t) = probability`. [`_cfar_decide`](@ref)
uses it as a **small-sample penalty** on its threshold, not because its z-score
is exactly Student-t: the per-bin values are χ² energies and the competing
hypotheses correlated, so the nominal `dof = peak_bin_count − 1` is a heuristic
and the realised false-alarm rate only approximately the nominal one. The shape
is what matters: steep as `dof → 1`, relaxing to the normal quantile as
`dof → ∞`, so mature many-bin locks are unaffected.

# Implementation

Hill's algorithm (Hill, G. W. (1970), *Algorithm 396: Student's t-quantiles*,
Comm. ACM 13(10), 619–620) on the two-tailed probability `2·(1 − probability)`:
`dof = 1` and `dof = 2` are inverted exactly in closed form, above that a series
in either the normal deviate [`_norm_quantile`](@ref) (near-normal branch) or
the two-tailed probability (deep-tail branch). The lower tail is reflected by
symmetry; `t(0.5) = 0` is returned directly.

Accuracy is ≲2e-4 relative (worst at `dof ≈ 3–10`), far below the nominal-d.o.f.
modelling error above. It is used instead of an exact
`SpecialFunctions.beta_inc_inv` inversion because that allocates over part of
the `dof` range — this runs once per block for every unsynced satellite, and
`test/track_in_place.jl`'s pre-sync guard requires zero allocations — and is
10–170× slower.
"""
function _t_quantile(probability::Float64, dof::Real)
    probability == 0.5 && return 0.0
    probability < 0.5 && return -_t_quantile(1 - probability, dof)
    two_tailed_probability = 2 * (1 - probability)
    # `dof = 1` is the standard Cauchy, `quantile(p) = cot(π·two_tailed/2)`
    # (the cotangent form keeps its accuracy in the far tail, where
    # `tan(π(probability − ½))` cancels). At `dof = 2` the CDF
    # `½·(1 + t/√(2 + t²))` inverts in closed form. The series below is weakest
    # exactly there.
    dof == 1 &&
        return cos(two_tailed_probability * pi / 2) / sin(two_tailed_probability * pi / 2)
    dof == 2 && return sqrt(2 / (two_tailed_probability * (2 - two_tailed_probability)) - 2)
    a = 1 / (dof - 0.5)
    b = 48 / a^2
    c = ((20700 * a / b - 98) * a - 16) * a + 96.36
    d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * pi / 2) * dof
    y = (d * two_tailed_probability)^(2 / dof)
    if y > 0.05 + a
        # Near-normal branch: correct the normal deviate (the `dof → ∞` limit)
        # by an asymptotic series in `1/dof`. Mature large-dof locks land here.
        x = _norm_quantile(probability)
        y = x^2
        dof < 5 && (c += 0.3 * (dof - 4.5) * (x + 0.6))
        c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c
        y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x
        y = expm1(a * y^2)
    else
        # Deep-tail branch: a series in the two-tailed probability itself,
        # where the normal deviate is a poor starting point.
        y =
            (
                (
                    1 / (((dof + 6) / (dof * y) - 0.089 * d - 0.822) * (dof + 2) * 3) +
                    0.5 / (dof + 4)
                ) * y - 1
            ) * (dof + 1) / (dof + 2) + 1 / y
    end
    sqrt(dof * y)
end

"""
$(SIGNATURES)

Per-hypothesis bin statistics for the soft, maximum-energy CFAR sync
detectors — one entry per candidate timing hypothesis. It backs both
[`_detect_bit_edge_cfar`](@ref) (hypothesis = bit-edge phase, bin = one
navigation bit, updated by [`_update_phase_accumulators!`](@ref)) and
[`_detect_secondary_code_cfar`](@ref) (hypothesis = overlay rotation, bin = one
overlay-wiped secondary-code period, updated by
[`_update_secondary_accumulators!`](@ref)); a signal uses at most one. Advanced
one primary-code block at a time, so detection is O(hypotheses) per block with no
growing history. Below, `period` is `blocks_per_bit` or the secondary-code
length `N`.

The vectors are mutated in place across the immutable [`BitBuffer`](@ref)
reconstructions (like `soft_bits`): an immutable representation would either be
copied (~660 B) into every per-block `BitBuffer` or boxed each block, both of
which allocate on the pre-sync hot path.

Fields (all length `period` once seeded; empty before the first block):

  - `open_bin_sum` — coherent sum of the hypothesis's *currently open* bin
    (overlay-wiped for the secondary-code detector).
  - `mean_bin_energy` / `bin_energy_sum_of_squared_deviations` — Welford running mean and
    sum of squared deviations (`M₂ = Σ(energyᵢ − mean)²`) of the hypothesis's
    *completed*-bin energies, for a numerically stable variance. The bin
    count is not stored — it is `div(num_blocks - hypothesis, period)`.
  - `last_bin_polarity` — sign (`±1`, `0` before the first bin) of the most
    recently completed bin's real part, i.e. the lock polarity.
"""
struct PhaseAccumulators
    open_bin_sum::Vector{ComplexF64}
    mean_bin_energy::Vector{Float64}
    bin_energy_sum_of_squared_deviations::Vector{Float64}
    last_bin_polarity::Vector{Int8}
end

PhaseAccumulators() = PhaseAccumulators(ComplexF64[], Float64[], Float64[], Int8[])

# Are the accumulators seeded for a `blocks_per_bit`-phase search yet?
@inline _is_seeded(accumulators::PhaseAccumulators, blocks_per_bit::Int) =
    length(accumulators.mean_bin_energy) == blocks_per_bit

# Size the (initially empty) accumulator vectors to `blocks_per_bit` phases
# and zero them.
function _seed_phase_accumulators!(accumulators::PhaseAccumulators, blocks_per_bit::Int)
    for vector in (
        accumulators.open_bin_sum,
        accumulators.mean_bin_energy,
        accumulators.bin_energy_sum_of_squared_deviations,
        accumulators.last_bin_polarity,
    )
        resize!(vector, blocks_per_bit)
    end
    _reset_phase_accumulators!(accumulators)
end

# Zero the accumulators in place, keeping their current length. Needed by
# `_resync_bit_buffer`: `_is_seeded` is still true there, so without this the next
# block would fold into the old lock's statistics. No-op on the hard path (empty).
function _reset_phase_accumulators!(accumulators::PhaseAccumulators)
    fill!(accumulators.open_bin_sum, zero(ComplexF64))
    fill!(accumulators.mean_bin_energy, 0.0)
    fill!(accumulators.bin_energy_sum_of_squared_deviations, 0.0)
    fill!(accumulators.last_bin_polarity, Int8(0))
    accumulators
end

"""
$(SIGNATURES)

Fold the `prompt` of the primary-code block at 0-based `block_index` into
the [`PhaseAccumulators`](@ref) (in place) for a `blocks_per_bit`-phase
bit-edge search. Each `phase`'s bins start at index `phase` and span
`blocks_per_bit` blocks; `prompt` is added to every phase whose bin is open
at `block_index` (`block_index ≥ phase`), and the single phase whose bin
ends at `block_index` has that completed bin's energy `|bin_sum|²` folded
into its Welford mean / sum-of-squared-deviations and its polarity recorded.
O(blocks_per_bit) per call, allocation-free.
"""
function _update_phase_accumulators!(
    accumulators::PhaseAccumulators,
    prompt::ComplexF64,
    block_index::Int,
    blocks_per_bit::Int,
)
    @inbounds for phase = 0:(blocks_per_bit-1)
        block_index < phase && continue
        accumulators.open_bin_sum[phase+1] += prompt
        if (block_index - phase) % blocks_per_bit == blocks_per_bit - 1
            bin_sum = accumulators.open_bin_sum[phase+1]
            bin_energy = abs2(bin_sum)
            # Welford update of this phase's completed-bin energy mean / M₂.
            completed_bin_count = div(block_index + 1 - phase, blocks_per_bit)
            energy_delta = bin_energy - accumulators.mean_bin_energy[phase+1]
            accumulators.mean_bin_energy[phase+1] += energy_delta / completed_bin_count
            accumulators.bin_energy_sum_of_squared_deviations[phase+1] +=
                energy_delta * (bin_energy - accumulators.mean_bin_energy[phase+1])
            accumulators.last_bin_polarity[phase+1] = real(bin_sum) < 0 ? Int8(-1) : Int8(1)
            accumulators.open_bin_sum[phase+1] = zero(ComplexF64)
        end
    end
    accumulators
end

"""
$(SIGNATURES)

Shared CFAR (constant-false-alarm-rate) decision core for the soft,
maximum-energy sync detectors [`_detect_bit_edge_cfar`](@ref) and
[`_detect_secondary_code_cfar`](@ref). Given the per-hypothesis statistics in
[`PhaseAccumulators`](@ref) it finds the maximum-energy hypothesis and its
closest competitor and decides whether the peak is significant enough to lock.

Arguments:

  - `mean_bin_energy` — per-hypothesis mean completed-bin energy
    (`accumulators.mean_bin_energy`).
  - `bin_energy_sum_of_squared_deviations` — per-hypothesis Welford `M₂`
    (`accumulators.bin_energy_sum_of_squared_deviations`), used only for the
    peak.
  - `num_blocks` — total primary-code blocks folded so far.
  - `period` — number of primary-code blocks per bin, which is also the number
    of competing hypotheses: `blocks_per_bit` (20 for L1 C/A) or the
    secondary-code length `N`. A hypothesis `h ∈ 0:period-1` has completed-bin
    count `div(num_blocks - h, period)`.
  - `confidence` — target `1 − P(false lock)`.

Returns `(accepted, peak_index, peak_bin_count)`: `accepted` is whether the
peak's energy gap over its runner-up is statistically significant (the caller
still applies its own bin-boundary gate before firing); `peak_index` is the
0-based winning hypothesis (`-1` if none qualifies) and `peak_bin_count` its
completed-bin count.

# Statistic

Each hypothesis coherently sums its `period`-block bins and the per-bin energy
`|Σ|²` is folded into a running mean; the true hypothesis keeps full coherent
gain on every bin while wrong ones lose energy (L1 C/A: wrong phases straddle a
data-bit transition; secondary code: wrong rotations fail to wipe the overlay).
The per-hypothesis mean bin energy is therefore the maximum-likelihood timing
statistic.

# CFAR confidence

The noise scale is the *bin-to-bin* sample variance of the winning hypothesis's
own completed-bin energies (Welford). The true hypothesis's bins vary only with
thermal noise and slow drift — exactly the spread the gap must be compared
against. The peak is accepted only when

    z_score = energy_gap / standard_error   ≥   t⁻¹(1 - false_alarm_probability/(period - 1);  ν = peak_bin_count − 1)

where the standard error combines the peak's per-bin energy variance over the
peak and runner-up bin counts, `false_alarm_probability = 1 - confidence` is
Bonferroni-split over the `period - 1` competitors, and the Student-t quantile
[`_t_quantile`](@ref) is a small-sample penalty for a variance estimated over
few bins. A real peak's gap grows `z_score` like the square root of the bin
count, so the detector self-paces with C/N₀, while a drift-only asymmetry keeps
`z_score` bounded and never locks.
"""
@inline function _cfar_decide(
    mean_bin_energy::AbstractVector{Float64},
    bin_energy_sum_of_squared_deviations::AbstractVector{Float64},
    num_blocks::Int,
    period::Int,
    confidence::Float64,
)
    # Need at least two complete bins on some hypothesis before any hypothesis
    # can serve as a runner-up; below that a transition cannot have been
    # localized yet.
    num_blocks < 2 * period && return (false, -1, 0)

    # Single O(period) pass: the peak (highest mean bin energy among hypotheses
    # with ≥ 2 complete bins) and the two highest hypotheses overall. The
    # runner-up is the higher of those two that isn't the peak.
    peak_index = -1
    peak_energy = -1.0
    peak_bin_count = 0          # best, ≥ 2 bins
    best_index = -1
    best_energy = -1.0
    best_bin_count = 0          # highest, ≥ 1 bin
    second_best_energy = -1.0
    second_best_bin_count = 0             # 2nd highest, ≥ 1 bin
    @inbounds for h = 0:(period-1)
        bin_count = div(num_blocks - h, period)
        bin_count < 1 && continue
        energy = mean_bin_energy[h+1]
        if energy > best_energy
            second_best_energy = best_energy
            second_best_bin_count = best_bin_count
            best_energy = energy
            best_index = h
            best_bin_count = bin_count
        elseif energy > second_best_energy
            second_best_energy = energy
            second_best_bin_count = bin_count
        end
        if bin_count >= 2 && energy > peak_energy
            peak_energy = energy
            peak_index = h
            peak_bin_count = bin_count
        end
    end
    # `num_blocks >= 2 * period` guarantees hypothesis 0 has ≥ 2 bins, so
    # `peak_index` is always set.
    peak_index < 0 && return (false, -1, 0)

    # Runner-up: the highest-energy hypothesis that isn't the peak (any bin count).
    runner_up_energy, runner_up_bin_count =
        best_index == peak_index ? (second_best_energy, second_best_bin_count) :
        (best_energy, best_bin_count)
    runner_up_bin_count < 1 && return (false, peak_index, peak_bin_count)

    energy_gap = peak_energy - runner_up_energy
    energy_gap <= 0 && return (false, peak_index, peak_bin_count)

    # Noise scale: the peak's own bin-to-bin energy variance (see docstring). A
    # within-bin residual would miss slow drift and mistake a tiny systematic
    # asymmetry for a peak. The runner-up is assumed to share this variance (its
    # straddling bins can only make it larger, so this is the fastest-locking
    # choice that still rejects drift-only asymmetries).
    @inbounds bin_energy_variance =
        bin_energy_sum_of_squared_deviations[peak_index+1] / (peak_bin_count - 1)
    standard_error =
        sqrt(bin_energy_variance * (1 / peak_bin_count + 1 / runner_up_bin_count))
    z_score =
        standard_error > 0 ? energy_gap / standard_error : (energy_gap > 0 ? Inf : 0.0)

    # Clamp the quantile argument to the open interval (0, 1): `confidence = 1.0`
    # (or rounding up to it) would give a NaN threshold, and since NaN comparisons
    # are false the gate below would pass — "maximum confidence" would lock
    # immediately.
    false_alarm_probability = 1 - confidence
    quantile_argument =
        clamp(1 - false_alarm_probability / (period - 1), nextfloat(0.0), prevfloat(1.0))
    # Student-t threshold, not normal; see `_t_quantile`.
    z_threshold = _t_quantile(quantile_argument, peak_bin_count - 1)
    z_score < z_threshold && return (false, peak_index, peak_bin_count)

    (true, peak_index, peak_bin_count)
end

"""
$(SIGNATURES)

Soft-decision, CFAR bit-edge detector for signals whose navigation bit spans
more than one primary-code period with no secondary code (selected by
[`uses_soft_bit_edge_detection`](@ref); currently GPS L1 C/A, 20 blocks per bit).

`accumulators` holds the per-phase bin statistics (advanced by
[`_update_phase_accumulators!`](@ref)), `num_blocks` is the total number of
prompts seen, and the hypothesis test is the shared [`_cfar_decide`](@ref) over
the edge phases `0:blocks_per_bit-1`.

# Boundary firing

`found = true` is reported only when the most recent block also *ends* the
winning phase's bit (`num_blocks % blocks_per_bit == peak_phase`), so the
upcoming integration starts a fresh bit and the post-sync integration is aligned
to the true bit grid. Firing only at the winner's own boundary makes the
off-by-one lock of issue #124 impossible. `polarity` is the sign of the most
recently completed bin's coherent sum.
"""
function _detect_bit_edge_cfar(
    accumulators::PhaseAccumulators,
    blocks_per_bit::Int,
    confidence::Float64,
    num_blocks::Int,
)
    accepted, peak_phase, _ = _cfar_decide(
        accumulators.mean_bin_energy,
        accumulators.bin_energy_sum_of_squared_deviations,
        num_blocks,
        blocks_per_bit,
        confidence,
    )
    accepted || return SyncResult(false, 0, Int8(0))

    # Fire only at the winning phase's own bit boundary so the upcoming
    # integration starts a new bit.
    num_blocks % blocks_per_bit != peak_phase && return SyncResult(false, 0, Int8(0))

    @inbounds polarity =
        accumulators.last_bin_polarity[peak_phase+1] < 0 ? Int8(-1) : Int8(+1)
    SyncResult(true, 0, polarity)
end

"""
$(SIGNATURES)

Fold the `prompt` of the primary-code block at 0-based `block_index` into the
[`PhaseAccumulators`](@ref) (in place) for a soft, maximum-energy secondary-code
rotation search of period `N = secondary_code_length` — the secondary-code
analog of [`_update_phase_accumulators!`](@ref). Each rotation hypothesis
`d ∈ 0:N-1` has bins of `N` blocks starting at index `d`, and each block is
multiplied by the *known* secondary chip before being summed, so the correct
rotation wipes the overlay and adds coherently while wrong ones lose energy.

Under hypothesis `d`, block `i` carries secondary chip `mod(i - d, N)` of
`signal`'s code for `prn`, so chip 0 is the physical secondary chip 0 — the same
overlay the post-sync replica applies. O(N) per call, allocation-free.
"""
function _update_secondary_accumulators!(
    accumulators::PhaseAccumulators,
    prompt::ComplexF64,
    block_index::Int,
    secondary_code_length::Int,
    signal::AbstractGNSSSignal,
    prn::Integer,
)
    N = secondary_code_length
    secondary_code = get_secondary_code(signal)
    @inbounds for d = 0:(N-1)
        block_index < d && continue
        chip = (block_index - d) % N
        overlay = GNSSSignals.secondary_value(secondary_code, prn, chip)
        accumulators.open_bin_sum[d+1] += prompt * overlay
        if chip == N - 1
            bin_sum = accumulators.open_bin_sum[d+1]
            bin_energy = abs2(bin_sum)
            # Welford update of this rotation's completed-bin energy mean / M₂.
            completed_bin_count = div(block_index + 1 - d, N)
            energy_delta = bin_energy - accumulators.mean_bin_energy[d+1]
            accumulators.mean_bin_energy[d+1] += energy_delta / completed_bin_count
            accumulators.bin_energy_sum_of_squared_deviations[d+1] +=
                energy_delta * (bin_energy - accumulators.mean_bin_energy[d+1])
            accumulators.last_bin_polarity[d+1] = real(bin_sum) < 0 ? Int8(-1) : Int8(1)
            accumulators.open_bin_sum[d+1] = zero(ComplexF64)
        end
    end
    accumulators
end

"""
$(SIGNATURES)

Soft-decision, CFAR secondary-code sync detector — the analog of
[`_detect_bit_edge_cfar`](@ref) for signals with a short periodic overlay code
(selected by [`uses_soft_secondary_code_detection`](@ref)). `accumulators` holds
the per-rotation bin statistics (advanced by
[`_update_secondary_accumulators!`](@ref)); `secondary_code_length` (`N`) is both
the bin length and the number of hypotheses, and the test is the shared
[`_cfar_decide`](@ref). Unlike the hard [`_secondary_code_search`](@ref) it uses
the soft prompt and a CFAR confidence, so it rejects the noise-driven false locks
a short hard template match is prone to.

# Boundary firing

`found = true` is reported only when the most recent block also *ends* the
winning rotation's period (`num_blocks % N == peak_rotation`), so the upcoming
integration starts at physical secondary chip 0 and `SyncResult.phase` is always
`0` (the anchor [`_snap_code_phase_from_synced_signal`](@ref) uses). `polarity` is
the sign of the winning period's overlay-wiped sum.
"""
function _detect_secondary_code_cfar(
    accumulators::PhaseAccumulators,
    secondary_code_length::Int,
    confidence::Float64,
    num_blocks::Int,
)
    accepted, peak_rotation, _ = _cfar_decide(
        accumulators.mean_bin_energy,
        accumulators.bin_energy_sum_of_squared_deviations,
        num_blocks,
        secondary_code_length,
        confidence,
    )
    accepted || return SyncResult(false, 0, Int8(0))

    # Fire only at the winning rotation's own period boundary so the upcoming
    # integration starts at secondary chip 0.
    num_blocks % secondary_code_length != peak_rotation &&
        return SyncResult(false, 0, Int8(0))

    @inbounds polarity =
        accumulators.last_bin_polarity[peak_rotation+1] < 0 ? Int8(-1) : Int8(+1)
    SyncResult(true, 0, polarity)
end

"""
$(SIGNATURES)

**Hard-decision** secondary-code rotation search, used by signals not on a soft
detector — currently the 1800-chip overlay pilots GPS L1C-P and BeiDou B1C-P. It
locks after a single secondary-code period in the worst case and recovers the
secondary-code phase.

`received` is the prompt-sign sliding window with the newest block in bit 0.
`reference` is the secondary code packed in the same newest-first order — bit
`i` holds chip `N - 1 - i` — so that when the most recent `N` blocks span one
period ending on its last chip, `received & mask == reference`.

The search rotates the low `N` bits of `received` left by `d ∈ 0:N-1`, tracking
the best positive- and negated-polarity Hamming match in one pass. Rotation `d`
maps to the secondary chip of the **upcoming** integration as
`phase = mod(N - d, N)`, the anchor of
[`_snap_code_phase_from_synced_signal`](@ref). Returns `SyncResult(false, 0, 0)`
when the best distance exceeds `max_errors`. Inlined so the per-signal constants
fold at the call site.
"""
@inline function _secondary_code_search(
    received::B,
    reference::B,
    secondary_code_length::Int,
    max_errors::Int,
) where {B<:Unsigned}
    N = secondary_code_length
    # `reference` lives in the low `N` bits; mask the search to that window.
    # `one(B) << N` is undefined when `N` equals the full bit width
    # (e.g. UInt1800 for L1C-P), so special-case the exact-width buffer.
    mask = N == 8 * sizeof(B) ? ~zero(B) : (one(B) << N) - one(B)
    masked = received & mask
    best_d = 0
    best_dist = N + 1
    best_pol = Int8(0)
    @inbounds for d = 0:(N-1)
        # Rotate-left by `d` within the N-bit window. `d == 0` is special-
        # cased because `masked >> N` is undefined for an exact-width buffer.
        shifted = d == 0 ? masked : ((masked << d) | (masked >> (N - d))) & mask
        dist_pos = count_ones(shifted ⊻ reference)
        # Both operands occupy only the low `N` bits, so the negated-polarity
        # distance is the complement within that window.
        dist_neg = N - dist_pos
        if dist_pos < best_dist
            best_dist = dist_pos
            best_d = d
            best_pol = Int8(+1)
        end
        if dist_neg < best_dist
            best_dist = dist_neg
            best_d = d
            best_pol = Int8(-1)
        end
    end
    best_dist > max_errors && return SyncResult(false, 0, Int8(0))
    SyncResult(true, mod(N - best_d, N), best_pol)
end

"""
$(SIGNATURES)

Shared detector body for signals with no sub-block boundary to find, which
therefore report `found = true` from the very first integration.

Two shapes qualify. Signals that broadcast one channel symbol per primary code
period (GPS L1C-D, GPS L2CM, Galileo E1B / E6-B, BeiDou B2b-I / B1C-D): the
primary-block signs are the symbol stream, and downstream consumers
(GNSSDecoder.jl) resolve the residual ±1 polarity via the navigation preamble.
And Galileo E5a-QP, which has neither data nor an overlay, so every block
boundary is equivalent and the lock only gates its switch to whole-code-cycle
integration.
"""
@inline function _detect_symbol_is_code_block_sync(
    ::AbstractGNSSSignal,
    ::Integer,           # PRN — ignored
    ::Unsigned,
    ::Integer,
)
    SyncResult(true, 0, Int8(+1))
end

"""
$(SIGNATURES)

Shared hard-path detector body for signals with a periodic secondary / overlay
code: wait until the sliding `code_block_bits` window covers one full
secondary-code period, then run [`_secondary_code_search`](@ref) against
[`_packed_secondary_code`](@ref), allowing `floor(tolerance × N)` errors (see
[`get_bit_edge_or_secondary_code_tolerance`](@ref)). Returns
`SyncResult(false, 0, 0)` until `N` blocks are buffered.

Every secondary-coded signal's `detect_bit_or_secondary_code_sync` delegates
here; for those with `uses_soft_secondary_code_detection` it is reached only if
that trait is overridden to `false`. A new secondary-coded signal needs
`get_secondary_code`, a `get_code_block_buffer_type` at least `N` bits wide, and a
`detect_bit_or_secondary_code_sync` method delegating here; specialize
[`_packed_secondary_code`](@ref) only if its overlay isn't reachable through
`get_secondary_code`.
"""
@inline function _detect_secondary_code_sync(
    signal::AbstractGNSSSignal,
    prn::Integer,
    code_block_bits::B,
    num_code_blocks::Integer,
) where {B<:Unsigned}
    secondary_code_length = get_secondary_code_length(signal)
    num_code_blocks < secondary_code_length && return SyncResult(false, 0, Int8(0))
    max_errors =
        floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * secondary_code_length)
    _secondary_code_search(
        code_block_bits,
        _packed_secondary_code(B, signal, prn),
        secondary_code_length,
        max_errors,
    )
end

"""
$(SIGNATURES)

Reference for [`_detect_secondary_code_sync`](@ref): the signal's secondary code
for `prn`, packed into `B` in the newest-first order of
[`_secondary_code_search`](@ref) (bit `i` holds chip `N - 1 - i`).

The generic method below covers every signal GNSSSignals defines; specialize it
only for a signal whose overlay is not reachable through `get_secondary_code`.
"""
function _packed_secondary_code end

# Sets bit `N - 1 - k` iff secondary chip `k` of `get_secondary_code` is positive,
# for both `SharedSecondaryCode` and `PerPRNSecondaryCode`. Because every
# signal uses this one convention, `polarity = +1` means the same thing for all:
# the prompt signs follow `get_secondary_code`, as the post-sync replica does.
@inline function _packed_secondary_code(
    ::Type{B},
    signal::AbstractGNSSSignal,
    prn::Integer,
) where {B<:Unsigned}
    secondary_code = get_secondary_code(signal)
    N = get_secondary_code_length(signal)
    packed = zero(B)
    @inbounds for k = 0:(N-1)
        if GNSSSignals.secondary_value(secondary_code, prn, k) > 0
            packed |= one(B) << (N - 1 - k)
        end
    end
    packed
end

"""
$(SIGNATURES)

BitBuffer to buffer bits.

`code_block_buffer` is the pre-sync sliding window of prompt signs, of width `B`
([`get_code_block_buffer_type`](@ref)); it is dead state after sync. Decoded
navigation bits accumulate in `soft_bits` as **soft bits** (one
polarity-corrected coherent prompt sum per bit; the hard bit is `soft_bit > 0`),
the only bit store, unbounded between resets.

`phase_acc` holds the [`PhaseAccumulators`](@ref) of the signal's soft CFAR sync
detector, if any ([`uses_soft_bit_edge_detection`](@ref) /
[`uses_soft_secondary_code_detection`](@ref)); it stays empty on the hard path.
"""
struct BitBuffer{B<:Unsigned}
    code_block_buffer::B
    code_block_buffer_length::Int
    found::Bool
    secondary_phase::Int      # 0 until found; secondary-chip offset post-sync
    polarity::Int8            # +1 or -1 once found; 0 before sync
    prompt_accumulator::ComplexF64
    prompt_accumulator_integrated_code_blocks::Int
    soft_bits::Vector{Float32}
    phase_acc::PhaseAccumulators
end

# Untyped default: a `UInt128` search buffer. `TrackedSignal` uses the typed
# constructor below with the width from `get_code_block_buffer_type`.
function BitBuffer()
    BitBuffer{UInt128}(
        zero(UInt128),
        0,
        false,
        0,
        Int8(0),
        complex(0.0, 0.0),
        0,
        Float32[],
        PhaseAccumulators(),
    )
end

# Typed empty constructor used by the per-signal `TrackedSignal` path.
function BitBuffer{B}() where {B<:Unsigned}
    BitBuffer{B}(
        zero(B),
        0,
        false,
        0,
        Int8(0),
        complex(0.0, 0.0),
        0,
        Float32[],
        PhaseAccumulators(),
    )
end

# Convenience constructor without phase / polarity (both zero) and with empty
# accumulators, for tests and benchmarks. `soft_bits` is aliased, not copied, so
# the caller can observe in-place pushes.
function BitBuffer(
    code_block_buffer::B,
    code_block_buffer_length::Integer,
    found::Bool,
    prompt_accumulator::Complex,
    prompt_accumulator_integrated_code_blocks::Integer,
    soft_bits::Vector{Float32} = Float32[],
) where {B<:Unsigned}
    BitBuffer{B}(
        code_block_buffer,
        Int(code_block_buffer_length),
        found,
        0,
        Int8(0),
        ComplexF64(prompt_accumulator),
        Int(prompt_accumulator_integrated_code_blocks),
        soft_bits,
        PhaseAccumulators(),
    )
end

@inline length(bit_buffer::BitBuffer) = Base.length(bit_buffer.soft_bits)
@inline has_bit_or_secondary_code_been_found(bit_buffer::BitBuffer) = bit_buffer.found

# The soft bits: the summed filtered prompt of each completed bit (hard decision
# `soft_bit > 0`), emptied at the start of each `track` call. A plain comment, not
# a docstring, so `checkdocs = :exports` doesn't require it in the manual.
@inline get_soft_bits(bit_buffer::BitBuffer) = bit_buffer.soft_bits

"""
$(SIGNATURES)

Width `B` of the packed prompt-sign buffer (`BitBuffer.code_block_buffer`) for
`signal`, returned as a concrete `Unsigned` subtype.

Only the hard [`_secondary_code_search`](@ref) reads that buffer, so a hard-path
signal needs at least `N` bits (the search masks to the low `N`; the 1800-chip
overlay pilots use the exact-width `UInt1800`). Soft-detector signals keep a
width that holds their horizon anyway, so the hard path stays available if the
trait is overridden. Signals with nothing to search use `UInt8`. Per-signal
values are tabulated in docs/src/bit_sync.md.

The default is `UInt64`. The width is a type parameter of `BitBuffer{B}` and
`TrackedSignal`, keeping construction type-stable.
"""
@inline get_code_block_buffer_type(::AbstractGNSSSignal) = UInt64

"""
$(SIGNATURES)

Per-signal Hamming tolerance of the **hard-decision** sweep
[`_secondary_code_search`](@ref): the largest fraction of bit-flips accepted
before reporting `found = true`, converted at the call site to
`max_errors = floor(Int, tolerance × window_size)`.

Only hard-path signals read it — currently the 1800-chip overlay pilots GPS
L1C-P and BeiDou B1C-P. The soft CFAR detectors are tuned by
[`get_bit_edge_detection_confidence`](@ref) instead, and the
one-symbol-per-code-period detectors ignore it.

Default is `0.025` (2.5 %), i.e. `max_errors = 45` over 1800 chips.

# Overriding

To loosen the tolerance for low-C/N₀ work, dispatch the trait on the (hard-path)
signal type in your own module:

```julia
Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL1C_P) = 0.05
```

The override takes effect at the next detector call — no TrackState rebuild
needed. It has no effect on a soft-detector signal.
"""
@inline get_bit_edge_or_secondary_code_tolerance(::AbstractGNSSSignal) = 0.025

"""
$(SIGNATURES)

Whether `signal`'s bit edge is located with the soft-decision,
maximum-energy CFAR detector [`_detect_bit_edge_cfar`](@ref) rather than the
hard-decision `detect_bit_or_secondary_code_sync` path.

The detector is parameterised only by the blocks per navigation bit
([`_calc_num_code_blocks_that_form_a_bit`](@ref)), so the default enables it for
any signal whose bit spans **more than one** primary-code period **and** which
has no secondary code (currently only GPS L1 C/A). BeiDou B1I/B3I do not
qualify even on GEO PRNs, which carry no overlay, because their signal type
reports a 20-chip secondary code (see docs/src/bit_sync.md).

Override per signal type to force the choice, e.g. to disable it:

```julia
Tracking.uses_soft_bit_edge_detection(::SomeSignal) = false
```

Constant-folded per signal type, so the branch in [`_buffer_find_bit`](@ref)
compiles away.
"""
@inline uses_soft_bit_edge_detection(signal::AbstractGNSSSignal) =
    _calc_num_code_blocks_that_form_a_bit(signal) > 1 &&
    get_secondary_code_length(signal) == 1

"""
$(SIGNATURES)

Whether `signal`'s secondary/overlay code is located with the soft-decision,
maximum-energy CFAR detector [`_detect_secondary_code_cfar`](@ref) rather than
the hard-decision sweep [`_secondary_code_search`](@ref).

The soft detector coherently integrates one full secondary-code period per bin,
so the default enables it only for **short** codes,
`1 < get_secondary_code_length(signal) ≤ 100` (periods up to 100 ms). The
1800-chip overlay pilots (GPS L1C-P, BeiDou B1C-P; 18 s) are far too long to
integrate coherently — and too long to be false-lock-prone — so they keep the
hard sweep. Mutually exclusive with [`uses_soft_bit_edge_detection`](@ref),
which requires *no* secondary code.

Override per signal type to force the choice, e.g. to disable it:

```julia
Tracking.uses_soft_secondary_code_detection(::GPSL5I) = false
```

Constant-folded per signal type, so the branch in [`_buffer_find_bit`](@ref)
compiles away.
"""
@inline uses_soft_secondary_code_detection(signal::AbstractGNSSSignal) =
    1 < get_secondary_code_length(signal) <= 100

"""
$(SIGNATURES)

Target confidence (one minus the probability of a false lock) for the
soft-decision CFAR sync detectors [`_detect_bit_edge_cfar`](@ref) and
[`_detect_secondary_code_cfar`](@ref).

Default `0.999`: the detector integrates until the maximum-energy hypothesis
beats its closest competitor with this confidence — as little as two bins for a
clean signal, longer in noise. Lower it to lock faster at the cost of more false
locks.

# Overriding

```julia
Tracking.get_bit_edge_detection_confidence(::GPSL1CA) = 0.9999
```

Takes effect at the next detector call — no TrackState rebuild needed.
"""
@inline get_bit_edge_detection_confidence(::AbstractGNSSSignal) = 0.999

# Number of primary-code blocks that form one navigation bit.
# Returns 0 for pilot signals (`data_frequency = 0`), where the concept is
# undefined; callers must guard for that case before using the result.
@inline function _calc_num_code_blocks_that_form_a_bit(signal::AbstractGNSSSignal)
    data_freq = get_data_frequency(signal)
    iszero(data_freq) && return 0
    Int(get_code_frequency(signal) / (get_code_length(signal) * data_freq))
end

"""
$(SIGNATURES)

Buffer data bits based on the prompt accumulation and the current prompt value.

Post-sync, `integrated_code_blocks` is added to a running count and one soft bit
is emitted each time that count reaches the signal's blocks-per-bit. A record
that carries the count *past* the boundary (only possible from an external,
non-bit-aligned `CorrelatorOutput` producer) drops bit sync and warns instead;
see issue #238 and docs/src/bit_sync.md.
"""
function buffer(
    signal::AbstractGNSSSignal,
    prn::Integer,
    bit_buffer::BitBuffer{B},
    integrated_code_blocks,
    prompt,
) where {B<:Unsigned}
    # 0 for pilots (`get_data_frequency = 0`), see the helper.
    num_code_blocks_that_form_a_bit = _calc_num_code_blocks_that_form_a_bit(signal)

    if (bit_buffer.found == false)
        return _buffer_find_bit(
            signal,
            prn,
            bit_buffer,
            num_code_blocks_that_form_a_bit,
            integrated_code_blocks,
            prompt,
        )
    end

    # Pilots carry no data: nothing to decode post-sync. The lock state stays
    # for code-phase anchoring and longer coherent integration.
    num_code_blocks_that_form_a_bit == 0 && return bit_buffer

    prompt_accumulator = bit_buffer.prompt_accumulator + prompt
    prompt_accumulator_integrated_code_blocks =
        bit_buffer.prompt_accumulator_integrated_code_blocks + integrated_code_blocks

    if prompt_accumulator_integrated_code_blocks > num_code_blocks_that_form_a_bit
        # Overshot the bit boundary (e.g. 18 + a 3-block record = 21 for L1 C/A).
        # The equality test below could never be met again, so no further bit
        # would be pushed (issue #238). Only an external `CorrelatorOutput`
        # producer gets here. Resync rather than emit at `>=`: the straddling
        # record's energy belongs to two bits and would silently corrupt both.
        _warn_bit_boundary_overshoot(
            get_signal_id(signal),
            prn,
            prompt_accumulator_integrated_code_blocks,
            num_code_blocks_that_form_a_bit,
        )
        return _resync_bit_buffer(bit_buffer)
    elseif prompt_accumulator_integrated_code_blocks == num_code_blocks_that_form_a_bit
        # Flip the decoded bit if the detector locked at negative polarity.
        bit_acc =
            bit_buffer.polarity < 0 ? -real(prompt_accumulator) : real(prompt_accumulator)
        # Sign = hard decision, magnitude = confidence. Unbounded until `reset`.
        push!(bit_buffer.soft_bits, Float32(bit_acc))
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_length,
            true,
            bit_buffer.secondary_phase,
            bit_buffer.polarity,
            zero(prompt_accumulator),
            0,
            bit_buffer.soft_bits,
            bit_buffer.phase_acc,
        )
    else
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_length,
            true,
            bit_buffer.secondary_phase,
            bit_buffer.polarity,
            prompt_accumulator,
            prompt_accumulator_integrated_code_blocks,
            bit_buffer.soft_bits,
            bit_buffer.phase_acc,
        )
    end
end

# Return `bit_buffer` to its pre-sync state after a bit-boundary overshoot,
# zeroing the per-hypothesis statistics. `soft_bits` (complete, correct bits) is
# kept, including its vector identity, which the caller may be holding.
@inline function _resync_bit_buffer(bit_buffer::BitBuffer{B}) where {B<:Unsigned}
    BitBuffer{B}(
        zero(B),
        0,
        false,
        0,
        Int8(0),
        complex(0.0, 0.0),
        0,
        bit_buffer.soft_bits,
        _reset_phase_accumulators!(bit_buffer.phase_acc),
    )
end

# `_id` is signal- and PRN-specific so `maxlog = 1` doesn't hide how many
# satellites a producer's record-sizing bug hits. Warn rather than throw: it
# should cost a bit-sync, not the receiver.
@noinline function _warn_bit_boundary_overshoot(
    signal_id::Symbol,
    prn::Integer,
    accumulated_blocks::Integer,
    blocks_per_bit::Integer,
)
    @warn """
          Signal `:$signal_id` PRN $prn: a correlator record carried the bit \
          accumulator past the navigation-bit boundary ($accumulated_blocks of \
          $blocks_per_bit code blocks), which only a record not aligned to that \
          boundary can do. Dropping bit sync and re-running the detector; the \
          partial bit is discarded, the bits decoded so far are kept. Check the \
          record sizing of the external `CorrelatorOutput` producer — it must \
          not straddle a bit boundary.""" _id =
        Symbol(:bit_boundary_overshoot_, signal_id, :_, prn) maxlog = 1
    nothing
end

function _buffer_find_bit(
    signal,
    prn::Integer,
    bit_buffer::BitBuffer{B},
    num_code_blocks_that_form_a_bit,
    integrated_code_blocks,
    prompt,
) where {B<:Unsigned}
    if (integrated_code_blocks != 1)
        error(
            "The number code blocks must be equal to 1 if bit or secondary code hasn't been found yet.",
        )
    end
    code_block_buffer = (bit_buffer.code_block_buffer << 1) + B(real(prompt) > 0)
    code_block_buffer_length = bit_buffer.code_block_buffer_length + 1

    # Soft CFAR detectors (bit-edge or secondary-code; a signal uses at most
    # one) share `phase_acc`; everything else takes the hard sliding-window
    # path. The branch folds at compile time per signal type.
    phase_acc = bit_buffer.phase_acc
    if uses_soft_bit_edge_detection(signal)
        blocks_per_bit = num_code_blocks_that_form_a_bit
        _is_seeded(phase_acc, blocks_per_bit) ||
            _seed_phase_accumulators!(phase_acc, blocks_per_bit)
        _update_phase_accumulators!(
            phase_acc,
            ComplexF64(prompt),
            code_block_buffer_length - 1,
            blocks_per_bit,
        )
        sync = _detect_bit_edge_cfar(
            phase_acc,
            blocks_per_bit,
            get_bit_edge_detection_confidence(signal),
            code_block_buffer_length,
        )
    elseif uses_soft_secondary_code_detection(signal)
        secondary_code_length = get_secondary_code_length(signal)
        _is_seeded(phase_acc, secondary_code_length) ||
            _seed_phase_accumulators!(phase_acc, secondary_code_length)
        _update_secondary_accumulators!(
            phase_acc,
            ComplexF64(prompt),
            code_block_buffer_length - 1,
            secondary_code_length,
            signal,
            prn,
        )
        sync = _detect_secondary_code_cfar(
            phase_acc,
            secondary_code_length,
            get_bit_edge_detection_confidence(signal),
            code_block_buffer_length,
        )
    else
        sync = detect_bit_or_secondary_code_sync(
            signal,
            prn,
            code_block_buffer,
            code_block_buffer_length,
        )
    end
    if !sync.found
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_length,
            false,
            0,
            Int8(0),
            complex(0.0, 0.0),
            0,
            bit_buffer.soft_bits,
            phase_acc,
        )
    end
    if get_secondary_code_length(signal) > 1
        # Secondary-code signals: the pre-sync signs carry the overlay, not data,
        # so no bits are recovered. The hard search can lock at any secondary
        # chip, so seed the block count with `sync.phase` to make the first bit
        # end on the data-bit boundary (issue #125). Pilots never read the seed.
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_length,
            true,
            sync.phase,
            sync.polarity,
            complex(0.0, 0.0),
            sync.phase,
            bit_buffer.soft_bits,
            phase_acc,
        )
    end
    if num_code_blocks_that_form_a_bit == 0
        # Dataless and no overlay (Galileo E5a-QP): no bits to recover, and the
        # divisions below would be by zero. Nothing reads `secondary_phase`.
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_length,
            true,
            0,
            sync.polarity,
            complex(0.0, 0.0),
            0,
            bit_buffer.soft_bits,
            phase_acc,
        )
    end
    num_bits = min(
        div(code_block_buffer_length, num_code_blocks_that_form_a_bit),
        div(sizeof(code_block_buffer) * 8, num_code_blocks_that_form_a_bit),
    )
    # Hoisted so the closure below doesn't capture `sync`, which is assigned in
    # several branches and would be boxed — a per-block allocation on the
    # pre-sync hot path.
    sync_polarity = Int(sync.polarity)
    for bit_index = num_bits:-1:1     # oldest recovered bit first
        # Apply the lock polarity to the recovered pre-sync bits too, as
        # `buffer` does post-sync (issue #127).
        bit_sum =
            sum(0:(num_code_blocks_that_form_a_bit-1)) do code_block_index
                buffer_code_block_index =
                    (bit_index - 1) * num_code_blocks_that_form_a_bit + code_block_index
                ((code_block_buffer & (one(B) << buffer_code_block_index)) > 0) * 2 - 1
            end * sync_polarity
        # `bit_sum` is a ±1 sign-vote count; scale it by the sync-time prompt
        # magnitude so it is in the same units as the post-sync soft bits.
        push!(bit_buffer.soft_bits, Float32(bit_sum * abs(prompt)))
    end
    return BitBuffer{B}(
        code_block_buffer,
        code_block_buffer_length,
        true,
        sync.phase,
        sync.polarity,
        complex(0, 0),
        0,
        bit_buffer.soft_bits,
        phase_acc,
    )
end

"""
$(SIGNATURES)

Walk a just-synced `bit_buffer`'s `secondary_phase` forward by
`num_code_blocks` primary-code blocks (modulo the secondary-code length). No-op
for signals without a secondary code, where the field is unused.

`secondary_phase` is read once, by
[`_snap_code_phase_from_synced_signal`](@ref) after the whole chunk is folded,
but the detector reports it for the block right after the syncing record. Every
further record in the same chunk (those `_apply_correlator_output` marks
`correlated_pre_sync`) must move it along, or the post-sync replica applies the
wrong overlay chip — the secondary-code sibling of issue #219.
"""
@inline function _advance_secondary_phase(
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer{B},
    num_code_blocks::Integer,
) where {B<:Unsigned}
    secondary_code_length = get_secondary_code_length(signal)
    secondary_code_length > 1 || return bit_buffer
    BitBuffer{B}(
        bit_buffer.code_block_buffer,
        bit_buffer.code_block_buffer_length,
        bit_buffer.found,
        mod(bit_buffer.secondary_phase + Int(num_code_blocks), secondary_code_length),
        bit_buffer.polarity,
        bit_buffer.prompt_accumulator,
        bit_buffer.prompt_accumulator_integrated_code_blocks,
        bit_buffer.soft_bits,
        bit_buffer.phase_acc,
    )
end

function reset(bit_buffer::BitBuffer{B}) where {B<:Unsigned}
    empty!(bit_buffer.soft_bits)
    BitBuffer{B}(
        bit_buffer.code_block_buffer,
        bit_buffer.code_block_buffer_length,
        bit_buffer.found,
        bit_buffer.secondary_phase,
        bit_buffer.polarity,
        bit_buffer.prompt_accumulator,
        bit_buffer.prompt_accumulator_integrated_code_blocks,
        bit_buffer.soft_bits,
        bit_buffer.phase_acc,
    )
end

# Whether the replica wipes every sign modulation off the signal's prompt, so
# that consecutive prompts share their sign: a dataless signal, synced to its
# secondary code where it has one. `synced` must hold for the record's
# correlation, i.e. before the fold.
@inline _is_wiped_off(signal::AbstractGNSSSignal, synced::Bool) =
    iszero(get_data_frequency(signal)) && (get_secondary_code_length(signal) == 1 || synced)

# The polarity of a pilot's prompt according to its secondary-code sync, which
# the four-quadrant PLL reads the prompt with, or 0 for the Costas one:
# a data signal, a pilot before its sync, and a pilot without a secondary code
# (GPS L2 CL, Galileo E5a-QP), whose sync reads no sign. The sync reads it off
# one whole secondary period of the prompt summed with the code wiped off. That
# sum was correlated with a pre-sync replica, which carries secondary chip 0 on
# every block, so the post-sync prompt has the sync polarity times chip 0. The
# switch may be a half-cycle jump; for both this and the pilots without a
# secondary code, see "Carrier loop staging" in docs/src/loop_filter.md.
@inline function _sync_polarity(
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer,
    prn::Integer,
)
    _is_wiped_off(signal, bit_buffer.found) && get_secondary_code_length(signal) > 1 ||
        return Int8(0)
    chip0 = GNSSSignals.secondary_value(get_secondary_code(signal), prn, 0)
    Int8(bit_buffer.polarity * sign(chip0))
end
