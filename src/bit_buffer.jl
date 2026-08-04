"""
    SyncResult

Outcome of a per-signal bit-sync / secondary-code-sync detector call.

Fields:

  - `found::Bool` — whether the detector locked on this update.
  - `phase::Int` — when `found = true`, the secondary-code chip the
    *upcoming* integration aligns to, in `0:secondary_code_length-1`
    (recovered by the hard rotation search in [`_secondary_code_search`](@ref)).
    Zero for signals without a secondary code (the L1 C/A bit-edge case
    fires at the data-bit boundary, where the upcoming integration starts a
    new bit, not at a secondary-chip offset) and also for the soft
    [`_detect_secondary_code_cfar`](@ref), which only fires at the winning
    rotation's own period boundary, so the upcoming integration always starts
    at secondary chip 0.
  - `polarity::Int8` — `+1` or `-1`; which match orientation the detector
    locked. Carries through to the post-sync prompt accumulator so that a
    negative-polarity lock doesn't trip the downstream bit decoder.
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
freedom, i.e. `t` such that `P(T ≤ t) = probability`. Used by
[`_detect_bit_edge_cfar`](@ref) in place of a standard-normal quantile as a
**small-sample penalty**, not because the detector's z-score is exactly
Student-t distributed. The z-score there divides the energy gap by a standard
error built from a variance *estimated* over `peak_bin_count` bins; a normal
threshold treats that estimate as exact and is badly overconfident when the
count is tiny (its ~3.9 threshold is a ~10% per-phase false alarm at one d.o.f.,
not the intended ~1e-4), which let a variance-collapsed 1–2-bin fluctuation lock
(JuliaGNSS/Tracking#124 fixed a 1-block hard-detector miss; this closes the soft
detector's own few-bin false-lock). The t-quantile is the right *shape* of
correction — it grows steeply as `dof → 1` and relaxes to the normal quantile as
`dof → ∞` (so mature many-bin locks are unaffected, no cutoff needed) — but its
nominal `dof = peak_bin_count − 1` is a heuristic, **not** the true degrees of
freedom: the per-bin values are χ² (non-central) energies rather than Gaussian,
and the competing phases share blocks (correlated), so the exact sampling
distribution is neither normal nor Student-t and the realised false-alarm rate
is only approximately the nominal one. It is used as a conservative,
integration-forcing penalty, not an exact calibration.

# Implementation

Hill's algorithm (Hill, G. W. (1970), *Algorithm 396: Student's t-quantiles*,
Comm. ACM 13(10), 619–620), driven by the two-tailed probability
`2·(1 − probability)`: `dof = 1` (Cauchy) and `dof = 2` invert in elementary
functions and are returned exactly, and above that a series in either the
normal deviate [`_norm_quantile`](@ref) (the near-normal branch, `y > 0.05 + a`)
or the two-tailed probability itself (the deep-tail branch, `y ≤ 0.05 + a`) is
used, where `y = (d·two_tailed)^(2/dof)`. The lower tail is reflected by
symmetry, and the median `t(0.5) = 0` is returned directly — which also avoids
a `1 − 0.5` self-recursion. The sole caller passes `probability > 0.5` (and
`dof ≥ 1`, since `peak_bin_count ≥ 2` at the call site).

Hill's series is accurate to ≲2e-4 *relative* — worst around `dof ≈ 3–10`,
better than 1e-5 by `dof ≈ 50` — which is orders of magnitude finer than the
nominal-d.o.f. modelling error above, and far too fine to move which block the
detector locks on. That error buys both of the properties this call site needs,
against an exact inversion through `SpecialFunctions.beta_inc_inv`:

  - **Allocation-free.** `beta_inc_inv` allocates internal scratch arrays
    (`zeros(22)` / `zeros(31)` inside the incomplete-beta asymptotic
    expansions) over part of the `dof` range the detector sweeps. This
    function runs once per primary-code block for every unsynced satellite, so
    an allocation here makes `track!` allocate in proportion to the signal
    length across the whole acquisition — the regression
    `test/track_in_place.jl`'s pre-sync guard exists to catch.
  - **Cheap.** `beta_inc_inv` is an iterative Newton solve that evaluates a
    full regularized incomplete beta per step (~0.2–2.8 µs here); Hill's is a
    single closed-form pass (~3–78 ns), 10–170× faster over the whole `dof`
    range and never slower. End to end that is a pre-sync `track!` at 1.7× (200
    blocks) to 1.3× (900 blocks).
"""
function _t_quantile(probability::Float64, dof::Real)
    probability == 0.5 && return 0.0
    probability < 0.5 && return -_t_quantile(1 - probability, dof)
    # Hill's algorithm is written in terms of the two-tailed probability.
    two_tailed_probability = 2 * (1 - probability)
    # `dof = 1` is the standard Cauchy, `quantile(p) = cot(π·two_tailed/2)`
    # (the cotangent form keeps its accuracy in the far tail, where the
    # equivalent `tan(π(probability − ½))` cancels). At `dof = 2` the CDF
    # `½·(1 + t/√(2 + t²))` inverts in closed form. Both are exact, and both
    # are where the series below is weakest.
    dof == 1 &&
        return cos(two_tailed_probability * pi / 2) / sin(two_tailed_probability * pi / 2)
    dof == 2 && return sqrt(2 / (two_tailed_probability * (2 - two_tailed_probability)) - 2)
    a = 1 / (dof - 0.5)
    b = 48 / a^2
    c = ((20700 * a / b - 98) * a - 16) * a + 96.36
    d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * pi / 2) * dof
    y = (d * two_tailed_probability)^(2 / dof)
    if y > 0.05 + a
        # Near-normal branch (large `y`): correct the normal deviate (the
        # `dof → ∞` limit) by an asymptotic series in `1/dof`. This is where
        # mature (large-dof) locks land — `y → 1` as the bin count grows — so
        # the threshold relaxes onto `Φ⁻¹` with no cutoff.
        x = _norm_quantile(probability)
        y = x^2
        dof < 5 && (c += 0.3 * (dof - 4.5) * (x + 0.6))
        c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c
        y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x
        y = expm1(a * y^2)
    else
        # Deep-tail branch (small `y`): a series in the two-tailed probability
        # itself, used at few-bin / small-dof arguments where the normal deviate
        # is a poor starting point.
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
detectors — one entry per candidate timing hypothesis. It backs **both** soft
detectors (a signal uses at most one): the GPS L1 C/A bit-edge detector
[`_detect_bit_edge_cfar`](@ref), where a hypothesis is an edge phase
`phase ∈ 0:blocks_per_bit-1` and the bin is one navigation bit, updated by
[`_update_phase_accumulators!`](@ref); and the secondary-code detector
[`_detect_secondary_code_cfar`](@ref), where a hypothesis is an overlay rotation
`d ∈ 0:N-1` and the bin is one (overlay-wiped) secondary-code period, updated by
[`_update_secondary_accumulators!`](@ref). In both cases it is advanced one
primary-code block at a time so detection stays O(hypotheses) per block with no
growing pre-sync history. Below, `period` is `blocks_per_bit` (L1 C/A) or the
secondary-code length `N`.

The vectors are mutated in place across the immutable [`BitBuffer`](@ref)
reconstructions (the same pattern as `soft_bits`). This is deliberate: the
accumulators are streaming state advanced every block, and an immutable
representation would either inline ~660 B into every per-block `BitBuffer`
copy (≈6× the heap traffic, dead weight post-sync) or box a fresh value
each block — both reintroduce the per-block allocation this design exists
to avoid. A single shared, in-place-updated buffer per satellite is the
allocation-free choice.

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
maximum-energy sync detectors — both the GPS L1 C/A bit-edge detector
[`_detect_bit_edge_cfar`](@ref) and the secondary-code detector
[`_detect_secondary_code_cfar`](@ref) route through here. Given the
per-hypothesis running statistics carried in [`PhaseAccumulators`](@ref) it
identifies the maximum-energy hypothesis, its closest competitor, and decides
whether the peak is significant enough to lock. The two detectors differ only in
what a "hypothesis" is (a bit-edge phase vs. a secondary-code rotation) and how
the accumulators are fed; the decision below is identical.

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
own completed-bin energies (Welford, numerically stable at any bin count). The
true hypothesis's bins vary only with thermal noise and slow drift — exactly the
run-to-run spread the test must compare the gap against. The peak is accepted
only when it beats the runner-up by a margin significant under that spread:

    z_score = energy_gap / standard_error   ≥   t⁻¹(1 - false_alarm_probability/(period - 1);  ν = peak_bin_count − 1)

where the standard error combines the peak's per-bin energy variance over the
peak and runner-up bin counts, `false_alarm_probability = 1 - confidence` is
Bonferroni-split over the `period - 1` competing hypotheses, and the quantile is
the Student-t inverse-CDF [`_t_quantile`](@ref) at a nominal `ν = peak_bin_count − 1` d.o.f. — a small-sample penalty for dividing by a variance estimated over
that few bins, not a claim that `z_score` is exactly Student-t (the per-bin
energies are χ² and the hypotheses correlated; see [`_t_quantile`](@ref)). A real
peak has a structural gap that dwarfs the thermal bin-to-bin spread, so `z_score`
grows like the square root of the bin count and crosses the threshold sooner at
high C/N₀ and later in noise — the detector self-paces — while a drift-only
asymmetry keeps `z_score` bounded and never locks.
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
    # runner-up is the higher of those two that isn't the peak. A hypothesis's
    # completed-bin count is `div(num_blocks - h, period)`.
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

    # Noise scale: the bin-to-bin variance of the peak hypothesis's own
    # completed-bin energies, maintained by Welford so it is numerically stable
    # at any bin count. The true hypothesis has pure (non-straddling) bins, so
    # this is the genuine run-to-run spread — thermal noise *and* slow drift /
    # quantization — which a within-bin residual would miss and then mistake a
    # tiny systematic cross-hypothesis asymmetry for a significant peak. The
    # runner-up is assumed to share this per-bin variance (it can only be larger,
    # from its straddling bins, so this is the optimistic — fastest-locking —
    # choice that still rejects drift-only asymmetries).
    @inbounds bin_energy_variance =
        bin_energy_sum_of_squared_deviations[peak_index+1] / (peak_bin_count - 1)
    standard_error =
        sqrt(bin_energy_variance * (1 / peak_bin_count + 1 / runner_up_bin_count))
    z_score =
        standard_error > 0 ? energy_gap / standard_error : (energy_gap > 0 ? Inf : 0.0)

    # `confidence` is the (1 - false-sync probability) target, Bonferroni-split
    # over the `period - 1` competing hypotheses. Clamp the quantile argument to
    # the open interval (0, 1): otherwise `confidence = 1.0` (or the rounding of
    # `1 - tiny/(period-1)` up to 1.0) yields a `quantile(1) = NaN`, and the
    # `z_score < threshold` gate would silently pass (NaN comparisons are false),
    # turning "maximum confidence" into "lock immediately".
    false_alarm_probability = 1 - confidence
    quantile_argument =
        clamp(1 - false_alarm_probability / (period - 1), nextfloat(0.0), prevfloat(1.0))
    # `standard_error` divides by a variance ESTIMATED from `peak_bin_count`
    # bins, so a normal threshold is badly overconfident at 1-2 bins (its ~3.9
    # value is a ~10% per-hypothesis false alarm at 1 d.o.f., not the intended
    # 1e-4), which let a variance-collapsed 2-bin fluctuation lock. Apply the
    # Student-t quantile at a nominal `peak_bin_count - 1` d.o.f. as a
    # small-sample penalty (NOT because `z_score` is exactly Student-t: the
    # per-bin energies are χ² and the hypotheses correlated, so the nominal
    # d.o.f. and false-alarm rate are approximate — see `_t_quantile`). It
    # demands z ≈ 6000 at 1 d.o.f. and relaxes to the normal threshold as bins
    # accumulate, forcing enough integration before a lock (which also lets the
    # carrier settle) while leaving mature clean locks (many bins ⇒ t ≈ normal)
    # unchanged.
    z_threshold = _t_quantile(quantile_argument, peak_bin_count - 1)
    z_score < z_threshold && return (false, peak_index, peak_bin_count)

    (true, peak_index, peak_bin_count)
end

"""
$(SIGNATURES)

Soft-decision, CFAR bit-edge detector — a signal-agnostic maximum-energy
timing synchronizer for any signal whose navigation bit spans more than one
primary-code period with no secondary code (selected by
[`uses_soft_bit_edge_detection`](@ref); GPS L1 C/A is the current example
with `blocks_per_bit = 20`). Signals with a periodic secondary code (GPS
L5I, GPS L1C-P) instead use the soft [`_detect_secondary_code_cfar`](@ref) or
the hard [`_secondary_code_search`](@ref), which correlate against a known
overlay.

The per-phase bin statistics are carried in `accumulators`
([`PhaseAccumulators`](@ref), advanced in place by
[`_update_phase_accumulators!`](@ref)); `blocks_per_bit` is the number of
primary-code blocks per navigation bit (20 for L1 C/A) and `num_blocks` is
the total number of prompts seen. The maximum-energy hypothesis test — peak
vs. runner-up, Welford noise scale, and Student-t CFAR threshold — is the
shared [`_cfar_decide`](@ref) core (the hypothesis here is the bit-edge
`phase ∈ 0:blocks_per_bit-1`). The detector is O(blocks_per_bit) per call — no
rescan of history — so it stays cheap on the pre-sync hot path no matter how
long sync takes.

# Boundary firing

Detection is reported (`found = true`) only when the most recent block also
*ends* the winning phase's bit (`num_blocks % blocks_per_bit == peak_phase`),
so the upcoming integration starts a fresh navigation bit. This keeps the
`phase = 0` contract — and, crucially, the post-sync coherent
`blocks_per_bit`-block integration — aligned to the true bit grid, which is
what makes the off-by-one lock of issue #124 structurally impossible: the
detector can only ever fire at the energy-maximizing phase's own boundary,
never at a neighbour's.

`polarity` is the sign of the most recently completed bin's coherent sum.
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
rotation search of period `N = secondary_code_length`. This is the
secondary-code analog of [`_update_phase_accumulators!`](@ref): each rotation
hypothesis `d ∈ 0:N-1` is a candidate secondary-chip alignment whose bins span
`N` blocks starting at index `d`. Unlike the bit-edge case, each block is
multiplied by the *known* secondary chip before being summed, so the correct
rotation wipes the overlay and every block in the bin adds coherently
(`|Σ p·s|²` large) while wrong rotations straddle the fixed overlay transitions
and lose energy.

`signal` / `prn` select the secondary code via `secondary_value`; under
hypothesis `d`, block `i` carries secondary chip `mod(i - d, N)` (`±1`), so the
bin completes when `chip == N - 1` and chip 0 is anchored to the physical
secondary-code chip 0 — the same overlay the replica applies post-sync, and the
same rotation the hard [`_secondary_code_search`](@ref) locks (its bespoke packed
reference differs only by an overall polarity, which the energy statistic here is
invariant to). Each completed bin's energy `|bin_sum|²` is folded into that
rotation's Welford mean / sum-of-squared-deviations and its polarity recorded.
O(N) per call, allocation-free.
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

Soft-decision, CFAR secondary-code sync detector — the maximum-energy analog of
[`_detect_bit_edge_cfar`](@ref) for signals carrying a short periodic secondary /
overlay code (selected by [`uses_soft_secondary_code_detection`](@ref); GPS L5I,
GPS L5Q, Galileo E1C / E5aI / E5aQ). The per-rotation bin statistics are carried
in `accumulators` ([`PhaseAccumulators`](@ref), advanced in place by
[`_update_secondary_accumulators!`](@ref)); `secondary_code_length` (`N`) is both
the bin length and the number of rotation hypotheses, and `num_blocks` is the
total number of prompts seen. The peak-vs-runner-up hypothesis test is the shared
[`_cfar_decide`](@ref) core.

Compared with the hard-decision rotation sweep [`_secondary_code_search`](@ref),
this uses the *soft* prompt magnitude and a CFAR confidence, so it rejects the
noise-driven false locks a short (e.g. 10-chip NH10) hard template match is prone
to and self-paces with C/N₀.

# Boundary firing

Detection is reported only when the most recent block also *ends* the winning
rotation's known-code period (`num_blocks % N == peak_rotation`), so the upcoming
integration starts at secondary chip 0 — the true chip-0 boundary, since each
rotation is anchored to the physical secondary chip 0 (see
[`_update_secondary_accumulators!`](@ref)). The reported `SyncResult.phase` is
therefore always `0`, and downstream code-phase snapping
([`_snap_code_phase_from_synced_signal`](@ref)) anchors on that boundary.
`polarity` is the sign of the winning period's coherent (overlay-wiped) sum
(resolved to the data-bit / carrier sign downstream by the navigation preamble).
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

Generic **hard-decision** secondary-code sync detector for signals on the
hard path — among the currently implemented signals only GPS L1C-P (its
1800-chip overlay) routes through here. The short-secondary-code signals
(GPS L5I/L5Q, Galileo E1C/E5aI/E5aQ) were moved to the soft, maximum-energy
[`_detect_secondary_code_cfar`](@ref) (see
[`uses_soft_secondary_code_detection`](@ref)); GPS L1 C/A uses the soft
bit-edge [`_detect_bit_edge_cfar`](@ref). Unlike those, this runs a full
rotation search against a *known* overlay code, so it locks after a single
secondary-code period in the worst case and recovers the true
secondary-code phase.

`received` is the prompt-sign sliding window with the newest block in
bit 0. `reference` is the secondary code packed in that **same**
newest-first order — bit `i` holds secondary chip `(N - 1 - i)` — so that
when the most recent `N` blocks span exactly one period ending on its
last chip, `received & mask == reference`. `N` is the secondary-code
length.

The search rotates the low `N` bits of `received` left by `d ∈ 0:N-1`,
tracking the best positive- and negated-polarity Hamming match in one
pass. The winning rotation `d` is how far the buffer leads the reference,
which maps to the secondary-chip offset of the **upcoming** integration
as `phase = mod(N - d, N)` — exactly the value the post-sync `code_phase`
snap ([`_snap_code_phase_from_synced_signal`](@ref)) anchors on. Returns
`SyncResult(false, 0, 0)` when the best distance exceeds `max_errors`.

Inlined so the per-signal `reference` / `N` constants fold at the call
site (the L1C-P window is 1800).
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

Shared detector body for signals that broadcast one channel symbol per
primary code period (GPS L1C-D, Galileo E1B): the buffer of primary-block
signs is itself the symbol stream — there is no sub-symbol boundary to
find — so the detector reports `found = true` from the very first
integration, leaving downstream consumers (GNSSDecoder.jl) to resolve the
residual ±1 polarity ambiguity via the navigation preamble.

A signal of this shape delegates its `detect_bit_or_secondary_code_sync`
method here.
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

Shared detector body for signals that lock onto a periodic secondary /
overlay code (GPS L5I's NH10, GPS L1C-P's 1800-chip overlay): wait until
the sliding `code_block_bits` window covers one full secondary-code
period, then run the [`_secondary_code_search`](@ref) rotation sweep
against the signal's packed reference ([`_packed_secondary_code`](@ref)).
The tolerance is a percentage of the secondary-code window, discretized
per signal as `floor(tolerance × N)` — see
[`get_bit_edge_or_secondary_code_tolerance`](@ref).

A new secondary-coded signal only needs a [`_packed_secondary_code`](@ref)
method (plus `get_code_block_buffer_type`) and a
`detect_bit_or_secondary_code_sync` method delegating here.
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

Hook for [`_detect_secondary_code_sync`](@ref): return the signal's
secondary / overlay code for `prn`, packed into the buffer type `B` in
the same newest-first order the prompt buffer fills — bit `i` holds
secondary chip `N - 1 - i`, so that when the most recent `N` blocks span
exactly one period ending on its last chip, `received & mask == reference` (see [`_secondary_code_search`](@ref)).
"""
function _packed_secondary_code end

# Generic packer: derive the reference directly from the signal's
# `SecondaryCode` (`get_secondary_code`), setting bit `N - 1 - k` iff
# secondary chip `k` is positive. Works for both `SharedSecondaryCode`
# (GPS L5Q's NH20, Galileo E1C's CS25 / E5aI's CS20) and
# `PerPRNSecondaryCode` (Galileo E5aQ's per-PRN CS100), so a new
# secondary-coded signal needs no bespoke method — it only picks a wide
# enough `get_code_block_buffer_type`. The rotation search in
# [`_secondary_code_search`](@ref) tries both polarities, so the absolute
# ±1 → bit convention here only sets the reported `polarity` sign, not
# whether the lock is found. GPS L5I and L1C-P keep their own (more
# specific) methods, which win dispatch.
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

The `code_block_buffer` field is the sync-search sliding window — its
width `B` is chosen per signal by [`get_code_block_buffer_type`](@ref) so
that a single integer can hold the entire pre-sync search horizon (one
NH10 period for GPS L5I, 40 primary blocks for GPS L1 C/A, 1800 chips
for the GPS L1C-P overlay, etc.). After sync the field is dead state and
the decoded navigation bits accumulate as **soft bits** in `soft_bits`
(one polarity-corrected coherent prompt sum per bit, sign = the hard
bit). Soft bits are the only bit store — decoders take them directly
(soft-decision Viterbi for Galileo E1B, LDPC for GPS L1C-D, confidence
weighting everywhere else), and a hard decision is just `soft_bit > 0`.
Because the store is a growable vector rather than the fixed `UInt128`
it used to be, there is no limit on how many bits may accumulate between
resets.

The `phase_acc` field holds the incremental per-hypothesis bin statistics
([`PhaseAccumulators`](@ref)) consumed by whichever soft CFAR sync detector the
signal uses — [`_detect_bit_edge_cfar`](@ref) for signals whose
[`uses_soft_bit_edge_detection`](@ref) is `true` (GPS L1 C/A), or
[`_detect_secondary_code_cfar`](@ref) for those whose
[`uses_soft_secondary_code_detection`](@ref) is `true` (GPS L5I/L5Q, Galileo
E1C/E5aI/E5aQ). It is seeded and updated only for those signals; for all others
(hard-decision path) it stays empty. Its size is bounded (one entry per
hypothesis), so there is no growing pre-sync history.
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

# Default constructor preserves the pre-refactor `UInt128`-backed search
# buffer. Once `get_code_block_buffer_type` lands (Step 2) the per-signal
# `TrackedSignal` constructor picks the right width instead.
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

# Convenience outer constructor for the pre-sync / post-sync search state
# without the phase / polarity arguments (assumed zero). Used by test and
# benchmark code that builds a `BitBuffer` from raw integer / Complex{Int}
# literals. `soft_bits` seeds the decoded-bit store; it is aliased, not
# copied, so the caller can observe in-place pushes.
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

# Get the soft bits, i.e. the accumulated (summed) filtered prompt of each
# completed bit. This is the only decoded-bit store; a hard decision is the
# sign, `soft_bit > 0`. Bits recovered from the pre-sync sign window at
# bit-sync time only have ±1 prompt signs available; their sign-vote sum is
# scaled by the sync-time prompt magnitude so the magnitudes stay comparable
# (as reliabilities) with the coherently accumulated post-sync bits. The
# buffer is reset to length 0 at the start of each `track` call. Kept as a
# plain comment (not a docstring) to match the sibling accessor
# `get_num_bits`, which `checkdocs = :exports` would otherwise require to
# appear in the manual.
@inline get_soft_bits(bit_buffer::BitBuffer) = bit_buffer.soft_bits

"""
$(SIGNATURES)

Width `B` of the packed prompt-sign buffer (`BitBuffer.code_block_buffer`) for
`signal`, returned as a concrete `Unsigned` subtype.

That packed buffer is the sliding-window search horizon **only for the
hard-decision path** — the rotation/Hamming sweep [`_secondary_code_search`](@ref),
which among the currently implemented signals is used by GPS L1C-P alone. The
soft-decision CFAR detectors ([`_detect_bit_edge_cfar`](@ref) for GPS L1 C/A,
[`_detect_secondary_code_cfar`](@ref) for the short-secondary-code signals) read
the incremental [`PhaseAccumulators`](@ref) instead, so for those signals the
packed buffer of this width is built but not consulted for detection — it is
vestigial (the width could be `UInt8`; it is left at the horizon width below for
uniformity and so the hard path stays available). The table gives, per signal,
the returned width and — for the hard path — the horizon it must hold:

| Signal        | Returns    | Detector / buffer role                              |
|:------------- |:---------- |:--------------------------------------------------- |
| GPS L1 C/A    | `UInt64`   | soft bit-edge CFAR — packed buffer vestigial        |
| Galileo E1B   | `UInt8`    | symbol = primary period, buffer unused              |
| GPS L5I       | `UInt32`   | soft secondary CFAR — packed buffer vestigial       |
| GPS L5Q       | `UInt32`   | soft secondary CFAR — packed buffer vestigial       |
| GPS L1C-D     | `UInt8`    | symbol = primary period, buffer unused              |
| GPS L1C-P     | `UInt1800` | **hard** rotation sweep — 1800-chip overlay horizon |
| GPS L2CM      | `UInt8`    | symbol = primary period, buffer unused              |
| GPS L2CL      | `UInt8`    | dataless pilot, no sync, buffer unused              |
| Galileo E1C   | `UInt32`   | soft secondary CFAR — packed buffer vestigial       |
| Galileo E5a-I | `UInt32`   | soft secondary CFAR — packed buffer vestigial       |
| Galileo E5a-Q | `UInt128`  | soft secondary CFAR — packed buffer vestigial       |

The default for any signal not specialized below is `UInt64`. The width
flows through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF, CN0}` so the
parameter chain stays type-stable at construction.
"""
@inline get_code_block_buffer_type(::AbstractGNSSSignal) = UInt64

"""
$(SIGNATURES)

Per-signal Hamming tolerance used by the **hard-decision** rotation/Hamming
sweep [`_secondary_code_search`](@ref), expressed as a fraction of the search
window.

Returns the largest **fraction** of bit-flips the per-signal
`detect_bit_or_secondary_code_sync` accepts before reporting
`found = true`. Each detector converts this to an integer error budget
at its call site: `max_errors = floor(Int, tolerance × window_size)`.

Among the currently implemented signals **only GPS L1C-P reads this trait** — it
is the sole signal still on the hard path. The short-secondary-code signals (GPS
L5I/L5Q, Galileo E1C/E5aI/E5aQ) were moved to the soft, confidence-driven
[`_detect_secondary_code_cfar`](@ref) (selected by
[`uses_soft_secondary_code_detection`](@ref)) and no longer consult it; likewise
GPS L1 C/A uses [`_detect_bit_edge_cfar`](@ref). Both soft detectors are tuned by
[`get_bit_edge_detection_confidence`](@ref) instead. Galileo E1B and GPS L1C-D
broadcast one channel symbol per primary code period, so their detectors return
`SyncResult(true, 0, +1)` unconditionally — the trait default applies but the
value is ignored — and GPS L2CL is a dataless pilot with no sync.

Default is `0.025` (2.5 %). At L1C-P's 1800-chip window that discretizes to
`max_errors = 45`. (For reference, the same 2.5 % at the now-soft signals' short
windows would floor to 0–2 errors — an exact or near-exact match — which is
exactly the noise-driven false-lock exposure the move to the soft detector
removed.)

# Overriding

To loosen the tolerance for low-C/N₀ work, dispatch the trait on the (hard-path)
signal type in your own module:

```julia
Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL1C_P) = 0.05
```

The override takes effect at the next call to
`detect_bit_or_secondary_code_sync` — there is no need to rebuild any
TrackState. The trait is `@inline`'d so the override folds at the
detector's call site. (Overriding it for a soft-detector signal has no effect;
tune [`get_bit_edge_detection_confidence`](@ref) there instead.)
"""
@inline get_bit_edge_or_secondary_code_tolerance(::AbstractGNSSSignal) = 0.025

"""
$(SIGNATURES)

Whether `signal`'s bit edge is located with the soft-decision,
maximum-energy CFAR detector [`_detect_bit_edge_cfar`](@ref) (which reads
the incremental [`PhaseAccumulators`](@ref)) rather than the hard-decision
`detect_bit_or_secondary_code_sync` path.

This is **signal-agnostic**: the detector and accumulators are
parameterised by the number of primary-code blocks per navigation bit
(`L`, from [`_calc_num_code_blocks_that_form_a_bit`](@ref)), with no
per-signal constants. The default enables it for any signal whose
navigation bit spans **more than one** primary-code period **and** which
carries no secondary/overlay code — i.e. the bit edge is a sub-bit timing
offset to be found, not a symbol boundary that is already aligned (Galileo
E1B, GPS L1C-D: one symbol per primary period) and not a periodic overlay
(GPS L5I, L1C-P: located by [`_secondary_code_search`](@ref)). Among the
currently implemented signals only GPS L1 C/A (20 blocks/bit) qualifies,
but a newly added signal with the same structure is picked up
automatically.

Override per signal type to force the choice, e.g. to disable it:

```julia
Tracking.uses_soft_bit_edge_detection(::SomeSignal) = false
```

The result is constant-folded per signal type, so the branch in
[`_buffer_find_bit`](@ref) compiles away and signals that don't use it
never seed or update `phase_acc`.
"""
@inline uses_soft_bit_edge_detection(signal::AbstractGNSSSignal) =
    _calc_num_code_blocks_that_form_a_bit(signal) > 1 &&
    get_secondary_code_length(signal) == 1

"""
$(SIGNATURES)

Whether `signal`'s secondary/overlay code is located with the soft-decision,
maximum-energy CFAR detector [`_detect_secondary_code_cfar`](@ref) (which reads
the incremental [`PhaseAccumulators`](@ref)) rather than the hard-decision
rotation/Hamming sweep [`_secondary_code_search`](@ref).

The soft detector coherently integrates one *full secondary-code period* per bin,
so it only makes sense while a whole period is a phase-coherent integration
length. The default therefore enables it for signals with a **short** secondary
code — `1 < get_secondary_code_length(signal) ≤ 100` — which covers GPS L5I
(NH10, 10), GPS L5Q (NH20, 20), Galileo E1C (CS25, 25), E5aI (CS20, 20) and E5aQ
(CS100, 100), whose periods are ≤ 100 ms. GPS L1C-P is deliberately excluded: its
1800-chip overlay is an 18 s period, far too long to integrate coherently (and
its long code is not false-lock-prone), so it keeps the hard
[`_secondary_code_search`](@ref).

Because it locates a periodic overlay, this is mutually exclusive with
[`uses_soft_bit_edge_detection`](@ref) (which requires *no* secondary code); a
signal routes to at most one soft detector.

Override per signal type to force the choice, e.g. to disable it:

```julia
Tracking.uses_soft_secondary_code_detection(::GPSL5I) = false
```

The result is constant-folded per signal type, so the branch in
[`_buffer_find_bit`](@ref) compiles away and signals that don't use it never
seed or update `phase_acc`.
"""
@inline uses_soft_secondary_code_detection(signal::AbstractGNSSSignal) =
    1 < get_secondary_code_length(signal) <= 100

"""
$(SIGNATURES)

Target confidence (one minus the probability of a false lock) for the
soft-decision CFAR sync detectors — the GPS L1 C/A bit-edge detector
[`_detect_bit_edge_cfar`](@ref) and the secondary-code detector
[`_detect_secondary_code_cfar`](@ref) both read it.

Default `0.999`: the detector keeps integrating primary-code blocks until
the maximum-energy hypothesis beats its closest competitor with this
confidence, so a clean signal locks in as little as two bins
while a noisy one self-paces to as long as it takes. Lower it to lock
faster at the cost of more false locks; raise it to be more conservative.

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
"""
function buffer(
    signal::AbstractGNSSSignal,
    prn::Integer,
    bit_buffer::BitBuffer{B},
    integrated_code_blocks,
    prompt,
) where {B<:Unsigned}
    # The divide is deferred to the helper — pilot signals
    # (`get_data_frequency = 0`) would otherwise blow up here with `Int(Inf)`.
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

    # Pilot signals (e.g. GPS L1C-P) carry no navigation data; once their
    # secondary code is synced there is nothing to decode, so the buffer is
    # left untouched (the `found` / `secondary_phase` / `polarity` state stays
    # for code-phase anchoring and longer coherent integration).
    num_code_blocks_that_form_a_bit == 0 && return bit_buffer

    prompt_accumulator = bit_buffer.prompt_accumulator + prompt
    prompt_accumulator_integrated_code_blocks =
        bit_buffer.prompt_accumulator_integrated_code_blocks + integrated_code_blocks

    if prompt_accumulator_integrated_code_blocks == num_code_blocks_that_form_a_bit
        # Flip the decoded bit if the detector locked at negative polarity:
        # the prompt accumulator's real-part sign is then inverted relative
        # to the data symbol's "0/1" convention.
        bit_acc =
            bit_buffer.polarity < 0 ? -real(prompt_accumulator) : real(prompt_accumulator)
        # The polarity-corrected accumulation IS the bit: its sign is the hard
        # decision, its magnitude the confidence. The vector grows without a
        # ceiling — the caller decides how often to read bits out and `reset`.
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

    # Signals that detect the sync point from soft prompts fold each block into
    # the per-hypothesis accumulators and run a maximum-energy CFAR detector:
    # GPS L1 C/A over bit-edge phases, and short-secondary-code signals (GPS
    # L5I/L5Q, Galileo E1C/E5aI/E5aQ) over overlay rotations. Everything else
    # stays on the hard-decision sliding-window path. Both soft branches share
    # the single `phase_acc` buffer (a signal is at most one of them). The
    # branch folds at compile time per signal type, so non-soft signals never
    # touch `phase_acc`.
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
        # Secondary-code signals (GPS L5I, GPS L1C-P): the buffered pre-sync
        # prompt signs are modulated by the secondary code, not the navigation
        # data, so there are no data bits to recover from them. Data-bit
        # decoding starts fresh post-sync. The rotation search locks at *any*
        # secondary chip, so the upcoming integration starts `sync.phase` blocks
        # into the current data bit — not at its boundary. Seed the accumulator
        # block count with `sync.phase` so the first emitted bit completes
        # exactly `num_code_blocks_that_form_a_bit − sync.phase` blocks later, on
        # the data-bit boundary, instead of `num_code_blocks_that_form_a_bit`
        # blocks after the lock instant (issue #125). For pilots
        # (`num_code_blocks_that_form_a_bit == 0`) the post-sync `buffer` path
        # returns early and never consults this seed, so it is harmless there.
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
    num_bits = min(
        div(code_block_buffer_length, num_code_blocks_that_form_a_bit),
        div(sizeof(code_block_buffer) * 8, num_code_blocks_that_form_a_bit),
    )
    # Hoist the one field the closure needs (`sync.polarity`) into a plain
    # local so the closure below does not capture `sync` itself. `sync` is
    # assigned in two branches above (the soft vs. hard detector), and
    # capturing a variable that is assigned in more than one place forces
    # Julia to box it: a per-call `Core.Box` plus heap-boxing of each
    # `SyncResult`. That fires on every code block until sync, so it shows up
    # as a per-code-block allocation on the pre-sync hot path (~80 B/block),
    # making `track!` allocate in proportion to the signal length instead of
    # staying allocation-free after warmup. Capturing the plain `Int` keeps
    # `sync` unboxed.
    sync_polarity = Int(sync.polarity)
    for bit_index = num_bits:-1:1     # oldest recovered bit first
        # Apply the lock polarity to the buffered pre-sync bits as well, so
        # they map symbol levels to bit values the same way as every
        # post-sync bit (which is sign-flipped via `bit_buffer.polarity` in
        # `buffer`). Otherwise a negative-polarity lock would invert only
        # the post-sync part of the bit stream (issue #127).
        bit_sum =
            sum(0:(num_code_blocks_that_form_a_bit-1)) do code_block_index
                buffer_code_block_index =
                    (bit_index - 1) * num_code_blocks_that_form_a_bit + code_block_index
                ((code_block_buffer & (one(B) << buffer_code_block_index)) > 0) * 2 - 1
            end * sync_polarity
        # The pre-sync window only stores prompt signs, so `bit_sum` is a
        # ±1-per-block vote count (already polarity-corrected above). Scale it
        # by the sync-time prompt magnitude (the best available amplitude
        # estimate) so these recovered soft bits live in the same
        # coherent-amplitude-sum units as the bits accumulated post-sync.
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

`secondary_phase` is the secondary chip the **upcoming** integration aligns to,
and it is read exactly once: by the code-phase snap
([`_snap_code_phase_from_synced_signal`](@ref)), which runs after the whole
chunk has been folded. The detector reports it for the block right after the
syncing record, so every further record folded in the same chunk — the ones
`_apply_correlator_output` marks `correlated_pre_sync`, correlated with the
pre-sync (un-wiped) replica — moves that anchor along by the blocks it covers.
Without this the snap anchors the replica `num_code_blocks` chips early, the
post-sync replica bakes the wrong overlay chip into every block, and the
coherent bit sum collapses to a fraction of its length — the secondary-code
sibling of the bit-grid slip in issue #219.
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
