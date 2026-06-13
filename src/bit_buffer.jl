"""
    SyncResult

Outcome of a per-signal bit-sync / secondary-code-sync detector call.

Fields:

- `found::Bool` — whether the detector locked on this update.
- `phase::Int` — when `found = true`, the secondary-code chip the
  *upcoming* integration aligns to, in `0:secondary_code_length-1`
  (recovered by the rotation search in [`_secondary_code_search`](@ref)).
  Zero for signals without a secondary code (the L1 C/A bit-edge case
  fires at the data-bit boundary, where the upcoming integration starts a
  new bit, not at a secondary-chip offset).
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
SpecialFunctions.jl). Used by [`_detect_bit_edge_cfar`](@ref) to turn a
false-sync probability into a detection z-score threshold. Returns `±Inf`
at `probability = 1` / `0`; callers keep the argument in the open interval.
"""
@inline _norm_quantile(probability::Float64) = sqrt(2.0) * erfinv(2 * probability - 1)

"""
$(SIGNATURES)

Per-phase bin statistics for the soft, maximum-energy bit-edge detector
[`_detect_bit_edge_cfar`](@ref) — one entry per candidate edge phase
`phase ∈ 0:blocks_per_bit-1`, where `blocks_per_bit` is the (per-signal)
number of primary-code blocks per navigation bit. Updated one primary-code
block at a time by [`_update_phase_accumulators!`](@ref) so detection stays
O(blocks_per_bit) per block with no growing pre-sync history.

The vectors are mutated in place across the immutable [`BitBuffer`](@ref)
reconstructions (the same pattern as `soft_bits`). This is deliberate: the
accumulators are streaming state advanced every block, and an immutable
representation would either inline ~660 B into every per-block `BitBuffer`
copy (≈6× the heap traffic, dead weight post-sync) or box a fresh value
each block — both reintroduce the per-block allocation this design exists
to avoid. A single shared, in-place-updated buffer per satellite is the
allocation-free choice.

Fields (all length `blocks_per_bit` once seeded; empty before the first block):

- `open_bin_sum` — coherent sum of the phase's *currently open* bin.
- `mean_bin_energy` / `bin_energy_sum_of_squared_deviations` — Welford running mean and
  sum of squared deviations (`M₂ = Σ(energyᵢ − mean)²`) of the phase's
  *completed*-bin energies, for a numerically stable variance. The bin
  count is not stored — it is `div(num_blocks - phase, blocks_per_bit)`.
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
    accumulators::PhaseAccumulators, prompt::ComplexF64, block_index::Int, blocks_per_bit::Int,
)
    @inbounds for phase in 0:(blocks_per_bit-1)
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

Soft-decision, CFAR bit-edge detector — a signal-agnostic maximum-energy
timing synchronizer for any signal whose navigation bit spans more than one
primary-code period with no secondary code (selected by
[`uses_soft_bit_edge_detection`](@ref); GPS L1 C/A is the current example
with `blocks_per_bit = 20`). Signals with a periodic secondary code (GPS
L5I, GPS L1C-P) instead use [`_secondary_code_search`](@ref), which
correlates against a known overlay.

The per-phase bin statistics are carried in `accumulators`
([`PhaseAccumulators`](@ref), advanced in place by
[`_update_phase_accumulators!`](@ref)); `blocks_per_bit` is the number of
primary-code blocks per navigation bit (20 for L1 C/A) and `num_blocks` is
the total number of prompts seen. The detector is O(blocks_per_bit) per
call — no rescan of history — so it stays cheap on the pre-sync hot path no
matter how long sync takes.

# Statistic

Under the hypothesis that the bit edge sits at `phase ∈ 0:blocks_per_bit-1`,
the prompts partition into `blocks_per_bit`-block bins aligned to `phase`.
Each complete bin is coherently summed (`Sₖₘ = Σ pᵢ`) and its energy
`|Sₖₘ|²` folded into the phase's running `mean_bin_energy`; the true phase
never straddles a data-bit transition, so it keeps full coherent gain on
every bin while wrong phases lose energy on transition bins. The per-phase
mean bin energy is therefore the maximum-likelihood timing statistic for
unknown i.i.d. data in AWGN.

# CFAR confidence

The noise scale is the *bin-to-bin* sample variance of the winning
phase's own completed-bin energies. The true edge's bins never straddle a
transition, so their energies vary only with thermal noise and slow
drift — exactly the run-to-run spread the test must compare the gap
against. (A within-bin residual would read ~0 on a clean but drifting
signal and then mistake a tiny systematic cross-phase asymmetry for a real
edge.) The peak is accepted only when it beats the runner-up by a margin
significant under that spread:

    z_score = energy_gap / standard_error   ≥   Φ⁻¹(1 - false_alarm_probability/(blocks_per_bit - 1))

where the standard error combines the peak phase's per-bin energy variance
over the peak and runner-up bin counts, `false_alarm_probability =
1 - confidence` is Bonferroni-split over the `blocks_per_bit - 1` competing
phases, and `Φ⁻¹` is [`_norm_quantile`](@ref). A real edge has a structural
gap that dwarfs the thermal bin-to-bin spread, so `z_score` grows like the
square root of the bin count and crosses the threshold sooner at high C/N₀
and later in noise — the detector self-paces — while a drift-only asymmetry
keeps `z_score` bounded and never locks.

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
    accumulators::PhaseAccumulators, blocks_per_bit::Int, confidence::Float64, num_blocks::Int,
)
    # Need at least two complete bins on some phase before any phase can
    # serve as a runner-up; below that a transition cannot have been
    # localized yet.
    num_blocks < 2 * blocks_per_bit && return SyncResult(false, 0, Int8(0))

    mean_bin_energy = accumulators.mean_bin_energy
    # Single O(blocks_per_bit) pass: the peak (highest mean bin energy among
    # phases with ≥ 2 complete bins) and the two highest phases overall. The
    # runner-up is the higher of those two that isn't the peak. A phase's
    # completed-bin count is `div(num_blocks - phase, blocks_per_bit)`.
    peak_phase = -1; peak_energy = -1.0; peak_bin_count = 0          # best, ≥ 2 bins
    best_phase = -1; best_energy = -1.0; best_bin_count = 0          # highest, ≥ 1 bin
    second_best_energy = -1.0; second_best_bin_count = 0             # 2nd highest, ≥ 1 bin
    @inbounds for phase in 0:(blocks_per_bit-1)
        bin_count = div(num_blocks - phase, blocks_per_bit)
        bin_count < 1 && continue
        energy = mean_bin_energy[phase+1]
        if energy > best_energy
            second_best_energy = best_energy
            second_best_bin_count = best_bin_count
            best_energy = energy
            best_phase = phase
            best_bin_count = bin_count
        elseif energy > second_best_energy
            second_best_energy = energy
            second_best_bin_count = bin_count
        end
        if bin_count >= 2 && energy > peak_energy
            peak_energy = energy
            peak_phase = phase
            peak_bin_count = bin_count
        end
    end
    # `num_blocks >= 2 * blocks_per_bit` guarantees phase 0 has ≥ 2 bins, so
    # `peak_phase` is always set.
    peak_phase < 0 && return SyncResult(false, 0, Int8(0))

    # Runner-up: the highest-energy phase that isn't the peak (any bin count).
    runner_up_energy, runner_up_bin_count =
        best_phase == peak_phase ? (second_best_energy, second_best_bin_count) :
        (best_energy, best_bin_count)
    runner_up_bin_count < 1 && return SyncResult(false, 0, Int8(0))

    energy_gap = peak_energy - runner_up_energy
    energy_gap <= 0 && return SyncResult(false, 0, Int8(0))

    # Noise scale: the bin-to-bin variance of the peak phase's own
    # completed-bin energies, maintained by Welford so it is numerically
    # stable at any bin count. The true edge has pure (non-straddling) bins,
    # so this is the genuine run-to-run spread — thermal noise *and* slow
    # drift / quantization — which a within-bin residual would miss and then
    # mistake a tiny systematic cross-phase asymmetry for a significant edge.
    # The runner-up is assumed to share this per-bin variance (it can only be
    # larger, from its straddling bins, so this is the optimistic —
    # fastest-locking — choice that still rejects drift-only asymmetries).
    @inbounds bin_energy_variance =
        accumulators.bin_energy_sum_of_squared_deviations[peak_phase+1] / (peak_bin_count - 1)
    standard_error =
        sqrt(bin_energy_variance * (1 / peak_bin_count + 1 / runner_up_bin_count))
    z_score = standard_error > 0 ? energy_gap / standard_error : (energy_gap > 0 ? Inf : 0.0)

    # `confidence` is the (1 - false-sync probability) target, Bonferroni-
    # split over the `blocks_per_bit - 1` competing phases. Clamp the quantile
    # argument to the open interval (0, 1): otherwise `confidence = 1.0` (or
    # the rounding of `1 - tiny/(blocks_per_bit-1)` up to 1.0) yields
    # `Φ⁻¹(1) = NaN`, and the `z_score < threshold` gate would silently pass
    # (NaN comparisons are false), turning "maximum confidence" into "lock
    # immediately".
    false_alarm_probability = 1 - confidence
    quantile_argument = clamp(
        1 - false_alarm_probability / (blocks_per_bit - 1),
        nextfloat(0.0),
        prevfloat(1.0),
    )
    z_threshold = _norm_quantile(quantile_argument)
    z_score < z_threshold && return SyncResult(false, 0, Int8(0))

    # Fire only at the winning phase's own bit boundary so the upcoming
    # integration starts a new bit.
    num_blocks % blocks_per_bit != peak_phase && return SyncResult(false, 0, Int8(0))

    @inbounds polarity =
        accumulators.last_bin_polarity[peak_phase+1] < 0 ? Int8(-1) : Int8(+1)
    SyncResult(true, 0, polarity)
end

"""
$(SIGNATURES)

Generic secondary-code sync detector shared by every signal that locks
onto a periodic secondary / overlay code — GPS L5I's NH10 and GPS L1C-P's
1800-chip overlay both route through here. Unlike
[`_detect_bit_edge_cfar`](@ref) (the soft maximum-energy bit-edge
detector used for GPS L1 C/A), this runs a full rotation search against a
*known* overlay code, so it locks after a single secondary-code period in
the worst case and recovers the true secondary-code phase.

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
site (the L5I window is 10, the L1C-P window 1800).
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
    @inbounds for d in 0:(N-1)
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

BitBuffer to buffer bits.

The `code_block_buffer` field is the sync-search sliding window — its
width `B` is chosen per signal by [`get_code_block_buffer_type`](@ref) so
that a single integer can hold the entire pre-sync search horizon (one
NH10 period for GPS L5I, 40 primary blocks for GPS L1 C/A, 1800 chips
for the GPS L1C-P overlay, etc.). After sync the field is dead state and
the decoded navigation bits accumulate in the fixed-width
`buffer::UInt128` instead.

The `phase_acc` field holds the incremental per-phase bin statistics
([`PhaseAccumulators`](@ref)) consumed by the soft bit-edge detector
[`_detect_bit_edge_cfar`](@ref). It is seeded and updated only for signals
whose [`uses_soft_bit_edge_detection`](@ref) trait is `true` (currently
GPS L1 C/A); for all other signals it stays empty. Its size is bounded
(one entry per phase), so there is no growing pre-sync history.
"""
struct BitBuffer{B<:Unsigned}
    code_block_buffer::B
    code_block_buffer_lengh::Int
    found::Bool
    secondary_phase::Int      # 0 until found; secondary-chip offset post-sync
    polarity::Int8            # +1 or -1 once found; 0 before sync
    buffer::UInt128
    length::Int
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
        zero(UInt128), 0, false, 0, Int8(0), zero(UInt128), 0, complex(0.0, 0.0), 0,
        Float32[], PhaseAccumulators(),
    )
end

# Typed empty constructor used by the per-signal `TrackedSignal` path.
function BitBuffer{B}() where {B<:Unsigned}
    BitBuffer{B}(
        zero(B), 0, false, 0, Int8(0), zero(UInt128), 0, complex(0.0, 0.0), 0,
        Float32[], PhaseAccumulators(),
    )
end

# Convenience outer constructor matching the legacy 7-arg form (no phase /
# polarity arguments — assumed zero; soft bits default to empty). Used by test
# code that builds a `BitBuffer` from raw integer / Complex{Int} literals.
function BitBuffer(
    code_block_buffer::B,
    code_block_buffer_lengh::Integer,
    found::Bool,
    buffer::Integer,
    length::Integer,
    prompt_accumulator::Complex,
    prompt_accumulator_integrated_code_blocks::Integer,
) where {B<:Unsigned}
    BitBuffer{B}(
        code_block_buffer,
        Int(code_block_buffer_lengh),
        found,
        0,
        Int8(0),
        UInt128(buffer),
        Int(length),
        ComplexF64(prompt_accumulator),
        Int(prompt_accumulator_integrated_code_blocks),
        Float32[],
        PhaseAccumulators(),
    )
end

@inline get_bits(bit_buffer::BitBuffer) = bit_buffer.buffer
@inline length(bit_buffer::BitBuffer) = bit_buffer.length
@inline has_bit_or_secondary_code_been_found(bit_buffer::BitBuffer) = bit_buffer.found

# Get the soft bits, i.e. the accumulated (summed) filtered prompt of each
# completed bit. The sign of each soft bit corresponds to the respective hard
# bit returned by `get_bits`. The buffer is reset to length 0 at the start of
# each `track` call, mirroring the hard bit buffer. Kept as a plain comment (not
# a docstring) to match the sibling accessors `get_bits` / `get_num_bits`, which
# `checkdocs = :exports` would otherwise require to appear in the manual.
@inline get_soft_bits(bit_buffer::BitBuffer) = bit_buffer.soft_bits

"""
$(SIGNATURES)

Width of the sync-search sliding window for `signal`, returned as a
concrete `Unsigned` subtype.

Each signal's bit/secondary-code sync detector matches the buffered
primary-code-block signs against a per-signal template. This trait
picks the smallest integer width that holds one full search horizon for
that signal:

| Signal            | Returns   | Window |
|-------------------|-----------|--------|
| GPS L1 C/A        | `UInt64`  | 40 blocks (2 × 20 blocks/symbol) |
| Galileo E1B       | `UInt8`   | 8 blocks (2 × 4 blocks/symbol)   |
| GPS L5I           | `UInt32`  | 20 blocks (2 × NH10 length)      |
| GPS L1C-D         | `UInt8`   | n/a — symbol = primary period, buffer unused |
| GPS L1C-P         | `UInt1800`| 1800-chip overlay (added in step 4) |

The default for any signal not specialized below is `UInt64`. The width
flows through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF}` so the
parameter chain stays type-stable at construction.
"""
@inline get_code_block_buffer_type(::AbstractGNSSSignal) = UInt64

"""
$(SIGNATURES)

Per-signal Hamming tolerance used by the secondary-code sync-search
detector [`_secondary_code_search`](@ref), expressed as a fraction of the
search window.

Returns the largest **fraction** of bit-flips the per-signal
`detect_bit_or_secondary_code_sync` accepts before reporting
`found = true`. Each detector converts this to an integer error budget
at its call site: `max_errors = floor(Int, tolerance × window_size)`.

Default is `0.025` (2.5 %), which discretizes per-signal as:

| Signal      | Window (blocks) | Effective `max_errors` |
|-------------|-----------------|------------------------|
| GPS L5I     | 10              | 0 (exact match)        |
| GPS L1C-P   | 1800            | 45                     |
| Galileo E1B | n/a — trivial   | unused                 |
| GPS L1C-D   | n/a — trivial   | unused                 |

GPS L1 C/A does **not** use this trait: its bit-edge detector is the
soft-decision, confidence-driven [`_detect_bit_edge_cfar`](@ref), tuned
by [`get_bit_edge_detection_confidence`](@ref) instead. Galileo E1B and
GPS L1C-D broadcast one channel symbol per primary code period, so their
detectors return `SyncResult(true, 0, +1)` unconditionally — the trait
default applies but the value is ignored.

# Overriding

To loosen the tolerance for low-C/N₀ work, dispatch the trait on the
signal type in your own module:

```julia
Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL5I) = 0.05
```

The override takes effect at the next call to
`detect_bit_or_secondary_code_sync` — there is no need to rebuild any
TrackState. The trait is `@inline`'d so the override folds at the
detector's call site.
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

Target confidence (one minus the probability of a false bit-edge lock)
for the GPS L1 C/A soft bit-edge detector [`_detect_bit_edge_cfar`](@ref).

Default `0.999`: the detector keeps integrating primary-code blocks until
the maximum-energy phase beats its closest competitor with this
confidence, so a clean signal locks in as little as two bits (~40 ms)
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
function buffer(signal::AbstractGNSSSignal, prn::Integer, bit_buffer::BitBuffer{B}, integrated_code_blocks, prompt) where {B<:Unsigned}
    # The divide is deferred to the helper — pilot signals
    # (`get_data_frequency = 0`) would otherwise blow up here with `Int(Inf)`.
    num_code_blocks_that_form_a_bit = _calc_num_code_blocks_that_form_a_bit(signal)

    if (bit_buffer.found == false)
        return _buffer_find_bit(
            signal, prn, bit_buffer, num_code_blocks_that_form_a_bit,
            integrated_code_blocks, prompt,
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
        bit_acc = bit_buffer.polarity < 0 ?
                  -real(prompt_accumulator) : real(prompt_accumulator)
        bit = bit_acc > 0
        # Store the polarity-corrected accumulation so the soft bit's sign
        # matches the decoded hard bit.
        push!(bit_buffer.soft_bits, Float32(bit_acc))
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            bit_buffer.secondary_phase,
            bit_buffer.polarity,
            get_bits(bit_buffer) << 1 + UInt64(bit),
            length(bit_buffer) + 1,
            zero(prompt_accumulator),
            0,
            bit_buffer.soft_bits,
            bit_buffer.phase_acc,
        )
    else
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            bit_buffer.secondary_phase,
            bit_buffer.polarity,
            bit_buffer.buffer,
            bit_buffer.length,
            prompt_accumulator,
            prompt_accumulator_integrated_code_blocks,
            bit_buffer.soft_bits,
            bit_buffer.phase_acc,
        )
    end
end

function _buffer_find_bit(signal, prn::Integer, bit_buffer::BitBuffer{B}, num_code_blocks_that_form_a_bit,
                          integrated_code_blocks, prompt) where {B<:Unsigned}
    if (integrated_code_blocks != 1)
        error(
            "The number code blocks must be equal to 1 if bit or secondary code hasn't been found yet.",
        )
    end
    code_block_buffer = bit_buffer.code_block_buffer << 1 + B(real(prompt) > 0)
    code_block_buffer_lengh = bit_buffer.code_block_buffer_lengh + 1

    # Signals that detect the bit edge from soft prompts (GPS L1 C/A) fold
    # each block into the per-phase accumulators and run the maximum-energy
    # CFAR detector; everything else stays on the hard-decision sliding-window
    # path. The branch folds at compile time per signal type, so non-soft
    # signals never update `phase_acc`.
    phase_acc = bit_buffer.phase_acc
    if uses_soft_bit_edge_detection(signal)
        blocks_per_bit = num_code_blocks_that_form_a_bit
        _is_seeded(phase_acc, blocks_per_bit) ||
            _seed_phase_accumulators!(phase_acc, blocks_per_bit)
        _update_phase_accumulators!(
            phase_acc, ComplexF64(prompt), code_block_buffer_lengh - 1, blocks_per_bit,
        )
        sync = _detect_bit_edge_cfar(
            phase_acc,
            blocks_per_bit,
            get_bit_edge_detection_confidence(signal),
            code_block_buffer_lengh,
        )
    else
        sync = detect_bit_or_secondary_code_sync(
            signal,
            prn,
            code_block_buffer,
            code_block_buffer_lengh,
        )
    end
    if !sync.found
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_lengh,
            false,
            0,
            Int8(0),
            zero(UInt128),
            0,
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
        # decoding starts fresh post-sync — the integration cadence then aligns
        # each integration to the secondary-code (data-bit) boundary, so no
        # accumulator seeding is needed here. We only record the recovered
        # `secondary_phase` / `polarity`, which anchor the shared `code_phase`.
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_lengh,
            true,
            sync.phase,
            sync.polarity,
            zero(UInt128),
            0,
            complex(0.0, 0.0),
            0,
            bit_buffer.soft_bits,
            phase_acc,
        )
    end
    num_bits = min(
        div(code_block_buffer_lengh, num_code_blocks_that_form_a_bit),
        div(sizeof(code_block_buffer) * 8, num_code_blocks_that_form_a_bit),
    )
    bits = reduce(num_bits:-1:1; init = UInt128(0)) do bits, bit_index
        # Don't shadow the outer `integrated_code_blocks` argument here:
        # because this closure reassigns its own local of the same name,
        # Julia conservatively boxes the outer parameter even though this
        # closure path does not always run, leading to a per-call
        # Core.Box allocation.
        bit_sum = sum(0:num_code_blocks_that_form_a_bit-1) do code_block_index
            buffer_code_block_index =
                (bit_index - 1) * num_code_blocks_that_form_a_bit + code_block_index
            ((code_block_buffer & (one(B) << buffer_code_block_index)) > 0) * 2 - 1
        end
        push!(bit_buffer.soft_bits, Float32(bit_sum))
        bits << 1 + (bit_sum > 0)
    end
    return BitBuffer{B}(
        code_block_buffer,
        code_block_buffer_lengh,
        true,
        sync.phase,
        sync.polarity,
        bits,
        num_bits,
        complex(0, 0),
        0,
        bit_buffer.soft_bits,
        phase_acc,
    )
end

function reset(bit_buffer::BitBuffer{B}) where {B<:Unsigned}
    empty!(bit_buffer.soft_bits)
    BitBuffer{B}(
        bit_buffer.code_block_buffer,
        bit_buffer.code_block_buffer_lengh,
        bit_buffer.found,
        bit_buffer.secondary_phase,
        bit_buffer.polarity,
        zero(UInt128),
        0,
        bit_buffer.prompt_accumulator,
        bit_buffer.prompt_accumulator_integrated_code_blocks,
        bit_buffer.soft_bits,
        bit_buffer.phase_acc,
    )
end
