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

Standard-normal quantile (inverse CDF) `Φ⁻¹(p)` for `p ∈ (0, 1)`, via
Acklam's rational approximation (relative error < 1.2e-9 across the
range). Used by [`_detect_bit_edge_cfar`](@ref) to turn a false-sync
probability into a detection z-threshold without pulling in a
special-functions dependency.
"""
function _norm_quantile(p::Float64)
    # Acklam's algorithm. Coefficients for the central and tail regions.
    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00)
    plow = 0.02425
    phigh = 1 - plow
    if p < plow
        q = sqrt(-2 * log(p))
        return (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
               ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    elseif p <= phigh
        q = p - 0.5
        r = q * q
        return (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
               (((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1)
    else
        q = sqrt(-2 * log(1 - p))
        return -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
                ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    end
end

"""
$(SIGNATURES)

Soft-decision, CFAR bit-edge detector for GPS L1 C/A — the
maximum-energy timing synchronizer used by
`detect_bit_or_secondary_code_sync(::GPSL1CA, …)`. Signals with a periodic
secondary code (GPS L5I, GPS L1C-P) instead use
[`_secondary_code_search`](@ref), which correlates against a known overlay.

`soft_prompts` is the full pre-sync history of filtered prompt
correlator outputs (one per primary-code period, oldest first); `L` is
the number of primary-code blocks per navigation bit (20 for L1 C/A).

# Statistic

Under the hypothesis that the bit edge sits at phase `k ∈ 0:L-1`, the
prompts partition into `L`-block bins aligned to `k`. Each complete bin
`m` is coherently summed (`Sₖₘ = Σ pᵢ`) and its energy `|Sₖₘ|²`
accumulated; the true phase never straddles a data-bit transition, so it
keeps full coherent gain on every bin while wrong phases lose energy on
transition bins. The per-phase mean bin energy `Ēₖ` is therefore the
maximum-likelihood timing statistic for unknown i.i.d. data in AWGN.

# CFAR confidence

The thermal-noise floor is estimated *within* the winning phase's bins
(residual of each sample about its bin mean) rather than across phases —
the cross-phase spread is dominated by the deterministic transition
structure, not noise, so it would never shrink with integration time.
The peak phase is accepted only when it beats the runner-up by a margin
that is significant under that noise estimate:

    z = (Ē₍₁₎ - Ē₍₂₎) / √(Var Ē₍₁₎ + Var Ē₍₂₎))   ≥  Φ⁻¹(1 - P_fa/(L-1))

where `Var Ēₖ = (N₀² + 2|μₖ|² N₀) / Mₖ` is the non-central-χ² energy
variance (`N₀` the noise-bin energy, `|μₖ|²` the coherent signal energy,
`Mₖ` the bin count) and `P_fa = 1 - confidence`, Bonferroni-split over the
`L-1` competing phases. As integration grows the gap stays fixed while
the standard error shrinks like `1/√M`, so `z` crosses the threshold
sooner at high C/N₀ and later in noise — the detector self-paces.

# Boundary firing

Detection is reported (`found = true`) only when the most recent block
also *ends* the winning phase's bit (`length(soft_prompts) % L == k`), so
the upcoming integration starts a fresh navigation bit. This keeps the
`phase = 0` contract — and, crucially, the post-sync coherent 20-block
integration — aligned to the true bit grid, which is what makes the
off-by-one lock of issue #124 structurally impossible: the detector can
only ever fire at the energy-maximizing phase's own boundary, never at a
neighbour's.

`polarity` is the sign of the most recently completed bin's coherent sum.
"""
function _detect_bit_edge_cfar(
    soft_prompts::AbstractVector{<:Complex},
    L::Int,
    confidence::Float64,
)
    n = length(soft_prompts)
    # Need at least two complete bins on some phase before any phase can
    # serve as a runner-up; below that a transition cannot have been
    # localized yet.
    n < 2L && return SyncResult(false, 0, Int8(0))

    mean_energy = fill(-1.0, L)   # Ēₖ, -1 marks "too few bins"
    bin_count = zeros(Int, L)
    @inbounds for k in 0:(L-1)
        nb = div(n - k, L)
        nb < 1 && continue
        energy = 0.0
        for m in 0:(nb-1)
            base = k + m * L
            s = zero(ComplexF64)
            for i in 1:L
                s += soft_prompts[base+i]
            end
            energy += abs2(s)
        end
        mean_energy[k+1] = energy / nb
        bin_count[k+1] = nb
    end

    # Peak phase (needs ≥ 2 bins) and best competing phase (needs ≥ 1).
    kpeak = 0
    peak = -1.0
    @inbounds for k in 0:(L-1)
        if bin_count[k+1] >= 2 && mean_energy[k+1] > peak
            peak = mean_energy[k+1]
            kpeak = k
        end
    end
    peak < 0 && return SyncResult(false, 0, Int8(0))

    ksecond = -1
    second = -1.0
    @inbounds for k in 0:(L-1)
        if k != kpeak && bin_count[k+1] >= 1 && mean_energy[k+1] > second
            second = mean_energy[k+1]
            ksecond = k
        end
    end
    ksecond < 0 && return SyncResult(false, 0, Int8(0))

    gap = peak - second
    gap <= 0 && return SyncResult(false, 0, Int8(0))

    # Thermal-noise estimate: pooled within-bin residual variance of the
    # peak phase's bins. Pure bins (true edge) leave only thermal noise;
    # if a wrong phase transiently wins, its straddling bins inflate this
    # estimate, which only makes the test more conservative.
    nb_peak = bin_count[kpeak+1]
    resid = 0.0
    last_bin_sum = zero(ComplexF64)
    @inbounds for m in 0:(nb_peak-1)
        base = kpeak + m * L
        s = zero(ComplexF64)
        sq = 0.0
        for i in 1:L
            p = soft_prompts[base+i]
            s += p
            sq += abs2(p)
        end
        # `Σ|pᵢ|² - |Σpᵢ|²/L = Σ|pᵢ - p̄|²` is non-negative in exact
        # arithmetic, but catastrophic cancellation (near-constant high-SNR
        # bins) can round it slightly below zero; clamp so the variance and
        # its square root below stay real.
        resid += max(sq - abs2(s) / L, 0.0)   # Σ|pᵢ - p̄|² within this bin
        last_bin_sum = s
    end
    sigma2 = resid / (nb_peak * (L - 1))   # per-sample complex noise power
    noise_bin = L * sigma2                 # mean energy of a noise-only bin

    var_bin(e, m) = (noise_bin^2 + 2 * max(e - noise_bin, 0.0) * noise_bin) / m
    se = sqrt(max(var_bin(peak, nb_peak) + var_bin(second, bin_count[ksecond+1]), 0.0))
    z = se > 0 ? gap / se : Inf

    z_threshold = _norm_quantile(1 - (1 - confidence) / (L - 1))
    z < z_threshold && return SyncResult(false, 0, Int8(0))

    # Fire only at the winning phase's own bit boundary so the upcoming
    # integration starts a new bit.
    n % L != kpeak && return SyncResult(false, 0, Int8(0))

    polarity = real(last_bin_sum) < 0 ? Int8(-1) : Int8(+1)
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

The `soft_prompts` field is the pre-sync soft-decision history (the
filtered prompt of each primary-code block, oldest first) consumed by the
soft bit-edge detector [`_detect_bit_edge_cfar`](@ref). It is populated
only for signals whose [`uses_soft_bit_edge_detection`](@ref) trait is
`true` (currently GPS L1 C/A), capped to a bounded window
([`SOFT_PROMPT_HISTORY_CAP`](@ref) blocks), and emptied once sync is
found.
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
    soft_prompts::Vector{ComplexF64}
end

# Default constructor preserves the pre-refactor `UInt128`-backed search
# buffer. Once `get_code_block_buffer_type` lands (Step 2) the per-signal
# `TrackedSignal` constructor picks the right width instead.
function BitBuffer()
    BitBuffer{UInt128}(
        zero(UInt128), 0, false, 0, Int8(0), zero(UInt128), 0, complex(0.0, 0.0), 0,
        Float32[], ComplexF64[],
    )
end

# Typed empty constructor used by the per-signal `TrackedSignal` path.
function BitBuffer{B}() where {B<:Unsigned}
    BitBuffer{B}(
        zero(B), 0, false, 0, Int8(0), zero(UInt128), 0, complex(0.0, 0.0), 0,
        Float32[], ComplexF64[],
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
        ComplexF64[],
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

Whether `signal`'s bit-edge sync uses the soft-decision, maximum-energy
CFAR detector [`_detect_bit_edge_cfar`](@ref) (which reads the
`soft_prompts` history) instead of the hard-decision
`detect_bit_or_secondary_code_sync` path.

`true` only for GPS L1 C/A. Returned as a compile-time constant so the
branch in [`_buffer_find_bit`](@ref) folds away per signal type — signals
that return `false` never allocate or populate `soft_prompts`.
"""
@inline uses_soft_bit_edge_detection(::AbstractGNSSSignal) = false

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

"""
$(SIGNATURES)

Maximum number of soft prompts retained in the pre-sync `soft_prompts`
history before the oldest full bit's worth is dropped. Bounds memory and
per-call work if a satellite never locks (e.g. a stream with no data
transitions); 400 blocks (~400 ms for L1 C/A) is far longer than any
trackable cold start needs.
"""
const SOFT_PROMPT_HISTORY_CAP = 400

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
            bit_buffer.soft_prompts,
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
            bit_buffer.soft_prompts,
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

    # Signals that detect the bit edge from soft prompts (GPS L1 C/A) keep a
    # bounded soft-prompt history and run the maximum-energy CFAR detector;
    # everything else stays on the hard-decision sliding-window path. The
    # branch folds at compile time per signal type, so non-soft signals never
    # touch `soft_prompts`.
    soft_prompts = bit_buffer.soft_prompts
    if uses_soft_bit_edge_detection(signal)
        push!(soft_prompts, ComplexF64(prompt))
        # Drop whole bits' worth of the oldest history when over the cap, on a
        # bit boundary, so the window start stays aligned to the bit grid and
        # the detector's phase indexing is unaffected.
        if length(soft_prompts) > SOFT_PROMPT_HISTORY_CAP &&
           num_code_blocks_that_form_a_bit > 0 &&
           code_block_buffer_lengh % num_code_blocks_that_form_a_bit == 0
            deleteat!(soft_prompts, 1:num_code_blocks_that_form_a_bit)
        end
        sync = _detect_bit_edge_cfar(
            soft_prompts,
            num_code_blocks_that_form_a_bit,
            get_bit_edge_detection_confidence(signal),
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
            soft_prompts,
        )
    end
    # Sync found: the soft-prompt history is dead state from here on.
    empty!(soft_prompts)
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
            soft_prompts,
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
        soft_prompts,
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
        bit_buffer.soft_prompts,
    )
end
