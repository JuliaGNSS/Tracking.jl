module BitBufferTest

using Test: @test, @testset, @inferred, @test_throws
using Random: MersenneTwister
using GNSSSignals:
    GPSL1CA, GPSL5I, get_secondary_code, secondary_value, get_secondary_code_length
using Tracking:
    BitBuffer,
    buffer,
    reset,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    SyncResult,
    PhaseAccumulators,
    _norm_quantile,
    _t_quantile,
    _cfar_decide,
    _detect_bit_edge_cfar,
    _detect_secondary_code_cfar,
    _seed_phase_accumulators!,
    _update_phase_accumulators!,
    _update_secondary_accumulators!,
    _secondary_code_search,
    _packed_secondary_code

@testset "SyncResult" begin
    r = @inferred SyncResult(false, 0, Int8(0))
    @test r.found == false
    @test r.phase == 0
    @test r.polarity == 0
end

@testset "_norm_quantile" begin
    @test _norm_quantile(0.5) == 0.0
    @test _norm_quantile(0.975) ≈ 1.959964 atol = 1e-6     # tabulated Φ⁻¹
    @test _norm_quantile(0.999) ≈ 3.090232 atol = 1e-6
    @test _norm_quantile(0.025) ≈ -_norm_quantile(0.975) atol = 1e-12
    # The open-interval contract: `_t_quantile` clamps before calling in, so
    # the ±Inf endpoints are documented behaviour rather than an error.
    @test _norm_quantile(1.0) == Inf
    @test _norm_quantile(0.0) == -Inf
end

# `@allocated` at module scope picks up boxing from untyped global lookups, so
# the measurement goes through a typed function — the same reason
# `test/track_in_place.jl` measures inside helpers.
_measure_t_quantile_alloc(probability::Float64, dof::Int) =
    (_t_quantile(probability, dof); @allocated _t_quantile(probability, dof))

@testset "_t_quantile" begin
    # Median is exactly 0 — and, as a regression guard, `probability == 0.5`
    # must return directly rather than recurse into `_t_quantile(1 - 0.5, dof)`
    # (which previously self-recursed forever → stack overflow).
    @test _t_quantile(0.5, 1) === 0.0
    @test _t_quantile(0.5, 30) === 0.0
    # dof = 1 is the standard Cauchy, quantile(p) = tan(π(p − 0.5)).
    @test _t_quantile(0.75, 1) ≈ 1.0 atol = 1e-6            # tan(π/4)
    @test _t_quantile(0.9, 1) ≈ tan(pi * 0.4) atol = 1e-6   # ≈ 3.0777
    # Symmetric about 0.
    @test _t_quantile(0.25, 1) ≈ -_t_quantile(0.75, 1) atol = 1e-9
    @test _t_quantile(0.001, 7) ≈ -_t_quantile(0.999, 7) atol = 1e-6
    # Accuracy against standard Student-t tables. `dof = 1` (Cauchy) and
    # `dof = 2` are closed forms and land on the exact value; above that Hill's
    # series carries a small relative error (worst ~2e-4 around dof 3–10),
    # which is still orders of magnitude finer than the nominal-d.o.f.
    # modelling error the threshold is built on — and far too fine to move
    # which block the detector locks on.
    # (The two closed forms are pinned at the precision of the table literal
    # itself, 7 significant digits — they agree with the exact value to ~1e-16.)
    @test _t_quantile(0.975, 1) ≈ 12.706205 rtol = 1e-7
    @test _t_quantile(0.975, 2) ≈ 4.302653 rtol = 1e-7
    @test _t_quantile(0.995, 3) ≈ 5.840909 rtol = 1e-4
    @test _t_quantile(0.95, 5) ≈ 2.015048 rtol = 1e-4
    @test _t_quantile(0.975, 10) ≈ 2.228139 rtol = 1e-4
    @test _t_quantile(0.999, 20) ≈ 3.551808 rtol = 1e-4
    @test _t_quantile(0.975, 30) ≈ 2.042272 rtol = 1e-4
    @test _t_quantile(0.975, 120) ≈ 1.979930 rtol = 1e-4
    # Heavier tails than the normal, converging to it as dof → ∞ (no cutoff).
    p = 1 - (1 - 0.999) / (20 - 1)   # the detector's Bonferroni-split argument, L1CA
    normal_limit = 3.87813           # tabulated Φ⁻¹(p), the dof → ∞ limit
    @test _t_quantile(p, 1) > _t_quantile(p, 10) > _t_quantile(p, 1000) > normal_limit
    @test _t_quantile(p, 100_000) ≈ normal_limit atol = 1e-3
    # The behaviour the fix relies on: at 1 d.o.f. (2 bins) the threshold is huge,
    # so the observed z ≈ 40 reacquisition fluctuation is rejected; by ~70 bins a
    # clean z ≈ 3.9 clears it.
    @test _t_quantile(p, 1) > 1000
    @test _t_quantile(p, 69) < 4.2
    # Monotone in dof over the whole range the detector sweeps — no branch seam
    # between the closed forms (dof ≤ 2) and Hill's two series.
    @test issorted([_t_quantile(p, dof) for dof = 1:400]; rev = true)
    # Allocation-free at every dof. The detector calls this once per
    # primary-code block for every unsynced satellite, so an allocating
    # quantile makes `track!` allocate in proportion to the signal length
    # (the exact `SpecialFunctions.beta_inc_inv` inversion did, via internal
    # `zeros(31)` scratch, from ~12 d.o.f. upwards). `test/track_in_place.jl`
    # guards the same property end-to-end; this pins it per dof, where the
    # allocating stretches of the range are unmissable.
    @test all(_measure_t_quantile_alloc(p, dof) == 0 for dof = 1:400)
    @test _measure_t_quantile_alloc(0.25, 7) == 0    # the reflected branch too
    # The detector clamps its argument to the open interval; the endpoint must
    # still produce a finite (huge) threshold rather than NaN, or the
    # `z_score < threshold` gate would silently pass.
    @test isfinite(_t_quantile(prevfloat(1.0), 1))
    @test isfinite(_t_quantile(prevfloat(1.0), 40))
end

const L1CA_BLOCKS_PER_BIT = 20  # primary-code blocks per L1 C/A navigation bit

# Build a noiseless ±1 soft-prompt stream from a list of data bits, one bit =
# `L1CA_BLOCKS_PER_BIT` blocks, scaled by `amp`.
_bitstream(bits; amp = 1.0) =
    ComplexF64[amp * (b == 1 ? 1.0 : -1.0) for b in bits for _ = 1:L1CA_BLOCKS_PER_BIT]

# Fold a prompt stream into a fresh set of phase accumulators (mirrors what
# `_buffer_find_bit` does) and run the detector after the n-th block. With
# `upto = 0` it runs after every block and returns the (1-based) block index
# at which it first locks (0 if never), else it returns the SyncResult after
# exactly `upto` blocks.
function _detect_over(prompts, confidence; upto = 0)
    accumulators = PhaseAccumulators()
    _seed_phase_accumulators!(accumulators, L1CA_BLOCKS_PER_BIT)
    for (block_number, prompt) in enumerate(prompts)
        _update_phase_accumulators!(
            accumulators,
            ComplexF64(prompt),
            block_number - 1,
            L1CA_BLOCKS_PER_BIT,
        )
        result = _detect_bit_edge_cfar(
            accumulators,
            L1CA_BLOCKS_PER_BIT,
            confidence,
            block_number,
        )
        upto == 0 && result.found && return block_number
        upto == block_number && return result
    end
    return upto == 0 ? 0 :
           _detect_bit_edge_cfar(
        accumulators,
        L1CA_BLOCKS_PER_BIT,
        confidence,
        length(prompts),
    )
end

@testset "_detect_bit_edge_cfar" begin
    @testset "Needs at least two bins" begin
        # 39 blocks (< 2L) — never enough to nominate a runner-up.
        @test _detect_over(_bitstream([0, 1])[1:39], 0.999; upto = 39).found == false
    end

    @testset "Noiseless lock fires at the true bit boundary, not one early" begin
        # Data 0,0,1: first transition preceded by a repeated bit (the
        # exact issue-#124 trigger). The true edge is at block 60; the old
        # tolerant matcher fired at 59.
        @test _detect_over(_bitstream([0, 0, 1])[1:59], 0.999; upto = 59).found == false
        res = _detect_over(_bitstream([0, 0, 1]), 0.999; upto = 60)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == +1   # last completed bin is the "1"

        # Mirror pattern 1,1,0 → same boundary, negative polarity.
        res = _detect_over(_bitstream([1, 1, 0]), 0.999; upto = 60)
        @test res.found == true
        @test res.polarity == -1
    end

    @testset "No data transition never locks" begin
        # A constant data stream carries no edge information.
        @test _detect_over(_bitstream(fill(1, 5)), 0.999) == 0
    end

    @testset "Confidence drives lock latency in noise" begin
        # Alternating bits give a transition at every boundary — the
        # cleanest possible edge evidence. Add noise and feed the stream
        # block-by-block; a higher confidence target locks no earlier than a
        # lower one, and always at a true bit boundary.
        function lock_block(confidence, seed)
            rng = MersenneTwister(seed)
            clean = _bitstream([bit % 2 for bit = 0:39]; amp = 8.0)
            noisy = ComplexF64[prompt + complex(randn(rng), randn(rng)) for prompt in clean]
            _detect_over(noisy, confidence)
        end
        for seed = 1:5
            low_confidence_lock = lock_block(0.95, seed)
            high_confidence_lock = lock_block(0.99999, seed)
            # Both lock on this high-SNR stream within the window.
            @test low_confidence_lock > 0
            @test high_confidence_lock > 0
            # More confidence ⇒ never an earlier lock.
            @test high_confidence_lock >= low_confidence_lock
            # Whenever it does lock it is at a true bit boundary.
            @test low_confidence_lock % L1CA_BLOCKS_PER_BIT == 0
            @test high_confidence_lock % L1CA_BLOCKS_PER_BIT == 0
        end
    end

    @testset "confidence = 1.0 stays conservative, not an instant lock" begin
        # On a *noisy* signal (finite z) confidence 1.0 must be maximally
        # conservative, not lock at the first boundary the way the old
        # NaN-threshold path did. (A noiseless edge has z = Inf and still
        # locks — genuine certainty, exercised above.)
        rng = MersenneTwister(7)
        clean = _bitstream([bit % 2 for bit = 0:39]; amp = 3.0)
        noisy = ComplexF64[prompt + complex(randn(rng), randn(rng)) for prompt in clean]
        low_confidence_lock = _detect_over(noisy, 0.99)
        high_confidence_lock = _detect_over(noisy, 1.0)
        @test low_confidence_lock > 0                                      # locks at moderate confidence
        @test high_confidence_lock == 0 || high_confidence_lock > low_confidence_lock  # 1.0 never earlier
    end

    @testset "no false lock on a long near-constant run (Welford stability)" begin
        # High-magnitude, near-constant prompts with a tiny systematic drift
        # and no real bit edge. A naive Σe² − (Σe)²/M variance would lose
        # precision and could read zero (→ infinite confidence → false lock)
        # at large bin counts; Welford keeps it stable, so this never locks.
        accumulators = PhaseAccumulators()
        _seed_phase_accumulators!(accumulators, L1CA_BLOCKS_PER_BIT)
        found = false
        for block_number = 1:5000
            _update_phase_accumulators!(
                accumulators,
                complex(1e4 * (1 + 1e-7 * (block_number - 1)), 0.0),
                block_number - 1,
                L1CA_BLOCKS_PER_BIT,
            )
            if _detect_bit_edge_cfar(accumulators, L1CA_BLOCKS_PER_BIT, 0.999, block_number).found
                found = true
                break
            end
        end
        @test !found
    end
end

@testset "_cfar_decide" begin
    # The shared CFAR decision core: peak vs. runner-up, Welford noise scale,
    # Student-t threshold. `period` is the bin length = hypothesis count.
    period = 5

    @testset "Needs two bins on some hypothesis" begin
        # num_blocks < 2 * period ⇒ no runner-up can exist yet.
        mean = [10.0, 1.0, 1.0, 1.0, 1.0]
        m2 = zeros(5)
        @test _cfar_decide(mean, m2, 2 * period - 1, period, 0.999) == (false, -1, 0)
    end

    @testset "Noiseless separation locks at the peak" begin
        mean = [10.0, 1.0, 1.0, 1.0, 1.0]
        m2 = zeros(5)                              # zero variance ⇒ z = Inf
        accepted, peak, count = _cfar_decide(mean, m2, 2 * period, period, 0.999)
        @test accepted
        @test peak == 0                            # hypothesis 0 (0-based) is the peak
        @test count == 2                           # div(10 - 0, 5)
    end

    @testset "No separation never locks" begin
        mean = fill(5.0, 5)                         # every hypothesis equal
        m2 = zeros(5)
        accepted, _, _ = _cfar_decide(mean, m2, 2 * period, period, 0.999)
        @test !accepted
    end

    @testset "A large peak variance suppresses the lock" begin
        # Same energy gap, but the peak's own bin-to-bin spread is huge, so the
        # z-score stays below threshold and it does not lock.
        mean = [10.0, 1.0, 1.0, 1.0, 1.0]
        quiet = zeros(5)
        loud = [1.0e6, 0.0, 0.0, 0.0, 0.0]         # M₂ only on the peak
        @test _cfar_decide(mean, quiet, 2 * period, period, 0.999)[1]
        @test !_cfar_decide(mean, loud, 2 * period, period, 0.999)[1]
    end
end

# --- Soft secondary-code CFAR detector (GPS L5I NH10 and friends) ------------

# The ±1 secondary chips in time order for `signal`/`prn`.
_secondary_chips(signal, prn) = [
    Int(secondary_value(get_secondary_code(signal), prn, k)) for
    k = 0:(get_secondary_code_length(signal)-1)
]

# NH10 packed newest-first, as the hard detector references it.
_packed_l5i() = _packed_secondary_code(UInt32, GPSL5I(), 1)

# Build a prompt stream for a secondary-coded signal: block `i` (0-based) carries
# secondary chip `(i + start_chip) % N` times a per-period data symbol (constant
# over each period, flipped pseudo-randomly at period boundaries) times `amp`,
# plus optional complex Gaussian noise.
function _secondary_stream(
    signal,
    prn,
    nblocks;
    start_chip = 0,
    amp = 5.0,
    noise = 0.0,
    seed = 1,
)
    N = get_secondary_code_length(signal)
    chips = _secondary_chips(signal, prn)
    rng = MersenneTwister(seed)
    data = 1
    prompts = ComplexF64[]
    for i = 0:(nblocks-1)
        chip = mod(i + start_chip, N)
        chip == 0 && (data = rand(rng, (-1, 1)))
        p = ComplexF64(amp * data * chips[chip+1])
        noise > 0 && (p += noise * complex(randn(rng), randn(rng)))
        push!(prompts, p)
    end
    prompts
end

# Fold a prompt stream into fresh secondary accumulators (mirrors what
# `_buffer_find_bit` does) and run the detector after each block. `upto = 0`
# returns the 1-based block index of the first lock (0 if never); otherwise it
# returns the SyncResult after exactly `upto` blocks.
function _secondary_detect_over(prompts, signal, prn, confidence; upto = 0)
    N = get_secondary_code_length(signal)
    accumulators = PhaseAccumulators()
    _seed_phase_accumulators!(accumulators, N)
    local result
    for (block_number, prompt) in enumerate(prompts)
        _update_secondary_accumulators!(
            accumulators,
            ComplexF64(prompt),
            block_number - 1,
            N,
            signal,
            prn,
        )
        result = _detect_secondary_code_cfar(accumulators, N, confidence, block_number)
        upto == 0 && result.found && return block_number
        upto == block_number && return result
    end
    return upto == 0 ? 0 : result
end

@testset "_update_secondary_accumulators!" begin
    # Feed two clean NH10 periods aligned to chip 0. The correct rotation wipes
    # the overlay so its bin sums coherently (energy ≈ (N·amp)²); every other
    # rotation straddles the fixed NH transitions and loses energy.
    signal = GPSL5I()
    prn = 1
    N = get_secondary_code_length(signal)          # 10
    amp = 5.0
    accumulators = PhaseAccumulators()
    _seed_phase_accumulators!(accumulators, N)
    prompts = _secondary_stream(signal, prn, 2N; amp)   # data constant per period
    for (i, p) in enumerate(prompts)
        _update_secondary_accumulators!(accumulators, ComplexF64(p), i - 1, N, signal, prn)
    end
    energies = accumulators.mean_bin_energy
    peak = argmax(energies) - 1                    # 0-based winning rotation
    # Aligned to chip 0 ⇒ the winning rotation is 0.
    @test peak == 0
    @test energies[peak+1] ≈ (N * amp)^2
    # The peak dwarfs every competitor by a wide margin.
    runner_up = maximum(energies[2:end])
    @test energies[1] > 20 * runner_up
end

@testset "_detect_secondary_code_cfar" begin
    signal = GPSL5I()
    prn = 1
    N = get_secondary_code_length(signal)          # 10

    @testset "Needs at least two periods" begin
        prompts = _secondary_stream(signal, prn, 2N - 1; amp = 8.0)
        @test _secondary_detect_over(prompts, signal, prn, 0.999; upto = 2N - 1).found ==
              false
    end

    @testset "Clean lock fires at the period boundary with phase 0" begin
        for start_chip = 0:(N-1)
            prompts = _secondary_stream(signal, prn, 4N; start_chip, amp = 8.0)
            synced = _secondary_detect_over(prompts, signal, prn, 0.999)
            @test synced >= 2N                     # never before two periods
            # Fires only at a true NH10 boundary: the just-processed block is the
            # last (chip N-1) of the winning rotation's period.
            @test (start_chip + synced - 1) % N == N - 1
            res = _secondary_detect_over(prompts, signal, prn, 0.999; upto = synced)
            @test res.found
            @test res.phase == 0                   # upcoming integration at chip 0
        end
    end

    @testset "Polarity follows the winning period's sign" begin
        chips = _secondary_chips(signal, prn)
        # A period that is exactly +overlay locks at +1; the negated overlay at −1.
        pos = ComplexF64[8.0 * c for c in chips]
        neg = ComplexF64[-8.0 * c for c in chips]
        @test _secondary_detect_over(vcat(pos, pos, pos), signal, prn, 0.999; upto = 3N).polarity ==
              +1
        @test _secondary_detect_over(vcat(neg, neg, neg), signal, prn, 0.999; upto = 3N).polarity ==
              -1
    end

    @testset "No false lock on pure noise (the TEX-CUP PRN 8 symptom)" begin
        # The hard sign-template match (NH10, exact-match budget) accepts a random
        # 10-bit window whenever noise happens to align it with one of the ~20
        # rotation×polarity templates — a ~2 % per-window chance that accumulates
        # into false locks over a long pre-lock stretch. The soft CFAR detector,
        # feeding on the same noise, requires a significant energy gap over many
        # periods, so it does not lock on noise.
        hard_false_locks = 0
        soft_locks = 0
        for seed = 1:200
            rng = MersenneTwister(1000 + seed)
            # Soft: 40 blocks of pure complex-Gaussian noise (no signal).
            noise = ComplexF64[complex(randn(rng), randn(rng)) for _ = 1:4N]
            _secondary_detect_over(noise, signal, prn, 0.999) != 0 && (soft_locks += 1)
            # Hard: a random 10-bit sign window through the rotation sweep.
            window = rand(rng, UInt32) & UInt32((1 << N) - 1)
            _secondary_code_search(window, _packed_l5i(), N, 0).found &&
                (hard_false_locks += 1)
        end
        @test soft_locks == 0
        # Sanity: the hard exact-match detector really is noise-prone here.
        @test hard_false_locks > 0
    end
end

# Feed a noiseless ±1 prompt stream (one code block per call) through
# `buffer()` and return the 1-based block index at which the detector locked
# (0 if it never locked) together with the final bit buffer.
function _feed_prompts(prompts)
    signal = GPSL1CA()
    bit_buffer = BitBuffer{UInt64}()
    found_at = 0
    for (i, prompt) in enumerate(prompts)
        bit_buffer = buffer(signal, 1, bit_buffer, 1, prompt)
        if found_at == 0 && has_bit_or_secondary_code_been_found(bit_buffer)
            found_at = i
        end
    end
    found_at, bit_buffer
end

@testset "L1CA bit-edge lock is not one block early through buffer() (issue #124)" begin
    # Data bits 0, 0, 1 — the first transition is preceded by a repeated
    # bit. The buggy tolerant matcher fired at block 59 (one block before
    # the true edge); the soft edge-locked detector fires exactly at 60.
    found_at, bit_buffer = _feed_prompts([fill(-1.0 + 0.0im, 40); fill(1.0 + 0.0im, 20)])
    @test found_at == 60
    @test bit_buffer.polarity == +1

    # Mirror pattern 1, 1, 0 locks at the same block with negative polarity.
    found_at, bit_buffer = _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
    @test found_at == 60
    @test bit_buffer.polarity == -1

    # A constant prompt stream (no data transition) never locks.
    found_at, _ = _feed_prompts(fill(1.0 + 0.0im, 80))
    @test found_at == 0
end

@testset "Bit buffer" begin
    @testset "Initialize" begin
        bit_buffer = @inferred BitBuffer()
        @test bit_buffer.code_block_buffer == 0
        @test bit_buffer.code_block_buffer_length == 0
        @test has_bit_or_secondary_code_been_found(bit_buffer) == false
        @test bit_buffer.secondary_phase == 0
        @test bit_buffer.polarity == 0
        @test isempty(get_soft_bits(bit_buffer))
        @test length(bit_buffer) == 0
        @test bit_buffer.prompt_accumulator == complex(0, 0)
        @test bit_buffer.prompt_accumulator_integrated_code_blocks == 0
    end

    @testset "Throw error if bit hasn't been found yet and integrated code blocks is greater than 1" begin
        bit_buffer = @inferred BitBuffer()

        signal = GPSL1CA()
        next_bit_buffer = @test_throws "The number code blocks must be equal to 1" buffer(
            signal,
            1,
            bit_buffer,
            2,
            2 + 0im,
        )
    end

    @testset "Buffer" begin
        bit_buffer = BitBuffer()
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 1, 2 + 0im)
        @test length(next_bit_buffer) == 0
        @test isempty(get_soft_bits(next_bit_buffer))
        @test next_bit_buffer.code_block_buffer == 1
        @test next_bit_buffer.code_block_buffer_length == 1
    end

    @testset "Find bit start and buffer the pre-sync bits (negative polarity)" begin
        # Drive a clean data-1,1,0 stream (20 blocks each, +,+,-) through
        # `buffer()`. The soft detector locks at the true bit boundary
        # (block 60); the last completed bin is the "0", so the lock is at
        # negative polarity. The three pre-sync bits are decoded with that
        # polarity applied — the same sign-flip every post-sync bit gets — so
        # they stay consistent across the sync boundary (issue #127): block
        # signs +,+,- decode as bits 0,0,1, i.e. soft bits -,-,+.
        found_at, bit_buffer =
            _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
        @test found_at == 60
        @test bit_buffer.found == true
        @test bit_buffer.polarity == -1
        @test length(bit_buffer) == 3
        @test bit_buffer.code_block_buffer_length == 60
        soft = get_soft_bits(bit_buffer)
        @test length(soft) == 3
        @test soft[1] < 0 && soft[2] < 0 && soft[3] > 0
    end

    @testset "Find bit start and buffer the pre-sync bits (positive polarity)" begin
        # The whole-signal sign flip of the case above: data 0,0,1 (block
        # signs -,-,+). The last completed bin is the "1", so the lock is at
        # positive polarity. Because a global RF sign flip is exactly the
        # ambiguity polarity resolves, the recovered soft-bit signs are
        # identical to the negative-polarity case (issue #127): soft -,-,+.
        found_at, bit_buffer =
            _feed_prompts([fill(-1.0 + 0.0im, 40); fill(1.0 + 0.0im, 20)])
        @test found_at == 60
        @test bit_buffer.found == true
        @test bit_buffer.polarity == +1
        @test length(bit_buffer) == 3
        @test bit_buffer.code_block_buffer_length == 60
        soft = get_soft_bits(bit_buffer)
        @test length(soft) == 3
        @test soft[1] < 0 && soft[2] < 0 && soft[3] > 0
    end

    @testset "Buffer prompt when bit has been found" begin
        code_blocks_buffer = 0xfffffffffff0000
        code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
        bit_buffer = BitBuffer(
            code_blocks_buffer,
            code_blocks_buffer_length,
            true,
            complex(-1, 0),
            1,
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 1, -2 + 0im)
        @test isempty(get_soft_bits(next_bit_buffer))
        @test length(next_bit_buffer) == 0
        @test next_bit_buffer.prompt_accumulator == -3 + 0im
        @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 2
    end

    @testset "Buffer bit when bit has been found" begin
        code_blocks_buffer = 0xfffffffffff0000
        code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
        bit_buffer = BitBuffer(
            code_blocks_buffer,
            code_blocks_buffer_length,
            true,
            complex(-10, 2),
            19,
            Float32[-1.0, 1.0],
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 1, -2 + 0im)
        # Seeded soft bits -,+ plus the completed bit (-10+2im plus -2 → -12).
        @test sign.(get_soft_bits(next_bit_buffer)) == Float32[-1, 1, -1]
        @test length(next_bit_buffer) == 3
        @test next_bit_buffer.prompt_accumulator == 0 + 0im
        @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
    end

    @testset "Buffer bit when bit has been found" begin
        code_blocks_buffer = 0xfffffffffff0000
        code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
        bit_buffer = BitBuffer(
            code_blocks_buffer,
            code_blocks_buffer_length,
            true,
            complex(10, 2),
            10,
            Float32[1.0, 1.0],
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 10, 10 + 1im)
        # Seeded soft bits +,+ plus the completed bit (10+2im plus 10+1im → +20).
        @test sign.(get_soft_bits(next_bit_buffer)) == Float32[1, 1, 1]
        @test length(next_bit_buffer) == 3
        @test next_bit_buffer.prompt_accumulator == 0 + 0im
        @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
    end

    @testset "Soft bits" begin
        @testset "Initialized empty" begin
            bit_buffer = BitBuffer()
            @test get_soft_bits(bit_buffer) isa Vector{Float32}
            @test isempty(get_soft_bits(bit_buffer))
        end

        @testset "Accumulated sum is stored on completed bit" begin
            code_blocks_buffer = 0xfffffffffff0000
            code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
            bit_buffer = BitBuffer(
                code_blocks_buffer,
                code_blocks_buffer_length,
                true,
                complex(-10.0, 2.0),
                19,
                Float32[-1.0, 1.0],
            )
            signal = GPSL1CA()

            # 20th code block completes the bit; soft bit = real of the sum.
            # It lands after the two soft bits the buffer was seeded with.
            next_bit_buffer = buffer(signal, 1, bit_buffer, 1, -2 + 0im)
            @test get_soft_bits(next_bit_buffer) == Float32[-1.0, 1.0, -12.0]
            @test eltype(get_soft_bits(next_bit_buffer)) == Float32
        end

        @testset "No soft bit is stored before a bit completes" begin
            code_blocks_buffer = 0xfffffffffff0000
            code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
            bit_buffer = BitBuffer(
                code_blocks_buffer,
                code_blocks_buffer_length,
                true,
                complex(-1.0, 0.0),
                1,
            )
            signal = GPSL1CA()

            next_bit_buffer = buffer(signal, 1, bit_buffer, 1, -2 + 0im)
            @test isempty(get_soft_bits(next_bit_buffer))
        end

        @testset "Pre-sync recovered soft bits are amplitude-scaled (issue #134)" begin
            # Same data-1,1,0 stream as the negative-polarity lock test, but
            # at amplitude 2 rather than 1. At sync the pre-sync bits only
            # have ±1 prompt signs available, so each recovered bit's sign
            # vote (±20 over the 20 blocks/bit, polarity-corrected) is scaled
            # by the sync-time prompt magnitude — so these soft bits live in
            # the same coherent-amplitude-sum units as the post-sync bits
            # instead of being bare vote counts (issue #134).
            found_at, bit_buffer =
                _feed_prompts([fill(2.0 + 0.0im, 40); fill(-2.0 + 0.0im, 20)])
            @test found_at == 60
            @test length(bit_buffer) == 3
            soft = get_soft_bits(bit_buffer)
            # 20 sign votes × |prompt| (2) = magnitude 40; signs match the
            # decoded bits (0,0,1 → soft -,-,+, as in the polarity test).
            @test soft == Float32[-40.0, -40.0, 40.0]

            # Unit-amplitude prompts give magnitude 20 — i.e. the magnitude
            # tracks |prompt|, confirming the scaling (a raw vote count would
            # be ±20 in both runs).
            _, unit_buffer = _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
            @test get_soft_bits(unit_buffer) == Float32[-20.0, -20.0, 20.0]
        end

        @testset "Reset empties the soft bits but keeps the vector" begin
            bit_buffer = BitBuffer()
            push!(get_soft_bits(bit_buffer), 1.0f0, -2.0f0)
            soft_bits = get_soft_bits(bit_buffer)

            reset_bit_buffer = reset(bit_buffer)
            @test isempty(get_soft_bits(reset_bit_buffer))
            # Same vector is reused (non-allocating after the first track calls)
            @test get_soft_bits(reset_bit_buffer) === soft_bits
        end
    end

    @testset "Bit accumulation is unbounded" begin
        # Bits accumulate as soft bits with no ceiling (the fixed UInt128 that
        # used to overflow at 128 bits — issue #134 — is gone). A fold that
        # spans more than 128 bits must keep counting instead of throwing.
        signal = GPSL1CA()
        bit_buffer = BitBuffer(UInt128(0), 0, true, complex(0.0, 0.0), 0)
        for _ = 1:(200*20)
            bit_buffer = buffer(signal, 1, bit_buffer, 1, 1.0 + 0.0im)
        end
        @test length(bit_buffer) == 200
        @test length(get_soft_bits(bit_buffer)) == 200
        @test all(>(0), get_soft_bits(bit_buffer))
    end
end

end
