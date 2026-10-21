module BitBufferTest

using Test: @test, @testset, @inferred, @test_logs, @test_throws
using Random: MersenneTwister
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    get_secondary_code,
    secondary_value,
    get_secondary_code_length,
    get_signal_id
using Tracking: _warn_bit_boundary_overshoot
using TrackingLoops:
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
    # `_t_quantile` clamps first; the endpoints give ±Inf, not an error.
    @test _norm_quantile(1.0) == Inf
    @test _norm_quantile(0.0) == -Inf
end

# In a function: `@allocated` at module scope counts untyped-global boxing.
_measure_t_quantile_alloc(probability::Float64, dof::Int) =
    (_t_quantile(probability, dof); @allocated _t_quantile(probability, dof))

@testset "_t_quantile" begin
    # Exactly 0, returned directly: reflecting 0.5 would recurse forever.
    @test _t_quantile(0.5, 1) === 0.0
    @test _t_quantile(0.5, 30) === 0.0
    # dof = 1 is the standard Cauchy, quantile(p) = tan(π(p − 0.5)).
    @test _t_quantile(0.75, 1) ≈ 1.0 atol = 1e-6            # tan(π/4)
    @test _t_quantile(0.9, 1) ≈ tan(pi * 0.4) atol = 1e-6   # ≈ 3.0777
    # Symmetric about 0.
    @test _t_quantile(0.25, 1) ≈ -_t_quantile(0.75, 1) atol = 1e-9
    @test _t_quantile(0.001, 7) ≈ -_t_quantile(0.999, 7) atol = 1e-6
    # Against Student-t tables. `dof` 1 and 2 are closed forms, exact to the
    # table's 7 digits; Hill's series above them is off by up to ~2e-4 (dof
    # 3–10), far too little to move the lock block.
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
    # At 2 bins the threshold rejects the z ≈ 40 seen on reacquisition; by ~70
    # bins a clean z ≈ 3.9 clears it.
    @test _t_quantile(p, 1) > 1000
    @test _t_quantile(p, 69) < 4.2
    # Monotone over the swept range, with no seam between the branches.
    @test issorted([_t_quantile(p, dof) for dof = 1:400]; rev = true)
    # Allocation-free at every dof: it runs once per block per unsynced
    # satellite (end-to-end guard in `test/track_in_place.jl`).
    @test all(_measure_t_quantile_alloc(p, dof) == 0 for dof = 1:400)
    @test _measure_t_quantile_alloc(0.25, 7) == 0    # the reflected branch too
    # Finite, not NaN, at the clamped endpoint, or `z_score < threshold` would
    # silently pass.
    @test isfinite(_t_quantile(prevfloat(1.0), 1))
    @test isfinite(_t_quantile(prevfloat(1.0), 40))
end

const L1CA_BLOCKS_PER_BIT = 20  # primary-code blocks per L1 C/A navigation bit

# A noiseless ±`amp` prompt stream, `L1CA_BLOCKS_PER_BIT` blocks per data bit.
_bitstream(bits; amp = 1.0) =
    ComplexF64[amp * (b == 1 ? 1.0 : -1.0) for b in bits for _ = 1:L1CA_BLOCKS_PER_BIT]

# Fold prompts into fresh accumulators as `_buffer_find_bit` does, running the
# detector after each block. Returns the 1-based block of the first lock (0 if
# none) for `upto = 0`, else the SyncResult after `upto` blocks.
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
        # < 2L blocks: no runner-up yet.
        @test _detect_over(_bitstream([0, 1])[1:39], 0.999; upto = 39).found == false
    end

    @testset "Noiseless lock fires at the true bit boundary, not one early" begin
        # A repeated bit before the first transition (#124): the edge is at
        # block 60, not 59.
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
        @test _detect_over(_bitstream(fill(1, 5)), 0.999) == 0
    end

    @testset "Confidence drives lock latency in noise" begin
        # Alternating bits in noise: higher confidence never locks earlier,
        # and always on a bit boundary.
        function lock_block(confidence, seed)
            rng = MersenneTwister(seed)
            clean = _bitstream([bit % 2 for bit = 0:39]; amp = 8.0)
            noisy = ComplexF64[prompt + complex(randn(rng), randn(rng)) for prompt in clean]
            _detect_over(noisy, confidence)
        end
        for seed = 1:5
            low_confidence_lock = lock_block(0.95, seed)
            high_confidence_lock = lock_block(0.99999, seed)
            @test low_confidence_lock > 0
            @test high_confidence_lock > 0
            @test high_confidence_lock >= low_confidence_lock
            @test low_confidence_lock % L1CA_BLOCKS_PER_BIT == 0
            @test high_confidence_lock % L1CA_BLOCKS_PER_BIT == 0
        end
    end

    @testset "confidence = 1.0 stays conservative, not an instant lock" begin
        # With noise (finite z) confidence 1.0 is maximally conservative; a
        # noiseless edge (z = Inf) still locks.
        rng = MersenneTwister(7)
        clean = _bitstream([bit % 2 for bit = 0:39]; amp = 3.0)
        noisy = ComplexF64[prompt + complex(randn(rng), randn(rng)) for prompt in clean]
        low_confidence_lock = _detect_over(noisy, 0.99)
        high_confidence_lock = _detect_over(noisy, 1.0)
        @test low_confidence_lock > 0
        @test high_confidence_lock == 0 || high_confidence_lock > low_confidence_lock
    end

    @testset "no false lock on a long near-constant run (Welford stability)" begin
        # Large, slowly drifting prompts with no edge. A naive Σe² − (Σe)²/M
        # variance could cancel to zero (infinite z); Welford's does not.
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
    # `period` is both the bin length and the hypothesis count.
    period = 5

    @testset "Needs two bins on some hypothesis" begin
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
        # Same gap, but a huge spread on the peak keeps z below threshold.
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

# Block `i` (0-based) carries secondary chip `(i + start_chip) % N` times a
# random per-period data symbol times `amp`, plus optional noise.
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

# As `_detect_over`, for the secondary-code detector.
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
    # Two clean NH10 periods from chip 0: only rotation 0 wipes the overlay and
    # sums coherently to (N·amp)².
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
    @test peak == 0
    @test energies[peak+1] ≈ (N * amp)^2
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
            # The just-processed block is chip N-1.
            @test (start_chip + synced - 1) % N == N - 1
            res = _secondary_detect_over(prompts, signal, prn, 0.999; upto = synced)
            @test res.found
            @test res.phase == 0                   # upcoming integration at chip 0
        end
    end

    @testset "Polarity follows the winning period's sign" begin
        chips = _secondary_chips(signal, prn)
        pos = ComplexF64[8.0 * c for c in chips]
        neg = ComplexF64[-8.0 * c for c in chips]
        @test _secondary_detect_over(vcat(pos, pos, pos), signal, prn, 0.999; upto = 3N).polarity ==
              +1
        @test _secondary_detect_over(vcat(neg, neg, neg), signal, prn, 0.999; upto = 3N).polarity ==
              -1
    end

    @testset "No false lock on pure noise (the TEX-CUP PRN 8 symptom)" begin
        # The hard exact match accepts a random 10-bit window ~2 % of the time
        # (one of ~20 rotation×polarity templates), which adds up over a long
        # pre-lock stretch; the soft detector needs an energy gap over periods.
        hard_false_locks = 0
        soft_locks = 0
        for seed = 1:200
            rng = MersenneTwister(1000 + seed)
            noise = ComplexF64[complex(randn(rng), randn(rng)) for _ = 1:4N]
            _secondary_detect_over(noise, signal, prn, 0.999) != 0 && (soft_locks += 1)
            window = rand(rng, UInt32) & UInt32((1 << N) - 1)
            _secondary_code_search(window, _packed_l5i(), N, 0).found &&
                (hard_false_locks += 1)
        end
        @test soft_locks == 0
        @test hard_false_locks > 0
    end
end

# Feed prompts one block at a time through `buffer()`; return the 1-based lock
# block (0 if none) and the final bit buffer.
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
    found_at, bit_buffer = _feed_prompts([fill(-1.0 + 0.0im, 40); fill(1.0 + 0.0im, 20)])
    @test found_at == 60
    @test bit_buffer.polarity == +1

    found_at, bit_buffer = _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
    @test found_at == 60
    @test bit_buffer.polarity == -1

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
        # Block signs +,+,-: the last bin is negative, so the polarity is -1,
        # applied to the pre-sync bits as to every later one (#127).
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
        # The sign-flipped stream: polarity absorbs the flip, so the soft bits
        # match the case above (#127).
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
        # Seeded -,+ plus the completed bit, -12.
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
        # Seeded +,+ plus the completed bit, +20.
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

            # The 20th block completes the bit; its soft bit is the real sum.
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
            # Pre-sync only signs are kept, so each recovered bit's vote count
            # is scaled by the sync-time prompt magnitude, in the units of the
            # post-sync sums.
            found_at, bit_buffer =
                _feed_prompts([fill(2.0 + 0.0im, 40); fill(-2.0 + 0.0im, 20)])
            @test found_at == 60
            @test length(bit_buffer) == 3
            soft = get_soft_bits(bit_buffer)
            @test soft == Float32[-40.0, -40.0, 40.0]

            _, unit_buffer = _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
            @test get_soft_bits(unit_buffer) == Float32[-20.0, -20.0, 20.0]
        end

        @testset "Reset empties the soft bits but keeps the vector" begin
            bit_buffer = BitBuffer()
            push!(get_soft_bits(bit_buffer), 1.0f0, -2.0f0)
            soft_bits = get_soft_bits(bit_buffer)

            reset_bit_buffer = reset(bit_buffer)
            @test isempty(get_soft_bits(reset_bit_buffer))
            @test get_soft_bits(reset_bit_buffer) === soft_bits
        end
    end

    @testset "A record crossing the bit boundary resynchronises (issue #238)" begin
        # A record stepping the accumulator past the boundary (18 + 3 = 21)
        # would never meet it again and stall the bits silently.
        signal = GPSL1CA()

        @testset "Sync is dropped and the completed bits are kept" begin
            soft_bits = Float32[1.0, -1.0]
            bit_buffer = BitBuffer(UInt64(0xff), 8, true, complex(9.0, 1.0), 18, soft_bits)
            _seed_phase_accumulators!(bit_buffer.phase_acc, 20)
            bit_buffer.phase_acc.mean_bin_energy[3] = 42.0
            bit_buffer.phase_acc.last_bin_polarity[3] = Int8(-1)

            # Reported here, logged by the fold
            # (`Tracking._warn_bit_boundary_overshoot`).
            next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 3, 1.0 + 0.0im)

            @test has_bit_or_secondary_code_been_found(next_bit_buffer) == false
            @test next_bit_buffer.prompt_accumulator == complex(0.0, 0.0)
            @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
            @test next_bit_buffer.code_block_buffer == 0
            @test next_bit_buffer.code_block_buffer_length == 0
            @test next_bit_buffer.polarity == 0
            @test next_bit_buffer.secondary_phase == 0
            # The straddling record spans two bits and commits none.
            @test get_soft_bits(next_bit_buffer) === soft_bits
            @test get_soft_bits(next_bit_buffer) == Float32[1.0, -1.0]
            # The old lock's bin statistics must not bias the new search.
        @test all(iszero, next_bit_buffer.phase_acc.mean_bin_energy)
            @test all(iszero, next_bit_buffer.phase_acc.open_bin_sum)
            @test all(
                iszero,
                next_bit_buffer.phase_acc.bin_energy_sum_of_squared_deviations,
            )
            @test all(iszero, next_bit_buffer.phase_acc.last_bin_polarity)
        end

        @testset "The overshoot is reported to the caller, which warns" begin
            # What a user sees when external records cross a bit boundary. A
            # PRN of its own: the warning is `maxlog = 1` per signal and PRN.
            @test_logs (:warn, r"past the navigation-bit boundary") match_mode = :any _warn_bit_boundary_overshoot(
                get_signal_id(signal),
                99,
                21,
                20,
            )
        end

        @testset "Landing exactly on the boundary still commits a bit" begin
            # The guard fires on `>`; landing on the boundary is normal.
            bit_buffer = BitBuffer(UInt64(0xff), 8, true, complex(10.0, 0.0), 17)
            next_bit_buffer = @inferred buffer(signal, 3, bit_buffer, 3, 3.0 + 0.0im)
            @test has_bit_or_secondary_code_been_found(next_bit_buffer) == true
            @test get_soft_bits(next_bit_buffer) == Float32[13.0]
            @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
        end

        @testset "Secondary-code signals resynchronise too" begin
            # GPS L5I: the secondary-code phase goes with the lock.
            l5i = GPSL5I()
            bit_buffer = BitBuffer{UInt16}(
                UInt16(0x3ff),
                10,
                true,
                4,
                Int8(1),
                complex(5.0, 0.0),
                8,
                Float32[2.0],
                PhaseAccumulators(),
            )
            next_bit_buffer = @inferred buffer(l5i, 4, bit_buffer, 3, 1.0 + 0.0im)
            @test has_bit_or_secondary_code_been_found(next_bit_buffer) == false
            @test next_bit_buffer.secondary_phase == 0
            @test next_bit_buffer.polarity == 0
            @test next_bit_buffer.prompt_accumulator_integrated_code_blocks == 0
            @test get_soft_bits(next_bit_buffer) == Float32[2.0]
        end

        @testset "Bits keep coming after the straddling record" begin
            # Sync, decode, inject one 3-block record at 18 accumulated blocks,
            # and keep feeding: the detector re-locks on the true grid.
            rng = MersenneTwister(1234)
            data_bits = rand(rng, (-1.0, 1.0), 80)
            prompts = ComplexF64.(repeat(data_bits, inner = 20))

            bit_buffer = BitBuffer{UInt64}()
            bits_at_injection = -1
            bits_at_relock = -1
            max_accumulated_after = 0
            injected = false
            k = 1
            while k <= length(prompts)
                if !injected &&
                   has_bit_or_secondary_code_been_found(bit_buffer) &&
                   bit_buffer.prompt_accumulator_integrated_code_blocks == 18 &&
                   k + 2 <= length(prompts)
                    bits_at_injection = length(get_soft_bits(bit_buffer))
                    # One record spanning three primary blocks: 18 + 3 = 21.
                    bit_buffer = buffer(signal, 5, bit_buffer, 3, sum(prompts[k:(k+2)]))
                    @test has_bit_or_secondary_code_been_found(bit_buffer) == false
                    injected = true
                    k += 3
                    continue
                end
                bit_buffer = buffer(signal, 5, bit_buffer, 1, prompts[k])
                if injected
                    if bits_at_relock < 0 &&
                       has_bit_or_secondary_code_been_found(bit_buffer)
                        bits_at_relock = length(get_soft_bits(bit_buffer))
                    end
                    max_accumulated_after = max(
                        max_accumulated_after,
                        bit_buffer.prompt_accumulator_integrated_code_blocks,
                    )
                end
                k += 1
            end

            @test injected
            @test bits_at_injection > 0
            @test bits_at_relock >= 0
            @test has_bit_or_secondary_code_been_found(bit_buffer)
            @test length(get_soft_bits(bit_buffer)) > bits_at_relock
            @test max_accumulated_after < 20
            # Every later bit sums 20 equal-sign blocks.
            @test all(≈(20.0f0), abs.(get_soft_bits(bit_buffer)[(bits_at_relock+1):end]))
        end
    end

    @testset "Bit accumulation is unbounded" begin
        # More than 128 bits in one fold (#134).
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
