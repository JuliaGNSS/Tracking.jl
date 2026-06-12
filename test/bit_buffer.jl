module BitBufferTest

using Test: @test, @testset, @inferred, @test_throws
using Random: MersenneTwister
using GNSSSignals: GPSL1CA
using Tracking:
    BitBuffer,
    buffer,
    reset,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    SyncResult,
    PhaseAccumulators,
    _norm_quantile,
    _detect_bit_edge_cfar,
    _seed_phase_accumulators!,
    _update_phase_accumulators!

@testset "SyncResult" begin
    r = @inferred SyncResult(false, 0, Int8(0))
    @test r.found == false
    @test r.phase == 0
    @test r.polarity == 0
end

@testset "_norm_quantile" begin
    # Symmetry and a few tabulated standard-normal quantiles.
    @test _norm_quantile(0.5) ≈ 0.0 atol = 1e-9
    @test _norm_quantile(0.975) ≈ 1.959963985 atol = 1e-6
    @test _norm_quantile(0.025) ≈ -1.959963985 atol = 1e-6
    @test _norm_quantile(0.999) ≈ 3.090232306 atol = 1e-6
    @test _norm_quantile(1 - 0.999) ≈ -3.090232306 atol = 1e-6
    # Monotonically increasing.
    @test _norm_quantile(0.1) < _norm_quantile(0.2) < _norm_quantile(0.8)
end

const L1CA_BLOCKS_PER_BIT = 20  # primary-code blocks per L1 C/A navigation bit

# Build a noiseless ±1 soft-prompt stream from a list of data bits, one bit =
# `L1CA_BLOCKS_PER_BIT` blocks, scaled by `amp`.
_bitstream(bits; amp = 1.0) =
    ComplexF64[amp * (b == 1 ? 1.0 : -1.0) for b in bits for _ in 1:L1CA_BLOCKS_PER_BIT]

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
            accumulators, ComplexF64(prompt), block_number - 1, L1CA_BLOCKS_PER_BIT,
        )
        result = _detect_bit_edge_cfar(accumulators, L1CA_BLOCKS_PER_BIT, confidence, block_number)
        upto == 0 && result.found && return block_number
        upto == block_number && return result
    end
    return upto == 0 ? 0 :
           _detect_bit_edge_cfar(accumulators, L1CA_BLOCKS_PER_BIT, confidence, length(prompts))
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
            clean = _bitstream([bit % 2 for bit in 0:9]; amp = 8.0)
            noisy = ComplexF64[prompt + complex(randn(rng), randn(rng)) for prompt in clean]
            _detect_over(noisy, confidence)
        end
        for seed in 1:5
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
        clean = _bitstream([bit % 2 for bit in 0:9]; amp = 3.0)
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
        for block_number in 1:5000
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
        @test bit_buffer.buffer == 0
        @test bit_buffer.length == 0
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
        @test next_bit_buffer.length == 0
        @test next_bit_buffer.buffer == 0
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
        # signs +,+,- decode as bits 0,0,1 = 0b001 = 1, with soft bits -,-,+.
        found_at, bit_buffer =
            _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
        @test found_at == 60
        @test bit_buffer.found == true
        @test bit_buffer.polarity == -1
        @test bit_buffer.length == 3
        @test bit_buffer.buffer == 1
        @test bit_buffer.code_block_buffer_length == 60
        soft = get_soft_bits(bit_buffer)
        @test length(soft) == 3
        @test soft[1] < 0 && soft[2] < 0 && soft[3] > 0
    end

    @testset "Find bit start and buffer the pre-sync bits (positive polarity)" begin
        # The whole-signal sign flip of the case above: data 0,0,1 (block
        # signs -,-,+). The last completed bin is the "1", so the lock is at
        # positive polarity. Because a global RF sign flip is exactly the
        # ambiguity polarity resolves, the recovered hard bits and soft-bit
        # signs are identical to the negative-polarity case (issue #127):
        # buffer 0b001 = 1, soft -,-,+.
        found_at, bit_buffer =
            _feed_prompts([fill(-1.0 + 0.0im, 40); fill(1.0 + 0.0im, 20)])
        @test found_at == 60
        @test bit_buffer.found == true
        @test bit_buffer.polarity == +1
        @test bit_buffer.length == 3
        @test bit_buffer.buffer == 1
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
            0,
            0,
            complex(-1, 0),
            1,
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 1, -2 + 0im)
        @test next_bit_buffer.buffer == 0
        @test next_bit_buffer.length == 0
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
            1,
            2,
            complex(-10, 2),
            19,
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 1, -2 + 0im)
        @test next_bit_buffer.buffer == 2
        @test next_bit_buffer.length == 3
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
            3,
            2,
            complex(10, 2),
            10,
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 10, 10 + 1im)
        @test next_bit_buffer.buffer == 7
        @test next_bit_buffer.length == 3
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
                1,
                2,
                complex(-10.0, 2.0),
                19,
            )
            signal = GPSL1CA()

            # 20th code block completes the bit; soft bit = real of the sum
            next_bit_buffer = buffer(signal, 1, bit_buffer, 1, -2 + 0im)
            @test get_soft_bits(next_bit_buffer) == Float32[-12.0]
            @test eltype(get_soft_bits(next_bit_buffer)) == Float32
            # The hard bit sign matches the soft bit sign
            @test next_bit_buffer.buffer == 2
        end

        @testset "No soft bit is stored before a bit completes" begin
            code_blocks_buffer = 0xfffffffffff0000
            code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
            bit_buffer = BitBuffer(
                code_blocks_buffer,
                code_blocks_buffer_length,
                true,
                0,
                0,
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
            @test bit_buffer.length == 3
            soft = get_soft_bits(bit_buffer)
            # 20 sign votes × |prompt| (2) = magnitude 40; signs match the
            # decoded hard bits (0,0,1 → soft -,-,+, as in the polarity test).
            @test soft == Float32[-40.0, -40.0, 40.0]

            # Unit-amplitude prompts give magnitude 20 — i.e. the magnitude
            # tracks |prompt|, confirming the scaling (a raw vote count would
            # be ±20 in both runs).
            _, unit_buffer =
                _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
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

    @testset "Hard-bit buffer overflows loudly past 128 bits (issue #134)" begin
        # The hard bits live in a fixed UInt128. Pushing a 129th bit used to
        # silently shift the oldest bit out while `get_num_bits` kept
        # counting — now it must throw a descriptive error instead.
        signal = GPSL1CA()
        bit_buffer = BitBuffer(UInt128(0), 0, true, 0, 0, complex(0.0, 0.0), 0)
        for _ in 1:(128*20)
            bit_buffer = buffer(signal, 1, bit_buffer, 1, 1.0 + 0.0im)
        end
        @test bit_buffer.length == 128
        # 19 more blocks are fine; the block completing the 129th bit throws.
        for _ in 1:19
            bit_buffer = buffer(signal, 1, bit_buffer, 1, 1.0 + 0.0im)
        end
        @test_throws "hard-bit buffer is full" buffer(signal, 1, bit_buffer, 1, 1.0 + 0.0im)
    end
end

end
