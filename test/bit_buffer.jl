module BitBufferTest

using Test: @test, @testset, @inferred, @test_throws
using GNSSSignals: GPSL1CA
using Tracking:
    BitBuffer,
    buffer,
    reset,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    SyncResult,
    _try_match

@testset "SyncResult" begin
    r = @inferred SyncResult(false, 0, Int8(0))
    @test r.found == false
    @test r.phase == 0
    @test r.polarity == 0
end

@testset "_try_match" begin
    # Exact match at positive polarity.
    res = @inferred _try_match(UInt8(0x0f), UInt8(0x0f), UInt8(0xff), 0)
    @test res.found == true
    @test res.polarity == +1
    @test res.phase == 0

    # Negative polarity = template XOR mask.
    res = @inferred _try_match(UInt8(0xf0), UInt8(0x0f), UInt8(0xff), 0)
    @test res.found == true
    @test res.polarity == -1

    # Single bit-flip at tolerance 0 — rejected.
    res = @inferred _try_match(UInt8(0x0e), UInt8(0x0f), UInt8(0xff), 0)
    @test res.found == false

    # Same single bit-flip at tolerance 1 — accepted.
    res = @inferred _try_match(UInt8(0x0e), UInt8(0x0f), UInt8(0xff), 1)
    @test res.found == true
    @test res.polarity == +1

    # Out-of-band bits don't affect the match (the mask hides them).
    res = @inferred _try_match(UInt32(0xdeadbe0f), UInt32(0x0f), UInt32(0xff), 0)
    @test res.found == true
    @test res.polarity == +1

    # Width-generic: works on UInt128 buffers too.
    res = @inferred _try_match(UInt128(0x0f), UInt128(0x0f), UInt128(0xff), 0)
    @test res.found == true
    @test res.polarity == +1

    # Edge mask: an error inside the edge mask rejects even within budget…
    res = @inferred _try_match(UInt8(0x07), UInt8(0x0f), UInt8(0xff), 1, UInt8(0x18))
    @test res.found == false
    # …an error on the other edge-mask bit rejects too…
    res = @inferred _try_match(UInt8(0x1f), UInt8(0x0f), UInt8(0xff), 1, UInt8(0x18))
    @test res.found == false
    # …while the same budget spent outside the edge mask still matches.
    res = @inferred _try_match(UInt8(0x0e), UInt8(0x0f), UInt8(0xff), 1, UInt8(0x18))
    @test res.found == true
    @test res.polarity == +1
    # Negative polarity honors the edge mask as well.
    res = @inferred _try_match(UInt8(0xe0), UInt8(0x0f), UInt8(0xff), 1, UInt8(0x18))
    @test res.found == false
    res = @inferred _try_match(UInt8(0xf1), UInt8(0x0f), UInt8(0xff), 1, UInt8(0x18))
    @test res.found == true
    @test res.polarity == -1
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

@testset "L1CA bit-edge lock is not one block early (issue #124)" begin
    # Data bits 0, 0, 1 — the first transition is preceded by a repeated
    # bit. The buggy tolerant matcher fired at block 59 (one block before
    # the true edge); an edge-locked detector fires exactly at block 60.
    prompts = [fill(-1.0 + 0.0im, 40); fill(1.0 + 0.0im, 20)]
    found_at, bit_buffer = _feed_prompts(prompts)
    @test found_at == 60
    @test bit_buffer.polarity == +1

    # Mirror pattern 1, 1, 0 locks at the same block with negative polarity.
    found_at, bit_buffer = _feed_prompts([fill(1.0 + 0.0im, 40); fill(-1.0 + 0.0im, 20)])
    @test found_at == 60
    @test bit_buffer.polarity == -1
end

@testset "L1CA sync tolerance budget through buffer (issue #124)" begin
    # One flipped block (block 10) away from the edge — exactly at the
    # 1-error budget. Sync still locks at the true edge (block 40).
    prompts = [fill(-1.0 + 0.0im, 20); fill(1.0 + 0.0im, 20)]
    prompts[10] = 1.0 + 0.0im
    found_at, _ = _feed_prompts(prompts)
    @test found_at == 40

    # Two flipped blocks — budget + 1 — must not lock anywhere in the stream.
    prompts = [fill(-1.0 + 0.0im, 20); fill(1.0 + 0.0im, 40)]
    prompts[8] = 1.0 + 0.0im
    prompts[10] = 1.0 + 0.0im
    found_at, _ = _feed_prompts(prompts)
    @test found_at == 0
end

@testset "Bit buffer" begin
    @testset "Initialize" begin
        bit_buffer = @inferred BitBuffer()
        @test bit_buffer.code_block_buffer == 0
        @test bit_buffer.code_block_buffer_lengh == 0
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
        @test next_bit_buffer.code_block_buffer_lengh == 1
    end

    @testset "Find bit start and buffer bit" begin
        # `UInt128` literal so that the expected post-shift value
        # 0x1fffffffffff00000 (65 bits) fits without truncation. The
        # buffer width is now a per-`BitBuffer{B}` parameter, so the
        # test pins `B = UInt128` explicitly.
        code_blocks_buffer = UInt128(0xfffffffffff80000)
        code_blocks_buffer_length = ndigits(code_blocks_buffer; base = 2)
        bit_buffer = BitBuffer(
            code_blocks_buffer,
            code_blocks_buffer_length,
            false,
            0,
            0,
            complex(0, 0),
            0,
        )
        signal = GPSL1CA()

        next_bit_buffer = @inferred buffer(signal, 1, bit_buffer, 1, -2 + 0im)
        @test next_bit_buffer.found == true
        @test next_bit_buffer.length == 3
        @test next_bit_buffer.buffer == 6
        @test next_bit_buffer.code_block_buffer == 0x1fffffffffff00000
        @test next_bit_buffer.code_block_buffer_lengh == code_blocks_buffer_length + 1
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
end

end
