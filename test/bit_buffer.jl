module BitBufferTest

# The bit buffer's resync on a record that crosses the bit boundary, and the
# warning Tracking logs for it. The rest of `BitBuffer` and its sync detectors
# is tested in TrackingLoops (`test/bit_buffer.jl` there).

using Test: @test, @testset, @inferred, @test_logs
using Random: MersenneTwister
using GNSSSignals: GPSL1CA, GPSL5I, get_signal_id
using Tracking: _warn_bit_boundary_overshoot
using TrackingLoops:
    BitBuffer,
    buffer,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    PhaseAccumulators,
    _seed_phase_accumulators!

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
        @test all(iszero, next_bit_buffer.phase_acc.bin_energy_sum_of_squared_deviations)
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
                if bits_at_relock < 0 && has_bit_or_secondary_code_been_found(bit_buffer)
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

end
