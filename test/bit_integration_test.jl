module BitDetectionIntegrationTest

using Test: @test, @testset
using Unitful: Hz, s
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    GPSL5Q,
    GPSL1C_P,
    gen_code,
    get_code_frequency,
    get_code_length,
    get_secondary_code_length
using Tracking:
    TrackedSat,
    TrackState,
    track,
    get_code_phase,
    get_num_bits,
    get_filtered_prompts,
    get_last_fully_integrated_filtered_prompt,
    get_last_fully_integrated_num_code_blocks,
    set_preferred_num_code_blocks_to_integrate!
using TrackingLoops:
    get_soft_bits, has_bit_or_secondary_code_been_found, EarlyPromptLateCorrelator

@testset "Bit detection integration test" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)
    carrier_doppler = 0.0Hz

    # Both lock polarities. This near-noiseless signal locks at block 60, the
    # third bit boundary: the floating-point residual keeps the z-score finite,
    # so the Student-t threshold clears only with a third bin. It still fires on
    # a true boundary (#124). Pre-sync and post-sync bits must match the
    # transmitted ones up to one *global* inversion (#127).
    for data_bits in ([1, 0, 1, 1, 0, 0, 1], [0, 1, 0, 0, 1, 1, 0])
        track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, carrier_doppler)];)

        block_bits = repeat(data_bits, inner = 20)  # 20 code blocks per bit
        decoded_bits = Bool[]
        decoded_soft_bits = Float32[]
        foreach(enumerate(block_bits)) do (index, bit)
            code_phase = (index - 1) * num_samples * code_frequency / sampling_frequency
            carrier_phase =
                2π * (index - 1) * num_samples * carrier_doppler / sampling_frequency
            signal =
                (bit * 2 - 1) .* cis.(
                    2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency .+
                    carrier_phase,
                ) .* gen_code(
                    num_samples,
                    gpsl1,
                    1,
                    sampling_frequency,
                    code_frequency,
                    code_phase,
                )
            track_state = track(signal, track_state, sampling_frequency)
            @test has_bit_or_secondary_code_been_found(track_state) == (index >= 60)
            # Sync at block 60 recovers the 3 buffered bits (the UInt64 sign
            # window holds all 60), then one bit per 20 blocks.
            expected_num_bits = index == 60 ? 3 : (index > 60 && index % 20 == 0 ? 1 : 0)
            num_bits = get_num_bits(track_state)
            @test num_bits == expected_num_bits
            soft_bits = get_soft_bits(track_state)
            @test eltype(soft_bits) == Float32
            @test length(soft_bits) == num_bits
            append!(decoded_bits, soft_bits .> 0)
            append!(decoded_soft_bits, soft_bits)
        end

        @test length(decoded_bits) == length(data_bits)
        # A global inversion at most, never a mixed one (#127).
        @test decoded_bits == Bool.(data_bits) || decoded_bits == .!Bool.(data_bits)
        @test (decoded_soft_bits .> 0) == decoded_bits
    end
end

# Multi-block coherent integration: bits keep flowing with `preferred > 1`.
#
# The integration is one block until sync (block 60), then 4 blocks, five per
# bit; a length that does not divide the bit, like 3, is rejected (#128). The
# multi-block tracker must decode the same bits as a single-block reference.
@testset "Bit detection with multi-block coherent integration" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)

    multi_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])
    set_preferred_num_code_blocks_to_integrate!(multi_state, 1, 4)
    single_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])

    # The bit buffer is flushed on every `track` call.
    multi_bits = Bool[]
    single_bits = Bool[]
    multi_soft = Float32[]
    single_soft = Float32[]

    bits = vcat(ones(20), zeros(20), ones(20), zeros(20), ones(20))
    foreach(enumerate(bits)) do (index, bit)
        code_phase = (index - 1) * num_samples * code_frequency / sampling_frequency
        signal = ComplexF32.(
            (bit * 2 - 1) .* gen_code(
                num_samples,
                gpsl1,
                1,
                sampling_frequency,
                code_frequency,
                code_phase,
            ),
        )
        multi_state = track(signal, multi_state, sampling_frequency)
        single_state = track(signal, single_state, sampling_frequency)
        @test has_bit_or_secondary_code_been_found(multi_state) == (index >= 60)
        expected_num_bits = index == 60 ? 3 : (index > 60 && index % 20 == 0 ? 1 : 0)
        @test get_num_bits(multi_state) == expected_num_bits
        @test get_num_bits(single_state) == expected_num_bits
        append!(multi_bits, get_soft_bits(multi_state) .> 0)
        append!(single_bits, get_soft_bits(single_state) .> 0)
        append!(multi_soft, get_soft_bits(multi_state))
        append!(single_soft, get_soft_bits(single_state))
    end

    # The polarity between the replayed (1-3) and streamed (4-5) bits is not
    # asserted.
    @test length(multi_bits) == 5
    @test multi_bits == single_bits
    @test multi_bits[1] != multi_bits[2]
    @test multi_bits[3] != multi_bits[4]
    @test multi_bits[4] != multi_bits[5]

    # Post-sync a bit sums five 4-block prompts (~5) instead of twenty
    # single-block ones (~20).
    @test sign.(multi_soft) == sign.(single_soft)
    @test all(x -> abs(x) ≈ 5, multi_soft[4:end])
    @test all(x -> abs(x) ≈ 20, single_soft[4:end])
end

# Mid-fold bit sync: the accumulation window must stay on the bit grid.
#
# With a fold of several records, sync usually lands on a record that is not
# the last of its fold. The `k` records behind it still cover real code blocks;
# leaving them out of the block count shifts every later bit window `k` blocks
# off the grid (#219), which reads as a soft-bit magnitude of 20 − 2k, not 20.
@testset "mid-fold bit sync keeps the accumulation on the bit grid" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)
    num_blocks = 140

    signal = ComplexF32[]
    for index = 1:num_blocks
        # ±1 data bit, flipping every 20 code blocks.
        data_bit_sign = iseven(div(index - 1, 20)) ? 1.0 : -1.0
        code_phase = (index - 1) * num_samples * code_frequency / sampling_frequency
        append!(
            signal,
            ComplexF32.(
                data_bit_sign .* gen_code(
                    num_samples,
                    gpsl1,
                    1,
                    sampling_frequency,
                    code_frequency,
                    code_phase,
                ),
            ),
        )
    end

    # 1 ms is the baseline (one record per fold). The lock at block 60 is the
    # 4th of the 57..63 fold at 7 ms and the 8th of the 53..65 fold at 13 ms, so
    # 3 and 5 records trail it.
    for doppler_update_interval in (1e-3s, 3e-3s, 7e-3s, 13e-3s)
        track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])
        track_state =
            track(signal, track_state, sampling_frequency; doppler_update_interval)
        soft_bits = get_soft_bits(track_state, 1)
        # 3 replayed at the block-60 lock, then one per 20 blocks.
        @test length(soft_bits) == 3 + div(num_blocks - 60, 20)
        @test all(x -> abs(x) ≈ 20, soft_bits)
    end
end

# GPS L5I secondary-code (NH10) sync and phase recovery.
#
# Starting at every NH10 phase, the soft CFAR rotation search
# (`TrackingLoops._detect_secondary_code_cfar`) must not fire before two periods,
# must fire on a true NH10 boundary, and must anchor `code_phase` to the upcoming
# integration's chip 0.
@testset "GPS L5I secondary-code sync and phase recovery" begin
    gpsl5 = GPSL5I()
    prn = 1
    sampling_frequency = 25e6Hz   # must exceed L5I's 10.23 Mcps code rate
    code_frequency = get_code_frequency(gpsl5)
    primary_code_length = get_code_length(gpsl5)            # 10230 chips
    secondary_code_length = get_secondary_code_length(gpsl5)  # 10 (NH10)
    num_samples = round(Int, 25e6 / 1000)                  # 1 ms = one primary-code block

    for start_secondary_chip = 0:(secondary_code_length-1)
        track_state = TrackState(gpsl5, [TrackedSat(gpsl5, prn, 0.0, 0.0Hz)])
        synced_at_block = -1
        synced_code_phase = NaN
        for index = 1:40
            # Block `index` carries NH10 chip `(start_secondary_chip + index - 1) % 10`.
            gen_code_phase = (start_secondary_chip + (index - 1)) * primary_code_length
            signal = ComplexF32.(
                gen_code(
                    num_samples,
                    gpsl5,
                    prn,
                    sampling_frequency,
                    code_frequency,
                    gen_code_phase,
                ),
            )
            track_state = track(signal, track_state, sampling_frequency)
            if has_bit_or_secondary_code_been_found(track_state) && synced_at_block < 0
                synced_at_block = index
                synced_code_phase = get_code_phase(track_state)
            end
        end

        # A runner-up, and so a decision, needs two completed periods.
        @test synced_at_block >= 2 * secondary_code_length

        # The just-completed block is absolute chip 9: the true secondary-code
        # phase was recovered, not merely some lock.
        @test (start_secondary_chip + synced_at_block - 1) % secondary_code_length ==
              secondary_code_length - 1

        # The generator's fixed-point DDA (~2^-30 chip) leaves the phase within
        # ~1e-6 chip of the boundary; compare circularly.
        let wrap = primary_code_length * secondary_code_length,
            d = mod(synced_code_phase, wrap)

            @test min(d, wrap - d) < 1e-3
        end
    end

    # Start half a primary period into a secondary chip. The phase snap at sync
    # must keep the within-block phase (#117), so a 60-block buffer ends at
    # chip `k mod 10` plus the half block still in flight.
    for k in (0, 1, 5, 9, 13)
        start_phase = (k + 0.5) * primary_code_length
        num_blocks = 60
        signal = ComplexF32.(
            gen_code(
                num_blocks * num_samples,
                gpsl5,
                prn,
                sampling_frequency,
                code_frequency,
                start_phase,
            ),
        )
        track_state = TrackState(gpsl5, [TrackedSat(gpsl5, prn, start_phase, 0.0Hz)])
        track_state = track(signal, track_state, sampling_frequency)
        @test has_bit_or_secondary_code_been_found(track_state)
        # Up to a floating-point residual (~5e-7 chip) of the code loop.
        @test get_code_phase(track_state) ≈
              (k % secondary_code_length) * primary_code_length + 0.5 * primary_code_length atol =
            1e-4
    end
end

# GPS L5I post-sync data-bit decoding (#125).
#
# The decoder must wipe the NH chip from each prompt (else a 10-block bit sums to
# `data x Σ(NH signs)`, ≈2 instead of ≈10) and start its bit grid on the NH10 =
# data-bit boundary. Secondary-code lock leaves a ±1 polarity ambiguity (resolved
# downstream by the preamble), so bits are checked up to a global inversion.
@testset "GPS L5I post-sync data-bit decoding (issue #125)" begin
    gpsl5 = GPSL5I()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5)
    primary_code_length = get_code_length(gpsl5)              # 10230 chips
    secondary_code_length = get_secondary_code_length(gpsl5)  # 10 (NH10)
    num_samples = round(Int, 25e6 / 1000)                     # 1 ms = one primary block

    # One bit per NH10 period; the first two are equal so a lock window
    # straddling their boundary stays clean.
    data_bits = [1, 1, 0, 1, 0, 0, 1, 1, 0, 1]

    @testset "preferred = $preferred_blocks, start chip $start_secondary_chip" for preferred_blocks in
                                                                                   (1, 10),
        start_secondary_chip in (0, 3, 9)
        # One `track` call: fixed per-call chunks would drift against the loop's
        # sub-sample code Doppler and drop bits.
        total_blocks = secondary_code_length * length(data_bits) - start_secondary_chip
        signal = ComplexF32.(
            gen_code(
                total_blocks * num_samples,
                gpsl5,
                prn,
                sampling_frequency,
                code_frequency,
                start_secondary_chip * primary_code_length,
            ),
        )
        for b = 0:(total_blocks-1)
            abs_block = start_secondary_chip + b
            data_bit = data_bits[div(abs_block, secondary_code_length)+1]
            block_range = (b*num_samples+1):((b+1)*num_samples)
            @views signal[block_range] .*= ComplexF32(2 * data_bit - 1)
        end

        track_state = TrackState(
            gpsl5,
            [TrackedSat(gpsl5, prn, start_secondary_chip * primary_code_length, 0.0Hz)],
        )
        set_preferred_num_code_blocks_to_integrate!(track_state, prn, preferred_blocks)
        track_state = track(signal, track_state, sampling_frequency)

        @test has_bit_or_secondary_code_been_found(track_state)

        # Deterministic on this noiseless signal: the winning rotation is
        # `d* = mod(N − start, N)`, so the lock lands at the first block count
        # `≥ 2N` with `count % N == d*`, and decoding starts at data bit
        # `div(start + lock, N)`.
        N = secondary_code_length
        d_star = mod(N - start_secondary_chip, N)
        lock_block = let m = 2N
            while m % N != d_star
                m += 1
            end
            m
        end
        first_bit = div(start_secondary_chip + lock_block, N) + 1  # 1-based data-bit index
        expected_from_lock = data_bits[first_bit:end]              # decoded run, in order

        soft_bits = get_soft_bits(track_state)
        decoded_bits = [Int(s > 0) for s in soft_bits]
        n_decoded = length(decoded_bits)

        # The last bit may still be in flight when `track` returns.
        @test n_decoded in (length(expected_from_lock) - 1, length(expected_from_lock))
        expected = expected_from_lock[1:n_decoded]
        @test decoded_bits == expected || decoded_bits == 1 .- expected
        @test all((soft_bits .> 0) .== (decoded_bits .== 1))
        # NH wiped: every bit, the first included, sums one unit prompt per
        # integration.
        full_bit_magnitude = secondary_code_length / preferred_blocks
        @test all(>(0.9 * full_bit_magnitude), abs.(soft_bits))
    end
end

# Mid-fold secondary-code sync: the overlay anchor must travel with the fold.
#
# The secondary-code side of #219: the records trailing the syncing one in its
# fold also advance the secondary-chip anchor, since the code-phase snap runs
# after the fold. A stale anchor wipes the wrong NH chip from every block and the
# NH10 sum collapses from ~10 to ~0–2.
@testset "mid-fold secondary-code sync keeps the overlay anchor aligned" begin
    gpsl5 = GPSL5I()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5)
    primary_code_length = get_code_length(gpsl5)               # 10230 chips
    secondary_code_length = get_secondary_code_length(gpsl5)   # 10 (NH10)
    num_samples = round(Int, 25e6 / 1000)                      # 1 ms = one primary block
    data_bits = [1, 1, 0, 1, 0, 0, 1, 1, 0, 1]

    function decode(doppler_update_interval, start_secondary_chip)
        total_blocks = secondary_code_length * length(data_bits) - start_secondary_chip
        signal = ComplexF32.(
            gen_code(
                total_blocks * num_samples,
                gpsl5,
                prn,
                sampling_frequency,
                code_frequency,
                start_secondary_chip * primary_code_length,
            ),
        )
        for b = 0:(total_blocks-1)
            data_bit = data_bits[div(start_secondary_chip+b, secondary_code_length)+1]
            @views signal[(b*num_samples+1):((b+1)*num_samples)] .*=
                ComplexF32(2 * data_bit - 1)
        end
        track_state = TrackState(
            gpsl5,
            [TrackedSat(gpsl5, prn, start_secondary_chip * primary_code_length, 0.0Hz)],
        )
        get_soft_bits(
            track(signal, track_state, sampling_frequency; doppler_update_interval),
        )
    end

    @testset "start chip $start_secondary_chip" for start_secondary_chip in (0, 3)
        # Reference: one record per fold. The code-amplitude normalisation puts
        # magnitudes a few 1e-4 below nominal, hence the 0.99 floors.
        full_bit = 0.99 * secondary_code_length
        reference = copy(decode(1e-3s, start_secondary_chip))
        @test all(>(full_bit), abs.(reference))

        for doppler_update_interval in (3e-3s, 7e-3s)
            soft_bits = decode(doppler_update_interval, start_secondary_chip)
            @test length(soft_bits) == length(reference)
            @test sign.(soft_bits) == sign.(reference)
            # Only the sync bit is short, by the one dropped prompt.
            @test all(>(full_bit), abs.(soft_bits[2:end]))
            @test abs(soft_bits[1]) > 0.99 * (secondary_code_length - 1)
        end
    end
end

# Mid-fold secondary-code sync on a *pilot* (GPS L5Q, NH20).
#
# A pilot has no bit accumulator, but the stale overlay anchor of #219 hits it
# harder: the 20-block coherent sum cancels toward zero, the discriminators see
# 0/0 and the Doppler goes NaN, so this testset *errors* (an `InexactError` in
# the correlator) rather than fails when the anchor is wrong.
@testset "mid-fold pilot overlay anchor keeps coherent integration intact" begin
    gpsl5q = GPSL5Q()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5q)
    secondary_code_length = get_secondary_code_length(gpsl5q)   # 20 (NH20)
    num_samples = round(Int, 25e6 / 1000)                       # 1 ms = one primary block
    # Sync needs ≥ 2 NH20 periods; the rest are full-period integrations.
    num_blocks = 240
    signal = ComplexF32.(
        gen_code(
            num_blocks * num_samples,
            gpsl5q,
            prn,
            sampling_frequency,
            code_frequency,
            0.0,
        ),
    )

    for doppler_update_interval in (1e-3s, 3e-3s, 7e-3s)
        track_state = TrackState(gpsl5q, [TrackedSat(gpsl5q, prn, 0.0, 0.0Hz)])
        set_preferred_num_code_blocks_to_integrate!(track_state, prn, secondary_code_length)
        track_state =
            track(signal, track_state, sampling_frequency; doppler_update_interval)

        @test has_bit_or_secondary_code_been_found(track_state)
        @test get_last_fully_integrated_num_code_blocks(track_state, prn) ==
              secondary_code_length
        # The overlay is wiped across all 20 blocks, so the normalised prompt
        # keeps its full magnitude.
        @test abs(get_last_fully_integrated_filtered_prompt(track_state, prn)) > 0.99
        # The last five prompts all postdate sync.
        prompts = get_filtered_prompts(track_state, prn)
        @test all(>(0.99), abs.(prompts[(end-4):end]))
    end
end

# GPS L1C-P 1800-chip overlay sync — full 18 s end-to-end run.
#
# The rotation search over the per-PRN 1800-chip overlay must lock at the 1800th
# block, one full 18 s cycle, anchored to overlay chip 0.
@testset "GPS L1C-P overlay sync (18 s end-to-end)" begin
    gpsl1c_p = GPSL1C_P()
    prn = 1
    # TMBOC needs fs > ~12.28 MHz; 13 MHz keeps the 18 s run light.
    sampling_frequency = 13e6Hz
    code_frequency = get_code_frequency(gpsl1c_p)
    primary_code_length = get_code_length(gpsl1c_p)             # 10230 chips
    secondary_code_length = get_secondary_code_length(gpsl1c_p)  # 1800
    period_samples = round(Int, 13e6 / 100)                    # 10 ms = one primary period

    track_state = TrackState(
        gpsl1c_p,
        [
            TrackedSat(
                gpsl1c_p,
                prn,
                0.0,
                0.0Hz;
                correlator = EarlyPromptLateCorrelator(;
                    preferred_early_late_to_prompt_code_shift = 0.1,
                ),
            ),
        ],
    )
    synced_at_block = -1
    synced_code_phase = NaN
    for index = 1:(secondary_code_length+5)
        gen_code_phase = (index - 1) * primary_code_length
        signal = ComplexF32.(
            gen_code(
                period_samples,
                gpsl1c_p,
                prn,
                sampling_frequency,
                code_frequency,
                gen_code_phase,
            ),
        )
        track_state = track(signal, track_state, sampling_frequency)
        if has_bit_or_secondary_code_been_found(track_state) && synced_at_block < 0
            synced_at_block = index
            synced_code_phase = get_code_phase(track_state)
        end
    end

    @test synced_at_block == secondary_code_length
    # Chip 0 again, up to the residual of 1800 code-phase updates (~1e-7 chip);
    # the snap keeps the within-block phase (#117), it does not zero it.
    @test synced_code_phase ≈ 0.0 atol = 1e-4
end

end
