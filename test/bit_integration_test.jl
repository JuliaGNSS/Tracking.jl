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
    get_sat_state,
    get_code_phase,
    get_num_bits,
    get_soft_bits,
    get_filtered_prompts,
    get_last_fully_integrated_filtered_prompt,
    get_last_fully_integrated_num_code_blocks,
    has_bit_or_secondary_code_been_found,
    set_preferred_num_code_blocks_to_integrate!,
    EarlyPromptLateCorrelator

@testset "Bit detection integration test" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)
    carrier_doppler = 0.0Hz

    # Both lock polarities. This near-noiseless signal locks at block 60: its
    # finite z clears the Student-t threshold only at the third bin (issue #124:
    # never one block early). Pre- and post-sync bits must match the transmitted
    # ones up to a single global inversion (issue #127).
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
            # 3 buffered bits at sync (block 60), then one per 20 blocks.
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
        @test decoded_bits == Bool.(data_bits) || decoded_bits == .!Bool.(data_bits)
        @test (decoded_soft_bits .> 0) == decoded_bits
    end
end

# Multi-block coherent integration (`preferred = 4`, clamped to 1 until sync):
# post-sync, five 4-block integrations form a bit, and the decoded stream must
# match a single-block reference tracker. Exercises the widened code wrap, the
# multi-block boundary calculation and the loop-bandwidth caps end to end.
@testset "Bit detection with multi-block coherent integration" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)

    multi_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])
    set_preferred_num_code_blocks_to_integrate!(multi_state, 1, 4)
    single_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])

    # The bit buffer is flushed on every `track` call, so collect the bits
    # and soft bits across calls.
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
        # 3 buffered bits at sync (block 60), then one per 20 blocks.
        expected_num_bits = index == 60 ? 3 : (index > 60 && index % 20 == 0 ? 1 : 0)
        @test get_num_bits(multi_state) == expected_num_bits
        @test get_num_bits(single_state) == expected_num_bits
        append!(multi_bits, get_soft_bits(multi_state) .> 0)
        append!(single_bits, get_soft_bits(single_state) .> 0)
        append!(multi_soft, get_soft_bits(multi_state))
        append!(single_soft, get_soft_bits(single_state))
    end

    # All five bits match the reference; adjacent bits alternate (the pair 2/3
    # is not asserted).
    @test length(multi_bits) == 5
    @test multi_bits == single_bits
    @test multi_bits[1] != multi_bits[2]
    @test multi_bits[3] != multi_bits[4]
    @test multi_bits[4] != multi_bits[5]

    # Post-sync bits (from index 4) sum five 4-block prompts (≈5) vs. twenty
    # single-block prompts (≈20).
    @test sign.(multi_soft) == sign.(single_soft)
    @test all(x -> abs(x) ≈ 5, multi_soft[4:end])
    @test all(x -> abs(x) ≈ 20, single_soft[4:end])
end

# Mid-fold bit sync: the records trailing the syncing one inside its fold must
# still count towards the bit grid (issue #219). With the data flipping every 20
# blocks, a window `k` blocks off the grid has magnitude 20 − 2k instead of 20.
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

    # 1 ms: baseline, nothing trails the sync. 7 / 13 ms: the block-60 lock
    # has 3 / 5 records trailing it in its fold.
    for doppler_update_interval in (1e-3s, 3e-3s, 7e-3s, 13e-3s)
        track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])
        track_state =
            track(signal, track_state, sampling_frequency; doppler_update_interval)
        soft_bits = get_soft_bits(track_state, 1)
        # 3 replayed bits at the block-60 lock, then one per 20 blocks.
        @test length(soft_bits) == 3 + div(num_blocks - 60, 20)
        # Every bit sums a full, unstraddled 20 blocks.
        @test all(x -> abs(x) ≈ 20, soft_bits)
    end
end

# GPS L5I NH10 sync from every start chip (see `_detect_secondary_code_cfar`):
# no lock before two periods, lock on a true NH10 boundary, and `code_phase`
# anchored to the upcoming integration's chip 0.
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

        # No runner-up, hence no decision, before two NH10 periods.
        @test synced_at_block >= 2 * secondary_code_length

        # The syncing block must be absolute chip 9: checks the true secondary
        # phase was recovered, not merely that some lock happened.
        @test (start_secondary_chip + synced_at_block - 1) % secondary_code_length ==
              secondary_code_length - 1

        # `code_phase` at chip 0 (`SyncResult.phase == 0`), up to the generator's
        # fixed-point residual (~1e-6 chip); compare circular distance.
        let wrap = primary_code_length * secondary_code_length,
            d = mod(synced_code_phase, wrap)

            @test min(d, wrap - d) < 1e-3
        end
    end

    # Start half a primary period into chip k: the sync-time phase snap must
    # keep the within-block phase (issue #117), so the final `code_phase` is
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
        # Tolerance covers the sub-microchip code-Doppler drift over the run.
        @test get_code_phase(track_state) ≈
              (k % secondary_code_length) * primary_code_length + 0.5 * primary_code_length atol =
            1e-4
    end
end

# GPS L5I post-sync data-bit decoding (issue #125): the replica must wipe the NH
# chip per block and bits must start on the NH10 boundary. Checked at several
# start chips, 1- and 10-block integration, up to the inherent global ±1 polarity.
@testset "GPS L5I post-sync data-bit decoding (issue #125)" begin
    gpsl5 = GPSL5I()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5)
    primary_code_length = get_code_length(gpsl5)              # 10230 chips
    secondary_code_length = get_secondary_code_length(gpsl5)  # 10 (NH10)
    num_samples = round(Int, 25e6 / 1000)                     # 1 ms = one primary block

    # One data bit per NH10 period (10 primary blocks). The first two bits are
    # equal so the rotation search locks cleanly on a 10-block window that may
    # straddle their boundary.
    data_bits = [1, 1, 0, 1, 0, 0, 1, 1, 0, 1]

    @testset "preferred = $preferred_blocks, start chip $start_secondary_chip" for preferred_blocks in
                                                                                   (1, 10),
        start_secondary_chip in (0, 3, 9)
        # One continuous signal fed in a single `track` call, so the code-Doppler
        # drift cannot straddle fixed per-call chunks and drop bits.
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

        # Noiseless, so the lock is exact: the smallest block count ≥ 2N on the
        # winning rotation `d* = mod(N − start, N)`; decoding starts at data bit
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

        # The final bit may still be in flight when `track` returns.
        @test n_decoded in (length(expected_from_lock) - 1, length(expected_from_lock))
        expected = expected_from_lock[1:n_decoded]
        @test decoded_bits == expected || decoded_bits == 1 .- expected
        @test all((soft_bits .> 0) .== (decoded_bits .== 1))
        # The load-bearing check: every bit (none truncated) sums coherently to
        # ≈ N / preferred_blocks; with the NH code left in it would be ≈2.
        full_bit_magnitude = secondary_code_length / preferred_blocks
        @test all(>(0.9 * full_bit_magnitude), abs.(soft_bits))
    end
end

# Mid-fold secondary-code sync (issue #219, secondary-code path): the records
# trailing the syncing one must also advance the secondary-chip anchor, or the
# replica bakes the wrong NH chip and the NH10 sum collapses from ~10 to ~0–2.
@testset "mid-fold secondary-code sync keeps the overlay anchor aligned" begin
    gpsl5 = GPSL5I()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5)
    primary_code_length = get_code_length(gpsl5)               # 10230 chips
    secondary_code_length = get_secondary_code_length(gpsl5)   # 10 (NH10)
    num_samples = round(Int, 25e6 / 1000)                      # 1 ms = one primary block
    data_bits = [1, 1, 0, 1, 0, 0, 1, 1, 0, 1]

    # Decode the whole run in a single `track` call at the given chunk interval.
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
        # 1 ms reference (nothing trails the sync). 0.99 floors allow for the
        # code-amplitude normalisation.
        full_bit = 0.99 * secondary_code_length
        reference = copy(decode(1e-3s, start_secondary_chip))
        @test all(>(full_bit), abs.(reference))

        for doppler_update_interval in (3e-3s, 7e-3s)
            soft_bits = decode(doppler_update_interval, start_secondary_chip)
            # Same bits, on the same grid, as the one-record-per-fold reference.
            @test length(soft_bits) == length(reference)
            @test sign.(soft_bits) == sign.(reference)
            # Only the sync bit is short, by the one dropped block; the rest are full.
            @test all(>(full_bit), abs.(soft_bits[2:end]))
            @test abs(soft_bits[1]) > 0.99 * (secondary_code_length - 1)
        end
    end
end

# Mid-fold overlay anchor on a pilot (GPS L5Q, NH20; issue #219). A wrong anchor
# cancels the 20-block sum, the loop diverges to a NaN Doppler, and this testset
# errors (`InexactError`) rather than fails.
@testset "mid-fold pilot overlay anchor keeps coherent integration intact" begin
    gpsl5q = GPSL5Q()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5q)
    secondary_code_length = get_secondary_code_length(gpsl5q)   # 20 (NH20)
    num_samples = round(Int, 25e6 / 1000)                       # 1 ms = one primary block
    # ≥ 2 NH20 periods for the detector, plus many full-period integrations.
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

    # 1 ms: nothing trails the sync. 3 / 7 ms: records do.
    for doppler_update_interval in (1e-3s, 3e-3s, 7e-3s)
        track_state = TrackState(gpsl5q, [TrackedSat(gpsl5q, prn, 0.0, 0.0Hz)])
        # Coherently integrate one full NH20 period — the point of a pilot lock.
        set_preferred_num_code_blocks_to_integrate!(track_state, prn, secondary_code_length)
        track_state =
            track(signal, track_state, sampling_frequency; doppler_update_interval)

        @test has_bit_or_secondary_code_been_found(track_state)
        # The long integration is really in effect...
        @test get_last_fully_integrated_num_code_blocks(track_state, prn) ==
              secondary_code_length
        # ...and the overlay is wiped across all 20 blocks (normalised prompt ≈ 1).
        @test abs(get_last_fully_integrated_filtered_prompt(track_state, prn)) > 0.99
        # The last few integrations are all post-sync and must hold up too.
        prompts = get_filtered_prompts(track_state, prn)
        @test all(>(0.99), abs.(prompts[(end-4):end]))
    end
end

# GPS L1C-P 1800-chip overlay sync, full 18 s run: a clean aligned signal locks
# at block 1800 with `code_phase` anchored to overlay chip 0.
@testset "GPS L1C-P overlay sync (18 s end-to-end)" begin
    gpsl1c_p = GPSL1C_P()
    prn = 1
    # L1C-P's TMBOC modulation needs fs > ~12.28 MHz; 13 MHz keeps the run
    # (1800+ x 10 ms periods) as light as possible while staying valid.
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
        # Continuous code phase so `gen_code` lays down the correct overlay
        # chip on each successive primary-code period.
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

    # Locks after exactly one overlay cycle (1800 blocks).
    @test synced_at_block == secondary_code_length
    # Chip 0 again, up to floating-point residual (the snap keeps the
    # within-block phase, issue #117).
    @test synced_code_phase ≈ 0.0 atol = 1e-4
end

end
