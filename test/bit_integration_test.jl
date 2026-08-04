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
    get_bits,
    get_num_bits,
    get_soft_bits,
    get_filtered_prompts,
    get_last_fully_integrated_filtered_prompt,
    has_bit_or_secondary_code_been_found,
    set_preferred_num_code_blocks_to_integrate!,
    EarlyPromptLateCorrelator

@testset "Bit detection integration test" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)
    carrier_doppler = 0.0Hz

    # Run the same bit stream at both lock polarities. Under the Student-t
    # small-sample threshold this near-noiseless Float64 signal locks at block
    # 60 (the third data-bit boundary): a 2-bin lock needs an exactly-infinite
    # z-score (zero bin-to-bin variance), and the floating-point residual here
    # yields a large-but-finite z that clears the threshold only once dof grows
    # to a third bin. It still fires at a true bit boundary (the #124 property —
    # never one block early). The bits recovered from the pre-sync buffer and
    # the bits decoded after sync must be consistent: equal to the transmitted
    # bits up to a single *global* inversion across the whole stream (issue #127).
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
            # window still holds all 60); afterwards one bit completes every 20
            # blocks.
            expected_num_bits = index == 60 ? 3 : (index > 60 && index % 20 == 0 ? 1 : 0)
            num_bits = get_num_bits(track_state)
            @test num_bits == expected_num_bits
            bits_word = get_bits(track_state)
            for bit_index = num_bits:-1:1
                push!(decoded_bits, (bits_word >> (bit_index - 1)) & 1 == 1)
            end
            soft_bits = get_soft_bits(track_state)
            # There is exactly one soft bit (Float32 accumulation) per hard bit
            @test eltype(soft_bits) == Float32
            @test length(soft_bits) == num_bits
            append!(decoded_soft_bits, soft_bits)
        end

        @test length(decoded_bits) == length(data_bits)
        # Pre-sync (first 2) and post-sync bits must agree on the symbol
        # mapping: the decoded stream matches the transmitted one up to a
        # global inversion, never a mixed one (issue #127).
        @test decoded_bits == Bool.(data_bits) || decoded_bits == .!Bool.(data_bits)
        # The sign of each soft bit must match the respective hard bit —
        # including the soft bits recovered from the pre-sync buffer.
        @test (decoded_soft_bits .> 0) == decoded_bits
    end
end

# Multi-block coherent integration: bits keep flowing with `preferred > 1`.
#
# Until bit sync the integration length is clamped to one code block, so the
# preferred length can be set from the start. After sync at block 60 each
# integration spans 4 code blocks, five integrations form a bit, and a hard +
# soft bit must be
# emitted at every 20-block bit boundary — with `preferred = 3` (which does
# not divide the 20 blocks per bit and is therefore rejected since issue
# #128) the integrations would straddle bit boundaries and bit emission
# would stall forever. A single-block tracker runs alongside as the
# reference: the multi-block tracker must decode the identical bit stream.
# This exercises the post-sync multi-block path end to end: the widened
# replica-code wrap, the multi-block integration boundary calculation, and
# the `1/N` loop-bandwidth scaling.
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
        # Three bits are recovered from the buffered history at sync (block 60,
        # the near-noiseless lock latency under the t-quantile threshold);
        # afterwards exactly one bit must appear at every 20-block boundary.
        expected_num_bits = index == 60 ? 3 : (index > 60 && index % 20 == 0 ? 1 : 0)
        @test get_num_bits(multi_state) == expected_num_bits
        @test get_num_bits(single_state) == expected_num_bits
        num_bits = get_num_bits(multi_state)
        for bit_index = (num_bits-1):-1:0
            push!(multi_bits, (get_bits(multi_state) >> bit_index) & 1 == 1)
            push!(single_bits, (get_bits(single_state) >> bit_index) & 1 == 1)
        end
        append!(multi_soft, get_soft_bits(multi_state))
        append!(single_soft, get_soft_bits(single_state))
    end

    # All five data bits arrive and match the single-block reference bit for
    # bit. The three bits replayed from the buffered history at sync alternate,
    # and so do the two streamed post-sync bits; the relative polarity
    # between the replayed and the streamed portion is a property of the
    # bit-sync polarity convention and not asserted here.
    @test length(multi_bits) == 5
    @test multi_bits == single_bits
    @test multi_bits[1] != multi_bits[2]
    @test multi_bits[3] != multi_bits[4]
    @test multi_bits[4] != multi_bits[5]

    # The soft bits agree with the reference in sign; post-sync each bit is
    # the sum of five 4-block prompts (magnitude ~5) instead of twenty
    # single-block prompts. Sync at block 60 replays three pre-sync bits, so the
    # streamed post-sync bits start at index 4.
    @test sign.(multi_soft) == sign.(single_soft)
    @test all(x -> abs(x) ≈ 5, multi_soft[4:end])
    @test all(x -> abs(x) ≈ 20, single_soft[4:end])
end

# Mid-fold bit sync: the accumulation window must stay on the bit grid.
#
# With a `doppler_update_interval` longer than one code period a chunk's fold
# covers several records, so bit sync is generally detected on a record that is
# NOT the last of its fold. The records behind it were correlated with pre-sync
# replicas — their prompt may be unusable — but the code blocks they cover are
# real: dropping them from the accumulator's block count slides every following
# bit window `k` blocks off the navigation-bit grid, where `k` is the number of
# trailing records, permanently (issue #219).
#
# The data bit flips every 20 blocks here, so a window `k` blocks off the grid
# sums 20 − k blocks of its own bit against k of the neighbour's and lands at
# magnitude 20 − 2k instead of 20 — the misalignment is read straight off the
# soft bits, with no dependence on where the detector happened to lock.
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

    # 1 ms: one record per fold, so no record can ever trail the syncing one —
    # the baseline. 7 ms / 13 ms: this near-noiseless signal locks at block 60,
    # which is the 4th of the 57..63 fold and the 8th of the 53..65 fold, so 3
    # respectively 5 records trail it. Before the fix those runs lost a bit and
    # decoded the rest at magnitude 14 / 10.
    for doppler_update_interval in (1e-3s, 3e-3s, 7e-3s, 13e-3s)
        track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0, 0.0Hz)])
        track_state =
            track(signal, track_state, sampling_frequency; doppler_update_interval)
        soft_bits = get_soft_bits(track_state, 1)
        # 3 bits replayed from the pre-sync sign window at the block-60 lock
        # plus one per completed 20-block bit for the rest of the run.
        @test length(soft_bits) == 3 + div(num_blocks - 60, 20)
        # Every bit sums a full, unstraddled 20 blocks.
        @test all(x -> abs(x) ≈ 20, soft_bits)
    end
end

# GPS L5I secondary-code (NH10) sync and phase recovery.
#
# The L5I detector runs the soft, maximum-energy CFAR rotation search over NH10
# (`_detect_secondary_code_cfar`): it accumulates per-rotation overlay-wiped
# period energies, needs at least two NH10 periods before it can decide, and
# fires only at the winning rotation's own period boundary — recovering the true
# secondary-code phase. This test starts the signal at every NH10 phase and
# checks that:
#   * sync never fires before two NH10 periods (20 blocks),
#   * it fires on a true NH10 boundary (the just-completed block is chip 9), so
#     the recovered phase is the actual absolute alignment, and
#   * `code_phase` is anchored to the upcoming integration's chip 0.
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
            # Advance the *generated* signal's absolute code phase from an
            # initial offset of `start_secondary_chip` primary periods, so
            # block `index` carries NH10 chip `(start_secondary_chip + index - 1) % 10`.
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

        # The soft maximum-energy CFAR detector needs at least two completed NH10
        # periods before a runner-up (and thus a decision) exists, so it never
        # locks before the 20th block.
        @test synced_at_block >= 2 * secondary_code_length

        # It fires only at the winning rotation's own period boundary, i.e. when
        # the just-completed block is the last NH10 chip (absolute chip 9), so the
        # *upcoming* integration starts a fresh NH10 period (absolute chip 0). The
        # generated block `synced_at_block` carries absolute chip
        # `(start_secondary_chip + synced_at_block - 1) % 10`; boundary firing
        # demands that be chip 9 — this is the load-bearing check that the true
        # secondary-code phase was recovered (the absolute alignment), not merely
        # that *some* lock happened.
        @test (start_secondary_chip + synced_at_block - 1) % secondary_code_length ==
              secondary_code_length - 1

        # Because it fires at that true boundary, `code_phase` is anchored to the
        # upcoming integration's chip 0 (`SyncResult.phase == 0`). The embedded-LUT
        # generator's fixed-point DDA (~2^-30 chip) lands the phase within ~1e-6
        # chip of the integer boundary, so compare the circular distance within a
        # sub-sample tolerance.
        let wrap = primary_code_length * secondary_code_length,
            d = mod(synced_code_phase, wrap)

            @test min(d, wrap - d) < 1e-3
        end
    end

    # Sub-primary-block start phase: begin tracking half a primary-code period
    # into a secondary chip (`code_phase = (k + 0.5) x primary_code_length`).
    # The secondary-code phase snap runs once, at the sync transition, and
    # anchors `code_phase` to the right NH10 chip while *preserving* the
    # within-primary-block phase (issue #117): erasing it on every call would
    # discard a chunk-bounded partial integration and could wedge the
    # satellite into a state where no chunk ever completes a block. Here a
    # 60-block buffer fed in one `track` call leaves the loop mid-block with a
    # half-primary-period partial in flight, so the final `code_phase` lands
    # at secondary chip `k mod 10` plus that leftover half-block phase.
    # Confirms the rotation search + phase snap handle a non-block-aligned
    # start without dropping in-flight integration progress.
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
        # Post-sync the replica bakes the NH chip per block (issue #125), so the
        # loop sees a sign-consistent prompt and code Doppler nudges by floating-
        # point noise over the continuous 60-block run; the converged code phase
        # lands on chip `k mod 10` + the leftover half block up to a sub-microchip
        # residual (~5e-7 of 10230 chips). Same tolerance as the L1C-P run below.
        @test get_code_phase(track_state) ≈
              (k % secondary_code_length) * primary_code_length + 0.5 * primary_code_length atol =
            1e-4
    end
end

# GPS L5I post-sync data-bit decoding (issue #125).
#
# The soft NH10 CFAR search fires on an NH10 = data-bit boundary, so the
# post-sync bit decoder must (a) wipe the NH code from each prompt — at the
# default 1-block integration the replica covers a single primary period and has
# to bake the correct NH chip, otherwise a 10-block "bit" sum collapses to
# `data x Σ(NH signs)` and loses ~14 dB of decision margin — and (b) start its
# accumulation grid on that boundary, so bits are emitted on the NH10 / data-bit
# boundary. This decodes a known bit sequence at several starting chips and at
# both the default 1-block and a 10-block coherent integration.
#
# Note on polarity: locking onto a periodic secondary code carries an inherent
# ±1 data-polarity ambiguity (the rotation search cannot tell a data "1" period
# from a "0" period — both match the NH pattern in one orientation), resolved
# downstream by the nav-message preamble. The fix's job is a *clean, aligned,
# uncancelled* decode, so the sequence is checked up to a single global
# inversion, with the per-bit coherent magnitude asserted near the full NH10
# length (≈10, vs ≈2 if the NH code were left in).
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
        # One continuous signal: `gen_code` bakes the NH chip per block via the
        # advancing code phase, and each 1 ms block is scaled by its data bit.
        # Fed in a single `track` call so the inner loop owns all block-boundary
        # bookkeeping — robust to the sub-sample code-Doppler the loop develops,
        # which would otherwise straddle fixed per-call chunks and drop bits.
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

        # The soft CFAR detector is deterministic on this noiseless signal, so the
        # exact first decoded bit is derivable (no "matches at some offset"
        # slack). `_cfar_decide` needs ≥ 2 NH10 periods before a runner-up exists,
        # then `_detect_secondary_code_cfar` fires at the winning rotation's own
        # period boundary. The winning rotation is `d* = mod(N − start, N)` (the
        # one whose overlay wipe is coherent), so the lock lands at the smallest
        # block count `≥ 2N` with `count % N == d*`; the upcoming integration then
        # starts at absolute chip 0, i.e. data-bit `div(start + lock, N)`.
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

        # Every bit from the lock boundary to the end decodes, except at most the
        # final one (its integration can still be in flight when `track` returns).
        @test n_decoded in (length(expected_from_lock) - 1, length(expected_from_lock))
        # Exact, boundary-aligned decode at the derived offset — up to the single
        # global ±1 polarity ambiguity inherent to secondary-code lock.
        expected = expected_from_lock[1:n_decoded]
        @test decoded_bits == expected || decoded_bits == 1 .- expected
        # Soft-bit signs are internally consistent with the hard bits.
        @test all((soft_bits .> 0) .== (decoded_bits .== 1))
        # The NH code is wiped: every post-sync bit sums coherently. A full bit
        # spans `secondary_code_length` primary blocks; the soft bit sums one
        # unit-magnitude prompt per `preferred_blocks`-block integration, so its
        # magnitude is ≈ `secondary_code_length / preferred_blocks`. Because the
        # detector fires exactly on the data-bit boundary there is no truncated
        # first bit — *every* emitted bit is full. Left in, the NH signs would
        # nearly cancel (≈2 for NH10 at 1-block), so this near-nominal floor is
        # the load-bearing assertion of the fix.
        full_bit_magnitude = secondary_code_length / preferred_blocks
        @test all(>(0.9 * full_bit_magnitude), abs.(soft_bits))
    end
end

# Mid-fold secondary-code sync: the overlay anchor must travel with the fold.
#
# Same defect as the L1 C/A bit grid above (issue #219), on the secondary-code
# path. With a `doppler_update_interval` longer than one code period the records
# behind the syncing one inside its fold were correlated with the pre-sync (not
# overlay-wiped) replica: their prompt is dropped, but they still advance both
# the bit accumulator's block count and — because the code-phase snap runs after
# the fold — the secondary-chip anchor it aligns the upcoming integration to.
# Leaving the anchor where the detector reported it makes the post-sync replica
# bake the wrong NH chip into every block, and the coherent NH10 sum collapses
# from ~10 to ~0–2 (which is also a coin flip on every decoded bit).
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
        # One record per fold: no record can trail the syncing one — the
        # reference decode. (Magnitudes sit a few 1e-4 below the nominal NH10
        # length, from the code-amplitude normalisation, hence the 0.99 floors.)
        full_bit = 0.99 * secondary_code_length
        reference = copy(decode(1e-3s, start_secondary_chip))
        @test all(>(full_bit), abs.(reference))

        for doppler_update_interval in (3e-3s, 7e-3s)
            soft_bits = decode(doppler_update_interval, start_secondary_chip)
            # Same bits, on the same grid, as the one-record-per-fold reference.
            @test length(soft_bits) == length(reference)
            @test sign.(soft_bits) == sign.(reference)
            # Only the bit the sync lands in is short, and only by the prompts
            # that were dropped — one block each, one such record at both
            # intervals here. Every later bit is full. A misaligned overlay
            # would put all of them near 0–2 instead.
            @test all(>(full_bit), abs.(soft_bits[2:end]))
            @test abs(soft_bits[1]) > 0.99 * (secondary_code_length - 1)
        end
    end
end

# Mid-fold secondary-code sync on a *pilot* (GPS L5Q, NH20).
#
# A pilot carries no navigation data, so the accumulator part of issue #219 does
# not apply — `buffer` returns early for it. The overlay anchor does, and it hits
# harder: the overlay is the only thing a pilot locks, and its whole purpose is
# coherent integration over a full secondary-code period. Anchored to the chip
# the detector reported for the syncing record — instead of the one the records
# behind it moved on to — the replica wipes the wrong NH chip and a 20-block
# coherent sum cancels toward zero instead of reaching the nominal 1. That is not
# a degraded prompt: the discriminators then see 0/0 and the loop diverges into a
# NaN Doppler (an `InexactError` out of the correlator's shift arithmetic), so
# this testset *errors* rather than fails when the anchor is wrong.
@testset "mid-fold pilot overlay anchor keeps coherent integration intact" begin
    gpsl5q = GPSL5Q()
    prn = 1
    sampling_frequency = 25e6Hz
    code_frequency = get_code_frequency(gpsl5q)
    secondary_code_length = get_secondary_code_length(gpsl5q)   # 20 (NH20)
    num_samples = round(Int, 25e6 / 1000)                       # 1 ms = one primary block
    # Long enough for the NH20 detector (≥ 2 periods) plus a good run of
    # full-period coherent integrations afterwards.
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

    # 1 ms: one record per fold, nothing can trail the syncing record. 3 / 7 ms:
    # it can, and does.
    for doppler_update_interval in (1e-3s, 3e-3s, 7e-3s)
        track_state = TrackState(gpsl5q, [TrackedSat(gpsl5q, prn, 0.0, 0.0Hz)])
        # Coherently integrate one full NH20 period — the point of a pilot lock.
        set_preferred_num_code_blocks_to_integrate!(track_state, prn, secondary_code_length)
        track_state =
            track(signal, track_state, sampling_frequency; doppler_update_interval)

        @test has_bit_or_secondary_code_been_found(track_state)
        # The long integration is really in effect: until sync every record is
        # clamped to a single block, afterwards each spans a full NH20 period,
        # so this noiseless run produces ~40 + (240 − 40)/20 = 50 records — a
        # run stuck at one block per record would produce ~240. (v4 has no
        # `get_last_fully_integrated_num_code_blocks` to read the count
        # directly.)
        prompts = get_filtered_prompts(track_state, prn)
        @test length(prompts) < 60
        # The overlay is wiped across all 20 blocks of it, so the
        # sample-count-normalised prompt keeps its full magnitude. A one-chip
        # anchor error would put this near 1/20 of that.
        @test abs(get_last_fully_integrated_filtered_prompt(track_state, prn)) > 0.99
        # Not just the final record: every full-period integration after sync
        # holds up. Sync needs ≥ 2 NH20 periods, so the tail of the prompt
        # record covers post-sync integrations only.
        @test all(>(0.99), abs.(prompts[(end-4):end]))
    end
end

# GPS L1C-P 1800-chip overlay sync — full 18 s end-to-end run.
#
# L1C-P is a pilot (no navigation data); its only "bit"-like event is
# locking the per-PRN 1800-chip overlay, which takes one full 18 s overlay
# cycle. This drives a clean, perfectly-aligned, zero-Doppler signal for
# slightly more than 1800 primary-code periods and confirms the overlay
# rotation search locks at the 1800th block with the code phase anchored to
# overlay chip 0 (the upcoming chip after a full cycle started at chip 0).
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
    # Upcoming overlay chip after a full cycle from chip 0 is chip 0 again.
    # The phase snap preserves the within-primary-block phase (issue #117)
    # rather than zeroing it, so the block-aligned start lands at chip 0 up
    # to the floating-point residual accumulated over 1800 code-phase
    # updates (~1e-7 chips out of 10230).
    @test synced_code_phase ≈ 0.0 atol = 1e-4
end

end
