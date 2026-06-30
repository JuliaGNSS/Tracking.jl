module BitDetectionIntegrationTest

using Test: @test, @testset
using Unitful: Hz
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
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
    has_bit_or_secondary_code_been_found,
    set_preferred_num_code_blocks_to_integrate!

@testset "Bit detection integration test" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(gpsl1)
    carrier_doppler = 0.0Hz

    # Run the same bit stream at both lock polarities: the bit-edge detector
    # locks at block 40, where the newest 20 buffered blocks carry the second
    # data bit — `1, 0, …` yields a negative-polarity lock, `0, 1, …` a
    # positive one. The bits recovered from the pre-sync buffer and the bits
    # decoded after sync must be consistent: equal to the transmitted bits up
    # to a single *global* inversion across the whole stream (issue #127).
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
            @test has_bit_or_secondary_code_been_found(track_state) == (index >= 40)
            # Sync at block 40 recovers the 2 buffered bits; afterwards one
            # bit completes every 20 blocks.
            expected_num_bits = index == 40 ? 2 : (index > 40 && index % 20 == 0 ? 1 : 0)
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
# preferred length can be set from the start. After sync at block 40 (one
# full bit of each polarity plus the transition) each integration spans 4
# code blocks, five integrations form a bit, and a hard + soft bit must be
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
        @test has_bit_or_secondary_code_been_found(multi_state) == (index >= 40)
        # Two bits are recovered from the buffered history at sync (block 40);
        # afterwards exactly one bit must appear at every 20-block boundary.
        expected_num_bits = index == 40 ? 2 : (index > 40 && index % 20 == 0 ? 1 : 0)
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
    # bit. The two bits replayed from the buffered history at sync alternate,
    # and so do the three streamed post-sync bits; the relative polarity
    # between the replayed and the streamed portion is a property of the
    # bit-sync polarity convention and not asserted here.
    @test length(multi_bits) == 5
    @test multi_bits == single_bits
    @test multi_bits[1] != multi_bits[2]
    @test multi_bits[3] != multi_bits[4]
    @test multi_bits[4] != multi_bits[5]

    # The soft bits agree with the reference in sign; post-sync each bit is
    # the sum of five 4-block prompts (magnitude ~5) instead of twenty
    # single-block prompts.
    @test sign.(multi_soft) == sign.(single_soft)
    @test all(x -> abs(x) ≈ 5, multi_soft[3:end])
    @test all(x -> abs(x) ≈ 20, single_soft[3:end])
end

# GPS L5I secondary-code (NH10) sync and phase recovery.
#
# The L5I detector runs the generic rotation search over NH10, so it locks
# after a single NH10 period (10 blocks) in the worst case — regardless of
# where in the NH10 cycle tracking started — and recovers the true
# secondary-code phase. This test starts the signal at every NH10 phase and
# checks that:
#   * sync fires exactly at the 10th block (1x the secondary period), and
#   * `code_phase` is anchored to the *upcoming* integration's NH10 chip,
#     which — after exactly 10 blocks — equals the starting chip.
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

        # Sync locks after exactly one NH10 period (10 blocks), no matter the
        # starting phase — the rotation search does not wait for a boundary.
        @test synced_at_block == secondary_code_length

        # `code_phase` is anchored to the upcoming integration's NH10 chip.
        # After 10 blocks the upcoming chip equals the starting chip, so the
        # phase lands at `start_secondary_chip x primary_code_length` within
        # the widened (primary x secondary) wrap window. The embedded-LUT
        # generator's fixed-point DDA (~2^-30 chip) lands the phase within ~1e-6
        # chip of the integer boundary rather than exactly on it, so compare the
        # circular distance within a sub-sample tolerance.
        let wrap = primary_code_length * secondary_code_length,
            expected = start_secondary_chip * primary_code_length,
            d = mod(synced_code_phase - expected, wrap)

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
# The NH10 rotation search locks at *any* secondary chip, so the post-sync bit
# decoder must (a) wipe the NH code from each prompt — at the default 1-block
# integration the replica covers a single primary period and has to bake the
# correct NH chip, otherwise a 10-block "bit" sum collapses to
# `data x Σ(NH signs)` and loses ~14 dB of decision margin — and (b) start its
# accumulation grid at the recovered secondary phase rather than at the lock
# instant, so bits are emitted on the NH10 / data-bit boundary. This decodes a
# known bit sequence through boundary and non-boundary locks at both the default
# 1-block and a 10-block coherent integration.
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

        # The pre-sync rotation search consumes the first NH10 period (the lead
        # of the lock data bit); decoding then begins with the remainder of that
        # same data bit (`data_bits[2]`) and proceeds bit-by-bit. The trailing
        # data bit may be left mid-integration, so compare the bits that
        # completed.
        soft_bits = get_soft_bits(track_state)
        decoded_bits = [Int(s > 0) for s in soft_bits]
        n_decoded = length(decoded_bits)
        expected = data_bits[2:(1+n_decoded)]

        # At least every full data bit but the last must have decoded.
        @test n_decoded >= length(data_bits) - 2
        # Clean, boundary-aligned decode — up to the global polarity ambiguity.
        @test decoded_bits == expected || decoded_bits == 1 .- expected
        # Soft-bit signs are internally consistent with the hard bits.
        @test all((soft_bits .> 0) .== (decoded_bits .== 1))
        # The NH code is wiped: every *full* post-sync bit (all but a possibly
        # truncated first one, when the lock lands mid data bit) sums coherently.
        # A full bit spans `secondary_code_length` primary blocks; the soft bit
        # sums one unit-magnitude prompt per `preferred_blocks`-block integration,
        # so its magnitude is ≈ `secondary_code_length / preferred_blocks`. Left
        # in, the NH signs would nearly cancel (≈2 for NH10 at 1-block), so this
        # half-of-nominal floor is the load-bearing assertion of the fix.
        full_bit_magnitude = secondary_code_length / preferred_blocks
        @test all(>(0.5 * full_bit_magnitude), abs.(soft_bits[2:end]))
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

    track_state = TrackState(gpsl1c_p, [TrackedSat(gpsl1c_p, prn, 0.0, 0.0Hz)])
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
