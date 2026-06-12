module BitDetectionIntegrationTest

using Test: @test, @testset
using Unitful: Hz
using GNSSSignals: GPSL1CA, GPSL5I, GPSL1C_P, gen_code, get_code_frequency,
    get_code_length, get_secondary_code_length
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
                (bit * 2 - 1) .*
                cis.(
                    2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency .+
                    carrier_phase,
                ) .*
                gen_code(
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
            for bit_index in num_bits:-1:1
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

    for start_secondary_chip in 0:(secondary_code_length-1)
        track_state =
            TrackState(gpsl5, [TrackedSat(gpsl5, prn, 0.0, 0.0Hz)])
        synced_at_block = -1
        synced_code_phase = NaN
        for index in 1:40
            # Advance the *generated* signal's absolute code phase from an
            # initial offset of `start_secondary_chip` primary periods, so
            # block `index` carries NH10 chip `(start_secondary_chip + index - 1) % 10`.
            gen_code_phase =
                (start_secondary_chip + (index - 1)) * primary_code_length
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
        # the widened (primary x secondary) wrap window.
        @test synced_code_phase == start_secondary_chip * primary_code_length
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
              (k % secondary_code_length) * primary_code_length + 0.5 * primary_code_length atol = 1e-4
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

    @testset "preferred = $preferred_blocks, start chip $start_secondary_chip" for
            preferred_blocks in (1, 10), start_secondary_chip in (0, 3, 9)
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
        for b in 0:(total_blocks-1)
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
    for index in 1:(secondary_code_length + 5)
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
