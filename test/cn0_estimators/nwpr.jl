module NWPRCN0EstimatorTest

using Test: @test, @testset, @inferred, @test_throws
using Random: Random, Xoshiro
using Unitful: kHz, MHz, Hz, ms, dBHz, ustrip, uconvert
import Unitful
using GNSSSignals: GPSL1CA, GPSL1C_D, GPSL1C_P, GPSL5I, get_code, gen_code
using GNSSSignals: get_code_frequency
import Tracking
using Tracking:
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    CN0UpdateContext,
    BitBuffer,
    get_fallback_cn0_estimator,
    update,
    estimate_cn0,
    TrackedSat,
    TrackState,
    get_cn0_estimator,
    has_bit_or_secondary_code_been_found,
    track

@testset "NWPR's coherent window is sized from the signal's code period" begin
    # The window caps the coherent sum in *blocks*, so sizing it without knowing
    # the code period gets it wrong by the period's ratio: ~5 ms is 5 blocks of
    # L1 C/A's 1 ms code and 2 of L1C-P's 10 ms one. This is the form to reach
    # for on a correlator-ingest path, which is the one place NWPR is still the
    # estimator to choose over the default.
    @test NWPRCN0Estimator(GPSL1CA()).num_narrowband_code_blocks == 5
    @test NWPRCN0Estimator(GPSL1C_P()).num_narrowband_code_blocks == 2
    @test NWPRCN0Estimator(GPSL1CA(); num_records = 40).num_records == 40
    # An explicit window still wins.
    @test NWPRCN0Estimator(GPSL1C_P(); num_narrowband_code_blocks = 7).num_narrowband_code_blocks ==
          7
end

@testset "NWPR CN0 estimator on a bare prompt stream" begin
    @test_throws ArgumentError NWPRCN0Estimator(; num_records = 1)
    @test_throws ArgumentError NWPRCN0Estimator(; num_narrowband_code_blocks = 0)
    @test_throws ArgumentError NWPRCN0Estimator(; num_presync_narrowband_code_blocks = -1)

    # Bare stream: one block per record, back-to-back windows of
    # `num_narrowband_code_blocks` records. Nothing is reported before the first
    # window closes — the fallback's `0.0 dB-Hz` empty-buffer value shows through.
    estimator = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 20)
    @test @inferred(Base.length(estimator)) == 0
    @test @inferred(estimate_cn0(estimator, 1ms)) == 0.0dBHz
    for _ = 1:19
        estimator = @inferred update(estimator, 1.0 + 0.0im)
    end
    @test Base.length(estimator) == 0            # window still open
    estimator = update(estimator, 1.0 + 0.0im)   # 20th record closes it
    @test Base.length(estimator) == 1
    @test estimator.num_records_per_ratio == 20
    # 20 identical prompts: the coherent sum holds all the power, µ̂ = M, which
    # is the noise-free limit of the expression.
    @test @inferred(estimate_cn0(estimator, 1ms)) == Inf * dBHz
    # The ring holds `num_records ÷ M` ratios, so the memory stays ~100 records.
    for _ = 1:(20*10)
        estimator = update(estimator, 1.0 + 0.0im)
    end
    @test Base.length(estimator) == 5

    # A window of one record carries no information (NBP == WBP by construction),
    # so nothing is ever buffered and the fallback keeps reporting.
    single = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 1)
    for _ = 1:50
        single = update(single, 1.0 + 0.0im)
    end
    @test Base.length(single) == 0
    @test estimate_cn0(single, 1ms) == estimate_cn0(get_fallback_cn0_estimator(single), 1ms)
end

@testset "NWPR CN0 estimator noise floor beats the moment ratio (issue #217)" begin
    # The whole point of NWPR being the default: the moment ratio manufactures
    # signal power out of noise at a finite window, so it reports ~27.6 dB-Hz on
    # PURE NOISE at its default 100-prompt window and cannot separate a true
    # 20 dB-Hz signal from noise at all. Prompt model of the tests above and of
    # the docs: amplitude √(C/N₀·T) in unit-total-variance complex noise.
    T = 1ms
    db(x) = ustrip(uconvert(dBHz, x))
    fold(estimator, prompts) = foldl(update, prompts; init = estimator)
    prompts(cn0_db, n, seed) =
        (isnothing(cn0_db) ? 0.0 : sqrt(10^(cn0_db / 10) * ustrip(Unitful.s, T))) .+
        randn(Xoshiro(seed), ComplexF64, n)
    moments(cn0_db, seed) =
        db(estimate_cn0(fold(MomentsCN0Estimator(100), prompts(cn0_db, 100, seed)), T))
    nwpr(cn0_db, seed) = db(
        estimate_cn0(
            fold(
                NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 20),
                prompts(cn0_db, 100, seed),
            ),
            T,
        ),
    )
    median(x) = sort(collect(x))[div(length(x) + 1, 2)]
    seeds = 1:200

    # Pure noise: the moment ratio's median sits in the middle of the useful
    # range; NWPR's stays far below any threshold a receiver would use.
    @test median(moments(nothing, s) for s in seeds) > 25
    @test median(nwpr(nothing, s) for s in seeds) < 15
    # False alarms against a 25 dB-Hz code-lock threshold on pure noise.
    @test count(s -> moments(nothing, s) >= 25, seeds) / length(seeds) > 0.5
    @test count(s -> nwpr(nothing, s) >= 25, seeds) == 0

    # And it is an estimate, not just a detector: unbiased to ~1 dB from
    # 20 dB-Hz up, where the moment ratio is pinned to its noise floor.
    for cn0 in (20.0, 25.0, 30.0, 45.0)
        @test median(nwpr(cn0, s) for s in seeds) ≈ cn0 atol = 1.0
    end
    @test median(moments(20.0, s) for s in seeds) - 20 > 5
end

@testset "NWPR CN0 estimator follows the navigation-bit grid" begin
    # The context-driven form is what `track` calls. Feed a constant prompt and
    # walk the bit grid by hand, so only the window logic is under test.
    bit_buffer = BitBuffer{UInt64}()
    context(signal, blocks_per_bit, bit_block_index) =
        CN0UpdateContext(signal, 1, blocks_per_bit, bit_block_index, bit_buffer)
    fresh() = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 20)

    # GPS L1 C/A, synced: one window per 20-block navigation bit, and it only
    # starts on a bit boundary — the 10 blocks before the first boundary are
    # dropped rather than summed into a straddling window.
    estimator = fresh()
    for bit_block_index in [10:19; repeat(0:19, 2)]
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test Base.length(estimator) == 2
    @test estimator.num_records_per_ratio == 20

    # ... but the window is capped by `num_narrowband_code_blocks`, the loop's
    # coherence budget: shorter windows tile the bit from its start, so none
    # straddles a flip either, and four of them fit in a 20-block bit.
    capped = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 5)
    for bit_block_index in repeat(0:19, 2)
        capped = update(capped, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test Base.length(capped) == 8
    @test capped.num_records_per_ratio == 5
    # A window that would run past the end of the bit is not opened at all: with
    # a 6-block window three fit, and the bit's last two blocks are left out.
    ragged = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 6)
    for bit_block_index in repeat(0:19, 2)
        ragged = update(ragged, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test Base.length(ragged) == 6
    @test ragged.num_records_per_ratio == 6
    # A cap above the bit period cannot lengthen the window past the bit.
    wide = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 50)
    for bit_block_index in repeat(0:19, 2)
        wide = update(wide, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test Base.length(wide) == 2
    @test wide.num_records_per_ratio == 20

    # Not synced yet (bit grid unknown): a short unaligned window is used, since
    # the bit-edge detector needs seconds to lock at low C/N₀ and the moment
    # ratio's noise floor is useless there. Five blocks by default.
    estimator = fresh()
    for _ = 1:20
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, -1))
    end
    @test Base.length(estimator) == 4
    @test estimator.num_records_per_ratio == 5
    # An explicit pre-sync length is honoured.
    two_block = NWPRCN0Estimator(;
        num_records = 100,
        num_narrowband_code_blocks = 20,
        num_presync_narrowband_code_blocks = 2,
    )
    for _ = 1:8
        two_block = update(two_block, 1.0 + 0.0im, context(GPSL1CA(), 20, -1))
    end
    @test Base.length(two_block) == 4
    @test two_block.num_records_per_ratio == 2
    # Switching to the post-sync window restarts the ring: `M` enters the
    # estimate, so windows formed at the old record count cannot be combined in.
    for bit_block_index = 0:19
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test Base.length(estimator) == 1
    @test estimator.num_records_per_ratio == 20

    # Even at an unchanged window length, the buffered pre-sync windows are
    # dropped when the first bit-aligned one completes: they ran unaligned, so
    # some of them straddled a flip, and averaging them with clean windows would
    # drag the estimate down for a whole `num_records` after sync.
    same_length = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 5)
    for _ = 1:20
        same_length = update(same_length, 1.0 + 0.0im, context(GPSL1CA(), 20, -1))
    end
    @test Base.length(same_length) == 4
    @test same_length.num_records_per_ratio == 5
    @test !same_length.ratios_are_bit_aligned
    for bit_block_index = 0:4
        same_length =
            update(same_length, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test same_length.num_records_per_ratio == 5     # window length unchanged ...
    @test same_length.ratios_are_bit_aligned
    @test Base.length(same_length) == 1              # ... yet the ring restarted

    # A pre-sync window that is still open when sync arrives is dropped, not
    # carried into the bit-aligned window: it was opened at the free pre-sync
    # length and now sits at an unknown offset inside a bit, so continuing it
    # would sum across a data-bit transition. Nine unaligned records leave four
    # in an open five-block window; sync then reports offset 5.
    estimator = fresh()
    for _ = 1:9
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, -1))
    end
    @test estimator.num_accumulated_records == 4
    estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, 5))
    @test estimator.num_accumulated_records == 0     # dropped, not continued
    # ... and the next bit boundary opens the aligned window.
    for bit_block_index in [6:19; 0:19]
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test estimator.num_records_per_ratio == 20
    @test Base.length(estimator) == 1

    # A signal carrying a secondary code has no pre-sync coherence at all (the
    # unknown overlay flips sign per code block), so no window is opened and the
    # fallback estimator is what gets reported.
    estimator = fresh()
    for _ = 1:40
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL5I(), 10, -1))
    end
    @test Base.length(estimator) == 0
    @test estimate_cn0(estimator, 1ms) ==
          estimate_cn0(get_fallback_cn0_estimator(estimator), 1ms)
    # Once synced its bit is 10 blocks long, and the replica wipes the overlay.
    for bit_block_index in repeat(0:9, 3)
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL5I(), 10, bit_block_index))
    end
    @test Base.length(estimator) == 3
    @test estimator.num_records_per_ratio == 10

    # One symbol per code block (GPS L1C-D, Galileo E1B): every record may flip
    # sign, so no coherent window exists, synced or not.
    estimator = fresh()
    for _ = 1:40
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1C_D(), 1, 0))
    end
    @test Base.length(estimator) == 0

    # A pilot (no data bits at all) has no bit grid to respect post-sync, so
    # windows run back to back at the configured length.
    estimator = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 5)
    for _ = 1:20
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1C_P(), 0, 0))
    end
    @test Base.length(estimator) == 4
    @test estimator.num_records_per_ratio == 5

    # `CN0UpdateContext` derives the grid from the signal and the bit buffer, and
    # marks it unknown while sync has not been found.
    context_l1ca = @inferred CN0UpdateContext(GPSL1CA(), bit_buffer, 1)
    @test context_l1ca.num_code_blocks_per_bit == 20
    @test context_l1ca.bit_code_block_index == -1
    synced_buffer = BitBuffer{UInt64}(
        zero(UInt64),
        40,
        true,
        0,
        Int8(1),
        complex(0.0, 0.0),
        3,
        Float32[],
        Tracking.PhaseAccumulators(),
    )
    @test CN0UpdateContext(GPSL1CA(), synced_buffer, 1).bit_code_block_index == 3
    # Records of a fold that follow a sync detected in that same fold were
    # correlated with pre-sync replicas — their alignment is not trustworthy.
    @test CN0UpdateContext(GPSL1CA(), synced_buffer, 1, false).bit_code_block_index == -1
end

@testset "NWPR CN0 estimator when a record outgrows its own window" begin
    bit_buffer = BitBuffer{UInt64}()
    context(num_code_blocks, bit_block_index) =
        CN0UpdateContext(GPSL1CA(), num_code_blocks, 20, bit_block_index, bit_buffer)

    # `set_preferred_num_code_blocks_to_integrate!` at a whole navigation bit
    # closes every window on a single record, where `NBP == WBP` by construction
    # and the window carries no information. The windows buffered before the
    # switch have to go with it: keeping them would freeze the estimate on
    # windows formed at the old record length for good — `estimate_cn0` consults
    # the `fallback` only while the ring is empty — and report them against the
    # new record's integration time on top of that.
    rng = Xoshiro(1)
    # A true ~35 dB-Hz at T = 1 ms, so the ring holds a finite value that is
    # distinguishable from the fallback's.
    noisy() = sqrt(10^3.5 * 1e-3) + randn(rng, ComplexF64)
    estimator = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 5)
    for _ = 1:5, bit_block_index = 0:19
        estimator = update(estimator, noisy(), context(1, bit_block_index))
    end
    @test Base.length(estimator) == 20
    @test estimator.num_records_per_ratio == 5
    for _ = 1:20
        estimator = update(estimator, noisy(), context(20, 0))
    end
    @test Base.length(estimator) == 0
    @test estimate_cn0(estimator, 20ms) ==
          estimate_cn0(get_fallback_cn0_estimator(estimator), 20ms)
end

@testset "NWPR CN0 estimator skips a record spanning no whole code block" begin
    # The fractional record right after a sync phase-snap resets the accumulator
    # covers no whole code block, so it has no position on the bit grid: it
    # leaves the open window untouched rather than closing it one record over the
    # block count the inversion assumes. Its prompt still reaches the fallback.
    bit_buffer = BitBuffer{UInt64}()
    context(num_code_blocks, bit_block_index) =
        CN0UpdateContext(GPSL1CA(), num_code_blocks, 20, bit_block_index, bit_buffer)

    estimator = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 5)
    for bit_block_index = 0:2
        estimator = update(estimator, 1.0 + 0.0im, context(1, bit_block_index))
    end
    @test estimator.num_accumulated_records == 3
    estimator = update(estimator, 1.0 + 0.0im, context(0, 3))
    @test estimator.num_accumulated_records == 3
    @test estimator.num_accumulated_code_blocks == 3
    for bit_block_index = 3:4
        estimator = update(estimator, 1.0 + 0.0im, context(1, bit_block_index))
    end
    @test Base.length(estimator) == 1
    @test estimator.num_records_per_ratio == 5
    @test Base.length(get_fallback_cn0_estimator(estimator)) == 6
end

@testset "CN0 estimation integration test" begin
    Random.seed!(1234)
    carrier_doppler = 0Hz
    start_code_phase = 0
    code_frequency = 1023kHz
    sampling_frequency = 4MHz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    cn0_estimator = MomentsCN0Estimator(100)
    start_sample = 1
    gpsl1 = GPSL1CA()

    track_state = @inferred TrackState(
        gpsl1,
        [
            TrackedSat(
                gpsl1,
                prn,
                start_code_phase,
                carrier_doppler;
                cn0_estimator = NWPRCN0Estimator(gpsl1),
            ),
        ],
    )

    for i = 1:100
        signal =
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ start_code_phase,
                prn,
            ) .* 10^(45 / 20) .+ randn(ComplexF64, length(range)) .* sqrt(4e6)
        start_code_phase =
            code_frequency * length(range) ./ sampling_frequency + start_code_phase
        track_state = @inferred track(signal, track_state, sampling_frequency)
    end
    cn0_estimate = @inferred estimate_cn0(track_state)
    # 100 records is 20 windows at the default five-block pre-sync window, and a
    # 45 dB-Hz signal sits near NWPR's upper limit (`µ̂ → M`), where the dB
    # resolution per unit `µ̂` is poor — so a single realization spreads more than
    # the moment ratio did here. Measured over 300 seeds: median 44.4, extremes
    # 42.9 and 46.1, worst deviation 2.1 dB.
    @test cn0_estimate ≈ 45dBHz atol = 2.5dBHz
end

@testset "NWPR CN0 estimator in the loop (issue #217)" begin
    # End to end through `track`, on a data-modulated GPS L1 C/A signal: the
    # aligned post-sync window has to actually engage, and the estimate has to
    # land on the true C/N₀ — while pure noise must NOT read like a signal, which
    # is what the moment estimator does at any window length.
    gpsl1 = GPSL1CA()
    sampling_frequency = 4e6Hz
    num_samples = 4000
    code_frequency = get_code_frequency(gpsl1)

    # `fade_to` continues the same satellite at another C/N₀ for `fade_blocks`
    # more records, which is how the post-sync regime is reached at a C/N₀ where
    # the bit-edge detector would never have locked from cold.
    # NWPR is no longer the library default, so it is configured explicitly here
    # — sized for the signal, which is what `NWPRCN0Estimator(signal)` is for.
    function track_noisy(
        cn0_db,
        num_blocks;
        seed = 1,
        cn0_estimator = NWPRCN0Estimator(gpsl1),
        fade_to = nothing,
        fade_blocks = 0,
    )
        rng = Xoshiro(seed)
        amplitude(cn0) = isnothing(cn0) ? 0.0 : 10^(cn0 / 20)
        sat = TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator)
        track_state = TrackState(gpsl1, [sat])
        total_blocks = num_blocks + fade_blocks
        data_bits = rand(rng, (-1.0, 1.0), div(total_blocks, 20) + 2)
        for block_index = 1:total_blocks
            code_phase =
                (block_index - 1) * num_samples * code_frequency / sampling_frequency
            clean =
                amplitude(block_index > num_blocks ? fade_to : cn0_db) .*
                data_bits[div(block_index-1, 20)+1] .* gen_code(
                    num_samples,
                    gpsl1,
                    1,
                    sampling_frequency,
                    code_frequency,
                    ustrip(code_phase),
                )
            noise =
                randn(rng, ComplexF64, num_samples) .* sqrt(ustrip(Hz, sampling_frequency))
            track_state = track(clean .+ noise, track_state, sampling_frequency)
        end
        track_state
    end

    track_state = track_noisy(45.0, 900)
    estimator = get_cn0_estimator(track_state, 1)
    @test has_bit_or_secondary_code_been_found(track_state, 1)
    # Synced: four windows tile the 20-block navigation bit at the default
    # five-block coherence cap, and they follow the bit grid.
    @test estimator.num_records_per_ratio == 5
    @test estimator.ratios_are_bit_aligned
    @test Base.length(estimator) > 0
    # A median over seeds, not one run — it is an estimate of a statistical
    # quantity. Measured over 120 seeds: median 44.4, extremes 43.4 and 45.7, and
    # the worst five-seed median 0.96 dB off, so 1.5 dB is a real bound rather than
    # a generous one. (Before JuliaGNSS/Tracking.jl#219 was fixed this had a ~2 %
    # tail down to 36.5 dB-Hz, from a navigation-bit grid that had slipped a block
    # — which made every window that followed the grid straddle a flip.)
    median_of(xs) = sort(collect(xs))[div(length(xs) + 1, 2)]
    @test median_of(
        ustrip(uconvert(dBHz, estimate_cn0(track_noisy(45.0, 900; seed), 1))) for seed = 1:5
    ) ≈ 45 atol = 1.5

    # Pure noise, still pre-sync (which is where a code-lock detector has to
    # make its call): the moment estimator reports a signal that is not there,
    # NWPR does not. The bound is 25 dB-Hz rather than 20 on purpose — a short
    # pre-sync window has a heavy upper tail on noise (a few per cent of updates
    # clear 20 dB-Hz), which is what a lock threshold has to be set against.
    noise_state = track_noisy(nothing, 200)
    noise_cn0 = ustrip(uconvert(dBHz, estimate_cn0(noise_state, 1)))
    @test !has_bit_or_secondary_code_been_found(noise_state, 1)
    @test noise_cn0 < 25
    moments_noise_state =
        track_noisy(nothing, 200; cn0_estimator = MomentsCN0Estimator(100))
    @test ustrip(uconvert(dBHz, estimate_cn0(moments_noise_state, 1))) > 25

    # A signal that faded after bit sync was found is the case the coherence cap
    # exists for: the loop stops holding phase over a whole navigation bit long
    # before the signal becomes untrackable, so a bit-long coherent sum reads far
    # too low — often `-Inf`, i.e. "no signal", on a satellite that is being
    # tracked — while the capped window stays close to the truth. Locked in at
    # 45 dB-Hz (so bit sync is found in every run) and faded to 25 dB-Hz; a
    # median over seeds, since a single run of a low-C/N₀ estimate says little.
    fade = (; fade_to = 25.0, fade_blocks = 400)
    seeds = 1:16
    faded = [track_noisy(45.0, 1500; seed, fade...) for seed in seeds]
    full_bit = [
        track_noisy(
            45.0,
            1500;
            seed,
            cn0_estimator = NWPRCN0Estimator(; num_narrowband_code_blocks = 20),
            fade...,
        ) for seed in seeds
    ]
    @test all(state -> has_bit_or_secondary_code_been_found(state, 1), faded)
    @test get_cn0_estimator(faded[1], 1).num_records_per_ratio == 5
    @test get_cn0_estimator(full_bit[1], 1).num_records_per_ratio == 20
    capped_cn0 = [ustrip(uconvert(dBHz, estimate_cn0(state, 1))) for state in faded]
    full_bit_cn0 = [ustrip(uconvert(dBHz, estimate_cn0(state, 1))) for state in full_bit]
    # Margins from 96 seeds split into disjoint groups of 16: the capped median
    # never left 22.4–24.8 and the gap to the whole-bit window never fell below
    # 1.5 dB. Over all 96, capped reads a median 23.6 dB-Hz (p10 19.5) against the
    # whole-bit window's 21.0 (p10 12.0), and the whole-bit window reports `-Inf` —
    # "no signal" on a satellite that is being tracked — in 17 of 96 runs against
    # 1. Eight seeds are not enough for either bound: on some groups of eight the
    # whole-bit window happens to win.
    @test median_of(capped_cn0) ≈ 25 atol = 3.5
    @test median_of(capped_cn0) - median_of(full_bit_cn0) > 1
end

@testset "a window with no power is discarded, not stored as a ratio" begin
    # `NBP/WBP` is undefined for a window whose prompts are all zero — a component
    # the satellite does not transmit, or a correlator handed an empty buffer. The
    # guard drops such a window rather than parking a `(0, 0)` entry in the ring:
    # a stored zero contributes nothing to the ratio of sums but occupies a slot a
    # real measurement needed, survives `num_records` windows, and shifts every
    # later ratio's position in the ring.
    bit_buffer = BitBuffer{UInt64}()
    context(bit_block_index) =
        CN0UpdateContext(GPSL1CA(), 1, 20, bit_block_index, bit_buffer)
    estimator = NWPRCN0Estimator(; num_records = 100, num_narrowband_code_blocks = 20)

    # A full 20-block bit of zero prompts closes a window with no power at all.
    for i = 0:19
        estimator = update(estimator, complex(0.0, 0.0), context(i))
    end
    @test Base.length(estimator) == 0
    @test estimator.num_records_per_ratio == 0

    # So the next live bit lands in the *first* ring slot, not the second.
    for i = 0:19
        estimator = update(estimator, complex(1.0, 0.0), context(i))
    end
    @test Base.length(estimator) == 1
    @test Tracking.get_buffered_narrowband_powers(estimator)[1] == 400.0
    @test Tracking.get_buffered_wideband_powers(estimator)[1] == 20.0
end

end
