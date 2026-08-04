module CN0EstimationTest

using Test: @test, @testset, @inferred, @test_throws
using Random: Random, Xoshiro
using Unitful: kHz, MHz, Hz, ms, dBHz, ustrip, uconvert
import Unitful
using StaticArrays: SVector
using GNSSSignals: GPSL1CA, GPSL1C_D, GPSL1C_P, GPSL5I, get_code, gen_code
using GNSSSignals: get_code_frequency
import Tracking
using Tracking:
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    CN0UpdateContext,
    BitBuffer,
    default_cn0_estimator,
    get_fallback_cn0_estimator,
    get_prompt_buffer,
    get_current_index,
    update,
    EarlyPromptLateCorrelator,
    get_prompt,
    get_correlator_sample_shifts,
    estimate_cn0,
    TrackedSignal,
    TrackedSat,
    TrackState,
    add_satellite!,
    get_cn0_estimator,
    get_last_fully_integrated_integration_time,
    has_bit_or_secondary_code_been_found,
    track

@testset "Moments CN0 estimator" begin
    cn0_estimator = MomentsCN0Estimator(20)
    @test @inferred(get_prompt_buffer(cn0_estimator)) == zero(SVector{20,ComplexF64})
    @test @inferred(get_current_index(cn0_estimator)) == 0
    @test @inferred(Base.length(cn0_estimator)) == 0

    next_cn0_estimator = @inferred update(cn0_estimator, 1 + 2im)
    @test @inferred(get_prompt_buffer(next_cn0_estimator))[1] == 1 + 2im
    @test @inferred(get_current_index(next_cn0_estimator)) == 1
    @test @inferred(Base.length(next_cn0_estimator)) == 1

    cn0_estimator = MomentsCN0Estimator(ones(SVector{20,ComplexF64}), 20, 20)
    next_cn0_estimator = @inferred update(cn0_estimator, 1 + 2im)
    @test @inferred(get_prompt_buffer(next_cn0_estimator))[1] == 1 + 2im
    @test @inferred(get_current_index(next_cn0_estimator)) == 1
    @test @inferred(Base.length(next_cn0_estimator)) == 20

    cn0_estimator = MomentsCN0Estimator(ones(SVector{20,ComplexF64}), 19, 20)
    @test @inferred(get_current_index(cn0_estimator)) == 19
    @test @inferred(Base.length(cn0_estimator)) == 20
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

    # Not synced yet (bit grid unknown): a short unaligned window is used, since
    # the bit-edge detector needs seconds to lock at low C/N₀ and the moment
    # ratio's noise floor is useless there. Two blocks by default.
    estimator = fresh()
    for _ = 1:8
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, -1))
    end
    @test Base.length(estimator) == 4
    @test estimator.num_records_per_ratio == 2
    # Switching to the post-sync window restarts the ring: `M` enters the
    # estimate, so ratios formed at the old record count cannot be averaged in.
    for bit_block_index = 0:19
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, bit_block_index))
    end
    @test Base.length(estimator) == 1
    @test estimator.num_records_per_ratio == 20

    # A pre-sync window that is still open when sync arrives is dropped, not
    # carried into the bit-aligned window: it was opened at the free pre-sync
    # length and now sits at an unknown offset inside a bit, so continuing it
    # would sum across a data-bit transition. Nine unaligned records leave one
    # in an open two-block window; sync then reports offset 5.
    estimator = fresh()
    for _ = 1:9
        estimator = update(estimator, 1.0 + 0.0im, context(GPSL1CA(), 20, -1))
    end
    @test estimator.num_accumulated_records == 1
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

@testset "CN0 estimation" begin
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
    num_samples = 4000
    gpsl1 = GPSL1CA()

    for i = 1:100
        signal =
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ start_code_phase,
                prn,
            ) .* 10^(45 / 20) .+ randn(ComplexF64, length(range)) .* sqrt(4e6)
        code = get_code.(
            gpsl1,
            code_frequency .* (-2:4001) ./ sampling_frequency .+ start_code_phase,
            prn,
        )
        correlator = EarlyPromptLateCorrelator()
        sample_shifts =
            get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
        correlator = Tracking.downconvert_and_correlate_fused!(
            correlator,
            signal,
            code,
            sample_shifts,
            0.0Hz,
            sampling_frequency,
            0.0,
            start_sample,
            num_samples,
        )
        cn0_estimator = update(cn0_estimator, get_prompt(correlator))
    end
    @test @inferred(get_current_index(cn0_estimator)) == 100
    @test @inferred(Base.length(cn0_estimator)) == 100
    cn0_estimate = @inferred estimate_cn0(cn0_estimator, 1ms)

    @test cn0_estimate ≈ 45dBHz atol = 1.0dBHz
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
        [TrackedSat(gpsl1, prn, start_code_phase, carrier_doppler)];
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
    @test cn0_estimate ≈ 45dBHz atol = 1.0dBHz
end

# A minimal custom estimator: it only proves that the type is stored and its
# methods are the ones called, so it does no estimating at all.
struct CountingCN0Estimator <: Tracking.AbstractCN0Estimator
    num_prompts::Int
end
Tracking.update(estimator::CountingCN0Estimator, prompt) =
    CountingCN0Estimator(estimator.num_prompts + 1)
Tracking.estimate_cn0(estimator::CountingCN0Estimator, integration_time) =
    estimator.num_prompts * dBHz

@testset "custom CN0 estimator is pluggable (issue #217)" begin
    gpsl1 = GPSL1CA()
    # The estimator is a type parameter of `TrackedSignal`, so a custom one is
    # stored as is instead of being converted to the default estimator's type.
    tracked_signal = @inferred TrackedSignal(gpsl1; cn0_estimator = CountingCN0Estimator(0))
    @test get_cn0_estimator(tracked_signal) isa CountingCN0Estimator
    @test typeof(tracked_signal).parameters[5] === CountingCN0Estimator

    # ... including when swapped onto an existing signal via the kwarg-update
    # constructor, which is why that kwarg is not pinned to the old type.
    swapped = TrackedSignal(TrackedSignal(gpsl1); cn0_estimator = CountingCN0Estimator(7))
    @test get_cn0_estimator(swapped) isa CountingCN0Estimator
    @test estimate_cn0(swapped) == 7dBHz

    # ... and through both `TrackedSat` constructors into a `TrackState`.
    sat = TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = CountingCN0Estimator(3))
    @test estimate_cn0(sat) == 3dBHz
    track_state = TrackState(gpsl1, [sat])
    @test estimate_cn0(track_state, 1) == 3dBHz

    multi = TrackedSat(
        (GPSL1C_P(), GPSL1CA()),
        1,
        0.0,
        0.0Hz;
        cn0_estimator = (CountingCN0Estimator(1), MomentsCN0Estimator(10)),
    )
    @test estimate_cn0(multi, GPSL1C_P) == 1dBHz
    @test get_cn0_estimator(multi, GPSL1CA) isa MomentsCN0Estimator
    # Estimators buffer into a vector, so one instance may not be shared by two
    # signals of the same satellite — that would corrupt both.
    @test_throws ArgumentError TrackedSat(
        (GPSL1C_P(), GPSL1CA()),
        1,
        0.0,
        0.0Hz;
        cn0_estimator = MomentsCN0Estimator(10),
    )
    @test_throws ArgumentError TrackedSat(
        (GPSL1C_P(), GPSL1CA()),
        1,
        0.0,
        0.0Hz;
        cn0_estimator = (MomentsCN0Estimator(10),),
    )

    # The custom estimator is actually the one `track` feeds: 5 records in.
    signal =
        complex.(float.(gen_code(4000, gpsl1, 1, 4e6Hz, get_code_frequency(gpsl1), 0.0)))
    track_state = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = CountingCN0Estimator(0))],
    )
    for _ = 1:5
        track_state = track(signal, track_state, 4e6Hz)
    end
    @test estimate_cn0(track_state, 1) == 5dBHz

    # `num_prompts_for_cn0_estimation` still sizes the default estimator (both
    # its averaging span and its fallback's buffer).
    default = get_cn0_estimator(TrackedSignal(gpsl1; num_prompts_for_cn0_estimation = 40))
    @test default isa NWPRCN0Estimator
    @test default.num_records == 40
    @test Base.length(get_prompt_buffer(get_fallback_cn0_estimator(default))) == 40
    # The free-running window covers ~20 ms of code blocks: 20 for L1 C/A's 1 ms
    # code, 2 for L1C-P's 10 ms code.
    @test default_cn0_estimator(gpsl1, 100).num_narrowband_code_blocks == 20
    @test default_cn0_estimator(GPSL1C_P(), 100).num_narrowband_code_blocks == 2
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

    function track_noisy(cn0_db, num_blocks; seed = 1, cn0_estimator = nothing)
        rng = Xoshiro(seed)
        amplitude = isnothing(cn0_db) ? 0.0 : 10^(cn0_db / 20)
        sat = if isnothing(cn0_estimator)
            TrackedSat(gpsl1, 1, 0.0, 0.0Hz)
        else
            TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator)
        end
        track_state = TrackState(gpsl1, [sat])
        data_bits = rand(rng, (-1.0, 1.0), div(num_blocks, 20) + 2)
        for block_index = 1:num_blocks
            code_phase =
                (block_index - 1) * num_samples * code_frequency / sampling_frequency
            clean =
                amplitude .* data_bits[div(block_index-1, 20)+1] .* gen_code(
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
    # Synced: one window per navigation bit, i.e. 20 one-block records each.
    @test estimator.num_records_per_ratio == 20
    @test Base.length(estimator) > 0
    @test ustrip(uconvert(dBHz, estimate_cn0(track_state, 1))) ≈ 45 atol = 2

    # Pure noise, still pre-sync (which is where a code-lock detector has to
    # make its call): the moment estimator reports a signal that is not there,
    # NWPR does not.
    noise_state = track_noisy(nothing, 200)
    noise_cn0 = ustrip(uconvert(dBHz, estimate_cn0(noise_state, 1)))
    @test !has_bit_or_secondary_code_been_found(noise_state, 1)
    @test noise_cn0 < 20
    moments_noise_state =
        track_noisy(nothing, 200; cn0_estimator = MomentsCN0Estimator(100))
    @test ustrip(uconvert(dBHz, estimate_cn0(moments_noise_state, 1))) > 25
end

@testset "estimate_cn0 overloads on TrackState" begin
    # Multi-group + sat-id and single-group + sat-id variants. The
    # no-argument variant is already covered by the integration test
    # above; these two specialize on the group key / sat identifier
    # forwarding paths. With an unseeded CN0 estimator the value is
    # 0 dB-Hz — we only assert the methods dispatch and run.
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1, carrier_doppler = 0Hz)
    @test estimate_cn0(ts, :default, 1) == 0.0dBHz
    @test estimate_cn0(ts, 1) == 0.0dBHz
end

@testset "estimate_cn0 divides by the record's real integration time" begin
    # C/N₀ has units of Hz and belongs to the signal, not to how long the
    # correlate step chose to integrate. The estimator buffers *sample-normalized*
    # prompts, so a record spanning N code blocks arrives with N times the SNR of
    # a one-block record; `estimate_cn0` therefore has to divide by N × the code
    # period, not by the code period alone. It used to do the latter, which
    # over-reported by 10·log₁₀(N) — 13 dB at a full GPS L1 C/A bit — for anything
    # driven above one block by `set_preferred_num_code_blocks_to_integrate!` or
    # by an external correlator producer handing over longer records.
    gpsl1 = GPSL1CA()
    prn = 1

    # A prompt whose amplitude is fixed and whose noise shrinks as √N is exactly
    # what a longer coherent integration delivers after sample normalization, so
    # feeding the same shape at two different block counts must report the same
    # C/N₀ once the divisor is right.
    function cn0_at(num_blocks, prompts)
        tsig = Tracking.TrackedSignal(gpsl1)
        estimator = get_cn0_estimator(tsig)
        for p in prompts
            estimator = update(estimator, p)
        end
        tsig = Tracking.TrackedSignal(
            tsig;
            cn0_estimator = estimator,
            last_fully_integrated_num_code_blocks = num_blocks,
        )
        estimate_cn0(tsig)
    end

    Random.seed!(4321)
    amplitude = 1.0
    noise_1 = 0.2
    one_block = [amplitude + noise_1 * randn(ComplexF64) for _ = 1:100]
    # 20 blocks: same signal amplitude, noise down by √20.
    twenty_blocks = [amplitude + noise_1 / sqrt(20) * randn(ComplexF64) for _ = 1:100]

    cn0_1 = cn0_at(1, one_block)
    cn0_20 = cn0_at(20, twenty_blocks)
    # The 20-block record's SNR is 20× (13 dB) better and it spans 20× the time,
    # so the C/N₀ must come out the same.
    @test cn0_20 ≈ cn0_1 atol = 1.5dBHz

    # And the divisor is actually used: the same buffered prompts reported at a
    # different block count must move by exactly 10·log₁₀(N). Compared on the
    # linear values, since subtracting two logarithmic `Level`s is awkward.
    linear_cn0(x) = ustrip(Hz, Unitful.linear(x))
    same_prompts_1 = linear_cn0(cn0_at(1, one_block))
    same_prompts_20 = linear_cn0(cn0_at(20, one_block))
    @test 10 * log10(same_prompts_1 / same_prompts_20) ≈ 10 * log10(20) atol = 0.01

    # `estimate_cn0` and `get_last_fully_integrated_integration_time` are the two
    # halves of one contract: C/N₀ is per-Hz and says nothing on its own about
    # whether a record's peak clears the noise. Their product does — it is the
    # post-integration SNR — so a consumer gating on detectability multiplies
    # them. Same buffered prompts at 1 and 20 blocks: C/N₀ moves by 10·log₁₀(20)
    # and T moves by 20, so the SNR is unchanged.
    function snr_at(num_blocks, prompts)
        tsig = Tracking.TrackedSignal(gpsl1)
        estimator = get_cn0_estimator(tsig)
        for p in prompts
            estimator = update(estimator, p)
        end
        tsig = Tracking.TrackedSignal(
            tsig;
            cn0_estimator = estimator,
            last_fully_integrated_num_code_blocks = num_blocks,
        )
        cn0 = ustrip(Hz, Unitful.linear(estimate_cn0(tsig)))
        cn0 * ustrip(Unitful.s, get_last_fully_integrated_integration_time(tsig))
    end
    @test snr_at(20, one_block) ≈ snr_at(1, one_block) rtol = 1e-9
end

end
