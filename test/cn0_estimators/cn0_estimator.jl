module CN0EstimatorTest

using Test: @test, @testset, @inferred, @test_throws
using Random: Random
using Unitful: Hz, ms, dBHz, ustrip
import Unitful
using GNSSSignals: GPSL1CA, GPSL1C_P, gen_code, get_code_frequency
import Tracking
using Tracking:
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    NoiseRefCN0Estimator,
    default_cn0_estimator,
    get_fallback_cn0_estimator,
    get_prompt_buffer,
    requires_noise_density,
    update,
    estimate_cn0,
    TrackedSignal,
    TrackedSat,
    TrackState,
    add_satellite!,
    get_cn0_estimator,
    get_last_fully_integrated_integration_time,
    track

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

    # `num_prompts_for_cn0_estimation` still sizes the default estimator's
    # averaging span.
    default = get_cn0_estimator(TrackedSignal(gpsl1; num_prompts_for_cn0_estimation = 40))
    @test default isa NoiseRefCN0Estimator
    @test default.num_records == 40
    @test Base.length(default.buffered_cn0) == 40
    # The default reads a measured noise density, which is what makes the band
    # provision a source for it.
    @test requires_noise_density(default_cn0_estimator(gpsl1, 100))
    @test requires_noise_density(default_cn0_estimator(GPSL1C_P(), 100))
end

@testset "estimate_cn0 overloads on TrackState" begin
    # The group-key / sat-id forwarding paths. An unseeded estimator reports
    # `-Inf dB-Hz` (see `NoCN0Estimator`), so this only checks dispatch.
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1, carrier_doppler = 0Hz)
    @test estimate_cn0(ts, :default, 1) == -Inf * dBHz
    @test estimate_cn0(ts, 1) == -Inf * dBHz
end

@testset "estimate_cn0 divides by the record's real integration time" begin
    # A record spanning N code blocks carries N times the SNR of a one-block
    # record, so `estimate_cn0` must divide by N × the code period; dividing by
    # the code period alone over-reports by 10·log₁₀(N).
    gpsl1 = GPSL1CA()
    prn = 1

    # Fixed amplitude with noise shrinking as √N is what a longer coherent
    # integration delivers, so both block counts must report the same C/N₀.
    # Pinned to `MomentsCN0Estimator`, which takes one `T` at `estimate_cn0`
    # time; `NoiseRefCN0Estimator` applies `T` per record (tested in
    # `cn0_estimators/noise_ref.jl`).
    function cn0_at(num_blocks, prompts)
        tsig = Tracking.TrackedSignal(gpsl1; cn0_estimator = MomentsCN0Estimator(100))
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

    # C/N₀ × `get_last_fully_integrated_integration_time` is the
    # post-integration SNR, so it must not depend on the block count.
    function snr_at(num_blocks, prompts)
        tsig = Tracking.TrackedSignal(gpsl1; cn0_estimator = MomentsCN0Estimator(100))
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
