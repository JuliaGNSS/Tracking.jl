module CN0EstimationTest

using Test: @test, @testset, @inferred
using Random: Random
using Unitful: kHz, MHz, Hz, ms, dBHz, ustrip
import Unitful
using StaticArrays: SVector
using GNSSSignals: GPSL1CA, get_code
import Tracking
using Tracking:
    MomentsCN0Estimator,
    get_prompt_buffer,
    get_current_index,
    update,
    EarlyPromptLateCorrelator,
    get_prompt,
    get_correlator_sample_shifts,
    estimate_cn0,
    TrackedSat,
    TrackState,
    add_satellite!,
    get_cn0_estimator,
    get_last_fully_integrated_integration_time,
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
