module NWPRCN0EstimatorTest

using Test: @test, @testset, @inferred
using Random: Random, Xoshiro
using Unitful: kHz, MHz, Hz, dBHz, ustrip, uconvert
using GNSSSignals: GPSL1CA, get_code, gen_code
using GNSSSignals: get_code_frequency
import Tracking
using Tracking: TrackedSat, TrackState, get_cn0_estimator, track
using TrackingLoops:
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    estimate_cn0,
    has_bit_or_secondary_code_been_found

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
    # 20 pre-sync windows near NWPR's upper limit (`µ̂ → M`), where one
    # realisation spreads widely. Over 300 seeds: median 44.4, extremes 42.9 and
    # 46.1.
    @test cn0_estimate ≈ 45dBHz atol = 2.5dBHz
end

@testset "NWPR CN0 estimator in the loop (issue #217)" begin
    # End to end through `track` on data-modulated GPS L1 C/A: the aligned
    # post-sync window engages, the estimate lands on the truth, and pure noise
    # does not read like a signal.
    gpsl1 = GPSL1CA()
    sampling_frequency = 4e6Hz
    num_samples = 4000
    code_frequency = get_code_frequency(gpsl1)

    # `fade_to` continues at another C/N₀ for `fade_blocks` more records, to
    # reach post-sync at a C/N₀ where bit sync would never lock from cold.
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
    # Over 120 seeds: median 44.4, extremes 43.4 and 45.7, worst five-seed
    # median 0.96 dB off (issue #219 regression guard).
    median_of(xs) = sort(collect(xs))[div(length(xs) + 1, 2)]
    @test median_of(
        ustrip(uconvert(dBHz, estimate_cn0(track_noisy(45.0, 900; seed), 1))) for seed = 1:5
    ) ≈ 45 atol = 1.5

    # Pure noise, pre-sync: moments reports a signal, NWPR does not. 25 rather
    # than 20 dB-Hz because a short pre-sync window has a heavy upper tail.
    noise_state = track_noisy(nothing, 200)
    noise_cn0 = ustrip(uconvert(dBHz, estimate_cn0(noise_state, 1)))
    @test !has_bit_or_secondary_code_been_found(noise_state, 1)
    @test noise_cn0 < 25
    # Medians over seeds: one noise run is one draw of a wide distribution, and
    # which draw a seed gives differs between Julia versions. Over 100 seeds
    # moments read above 20 dB-Hz in 92 and above 25 in 64; its medians over
    # disjoint groups of 15 seeds lay in 24.8–28.4 dB-Hz, NWPR's were all `-Inf`.
    noise_seeds = 1:15
    moments_noise_cn0 = median_of(
        ustrip(
            uconvert(
                dBHz,
                estimate_cn0(
                    track_noisy(
                        nothing,
                        200;
                        seed,
                        cn0_estimator = MomentsCN0Estimator(100),
                    ),
                    1,
                ),
            ),
        ) for seed in noise_seeds
    )
    nwpr_noise_cn0 = median_of(
        ustrip(uconvert(dBHz, estimate_cn0(track_noisy(nothing, 200; seed), 1))) for
        seed in noise_seeds
    )
    @test moments_noise_cn0 > 20
    @test nwpr_noise_cn0 < 20

    # Faded after bit sync: the case the coherence cap exists for (see
    # `NWPRCN0Estimator`). Locked in at 45 dB-Hz, faded to 25 dB-Hz; a median over
    # seeds, since one low-C/N₀ run says little.
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
    # The staged loop holds phase over a whole bit at 25 dB-Hz, so the whole-bit
    # window does not read low either. Over 48 seeds in groups of 16: the capped
    # median stayed in 24.1–25.0 dB-Hz and 0.1–0.33 dB above the whole-bit
    # window's, and neither window reported `-Inf`.
    @test median_of(capped_cn0) ≈ 25 atol = 1.5
    @test median_of(capped_cn0) - median_of(full_bit_cn0) > -0.5
    @test !any(isinf, capped_cn0)
end

end
