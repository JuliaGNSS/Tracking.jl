module NoiseRefInTheLoopTest

# Why `NoiseReferencedCN0Estimator` is the default, in the loop (that it serves the
# signals `NWPRCN0Estimator` cannot, #217, is tested in TrackingLoops):
#
#   1. It beats NWPR where lock decisions are made, and is immune to the
#      phase-noise wash at long coherent integrations.
#   2. Its one bias, the reference's self-leakage, is the size the design says.

using Test: @test, @testset
using Random: Xoshiro, randn, rand
using Statistics: std
using Unitful: Hz, dBHz, ustrip, uconvert
using GNSSSignals: GPSL1CA, gen_code, get_chip_rate
import Tracking
using Tracking:
    CPUDownconvertAndCorrelator,
    TrackState,
    TrackedSat,
    reset_loop_filters!,
    set_preferred_num_code_periods_to_integrate!,
    track!
using TrackingLoops:
    CorrelatorNoiseEstimator,
    NWPRCN0Estimator,
    NoiseReferencedCN0Estimator,
    estimate_cn0,
    get_noise_density,
    is_bit_or_secondary_code_synced

_db(x) = ustrip(uconvert(dBHz, x))
_median(xs) = (v = sort(collect(xs)); v[div(length(v) + 1, 2)])

# One satellite through `track!`, data-modulated like a real L1 C/A signal.
# `code_periods > 1` lengthens the coherent integration once bit sync is in, which is the
# only point the cap allows it.
function _track_noisy(cn0_db, num_code_periods; seed = 1, cn0_estimator, code_periods = 1)
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    num_samples = 4000
    code_frequency = get_chip_rate(gpsl1)
    rng = Xoshiro(seed)
    amplitude = 10^(cn0_db / 20)
    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator)])
    dc = CPUDownconvertAndCorrelator()
    data_bits = rand(rng, (-1.0, 1.0), div(num_code_periods, 20) + 2)
    promoted = false
    for code_period_index = 1:num_code_periods
        code_phase = (code_period_index - 1) * num_samples * code_frequency / fs
        clean =
            amplitude .* data_bits[div(code_period_index-1, 20)+1] .*
            gen_code(num_samples, gpsl1, 1, fs, code_frequency, ustrip(code_phase))
        noise = randn(rng, ComplexF64, num_samples) .* sqrt(ustrip(Hz, fs))
        track!(clean .+ noise, track_state, fs; downconvert_and_correlator = dc)
        if code_periods > 1 && !promoted && is_bit_or_secondary_code_synced(track_state, 1)
            set_preferred_num_code_periods_to_integrate!(track_state, 1, code_periods)
            reset_loop_filters!(track_state, 1)
            promoted = true
        end
    end
    track_state
end

@testset "in the loop, against NWPR" begin
    # A quick version (5 seeds, 700 code periods) of the table in
    # `default_cn0_estimator`'s docstring. Medians over seeds.
    seeds = 1:5
    gpsl1 = GPSL1CA()
    run(cn0, est; code_periods = 1) = [
        _db(
            estimate_cn0(
                _track_noisy(cn0, 700; seed, cn0_estimator = est(), code_periods),
                1,
            ),
        ) for seed in seeds
    ]

    # 1-code-period integrations, 25 dB-Hz: the noise reference has no `1 < μ̂ < M` bound,
    # so it never degenerates. NWPR's degenerate count is not asserted: five
    # seeds cannot pin it and it varies with the Julia version's RNG stream.
    ref_low = run(25.0, () -> NoiseReferencedCN0Estimator())
    @test count(!isfinite, ref_low) == 0
    @test _median(ref_low) ≈ 25 atol = 3.0

    # 45 dB-Hz, 1-code-period integrations: both work, the reference is tighter and closer.
    nwpr_high = run(45.0, () -> NWPRCN0Estimator(gpsl1))
    ref_high = run(45.0, () -> NoiseReferencedCN0Estimator())
    @test _median(ref_high) ≈ 45 atol = 1.0
    @test abs(_median(ref_high) - 45) <= abs(_median(nwpr_high) - 45) + 0.1
    @test std(ref_high) < std(nwpr_high)

    # 20-code-period integrations (a whole L1 C/A bit): NWPR falls back for good, the
    # non-coherent reference is unaffected by the loop's phase noise.
    nwpr_long = run(45.0, () -> NWPRCN0Estimator(gpsl1); code_periods = 20)
    ref_long = run(45.0, () -> NoiseReferencedCN0Estimator(); code_periods = 20)
    @test _median(ref_long) ≈ 45 atol = 1.5
    @test _median(ref_long) - _median(nwpr_long) > 5
    # ... and the reference barely moves between the two integration lengths.
    @test abs(_median(ref_long) - _median(ref_high)) < 1.0
end

@testset "the reference's self-leakage is the size the design says" begin
    # `ε_self = C/f_chip`; see `NoiseReferencedCN0Estimator`. A paired comparison with
    # and without the satellite (same noise, same seeded draws), so the loop's σ
    # cancels.
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    num_samples = 40_000
    f_chip = ustrip(Hz, get_chip_rate(gpsl1))
    rng = Xoshiro(20260806)
    noise = sqrt(ustrip(Hz, fs)) .* randn(rng, ComplexF64, num_samples)

    function density(samples)
        # A freshly seeded generator per run, so both borrow the same codes in
        # the same order *and* despread them at the same phases and carriers.
        ts = TrackState(
            gpsl1,
            [TrackedSat(gpsl1, 1, 0.0, 0.0Hz)];
            noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(; rng = Xoshiro(1)),),
        )
        dc = CPUDownconvertAndCorrelator()
        for _ = 1:10
            track!(samples, ts, fs; downconvert_and_correlator = dc)
        end
        ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1CA))
    end

    for cn0 in (45.0, 50.0)
        carrier_power = 10^(cn0 / 10)
        satellite =
            10^(cn0 / 20) .* gen_code(num_samples, gpsl1, 1, fs, get_chip_rate(gpsl1), 0.0)
        ratio = density(ComplexF64.(satellite) .+ noise) / density(noise)
        # `N̂₀ → N₀ + C/f_chip`, i.e. a ratio of `1 + (C/N₀)/f_chip` against the
        # `N₀ = 1 Hz⁻¹` this noise is built to have: 1.031 at 45 dB-Hz, 1.098 at
        # 50. That is 0.13 dB and 0.40 dB of C/N₀ read low.
        @test ratio ≈ 1 + carrier_power / f_chip rtol = 0.25
        @test ratio > 1
    end
end

end
