module CorrelatorNoiseEstimatorTest

using Test: @test, @testset, @test_logs, collect_test_logs
using Logging: Warn
using Random: Xoshiro, randn
using StaticArrays: SMatrix
using Unitful: Hz, s, ms, dBHz, ustrip, uconvert
using GNSSSignals:
    GPSL1CA,
    GPSL1C_P,
    GPSL5I,
    GalileoE1B_BOC11,
    gen_code,
    get_code_frequency,
    get_code_length
import Tracking
using Tracking:
    BandMeasurement,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    CorrelatorNoiseEstimator,
    DefaultPostCorrFilter,
    EarlyPromptLateCorrelator,
    Int16DownconvertAndCorrelator,
    NWPRCN0Estimator,
    NoiseDensity,
    NoiseRefCN0Estimator,
    NumAnts,
    OneBitDownconvertAndCorrelator,
    TrackState,
    TrackedSat,
    TwoBitDownconvertAndCorrelator,
    add_satellite!,
    downconvert_and_correlate!,
    get_correlator_sample_shifts,
    estimate_cn0,
    get_noise_density,
    get_weights,
    noise_density_type,
    noise_window_looks,
    track!

const FS = 4e6Hz
const NUM_SAMPLES = 40_000                      # 10 ms at 4 MHz
const GPSL1 = GPSL1CA()

_db(x) = ustrip(uconvert(dBHz, x))
_median(xs) = (v = sort(collect(xs)); v[div(length(v) + 1, 2)])

# One satellite at `cn0_dbhz` in noise of per-sample variance `f_s`, so the true
# density is exactly `N₀ = 1 Hz⁻¹`.
function _sky(cn0_dbhz; seed = 1, num_samples = NUM_SAMPLES, fs = FS, prn = 1)
    rng = Xoshiro(seed)
    code = gen_code(num_samples, GPSL1, prn, fs, get_code_frequency(GPSL1), 0.0)
    amplitude = isfinite(cn0_dbhz) ? 10^(cn0_dbhz / 20) : 0.0
    ComplexF64.(amplitude .* code) .+
    sqrt(ustrip(Hz, fs)) .* randn(rng, ComplexF64, num_samples)
end

# A hardware-style source with no look count; the rank gate must not withhold it.
struct UncountedNoiseSource <: Tracking.AbstractNoiseEstimator end

_to_int16(signal) = Complex{Int16}.(
    round.(Int16, clamp.(real(signal) ./ 4, -2047, 2047)),
    round.(Int16, clamp.(imag(signal) ./ 4, -2047, 2047)),
)

# 20 calls: the window's relative error is `1/√(3K)`, so `K = 200` gives σ ≈ 0.22 dB
# and a 1 dB gate is ≈4.5σ for any draw (at `K = 50` it would be ≈0.6 dB).
function _tracked(signal, dc; fs = FS, num_calls = 20, kwargs...)
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())];
        kwargs...,
    )
    for _ = 1:num_calls
        track!((L1 = BandMeasurement(signal, fs),), ts; downconvert_and_correlator = dc)
    end
    ts
end

@testset "the software path measures the signal's noise density" begin
    # Empty sky, known floor; window relative error ≈1/√(3·200) ≈ 4 %.
    ts = _tracked(_sky(-Inf; seed = 7), CPUDownconvertAndCorrelator())
    estimator = ts.noise_estimators.GPSL1CA
    @test Base.length(estimator) == 200
    @test ustrip(Hz^-1, get_noise_density(estimator)) ≈ 1.0 rtol = 0.15
end

@testset "the software path needs no `append_noise_observation!` and never warns" begin
    # Counterpart of the never-fed warning in `test/cn0_estimators/noise_ref.jl`:
    # the sample-driven path measures before the fold reads, so no warning.
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    signal = _sky(45.0)
    logs, _ = collect_test_logs() do
        track!((L1 = BandMeasurement(signal, FS),), ts)
    end
    @test isempty(filter(r -> r.level == Warn, logs))
    @test Base.length(ts.noise_estimators.GPSL1CA) > 0
    @test isfinite(_db(estimate_cn0(ts, 1)))
end

@testset "every backend measures on its own kernel ($(nameof(DC)))" for (
    DC,
    build,
    quantise,
) in (
    (CPUDownconvertAndCorrelator, () -> CPUDownconvertAndCorrelator(), false),
    (
        CPUThreadedDownconvertAndCorrelator,
        () -> CPUThreadedDownconvertAndCorrelator(),
        false,
    ),
    (Int16DownconvertAndCorrelator, () -> Int16DownconvertAndCorrelator(NUM_SAMPLES), true),
    (OneBitDownconvertAndCorrelator, () -> OneBitDownconvertAndCorrelator(), true),
    (TwoBitDownconvertAndCorrelator, () -> TwoBitDownconvertAndCorrelator(), true),
)
    # The reference must use the prompt's kernel (see `CorrelatorNoiseEstimator`),
    # so each backend reports the C/N₀ its own quantiser delivers with no
    # correction: 1-bit ≈2 dB lower (plus ≈1 dB for the 1-bit carrier), 2-bit ≈1 dB.
    signal = _sky(45.0; seed = 3)
    ts = _tracked(quantise ? _to_int16(signal) : ComplexF32.(signal), build())
    cn0 = _db(estimate_cn0(ts, 1))
    @test Base.length(ts.noise_estimators.GPSL1CA) == 200
    if DC === OneBitDownconvertAndCorrelator
        @test 41.0 < cn0 < 44.5           # ≈2 dB of hard-limiting loss
    elseif DC === TwoBitDownconvertAndCorrelator
        @test 43.0 < cn0 < 45.5
    else
        @test cn0 ≈ 45.0 atol = 1.0
    end
end

@testset "observation cadence follows the chunk, not the buffer" begin
    # `num_sub = max(1, round(chunk_duration / code_period))` per chunk, and none
    # from the drain pass — gating on `samples_unchanged` would freeze the window
    # after chunk 0.
    signal = _sky(-Inf; seed = 11)
    for (interval, expected_num_sub) in ((1ms, 1), (4ms, 4), (20ms, 20))
        ts = TrackState(
            GPSL1,
            [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
        )
        # 40 ms of buffer: 40 / 20 / 2 chunks at the three intervals.
        long = repeat(signal, 4)
        num_chunks = ceil(Int, 40 / ustrip(ms, interval))
        track!((L1 = BandMeasurement(long, FS),), ts; doppler_update_interval = interval)
        @test Base.length(ts.noise_estimators.GPSL1CA) == num_chunks * expected_num_sub
    end

    # Another sampling frequency: the count follows duration, not sample count.
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    fast = _sky(-Inf; seed = 12, num_samples = 100_000, fs = 10e6Hz)
    track!((L1 = BandMeasurement(fast, 10e6Hz),), ts; doppler_update_interval = 4ms)
    # 10 ms at a 4 ms chunk: 4 + 4, then a 2 ms trailing chunk measured at its
    # true length.
    @test Base.length(ts.noise_estimators.GPSL1CA) == 4 + 4 + 2
end

@testset "a code period longer than the chunk still yields one observation" begin
    # The `max(1, …)`: GPS L1C-P (10 ms code) at 1 ms chunks, reachable when
    # another signal forces a shorter `doppler_update_interval`.
    fs = 20e6Hz
    l1cp = GPSL1C_P()
    num_samples = 200_000                       # 10 ms — one L1C-P code period
    signal = ComplexF32.(sqrt(ustrip(Hz, fs)) .* randn(Xoshiro(5), ComplexF64, num_samples))
    ts = TrackState(
        l1cp,
        [TrackedSat(l1cp, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    track!((L1 = BandMeasurement(signal, fs),), ts; doppler_update_interval = 1ms)
    # Ten 1 ms chunks, each a tenth of a code period — one slice apiece, not zero.
    @test Base.length(ts.noise_estimators.GPSL1C_P) == 10
    # 30 looks ⇒ σ ≈ 18 %: a sanity check, not a precision claim.
    @test ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1C_P)) ≈ 1.0 rtol = 0.55
end

@testset "the PRN rotates through the whole family, tracked codes included" begin
    # The position is carried in the window's observations (see `_next_noise_prn`).
    # Tracked PRNs are not skipped; the random phase makes that safe (see
    # `CorrelatorNoiseEstimator` and the two tests below).
    ts = TrackState(;
        signal = GPSL1,
        noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
    )
    add_satellite!(ts; prn = 3, carrier_doppler = 0.0Hz)
    add_satellite!(ts; prn = 4, carrier_doppler = 0.0Hz)
    track!((L1 = BandMeasurement(_sky(-Inf; seed = 2), FS),), ts)
    prns = [obs.prn for obs in ts.noise_estimators.GPSL1CA.buffered]
    @test length(prns) == 10
    @test allunique(prns)                       # it really does advance
    @test prns == Int16[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
end

@testset "the code phase and carrier are drawn per sub-integration" begin
    # The draws must reach the kernel and be reproducible from a supplied `rng`
    # (`cn0_estimators/noise_ref_in_the_loop.jl` depends on the latter).
    m = BandMeasurement(_sky(-Inf; seed = 4), FS)
    dc = CPUDownconvertAndCorrelator()
    function densities(seed)
        estimator =
            CorrelatorNoiseEstimator(; window_duration = 1000.0s, rng = Xoshiro(seed))
        Tracking.update_noise!(
            estimator,
            m,
            1,
            NUM_SAMPLES,
            Tracking.NoiseUpdateContext(GPSL1, 0, dc),
        )
        [ustrip(Hz^-1, o.noise_density) for o in estimator.buffered]
    end
    @test densities(1) == densities(1)          # seeded, so a run is reproducible
    @test densities(1) != densities(2)          # ... and the draws reach the kernel
end

@testset "a random code phase keeps a present-but-untracked signal off the reference" begin
    # Why the phase is randomised (see `CorrelatorNoiseEstimator`): `_sky` puts an
    # untracked PRN 7 at code phase 0 and zero Doppler, where a fixed-phase
    # reference would hit it on every PRN-7 observation (≈10.5·N₀ at 45 dB-Hz).
    m = BandMeasurement(_sky(45.0; seed = 3, prn = 7), FS)
    dc = CPUDownconvertAndCorrelator()
    # `carrier_dither = 0` isolates the code-phase draw.
    estimator = CorrelatorNoiseEstimator(;
        window_duration = 1000.0s,
        carrier_dither = 0.0Hz,
        rng = Xoshiro(7),
    )
    context = Tracking.NoiseUpdateContext(GPSL1, 0, dc)
    for _ = 1:40
        Tracking.update_noise!(estimator, m, 1, NUM_SAMPLES, context)
    end
    on_seven = [ustrip(Hz^-1, o.noise_density) for o in estimator.buffered if o.prn == 7]
    others = [ustrip(Hz^-1, o.noise_density) for o in estimator.buffered if o.prn != 7]
    @test length(on_seven) >= 10
    # Medians: a rare outlier is expected, a standing bias (ratio ≈11 at a fixed
    # phase) is not.
    @test _median(on_seven) / _median(others) ≈ 1.0 rtol = 0.6
end

@testset "the reference's taps are spaced wide enough to be independent" begin
    # Why ≥1 chip: see `tap_code_shift` in `CorrelatorNoiseEstimator`. This pins
    # that sample rounding at the band's `f_s` keeps 1.5 chips at ≥1 chip.
    estimator = CorrelatorNoiseEstimator()
    correlator = EarlyPromptLateCorrelator(;
        preferred_early_late_to_prompt_code_shift = estimator.tap_code_shift,
    )
    for fs in (2e6Hz, 4e6Hz, 10e6Hz, 20e6Hz)
        shifts = get_correlator_sample_shifts(correlator, fs, get_code_frequency(GPSL1))
        spacing_chips =
            (shifts[end] - shifts[1]) / 2 * ustrip(Hz, get_code_frequency(GPSL1)) /
            ustrip(Hz, fs)
        @test spacing_chips >= 1.0
    end
end

@testset "two bands keep separate windows at different sampling rates" begin
    l5 = GPSL5I()
    fs_l1, fs_l5 = 4e6Hz, 20e6Hz
    ts = TrackState(;
        signals = (l1 = (GPSL1,), l5 = (l5,)),
        noise_estimators = (
            GPSL1CA = CorrelatorNoiseEstimator(),
            GPSL5I = CorrelatorNoiseEstimator(),
        ),
    )
    add_satellite!(ts; prn = 1, group = :l1, carrier_doppler = 0.0Hz)
    add_satellite!(ts; prn = 1, group = :l5, carrier_doppler = 0.0Hz)
    # A four-times louder L5 front end, so the two cannot be confused.
    rng = Xoshiro(21)
    l1_samples = ComplexF32.(sqrt(ustrip(Hz, fs_l1)) .* randn(rng, ComplexF64, 40_000))
    l5_samples = ComplexF32.(2 * sqrt(ustrip(Hz, fs_l5)) .* randn(rng, ComplexF64, 200_000))
    measurements =
        (L1 = BandMeasurement(l1_samples, fs_l1), L5 = BandMeasurement(l5_samples, fs_l5))
    # 20 calls, see `_tracked`.
    for _ = 1:20
        track!(measurements, ts)
    end
    d_l1 = ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1CA))
    d_l5 = ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL5I))
    @test Base.length(ts.noise_estimators.GPSL1CA) == 200
    @test Base.length(ts.noise_estimators.GPSL5I) == 200
    @test d_l1 ≈ 1.0 rtol = 0.2
    @test d_l5 ≈ 4.0 rtol = 0.2
end

@testset "two signals on one band measure different floors under a CW tone" begin
    # Why the estimator is keyed by signal (see `AbstractNoiseEstimator`): the tone
    # moves BPSK(1) and BOC(1,1) floors in *opposite* directions, which no common
    # scale factor can. `GalileoE1B_BOC11` because CBOC needs f_s ≥ 12.276 MHz.
    fs = 4e6Hz
    n = 40_000
    e1b = GalileoE1B_BOC11()

    function densities(tone_frequency; seed = 17, jammer_to_noise_db = 13)
        rng = Xoshiro(seed)
        σ = sqrt(ustrip(Hz, fs))
        samples = σ .* randn(rng, ComplexF64, n)
        if !isnothing(tone_frequency)
            t = (0:(n-1)) ./ ustrip(Hz, fs)
            samples .+= (σ * 10^(jammer_to_noise_db / 20)) .* cis.(2π * tone_frequency .* t)
        end
        ts = TrackState(; signals = (gps = (GPSL1,), galileo = (e1b,)))
        add_satellite!(ts; prn = 1, group = :gps, carrier_doppler = 0.0Hz)
        add_satellite!(ts; prn = 1, group = :galileo, carrier_doppler = 0.0Hz)
        for _ = 1:5
            track!((L1 = BandMeasurement(ComplexF32.(samples), fs),), ts)
        end
        (
            l1ca = ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1CA)),
            e1b = ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GalileoE1B_BOC11)),
        )
    end

    # White noise: the two agree.
    white = densities(nothing)
    @test white.l1ca ≈ 1.0 rtol = 0.2
    @test white.e1b ≈ 1.0 rtol = 0.2
    @test white.e1b / white.l1ca ≈ 1.0 rtol = 0.25

    # At the chip rate (BPSK(1) null, BOC(1,1) peak): C/A ≈1.4× thermal, E1B ≈39×.
    at_chip_rate = densities(1.023e6)
    @test at_chip_rate.l1ca < 2 * white.l1ca
    @test at_chip_rate.e1b > 10 * at_chip_rate.l1ca

    # Inside C/A's main lobe the order reverses.
    in_main_lobe = densities(0.4e6)
    @test in_main_lobe.l1ca > 10 * white.l1ca
    @test in_main_lobe.l1ca > 1.15 * in_main_lobe.e1b
end

@testset "a signal carried by two groups is measured once" begin
    # Noise work is resolved per signal (`_noise_items`), not per group, so two
    # groups carrying one signal must not double its window.
    one_group = TrackState(;
        signals = (legacy_gps = (GPSL1,),),
        noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
    )
    two_groups = TrackState(;
        signals = (legacy_gps = (GPSL1,), also_l1 = (GPSL1,)),
        noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
    )
    add_satellite!(one_group; prn = 1, group = :legacy_gps, carrier_doppler = 0.0Hz)
    add_satellite!(two_groups; prn = 1, group = :legacy_gps, carrier_doppler = 0.0Hz)
    add_satellite!(two_groups; prn = 2, group = :also_l1, carrier_doppler = 0.0Hz)
    signal = _sky(-Inf; seed = 4)
    for ts in (one_group, two_groups)
        track!((L1 = BandMeasurement(signal, FS),), ts)
    end
    @test Base.length(one_group.noise_estimators.GPSL1CA) ==
          Base.length(two_groups.noise_estimators.GPSL1CA)
end

@testset "the group whose loop carries the despread may be empty" begin
    # The despread rides the first group's loop, which may hold no satellites;
    # an early return there would drop the measurement.
    ts = TrackState(;
        signals = (empty_first = (GPSL1,), with_sats = (GPSL1,)),
        noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
    )
    add_satellite!(ts; prn = 1, group = :with_sats, carrier_doppler = 0.0Hz)
    @test isempty(Tracking.get_sat_states(ts, :empty_first))
    track!((L1 = BandMeasurement(_sky(-Inf; seed = 11), FS),), ts)
    @test Base.length(ts.noise_estimators.GPSL1CA) == 10
    @test ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1CA)) ≈ 1.0 rtol = 0.3
end

@testset "the serial and threaded loops measure bit-identically" begin
    # One work item is the sole writer of its window and RNG, so the seeded draws
    # match and the windows agree bit for bit; catches shared state in `@batch`.
    signal = ComplexF32.(_sky(45.0; seed = 12))
    densities =
        map((CPUDownconvertAndCorrelator(), CPUThreadedDownconvertAndCorrelator()),) do dc
            ts = _tracked(signal, dc)
            (
                Base.length(ts.noise_estimators.GPSL1CA),
                get_noise_density(ts.noise_estimators.GPSL1CA),
            )
        end
    @test densities[1][1] == densities[2][1]
    @test densities[1][2] === densities[2][2]
end

@testset "an unchunked call measures once, and a re-pass says so" begin
    # `measure_noise` decides; `samples_unchanged` (a cache-reuse hint) must not.
    signal = _sky(-Inf; seed = 6)
    measurements = (L1 = BandMeasurement(signal, FS),)
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    dc = CPUDownconvertAndCorrelator()
    # Unchunked: 10 ms in 1 ms code-period slices.
    downconvert_and_correlate!(dc, measurements, ts)
    @test Base.length(ts.noise_estimators.GPSL1CA) == 10
    # A declared re-pass must not re-measure what has already been covered.
    downconvert_and_correlate!(dc, measurements, ts; measure_noise = false)
    @test Base.length(ts.noise_estimators.GPSL1CA) == 10
    # `samples_unchanged` alone suppresses nothing.
    downconvert_and_correlate!(dc, measurements, ts; samples_unchanged = true)
    @test Base.length(ts.noise_estimators.GPSL1CA) == 20
end

@testset "a band this call brought no samples for is skipped, not thrown on" begin
    # Advancing a multi-band `TrackState` one band at a time is legitimate; the
    # noise pass used to index the absent band and throw a `FieldError`.
    for measure_noise in (true, false)
        ts = TrackState(; signals = (gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
        add_satellite!(ts; prn = 1, group = :gps_l1, carrier_doppler = 0.0Hz)
        @test keys(ts.noise_estimators) == (:GPSL1CA, :GPSL5I)
        downconvert_and_correlate!(
            CPUDownconvertAndCorrelator(),
            (L1 = BandMeasurement(_sky(-Inf; seed = 21), FS),),
            ts;
            measure_noise,
        )
        # L5 brought no samples, so it measured nothing either way ...
        @test Base.length(ts.noise_estimators.GPSL5I) == 0
        # ... and L1 measured exactly as `measure_noise` says.
        @test (Base.length(ts.noise_estimators.GPSL1CA) > 0) == measure_noise
    end

    # `_noise_items` stays type-stable and allocation-free.
    ts = TrackState(; signals = (gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
    add_satellite!(ts; prn = 1, group = :gps_l1, carrier_doppler = 0.0Hz)
    dc = CPUDownconvertAndCorrelator()
    measurements = (L1 = BandMeasurement(_sky(-Inf; seed = 22), FS),)
    items = Tracking._noise_items(dc, ts, measurements, 0, nothing)
    # One descriptor, for L1 alone — its length is a property of the types.
    @test items isa Base.RefValue{<:Tuple{Any}}
    @test Base.return_types(
        Tracking._noise_items,
        (typeof(dc), typeof(ts), typeof(measurements), Int, Nothing),
    ) == [typeof(items)]
    resolve(dc, ts, m) = Tracking._noise_items(dc, ts, m, 0, nothing)
    resolve(dc, ts, measurements)
    @test @allocated(resolve(dc, ts, measurements)) == 0
end

@testset "a signal without a consumer is never despread" begin
    # Provisioning is gated on `requires_noise_density`: NWPR users pay nothing.
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator())],
    )
    @test ts.noise_estimators === NamedTuple()
    track!((L1 = BandMeasurement(_sky(45.0), FS),), ts)
    @test ts.noise_estimators === NamedTuple()
end

# Julia >= 1.11 only (as in `test/track_in_place.jl`): on 1.10 the measurement
# path leaves ~670 B per sub-integration (~5 kB per call), 0 B without the reference.
# Every backend is checked because the kernels differ; the one-/two-bit kernels once
# returned per-tap sums in a heap `Vector` here (now an `SVector{NC}`).
if VERSION >= v"1.11"
    @testset "the software measurement is allocation-free in steady state ($(nameof(DC)))" for (
        DC,
        build,
        quantise,
    ) in (
        (CPUDownconvertAndCorrelator, () -> CPUDownconvertAndCorrelator(), false),
        (
            Int16DownconvertAndCorrelator,
            () -> Int16DownconvertAndCorrelator(NUM_SAMPLES),
            true,
        ),
        (OneBitDownconvertAndCorrelator, () -> OneBitDownconvertAndCorrelator(), true),
        (TwoBitDownconvertAndCorrelator, () -> TwoBitDownconvertAndCorrelator(), true),
    )
        function measure(dc, measurements, ts)
            for _ = 1:8
                downconvert_and_correlate!(dc, measurements, ts; chunk_index = 0)
            end
            @allocated downconvert_and_correlate!(dc, measurements, ts; chunk_index = 0)
        end
        raw = _sky(45.0; seed = 8)
        signal = quantise ? _to_int16(raw) : ComplexF32.(raw)
        measurements = (L1 = BandMeasurement(signal, FS),)
        # Paired with a reference-free state, isolating the measurement.
        state(cn0) =
            TrackState(GPSL1, [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = cn0())])
        ts = state(NoiseRefCN0Estimator)
        plain_ts = state(() -> NWPRCN0Estimator(GPSL1))
        dc = build()
        plain_dc = build()
        track!(measurements, ts; downconvert_and_correlator = dc)
        track!(measurements, plain_ts; downconvert_and_correlator = plain_dc)
        @test measure(dc, measurements, ts) == 0
        @test measure(plain_dc, measurements, plain_ts) == 0
    end
end

@testset "a DC offset only bites on a partial code period" begin
    # A full code period rejects a DC offset (see `CorrelatorNoiseEstimator`); only
    # a code period longer than the chunk, at zero IF, is exposed.
    fs = 4e6Hz
    n = 40_000
    rng = Xoshiro(31)
    offset = 0.10 * sqrt(ustrip(Hz, fs))       # d/σ = 0.10
    noise = sqrt(ustrip(Hz, fs)) .* randn(rng, ComplexF64, n)
    with_offset = ComplexF32.(noise .+ offset)

    function density(signal, interval)
        ts = TrackState(;
            signal = GPSL1,
            noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
        )
        add_satellite!(ts; prn = 1, carrier_doppler = 0.0Hz)
        track!((L1 = BandMeasurement(signal, fs),), ts; doppler_update_interval = interval)
        ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1CA))
    end

    # Whole code periods: the offset is rejected, so the density is unmoved.
    @test density(with_offset, 1ms) ≈ density(ComplexF32.(noise), 1ms) rtol = 0.1
end

@testset "the window's running totals stay exact" begin
    # The cached sums (`NoiseWindowTotals`) must equal a walk of the window after
    # thousands of appends with mixed entry sizes, where drift would show.
    exact_span(e) = sum(o.duration for o in e.buffered)
    exact_looks(e) = sum(o.num_sub_integrations for o in e.buffered)
    exact_density(e) =
        sum(o.num_sub_integrations * o.noise_density for o in e.buffered) / exact_looks(e)

    rng = Xoshiro(17)
    window = 50.0ms
    estimator = CorrelatorNoiseEstimator(; window_duration = window)
    for i = 1:5000
        # 1-in-20 entries is a long pre-averaged dump among short ones.
        duration = uconvert(s, (rand(rng) < 0.05 ? 200.0 : 1.0) * rand(rng) * ms)
        Tracking.append_noise_observation!(
            estimator,
            Tracking.NoiseObservation(
                (1.0 + 5rand(rng)) * 1e-10 / 1.0Hz,
                rand(rng, 1:64),
                duration,
                Int16(rand(rng, 1:32)),
            ),
        )
        i % 500 == 0 || continue
        @test estimator.totals[].span ≈ exact_span(estimator) rtol = 1e-10
        @test estimator.totals[].looks == exact_looks(estimator)
        @test get_noise_density(estimator) ≈ exact_density(estimator) rtol = 1e-10
    end

    # The window is minimal but never short of the configured span.
    @test estimator.totals[].span >= window
    @test estimator.totals[].span - first(estimator.buffered).duration < window
end

@testset "reading the window does not walk it" begin
    # Reads and appends are O(1) via the cached totals; timing is flaky, so check
    # a full window is allocation-free instead.
    estimator = CorrelatorNoiseEstimator(; window_duration = 1.0s)
    observation = Tracking.NoiseObservation(2e-10 / 1.0Hz, 1, uconvert(s, 0.4ms), Int16(3))
    for _ = 1:6000
        Tracking.append_noise_observation!(estimator, observation)
    end
    @test Base.length(estimator) > 2000        # a full window, not a short one
    get_noise_density(estimator)
    @test (@allocated get_noise_density(estimator)) == 0
    @test (@allocated Tracking.append_noise_observation!(estimator, observation)) == 0
end

@testset "the reference owns no scratch of its own ($(nameof(DC)))" for (
    DC,
    build,
    quantise,
) in (
    (CPUDownconvertAndCorrelator, () -> CPUDownconvertAndCorrelator(), false),
    (Int16DownconvertAndCorrelator, () -> Int16DownconvertAndCorrelator(NUM_SAMPLES), true),
    (OneBitDownconvertAndCorrelator, () -> OneBitDownconvertAndCorrelator(), true),
    (TwoBitDownconvertAndCorrelator, () -> TwoBitDownconvertAndCorrelator(), true),
)
    # No replica buffer in the estimator: `_despread_one_signal!` borrows the
    # backend's scratch (a byte per sample per signal saved).
    @test :code_replica ∉ fieldnames(CorrelatorNoiseEstimator)
    signal = _sky(45.0; seed = 5)
    ts = _tracked(quantise ? _to_int16(signal) : ComplexF32.(signal), build())
    estimator = ts.noise_estimators.GPSL1CA
    @test get_noise_density(estimator) !== nothing      # it still measured
end

@testset "the reference and the satellites despread out of one scratch slot" begin
    # One shared per-thread slot, however many signals; safe because the noise
    # pass completes before the satellite loop reuses it.
    dc = CPUDownconvertAndCorrelator()
    @test isempty(Tracking._scratch_buffers(dc).code_replica)
    ts = _tracked(ComplexF32.(_sky(45.0; seed = 5)), dc)
    @test get_noise_density(ts.noise_estimators.GPSL1CA) !== nothing
    # Sized to the whole buffer: `gen_code_replica!` writes *at* `start_sample`.
    @test length(Tracking._scratch_buffers(dc).code_replica) >= NUM_SAMPLES
end

# --- antenna arrays ---------------------------------------------------------

# `_sky` as an `M`-column array: one shared satellite plus independent noise per
# column scaled by `noise_scales[m]`, so off-diagonals are near zero and a wrong
# column is distinguishable (unlike `repeat`-based fixtures).
function _sky_array(
    cn0_dbhz,
    noise_scales;
    seed = 1,
    num_samples = NUM_SAMPLES,
    fs = FS,
    prn = 1,
)
    rng = Xoshiro(seed)
    code = gen_code(num_samples, GPSL1, prn, fs, get_code_frequency(GPSL1), 0.0)
    amplitude = isfinite(cn0_dbhz) ? 10^(cn0_dbhz / 20) : 0.0
    signal = ComplexF64.(amplitude .* code)
    sigma = sqrt(ustrip(Hz, fs))
    reduce(
        hcat,
        (
            signal .+ (scale * sigma) .* randn(rng, ComplexF64, num_samples) for
            scale in noise_scales
        ),
    )
end

function _tracked_array(signal, num_ants; fs = FS, num_calls = 20, kwargs...)
    ts = TrackState(
        GPSL1,
        [
            TrackedSat(
                GPSL1,
                1,
                0.0,
                0.0Hz;
                num_ants = NumAnts(num_ants),
                cn0_estimator = NoiseRefCN0Estimator(),
            ),
        ];
        kwargs...,
    )
    for _ = 1:num_calls
        track!((L1 = BandMeasurement(signal, fs),), ts)
    end
    ts
end

@testset "the default estimator is provisioned for its group's antenna count" begin
    # `noise_density_type` must follow the group's antenna count.
    @test noise_density_type(CorrelatorNoiseEstimator()) === NoiseDensity
    @test noise_density_type(CorrelatorNoiseEstimator(; num_ants = NumAnts(1))) ===
          NoiseDensity

    single = TrackState(GPSL1, [TrackedSat(GPSL1, 1, 0.0, 0.0Hz)])
    @test noise_density_type(single.noise_estimators.GPSL1CA) === NoiseDensity

    array = TrackState(GPSL1, [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; num_ants = NumAnts(3))])
    D = noise_density_type(array.noise_estimators.GPSL1CA)
    @test D <: SMatrix{3,3}
    @test Tracking._num_ants_of_density_type(D) === NumAnts(3)
end

@testset "an explicitly-passed estimator must match its group's antenna count" begin
    # Caught at construction rather than as a shape error inside the fold.
    err = try
        TrackState(
            GPSL1,
            [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; num_ants = NumAnts(3))];
            noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
        )
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("GPSL1CA", err.msg)
    @test occursin("NumAnts(3)", err.msg)

    # Matching counts are accepted.
    @test TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; num_ants = NumAnts(3))];
        noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(; num_ants = NumAnts(3)),),
    ) isa TrackState
    @test TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz)];
        noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(),),
    ) isa TrackState
end

@testset "an antenna array measures a spatial covariance" begin
    # Each antenna's floor is `scaleₘ²` Hz⁻¹.
    scales = (1.0, 1.5, 2.0)
    ts = _tracked_array(_sky_array(-Inf, scales; seed = 11), 3)
    R = get_noise_density(ts.noise_estimators.GPSL1CA)

    @test R isa SMatrix{3,3}
    # Diagonal: each antenna's own N₀, at the scalar case's 15 % gate.
    for m = 1:3
        @test ustrip(Hz^-1, real(R[m, m])) ≈ scales[m]^2 rtol = 0.15
    end
    # Off-diagonals: small (not zero — finite looks) against their diagonals.
    for m = 1:3, n = 1:3
        m == n && continue
        @test abs(R[m, n]) < 0.25 * sqrt(real(R[m, m]) * real(R[n, n]))
    end
    # Hermitian, since `Σ b·bᴴ` is.
    @test R ≈ R'
end

@testset "a covariance is withheld until it spans its own dimensions" begin
    # The rank gate, see `_sufficient_looks`: `M` looks, i.e. `⌈M/3⌉` observations.
    # 1 ms buffers give one observation per call, the only shape that hits the gate.
    one_ms(M, seed) = ComplexF32.(
        reduce(
            hcat,
            (sqrt(ustrip(Hz, FS)) .* randn(Xoshiro(seed), ComplexF64, 4000) for _ = 1:M),
        ),
    )
    array_state(M) = TrackState(
        GPSL1,
        [
            TrackedSat(
                GPSL1,
                1,
                0.0,
                0.0Hz;
                num_ants = NumAnts(M),
                cn0_estimator = NoiseRefCN0Estimator(),
            ),
        ],
    )

    # No gate at or below the tap count.
    for (M, expected_calls) in ((1, 1), (3, 1), (4, 2), (8, 3))
        ts = array_state(M)
        ready_at = 0
        for call = 1:4
            track!((L1 = BandMeasurement(one_ms(M, call), FS),), ts)
            _, ready = Tracking._noise_density_and_ready(ts.noise_estimators.GPSL1CA)
            if ready && ready_at == 0
                ready_at = call
            end
        end
        @test ready_at == expected_calls
    end

    # The look count is the window's own running total, three per observation.
    ts = array_state(4)
    track!((L1 = BandMeasurement(one_ms(4, 1), FS),), ts)
    est = ts.noise_estimators.GPSL1CA
    @test noise_window_looks(est) == 3
    # Only the fold defers; the public reader reports what the window holds.
    @test !isnothing(get_noise_density(est))
    @test Tracking._noise_window_filling(est)
    @test estimate_cn0(ts, 1) == -Inf * dBHz

    # A filling window is a normal transient: no warning.
    quiet = array_state(4)
    logs, _ = collect_test_logs() do
        track!((L1 = BandMeasurement(one_ms(4, 5), FS),), quiet)
    end
    @test isempty(filter(r -> r.level == Warn, logs))

    # A dead input (zero floor) is not filling and must still warn.
    dead = array_state(4)
    dead_logs, _ = collect_test_logs() do
        track!((L1 = BandMeasurement(zeros(ComplexF32, 4000, 4), FS),), dead)
    end
    @test !Tracking._noise_window_filling(dead.noise_estimators.GPSL1CA)
    @test !isempty(filter(r -> r.level == Warn, dead_logs))

    # A source that reports no look count is taken at its word rather than gated.
    @test noise_window_looks(UncountedNoiseSource()) === nothing
    @test Tracking._sufficient_looks(zero(SMatrix{4,4,ComplexF64,16}), nothing)
    @test !Tracking._sufficient_looks(zero(SMatrix{4,4,ComplexF64,16}), 3)
    @test Tracking._sufficient_looks(zero(SMatrix{4,4,ComplexF64,16}), 4)
    @test Tracking._sufficient_looks(1.0 / 1.0Hz, 1)      # a scalar has no rank
end

@testset "the default filter's weights reduce the covariance to its own column" begin
    # Through `DefaultPostCorrFilter`'s weights the covariance is the last
    # antenna's floor, so array and single-antenna runs match exactly.
    scales = (1.0, 1.5, 2.0)
    signal = _sky_array(45.0, scales; seed = 12)
    ts_array = _tracked_array(signal, 3)
    R = get_noise_density(ts_array.noise_estimators.GPSL1CA)
    w = get_weights(DefaultPostCorrFilter(), NumAnts(3))
    @test Tracking._reduce_noise_density(R, w) === real(R[3, 3])

    ts_single = _tracked(view(signal, :, 3), CPUDownconvertAndCorrelator())
    scalar = get_noise_density(ts_single.noise_estimators.GPSL1CA)
    # Same samples and seeded draws, so equality, not a tolerance.
    @test Tracking._reduce_noise_density(R, w) === scalar
    @test estimate_cn0(ts_array, 1) === estimate_cn0(ts_single, 1)
end

# Gated to Julia >= 1.11 as the scalar version above. Pins that the covariance stays
# an `SMatrix` (isbits `NoiseWindowTotals`); a `Matrix` would pass every other test.
if VERSION >= v"1.11"
    @testset "an array's software measurement is allocation-free in steady state" begin
        function measure(dc, measurements, ts)
            for _ = 1:8
                downconvert_and_correlate!(dc, measurements, ts; chunk_index = 0)
            end
            @allocated downconvert_and_correlate!(dc, measurements, ts; chunk_index = 0)
        end
        signal = ComplexF32.(_sky_array(45.0, (1.0, 1.5, 2.0); seed = 13))
        measurements = (L1 = BandMeasurement(signal, FS),)
        state(cn0) = TrackState(
            GPSL1,
            [
                TrackedSat(
                    GPSL1,
                    1,
                    0.0,
                    0.0Hz;
                    num_ants = NumAnts(3),
                    cn0_estimator = cn0(),
                ),
            ],
        )
        ts = state(NoiseRefCN0Estimator)
        plain_ts = state(() -> NWPRCN0Estimator(GPSL1))
        @test keys(ts.noise_estimators) == (:GPSL1CA,)
        dc = CPUDownconvertAndCorrelator()
        plain_dc = CPUDownconvertAndCorrelator()
        track!(measurements, ts; downconvert_and_correlator = dc)
        track!(measurements, plain_ts; downconvert_and_correlator = plain_dc)
        @test measure(dc, measurements, ts) == 0
        @test measure(plain_dc, measurements, plain_ts) == 0
    end
end

end
