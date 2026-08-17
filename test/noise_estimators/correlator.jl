module CorrelatorNoiseEstimatorTest

using Test: @test, @testset, @test_logs, collect_test_logs
using Logging: Warn
using Random: Xoshiro, randn
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
    EarlyPromptLateCorrelator,
    Int16DownconvertAndCorrelator,
    NWPRCN0Estimator,
    NoiseRefCN0Estimator,
    OneBitDownconvertAndCorrelator,
    TrackState,
    TrackedSat,
    TwoBitDownconvertAndCorrelator,
    add_satellite!,
    downconvert_and_correlate!,
    get_correlator_sample_shifts,
    estimate_cn0,
    get_noise_density,
    track!

const FS = 4e6Hz
const NUM_SAMPLES = 40_000                      # 10 ms at 4 MHz
const GPSL1 = GPSL1CA()

_db(x) = ustrip(uconvert(dBHz, x))
_median(xs) = (v = sort(collect(xs)); v[div(length(v) + 1, 2)])

# One satellite at `cn0_dbhz` in noise of per-sample variance `f_s`, so the true
# noise density is exactly `N₀ = σ²/f_s = 1 Hz⁻¹` and the measured density can be
# compared against a number rather than against another measurement.
function _sky(cn0_dbhz; seed = 1, num_samples = NUM_SAMPLES, fs = FS, prn = 1)
    rng = Xoshiro(seed)
    code = gen_code(num_samples, GPSL1, prn, fs, get_code_frequency(GPSL1), 0.0)
    amplitude = isfinite(cn0_dbhz) ? 10^(cn0_dbhz / 20) : 0.0
    ComplexF64.(amplitude .* code) .+
    sqrt(ustrip(Hz, fs)) .* randn(rng, ComplexF64, num_samples)
end

_to_int16(signal) = Complex{Int16}.(
    round.(Int16, clamp.(real(signal) ./ 4, -2047, 2047)),
    round.(Int16, clamp.(imag(signal) ./ 4, -2047, 2047)),
)

# 20 calls, not 5. The window's own relative error is `1/√(3K)` over `K`
# observations, so a 50-entry window carries σ ≈ 0.6 dB — under half the `atol`
# the C/N₀ assertions below use. That was invisible while the reference despread
# at a fixed phase (one deterministic draw, which happened to land well); with the
# phase and carrier drawn per sub-integration the assertions have to be backed by
# a window long enough to support them. At `K = 200` σ is 0.22 dB, so a 1 dB gate
# is ≈4.5σ and the per-backend quantisation bands stay meaningful.
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
    # An empty sky measured against a known floor. The window holds
    # `num_calls · chunks · num_sub` observations of three taps each, so its own
    # relative error is ≈1/√(3·200) ≈ 4 %.
    ts = _tracked(_sky(-Inf; seed = 7), CPUDownconvertAndCorrelator())
    estimator = ts.noise_estimators.GPSL1CA
    @test Base.length(estimator) == 200
    @test ustrip(Hz^-1, get_noise_density(estimator)) ≈ 1.0 rtol = 0.15
end

@testset "the software path needs no `append_noise_observation!` and never warns" begin
    # The counterpart of the configured-but-never-fed regression in
    # `test/cn0_estimators/noise_ref.jl`: on the sample-driven path
    # `downconvert_and_correlate!` measures before the fold reads, so the window
    # is never empty at fold time and the diagnostic must stay quiet.
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
    # The reference has to traverse the *same* kernel as the prompt, and on the
    # bit-wise backends that is not a preference: their accumulators are popcount
    # counts rather than sample sums, so a float-kernel reference would put `|P|²`
    # and `N̂₀` on incompatible scales and the reported C/N₀ would be meaningless
    # rather than merely optimistic.
    #
    # What the numbers then say is the model-freeness itself: each backend reports
    # the C/N₀ *its own quantiser delivers*, with no per-backend correction
    # anywhere in Tracking. 1-bit costs ≈2 dB (plus ≈1 dB for the 1-bit carrier)
    # and 2-bit ≈1 dB, which is exactly what shows up.
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
    # `num_sub = max(1, round(chunk_duration / code_period))`, `C · num_sub`
    # observations per `track!` over `C` chunks, and **none** from the final drain
    # pass. Getting the drain gate wrong is the failure this pins: gating on
    # `samples_unchanged` alone would measure chunk 0 and then freeze the window
    # for the rest of the buffer.
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

    # A second sampling frequency, same chunk time: the count follows the chunk's
    # *duration* against the code period, not its sample count.
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    fast = _sky(-Inf; seed = 12, num_samples = 100_000, fs = 10e6Hz)
    track!((L1 = BandMeasurement(fast, 10e6Hz),), ts; doppler_update_interval = 4ms)
    # 10 ms of buffer at a 4 ms chunk: two whole chunks of four slices, then a
    # 2 ms remainder clamped to the buffer end, which yields two — the slice count
    # follows each chunk's own duration against the code period, so a short
    # trailing chunk is measured at its true length rather than padded.
    @test Base.length(ts.noise_estimators.GPSL1CA) == 4 + 4 + 2
end

@testset "a code period longer than the chunk still yields one observation" begin
    # `max(1, …)` is what stops a signal with a 10 ms code sampled with 1 ms
    # chunks from yielding zero observations forever. GPS L1C-P is that case, and
    # it is reachable whenever another signal forces a shorter
    # `doppler_update_interval`.
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
    # Ten observations of three taps is 30 looks, so the window's own relative σ is
    # `1/√30 ≈ 18 %` — this is a sanity check that the measurement is meaningful at
    # all, not a precision claim, and the tolerance says so. The precision claims
    # live in the tests above, on windows two orders of magnitude longer.
    @test ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1C_P)) ≈ 1.0 rtol = 0.55
end

@testset "the PRN rotates through the whole family, tracked codes included" begin
    # Advancing a PRN needs a position, and a position is scalar state with no
    # per-signal anchor — so it is carried in the window itself: each observation
    # records the code it was measured with, and the next one steps on from there.
    #
    # Tracked PRNs are deliberately **not** skipped. Skipping them was what made a
    # fixed code phase 0 safe; the phase is now drawn at random instead, which
    # both removes the need and removes the reference's last read of satellite
    # state. See the two tests below for the properties that replace the skip.
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
    # The draws have to actually reach the kernel, and they have to be
    # reproducible when the caller supplies a generator — the paired self-leakage
    # comparison in `cn0_estimators/noise_ref_in_the_loop.jl` depends on the
    # second half of that.
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
    # The property that replaces the untracked-PRN skip, and the reason the phase
    # is randomised at all. `_sky` puts its satellite at code phase 0 and zero
    # Doppler — exactly where a reference pinned to phase 0 would despread it. A
    # fixed phase therefore does not merely risk a hit, it lands one on *every*
    # observation measured with that PRN, worth `T·(C/N₀)/3 ≈ 10.5·N₀` at
    # 45 dB-Hz, and — because the replica is re-anchored on a nominal-rate grid —
    # it never drifts off again.
    #
    # Nothing tracks PRN 7 here, so the old skip would not have caught it either:
    # this is the untracked spoofed-or-unacquired signal, not the self-leakage
    # case.
    m = BandMeasurement(_sky(45.0; seed = 3, prn = 7), FS)
    dc = CPUDownconvertAndCorrelator()
    # `carrier_dither = 0` isolates the code-phase draw. With the dither on, a hit
    # needs both draws to land and the test would pass for the weaker reason.
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
    # Medians, not means: the randomisation converts the standing bias into a
    # ≈0.6 %-per-observation outlier, and pinning the median is what says the
    # *bias* is gone rather than that this seed drew no outlier at all. A fixed
    # phase 0 would put this ratio at ≈11.
    @test _median(on_seven) / _median(others) ≈ 1.0 rtol = 0.6
end

@testset "the reference's taps are spaced wide enough to be independent" begin
    # Taps are only worth having if they are statistically independent, and at the
    # *tracking* default of ±0.5 chip they are not: the sampled-code correlation
    # there is ≈0.5, so three taps are worth 2.25 looks instead of 3. This is a
    # regression on the rounding in `get_correlator_sample_shifts` at the band's
    # `f_s`, not on the code statistics — 1.5 chips must survive as ≥1 chip of
    # whole samples.
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
    # Same duration, different rates — and a four-times louder L5 front end, so
    # the two densities cannot be confused for one another.
    rng = Xoshiro(21)
    l1_samples = ComplexF32.(sqrt(ustrip(Hz, fs_l1)) .* randn(rng, ComplexF64, 40_000))
    l5_samples = ComplexF32.(2 * sqrt(ustrip(Hz, fs_l5)) .* randn(rng, ComplexF64, 200_000))
    measurements =
        (L1 = BandMeasurement(l1_samples, fs_l1), L5 = BandMeasurement(l5_samples, fs_l5))
    # 20 calls for the same reason `_tracked` uses 20 — see the note there.
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
    # This is the whole reason the estimator is keyed by signal. What a record
    # divides by is the *post-correlation* floor,
    # `N₀,eff = N₀ + ∫ S_I(f)·|G(f)|² df`, weighted by the despreading
    # modulation's own spectrum — so two signals sharing one band, one antenna and
    # one set of samples do not share a noise floor once the interference is
    # coloured. A per-band figure can only ever be right for one of them.
    #
    # BPSK(1) nulls at ±f_chip where BOC(1,1) peaks, and peaks inside the main
    # lobe where BOC(1,1) is weak, so a CW tone swept across the band moves the
    # two in *opposite* directions. That reversal is what makes this a spectral
    # result and not a scale artifact.
    #
    # `GalileoE1B_BOC11` rather than the CBOC `GalileoE1B` only because CBOC's
    # subchip factor needs f_s ≥ 12.276 MHz; the real E1B shows the same effect
    # about four times larger.
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

    # White noise: nothing to separate, and the two must agree — the property the
    # per-band design relied on and the only regime in which it held.
    white = densities(nothing)
    @test white.l1ca ≈ 1.0 rtol = 0.2
    @test white.e1b ≈ 1.0 rtol = 0.2
    @test white.e1b / white.l1ca ≈ 1.0 rtol = 0.25

    # A tone at the chip rate: BPSK(1)'s first null, BOC(1,1)'s peak. C/A barely
    # notices it (≈1.4× thermal); E1B's floor rises ≈39×, i.e. ≈28× C/A's.
    # Reporting C/A's figure to E1B here would overstate its C/N₀ by ≈14.5 dB.
    at_chip_rate = densities(1.023e6)
    @test at_chip_rate.l1ca < 2 * white.l1ca
    @test at_chip_rate.e1b > 10 * at_chip_rate.l1ca

    # ... and inside C/A's main lobe the order reverses, which no common scale
    # factor can produce.
    in_main_lobe = densities(0.4e6)
    @test in_main_lobe.l1ca > 10 * white.l1ca
    @test in_main_lobe.l1ca > 1.15 * in_main_lobe.e1b
end

@testset "a signal carried by two groups is measured once" begin
    # `_dc_one_group!` was rejected as the measurement site precisely because it
    # fires once per group: two groups may carry the same signal, and a per-group
    # site would double its window. Walking the `noise_estimators` NamedTuple —
    # one entry per signal — makes once-per-signal structural instead.
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

@testset "an unchunked call measures once, and a re-pass does not" begin
    # The gate is `chunk_duration === nothing && samples_unchanged`, and it is
    # reachable from a direct caller as well as from `track!`'s drain pass — which
    # is why `downconvert_and_correlate!`'s docstring says so.
    signal = _sky(-Inf; seed = 6)
    measurements = (L1 = BandMeasurement(signal, FS),)
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    dc = CPUDownconvertAndCorrelator()
    # Unchunked, samples not promised unchanged: the whole buffer in `num_sub`
    # slices — 10 ms against a 1 ms code period.
    downconvert_and_correlate!(dc, measurements, ts)
    @test Base.length(ts.noise_estimators.GPSL1CA) == 10
    # The same call promising the samples are unchanged is a re-pass, and must not
    # re-measure what has already been covered.
    downconvert_and_correlate!(dc, measurements, ts; samples_unchanged = true)
    @test Base.length(ts.noise_estimators.GPSL1CA) == 10
end

@testset "a signal without a consumer is never despread" begin
    # The whole reason provisioning is gated on `requires_noise_density`: someone
    # who stays on NWPR must not pay for a despread they never read.
    ts = TrackState(
        GPSL1,
        [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator())],
    )
    @test ts.noise_estimators === NamedTuple()
    track!((L1 = BandMeasurement(_sky(45.0), FS),), ts)
    @test ts.noise_estimators === NamedTuple()
end

# Gated to Julia >= 1.11, following the precedent in `test/track_in_place.jl`.
# The despread is allocation-free on 1.11+, but on 1.10 the compiler
# leaves ~670 B per sub-integration inside the correlate call — measured with the
# noise reference switched off it is 0 B on both versions, so it is the
# measurement path and not the walk that provisions it. 1.10 users get a working
# noise reference that allocates ~5 kB per `downconvert_and_correlate!` call;
# 1.11+ get zero. Asserting `== 0` unconditionally would either fail on a
# supported version or have to be weakened into something that no longer catches
# a real regression.
if VERSION >= v"1.11"
    @testset "the software measurement is allocation-free in steady state" begin
        function measure(dc, measurements, ts)
            for _ = 1:8
                downconvert_and_correlate!(dc, measurements, ts; chunk_index = 0)
            end
            @allocated downconvert_and_correlate!(dc, measurements, ts; chunk_index = 0)
        end
        signal = ComplexF32.(_sky(45.0; seed = 8))
        measurements = (L1 = BandMeasurement(signal, FS),)
        ts = TrackState(
            GPSL1,
            [TrackedSat(GPSL1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
        )
        dc = CPUDownconvertAndCorrelator()
        track!(measurements, ts; downconvert_and_correlator = dc)
        @test measure(dc, measurements, ts) == 0
    end
end

@testset "a DC offset only bites on a partial code period" begin
    # A full-period Gold code is balanced (`E[(Σc)²] ≈ 1`), so an ADC DC offset
    # averages out of a code-period despread; any sub-period window behaves like
    # iid signs and the offset adds a *positive*, non-averaging term to every
    # `|B|²`, inflating `N̂₀` and making C/N₀ read low. That is the one and only
    # configuration the risk applies to — a code period longer than the chunk —
    # and it is zero-IF only, since the reference downconverts.
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
    # `append_noise_observation!` and `get_noise_density` read cached sums rather
    # than walking the window, so the cache has to be provably equal to the walk
    # it replaces — including after thousands of appends have each added one entry
    # and subtracted another, and including for a producer that mixes entry sizes
    # (the case incremental floating-point arithmetic could drift on, and the one
    # the periodic exact refresh exists for).
    exact_span(e) = sum(o.duration for o in e.buffered)
    exact_looks(e) = sum(o.num_sub_integrations for o in e.buffered)
    exact_density(e) =
        sum(o.num_sub_integrations * o.noise_density for o in e.buffered) / exact_looks(e)

    rng = Xoshiro(17)
    window = 50.0ms
    estimator = CorrelatorNoiseEstimator(; window_duration = window)
    for i = 1:5000
        # 1-in-20 entries is a long pre-averaged dump among short ones, so the
        # running sums are not adding and removing like-sized numbers.
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

    # …and the window is still bounded in time: minimal, but never short of the
    # configured span.
    @test estimator.totals[].span >= window
    @test estimator.totals[].span - first(estimator.buffered).duration < window
end

@testset "reading the window does not walk it" begin
    # `get_noise_density` runs once per signal per chunk. It used to sum the whole
    # window on every call, which charged each chunk for the window length — i.e.
    # for exactly the `K` the estimator's accuracy is bought with, so configuring
    # a longer window made tracking slower. Timing is too flaky to assert on a
    # shared runner; that the cost does not grow with `K` is the property, so
    # assert the walk is gone by its allocation-free, K-independent shape instead.
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
    # The estimator holds no replica buffer at all: `_despread_one_signal!` draws
    # one from the backend's own scratch on the kernels that generate a replica and
    # ignores the size on the ones that pack the code sign plane themselves. That
    # removes a byte per sample of the whole measurement, per signal, on *every*
    # backend rather than on the three that never read it — half a megabyte per
    # signal on a 500 k-sample buffer.
    @test :code_replica ∉ fieldnames(CorrelatorNoiseEstimator)
    signal = _sky(45.0; seed = 5)
    ts = _tracked(quantise ? _to_int16(signal) : ComplexF32.(signal), build())
    estimator = ts.noise_estimators.GPSL1CA
    @test get_noise_density(estimator) !== nothing      # it still measured
end

@testset "the reference and the satellites despread out of one scratch slot" begin
    # The other half of that: the buffer the software reference writes its replica
    # into is the backend's shared per-thread slot — the same one the per-satellite
    # path uses — so there is exactly one, however many signals are measured. Safe
    # because the noise pass runs to completion before the per-group satellite loop
    # that reuses the slot, and the threaded backends index it by `threadid()`.
    dc = CPUDownconvertAndCorrelator()
    @test isempty(Tracking._scratch_buffers(dc).code_replica)
    ts = _tracked(ComplexF32.(_sky(45.0; seed = 5)), dc)
    @test get_noise_density(ts.noise_estimators.GPSL1CA) !== nothing
    # Sized against the whole buffer plus the tap spread, not the slice:
    # `gen_code_replica!` writes *at* `start_sample`, and a slice can start
    # anywhere in the measurement.
    @test length(Tracking._scratch_buffers(dc).code_replica) >= NUM_SAMPLES
end

end
