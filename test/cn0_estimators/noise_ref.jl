module NoiseRefCN0EstimatorTest

using Test: @test, @testset, @test_logs, collect_test_logs
using Logging: Warn
using Dictionaries: dictionary
using StaticArrays: SVector
using Unitful: Hz, s, dBHz, ustrip, uconvert
using GNSSSignals: GPSL1CA, GPSL5I
import Tracking
using Tracking:
    TrackState,
    TrackedSat,
    append_correlator_output!,
    estimate_dopplers_and_filter_prompt!,
    get_cn0_estimator
using TrackingLoops:
    CorrelatorNoiseEstimator,
    CorrelatorOutput,
    EarlyPromptLateCorrelator,
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    NoCN0Estimator,
    NoiseRefCN0Estimator,
    append_noise_observation!,
    estimate_cn0,
    noise_observation_from_samples,
    update_accumulator

# `ustrip` does not take a `Level` unit directly, so convert first — the idiom
# the other CN0 tests use.
_db(x) = ustrip(uconvert(dBHz, x))

# A pure **correlator-ingest** state (the hardware path): records come through
# `append_correlator_output!` with no samples, so the software noise measurement
# never runs — the only way to reach "configured but never fed".
#
# The records are noiseless: a raw accumulator `a = √power · N` normalises to a
# prompt of `√power`, so the reported C/N₀ follows exactly from the appended density.
function _ingested_state(
    cn0_estimator;
    fs = 4e6Hz,
    num_samples = 4000,
    num_records = 20,
    prompt_power = 1.0,
    signals = nothing,
)
    gpsl1 = GPSL1CA()
    sat =
        isnothing(signals) ? TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator) :
        TrackedSat(signals, 1, 0.0, 0.0Hz; cn0_estimator)
    ts = TrackState(gpsl1, [sat])
    a = sqrt(prompt_power) * num_samples
    correlator = update_accumulator(
        EarlyPromptLateCorrelator(),
        SVector(complex(0.5a, 0.0), complex(a, 0.0), complex(0.5a, 0.0)),
    )
    num_sigs = isnothing(signals) ? 1 : length(signals)
    for k = 1:num_records, i = 1:num_sigs
        output = CorrelatorOutput(correlator, num_samples, k * num_samples)
        if num_sigs == 1
            append_correlator_output!(ts, output, 1)
        else
            append_correlator_output!(ts, output, :default, 1, i)
        end
    end
    ts, fs
end

@testset "a configured-but-never-fed signal warns once and stays at -Inf" begin
    # The likeliest hardware integration mistake: a `CorrelatorNoiseEstimator` is
    # configured per the docs and `append_noise_observation!` is never called. The
    # static check does not fire (a source exists), so without this the runtime
    # skip would repeat in silence forever, leaving every satellite at
    # `-Inf dB-Hz` with no clue why.
    ts, fs = _ingested_state(NoiseRefCN0Estimator())
    @test keys(ts.noise_estimators) == (:GPSL1CA,)
    @test_logs (:warn,) estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
    @test estimate_cn0(ts, 1) == -Inf * dBHz
    @test Base.length(get_cn0_estimator(ts, 1)) == 0
end

@testset "two misconfigured signals both warn" begin
    # `maxlog` is keyed per callsite, so a bare `maxlog = 1` would let the first
    # signal to reach the skip silence the second — and since the message
    # interpolates the signal id, the second signal's name would never appear. The
    # per-signal `_id` is what prevents that, so pin the ids rather than the count
    # of messages the real logger would let through.
    gpsl1, gpsl5 = GPSL1CA(), GPSL5I()
    multi = TrackState((
        l1 = dictionary((
            1 => TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator()),
        )),
        l5 = dictionary((
            1 => TrackedSat(gpsl5, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator()),
        )),
    ))
    @test Set(keys(multi.noise_estimators)) == Set((:GPSL1CA, :GPSL5I))
    logs, _ = collect_test_logs() do
        estimate_dopplers_and_filter_prompt!(multi, (L1 = 4e6Hz, L5 = 4e6Hz))
    end
    warnings = filter(r -> r.level == Warn, logs)
    @test length(warnings) == 2
    @test Set(r.id for r in warnings) ==
          Set((:no_noise_density_GPSL1CA, :no_noise_density_GPSL5I))
    @test all(r -> occursin("append_noise_observation!", r.message), warnings)
end

@testset "the warm-up skip is per estimator, not per satellite" begin
    # See "the warm-up skip" in `requires_noise_density`. Two L1 C/A signals on
    # one satellite (synthetic, but the shortest setup with two different
    # estimators side by side); records are ingested, so the density stays
    # unavailable for the whole run.
    gpsl1 = GPSL1CA()
    signals = (gpsl1, gpsl1)

    function run(cn0_estimator)
        ts, fs = _ingested_state(cn0_estimator; signals)
        estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
        ts
    end

    # The mixed band: the noise-referenced signal is skipped throughout (nothing
    # ever fills its window), the NWPR one must not be.
    mixed = run((NoiseRefCN0Estimator(), NWPRCN0Estimator()))
    # The reference: the same NWPR signal on a satellite with no noise estimator
    # at all, so no skip can apply anywhere.
    reference = run((NoCN0Estimator(), NWPRCN0Estimator()))
    @test keys(mixed.noise_estimators) == (:GPSL1CA,)
    @test reference.noise_estimators === NamedTuple()

    skipped = get_cn0_estimator(mixed, :default, 1, 2)
    unskipped = get_cn0_estimator(reference, :default, 1, 2)
    @test skipped.num_records_per_ratio == unskipped.num_records_per_ratio
    @test skipped.ratios_are_bit_aligned == unskipped.ratios_are_bit_aligned
    @test Base.length(skipped) == Base.length(unskipped)
    @test Base.length(skipped) > 0
    @test estimate_cn0(mixed, :default, 1, 2) == estimate_cn0(reference, :default, 1, 2)
    # ... and the noise-referenced signal really did stay empty, so the skip was
    # exercised rather than bypassed.
    @test Base.length(get_cn0_estimator(mixed, :default, 1, 1)) == 0
end

@testset "the externally fed path reports a real C/N₀" begin
    # End to end on the hardware ingest path, with no sample buffer anywhere: a
    # producer that appends correlator outputs *and* its noise observations gets a
    # C/N₀ out. The records are noiseless with a normalised prompt power of 1, so
    # a density of `1e-6 Hz⁻¹` puts the answer at exactly
    # `10log10(1/1e-6 − 1/1e-3)`.
    ts, fs = _ingested_state(NoiseRefCN0Estimator())
    n = 4000
    append_noise_observation!(
        ts,
        noise_observation_from_samples(n * 1e-6 * ustrip(Hz, fs), n, fs),
        :GPSL1CA,
    )
    estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
    @test Base.length(get_cn0_estimator(ts, 1)) == 20
    @test _db(estimate_cn0(ts, 1)) ≈ 10log10(1 / 1e-6 - 1 / 1e-3) atol = 1e-6
end

@testset "a signal nobody asks a density of is not provisioned" begin
    gpsl1, gpsl5 = GPSL1CA(), GPSL5I()
    # Every shipped estimator but the noise-referenced one reads no density, so a
    # state that stays on one of them runs no despread at all.
    for estimator in (NWPRCN0Estimator(), MomentsCN0Estimator(100), NoCN0Estimator())
        ts =
            TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = estimator)])
        @test ts.noise_estimators === NamedTuple()
    end
    # The declared-signals constructor takes its slot type from
    # `default_cn0_estimator`, which *is* noise-referenced — so the default path
    # provisions, and needs no configuration at all.
    @test keys(TrackState(; signal = gpsl1).noise_estimators) == (:GPSL1CA,)

    # Mixed: one signal requires a density, the other does not — exactly one entry.
    mixed = TrackState((
        l1 = dictionary((
            1 => TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator()),
        )),
        l5 = dictionary((
            1 => TrackedSat(gpsl5, 1, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator()),
        )),
    ))
    @test keys(mixed.noise_estimators) == (:GPSL1CA,)
end

end
