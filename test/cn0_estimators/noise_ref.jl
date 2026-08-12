module NoiseRefCN0EstimatorTest

using Test: @test, @testset, @inferred, @test_throws, @test_logs, collect_test_logs
using Logging: Warn
using Random: Xoshiro, randn
using Statistics: mean, std
using Dictionaries: dictionary
using StaticArrays: SVector
using Unitful: Hz, s, ms, dBHz, ustrip, uconvert
using GNSSSignals: GPSL1CA, GPSL5I
import Tracking
using Tracking:
    BitBuffer,
    CN0UpdateContext,
    CorrelatorNoiseEstimator,
    CorrelatorOutput,
    EarlyPromptLateCorrelator,
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    NoCN0Estimator,
    NoiseRefCN0Estimator,
    TrackState,
    TrackedSat,
    append_correlator_output!,
    append_noise_observation!,
    estimate_cn0,
    estimate_dopplers_and_filter_prompt!,
    get_cn0_estimator,
    noise_observation_from_samples,
    requires_noise_density,
    update,
    update_accumulator

const T = 1ms
# `E|P|² = N₀/T`, so a unit-noise-power prompt stream is measured against
# `N₀ = T` — 1e-3 Hz⁻¹ at a 1 ms record.
const N₀ = uconvert(Hz^-1, T)

_context(; noise_density = N₀, integration_time = T) =
    CN0UpdateContext(GPSL1CA(), BitBuffer{UInt64}(), 1; noise_density, integration_time)

# The same post-correlation prompt model `test/cn0_estimator_comparison.jl` uses:
# `λ = (C/N₀)·T` is the per-record SNR, noise is CN(0,1), phase is perfect.
_prompts(λ, n, rng) = sqrt(λ) .+ (randn(rng, n) .+ im .* randn(rng, n)) ./ sqrt(2)

# `ustrip` does not take a `Level` unit directly, so convert first — the idiom
# the other CN0 tests use.
_db(x) = ustrip(uconvert(dBHz, x))

_fold(estimator, prompts) =
    foldl((e, p) -> update(e, p, _context()), prompts; init = estimator)

@testset "one record plus a density is already an estimate" begin
    estimator = NoiseRefCN0Estimator(; num_records = 4)
    @test @inferred(Base.length(estimator)) == 0
    # Empty ring: `-Inf dB-Hz`, the house convention for "not measured".
    @test @inferred(estimate_cn0(estimator, T)) == -Inf * dBHz

    # No warm-up, no window, no `M`: the very first record reports. A noiseless
    # prompt of power `p` against `N₀ = T` reads `(p − 1)/T`.
    one = @inferred update(estimator, complex(2.0, 0.0), _context())
    @test Base.length(one) == 1
    @test _db(estimate_cn0(one, T)) ≈ 10log10((4.0 - 1.0) / 1e-3)

    # `estimate_cn0`'s `integration_time` argument is ignored — `T` was applied
    # per record, where each record's own value was known.
    @test estimate_cn0(one, 20ms) == estimate_cn0(one, T)

    # The ring wraps at `num_records` and keeps only the newest.
    filled = _fold(NoiseRefCN0Estimator(; num_records = 4), fill(complex(2.0, 0.0), 10))
    @test Base.length(filled) == 4
    @test _db(estimate_cn0(filled, T)) ≈ 10log10((4.0 - 1.0) / 1e-3)

    @test requires_noise_density(estimator)
    @test !requires_noise_density(NWPRCN0Estimator())
    @test !requires_noise_density(MomentsCN0Estimator(10))
    @test !requires_noise_density(NoCN0Estimator())
end

@testset "records of different length are each divided by their own T" begin
    # `T` enters per record rather than once at `estimate_cn0`, which is what
    # lets one per-band density serve records of any length — and retires the
    # heterogeneous-`T` bug class structurally instead of guarding it. Against a
    # fixed `N₀` the sample-normalised prompt power of a record at C/N₀ = γ is
    # `N₀·(γ + 1/T)`, so a 20 ms record carrying the *same* γ has a visibly
    # different power: its share of the noise floor has shrunk twentyfold.
    γ = 3000.0                                   # ≈34.8 dB-Hz
    power(t_ms) = 1e-3 * (γ + 1 / (t_ms * 1e-3))
    short = update(NoiseRefCN0Estimator(), complex(sqrt(power(1)), 0.0), _context())
    long = update(
        NoiseRefCN0Estimator(),
        complex(sqrt(power(20)), 0.0),
        _context(; integration_time = 20ms),
    )
    @test _db(estimate_cn0(short, T)) ≈ 10log10(γ) atol = 1e-9
    @test _db(estimate_cn0(long, 20ms)) ≈ 10log10(γ) atol = 1e-9

    # ... and a ring holding both together still averages to the same γ, which is
    # the case `estimate_cn0(estimator, integration_time)` cannot express at all.
    both = update(short, complex(sqrt(power(20)), 0.0), _context(; integration_time = 20ms))
    @test _db(estimate_cn0(both, T)) ≈ 10log10(γ) atol = 1e-9
end

@testset "per-record terms are not clamped, so the mean is unbiased" begin
    # Individual terms go negative on noise, and that is exactly what keeps the
    # mean honest. Clamping per record would reintroduce a floor of the kind the
    # moment ratio has.
    rng = Xoshiro(20260806)
    noise_only = _fold(NoiseRefCN0Estimator(; num_records = 2000), _prompts(0.0, 2000, rng))
    terms = noise_only.buffered_cn0
    @test count(<(0), terms) > 500          # ≈63 % of a unit-mean exponential
    # ... and they cancel: the mean sits within a few σ of zero, where one σ of
    # the mean is `1/(√2000 · T)` ≈ 22 Hz.
    @test abs(mean(terms)) < 150
    # A whole ring of pure noise never reports a positive floor.
    @test _db(estimate_cn0(noise_only, T)) < 25
end

@testset "bias and σ match the pinned noise-reference columns" begin
    # `test/cn0_estimator_comparison.jl` pins what a noise-reference estimator
    # must reproduce; this is the shipped estimator measured the same way. Those
    # reference columns are the variance-free-reference limit, which is what an
    # exactly known `N₀` gives here.
    num_records = 100
    trials = 1500
    for (cn0_dbhz, σ_bound) in ((30.0, 0.9), (40.0, 0.30), (50.0, 0.10))
        λ = 10^(cn0_dbhz / 10) * 1e-3
        rng = Xoshiro(20260806 + round(Int, cn0_dbhz))
        estimates = map(1:trials) do _
            estimator = _fold(
                NoiseRefCN0Estimator(; num_records),
                _prompts(λ, num_records, rng),
            )
            # Back to λ̂, so the numbers are on the plan's scale.
            10^(_db(estimate_cn0(estimator, T)) / 10) * 1e-3
        end
        # Never degenerate — the first structural claim the plan rests on.
        @test all(isfinite, estimates)
        @test all(>(0), estimates)
        @test abs(10log10(mean(estimates) / λ)) < 0.05
        @test 4.342944819 * std(estimates) / λ < σ_bound
    end
end

@testset "no source configured is a loud, static error" begin
    # `Nothing` in the context's type parameter, so this is decided at compile
    # time and fires on the very first record. A silent substitution would hide a
    # backend that never populates the reference.
    estimator = NoiseRefCN0Estimator()
    @test_throws ArgumentError update(
        estimator,
        complex(1.0, 0.0),
        CN0UpdateContext(GPSL1CA(), BitBuffer{UInt64}(), 1),
    )
    # And there is no bare-prompt form at all: the density and `T` are exactly
    # what a prompt stream cannot carry.
    @test_throws ArgumentError update(estimator, complex(1.0, 0.0))
end

@testset "update is allocation-free and inferred" begin
    function fold_many(estimator, prompt, context, n)
        for _ = 1:n
            estimator = update(estimator, prompt, context)
        end
        estimator
    end
    estimator = NoiseRefCN0Estimator()
    context = _context()
    fold_many(estimator, complex(2.0, 0.0), context, 200)
    @test @allocated(fold_many(estimator, complex(2.0, 0.0), context, 100_000)) == 0
    @test @inferred(update(estimator, complex(2.0, 0.0), context)) isa NoiseRefCN0Estimator
    @test @inferred(estimate_cn0(estimator, T)) isa typeof(0.0dBHz)
end

# A pure **correlator-ingest** state: records are handed over with
# `append_correlator_output!` and no sample buffer is ever passed, so the software
# noise measurement cannot run. That is exactly the hardware/FPGA path, and it is
# the only way to reach the "configured but never fed" condition — on the
# sample-driven path `downconvert_and_correlate!` fills the window before the fold
# reads it.
#
# The records are noiseless: a raw accumulator of `a = √power · N` normalises to a
# prompt of `√power` over `N` samples, so the C/N₀ the fold reports follows exactly
# from the density that is (or is not) appended.
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

@testset "a configured-but-never-fed band warns once and stays at -Inf" begin
    # The likeliest hardware integration mistake: a `CorrelatorNoiseEstimator` is
    # configured per the docs and `append_noise_observation!` is never called. The
    # static check does not fire (a source exists), so without this the runtime
    # skip would repeat in silence forever, leaving every satellite at
    # `-Inf dB-Hz` with no clue why.
    ts, fs = _ingested_state(NoiseRefCN0Estimator())
    @test keys(ts.noise_estimators) == (:L1,)
    @test_logs (:warn,) estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
    @test estimate_cn0(ts, 1) == -Inf * dBHz
    @test Base.length(get_cn0_estimator(ts, 1)) == 0
end

@testset "two misconfigured bands both warn" begin
    # `maxlog` is keyed per callsite, so a bare `maxlog = 1` would let the first
    # band to reach the skip silence the second — and since the message
    # interpolates the band id, the second band's name would never appear. The
    # per-band `_id` is what prevents that, so pin the ids rather than the count
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
    @test Set(keys(multi.noise_estimators)) == Set((:L1, :L5))
    logs, _ = collect_test_logs() do
        estimate_dopplers_and_filter_prompt!(multi, (L1 = 4e6Hz, L5 = 4e6Hz))
    end
    warnings = filter(r -> r.level == Warn, logs)
    @test length(warnings) == 2
    @test Set(r.id for r in warnings) == Set((:no_noise_density_L1, :no_noise_density_L5))
    @test all(r -> occursin("append_noise_observation!", r.message), warnings)
end

@testset "the warm-up skip is per estimator, not per band" begin
    # A band carrying one `NoiseRefCN0Estimator` signal and one `NWPRCN0Estimator`
    # signal, driven through a warm-up where the density is unavailable. Skipping
    # the whole band's C/N₀ fold would drop NWPR's open narrowband window on every
    # skipped record — `_update_nwpr` drops it whenever the record is missing from
    # the bit grid — and silently demote NWPR to its fallback.
    # Two L1 C/A signals on one satellite — synthetic, but it is the shortest
    # configuration that puts two *different* C/N₀ estimators on one band, and it
    # covers both folds: the noise-referenced one drives the loop, the NWPR one
    # rides along as a passenger. Records are ingested rather than correlated, so
    # the density really is unavailable for the whole run.
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
    # The reference: the same NWPR signal on a band with no noise estimator at
    # all, so no skip can apply anywhere.
    reference = run((NoCN0Estimator(), NWPRCN0Estimator()))
    @test keys(mixed.noise_estimators) == (:L1,)
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
        :L1,
    )
    estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
    @test Base.length(get_cn0_estimator(ts, 1)) == 20
    @test _db(estimate_cn0(ts, 1)) ≈ 10log10(1 / 1e-6 - 1 / 1e-3) atol = 1e-6
end

@testset "a band nobody asks a density of is not provisioned" begin
    gpsl1, gpsl5 = GPSL1CA(), GPSL5I()
    # Every shipped estimator but the noise-referenced one reads no density, so a
    # state that stays on one of them runs no per-band despread at all.
    for estimator in (NWPRCN0Estimator(), MomentsCN0Estimator(100), NoCN0Estimator())
        ts =
            TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = estimator)])
        @test ts.noise_estimators === NamedTuple()
    end
    # The declared-signals constructor asks `default_cn0_estimator`, which is
    # still NWPR, so nothing is provisioned there either.
    @test TrackState(; signal = gpsl1).noise_estimators === NamedTuple()

    # Mixed: one band requires a density, the other does not — exactly one entry.
    mixed = TrackState((
        l1 = dictionary((
            1 => TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator()),
        )),
        l5 = dictionary((
            1 => TrackedSat(gpsl5, 1, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator()),
        )),
    ))
    @test keys(mixed.noise_estimators) == (:L1,)
end

end
