module CN0EstimatorComparisonTest

# Baseline for docs/plans/2026-08-06-per-band-noise-estimation-for-cn0.md: pins, by
# Monte Carlo over the post-correlation prompt model, the plan's structural claims —
# NWPR degenerates at low C/N₀, is biased, saturates in relative error, and loses to
# a coherent noise reference at every C/N₀.
#
# `NWPRCN0Estimator` is the shipped implementation (bare-prompt-stream `update`),
# given its best case: perfect phase, no bit transitions. The noise-reference columns
# are reference implementations with the noise power known exactly and no other
# satellites, so none of the cross-correlation / self-leakage bias the plan
# quantifies; "unbiased" below is about the estimator, not a full sky.

using Test: @test, @testset
using Random: Xoshiro, randn
using Statistics: mean, std
using Unitful: ms, s, dBHz, ustrip, uconvert
using Tracking: NWPRCN0Estimator, update, estimate_cn0

# Post-correlation prompts: `λ = (C/N₀)·T` is the per-record SNR; noise is CN(0,1),
# so λ is also the signal power. Perfect phase, signal on the real axis.
_prompts(λ, num_records, rng) =
    sqrt(λ) .+ (randn(rng, num_records) .+ im .* randn(rng, num_records)) ./ sqrt(2)

# Non-coherent noise reference: λ̂ = ⟨|P|²⟩ − 1 (known noise power).
_noise_ref_noncoherent(prompts) = mean(abs2, prompts) - 1

# Coherent noise reference: sum `M` prompts, subtract the sum's known noise power.
# Equivalent to `M`-times-longer records — NWPR's narrowband window is a longer
# coherent integration, not extra information.
function _noise_ref_coherent(prompts, M)
    # A remainder would shorten this column's observation and make it incomparable.
    @assert rem(length(prompts), M) == 0
    num_windows = div(length(prompts), M)
    acc = 0.0
    for w = 1:num_windows
        window = @view prompts[((w-1)*M+1):(w*M)]
        acc += (abs2(sum(window)) - M) / M^2
    end
    acc / num_windows
end

# The shipped estimator over the same prompts, returning λ̂ so all three share a
# scale. The two-argument `update` ignores the pre-sync window length, so it is left
# at its default.
function _nwpr(prompts, num_records, M, integration_time)
    estimator = NWPRCN0Estimator(; num_records, num_narrowband_code_blocks = M)
    for prompt in prompts
        estimator = update(estimator, prompt)
    end
    # λ = (C/N₀)·T. `estimate_cn0` returns dB-Hz, so undo the log and re-apply T.
    cn0_db = ustrip(uconvert(dBHz, estimate_cn0(estimator, integration_time)))
    10^(cn0_db / 10) * ustrip(s, integration_time)
end

# NWPR's pooled ratio outside `1 < µ̂ < M` gives `±Inf dB-Hz` (see `estimate_cn0`
# for `NWPRCN0Estimator`), i.e. λ̂ of exactly zero or Inf; so `degenerate` counts
# whole 100-record estimates, not discarded windows. A negative noise-reference λ̂
# is legitimate (that is what keeps it unbiased), so only zero and non-finite count.
_is_degenerate(λ̂) = !isfinite(λ̂) || iszero(λ̂)

# Relative σ in dB (delta method) and linear-mean bias over the non-degenerate
# estimates; excluding the degenerate ones is the generous reading for NWPR.
function _stats(estimates, λ)
    usable = filter(!_is_degenerate, estimates)
    degenerate = 1 - length(usable) / length(estimates)
    isempty(usable) && return (; σ_db = Inf, bias_db = Inf, degenerate)
    σ_db = 4.342944819 * std(usable) / λ
    bias_db = 10log10(max(mean(usable), eps()) / λ)
    (; σ_db, bias_db, degenerate)
end

const NUM_RECORDS = 100          # 100 ms of observation at 1 ms records
const M = 5                      # NWPR default window for a 1 ms code
const INTEGRATION_TIME = 1ms
const NUM_TRIALS = 4_000         # σ good to ~1 % per seed
const NUM_SEEDS = 9              # medianed over, see `_sweep`
const CN0S_DBHZ = 20:5:50

function _sweep_once(cn0_dbhz, seed)
    λ = 10^(cn0_dbhz / 10) * ustrip(ms, INTEGRATION_TIME) * 1e-3
    rng = Xoshiro(seed)
    nwpr = Vector{Float64}(undef, NUM_TRIALS)
    noncoh = Vector{Float64}(undef, NUM_TRIALS)
    coh = Vector{Float64}(undef, NUM_TRIALS)
    for t = 1:NUM_TRIALS
        prompts = _prompts(λ, NUM_RECORDS, rng)
        nwpr[t] = _nwpr(prompts, NUM_RECORDS, M, INTEGRATION_TIME)
        noncoh[t] = _noise_ref_noncoherent(prompts)
        coh[t] = _noise_ref_coherent(prompts, M)
    end
    (; λ, nwpr = _stats(nwpr, λ), noncoh = _stats(noncoh, λ), coh = _stats(coh, λ))
end

_median(xs) = (v = sort(collect(xs)); v[div(length(v) + 1, 2)])

_median_stats(per_seed, key) = (;
    σ_db = _median(getfield(s, key).σ_db for s in per_seed),
    bias_db = _median(getfield(s, key).bias_db for s in per_seed),
    degenerate = _median(getfield(s, key).degenerate for s in per_seed),
)

# Median over seeds so the thresholds clear the Monte-Carlo error, not just one
# seed: a single seed's `bias_db` has ≈0.014 dB sampling error at 30 dB-Hz against a
# 0.02 dB bound. Over 40 alternative seed bases at `NUM_SEEDS = 9` every assertion
# holds; the tightest are:
#
#   |noncoh.bias| < 0.02 @ 30 dB-Hz   worst 0.0120   (the binding one)
#   |coh.bias|    < 0.02 @ 30 dB-Hz   worst 0.0120
#   nwpr.bias     > 0.02 @ 40 dB-Hz   worst 0.0452
#   nwpr.σ(50) > 0.9·nwpr.σ(45)       worst ratio 0.969
#   nwpr.degenerate @ 20 dB-Hz        worst 0.0553 against a 0.01 bound
#
# (Five seeds left `|noncoh.bias|` at 0.0173.)
function _sweep(cn0_dbhz)
    per_seed = [_sweep_once(cn0_dbhz, 20260806 + cn0_dbhz + 1000k) for k = 1:NUM_SEEDS]
    (;
        λ = first(per_seed).λ,
        nwpr = _median_stats(per_seed, :nwpr),
        noncoh = _median_stats(per_seed, :noncoh),
        coh = _median_stats(per_seed, :coh),
    )
end

# One sweep, reused by every testset below so the Monte Carlo runs once.
const SWEEP = Dict(cn0 => _sweep(cn0) for cn0 in CN0S_DBHZ)

@testset "NWPR reports degenerate estimates at low C/N₀, a noise reference never does" begin
    # ≈5.6 % of whole estimates at 20 dB-Hz (see `_is_degenerate`); the threshold
    # pins the structural fact, not the rate.
    @test SWEEP[20].nwpr.degenerate > 0.01
    for cn0 in CN0S_DBHZ
        @test SWEEP[cn0].noncoh.degenerate == 0
        @test SWEEP[cn0].coh.degenerate == 0
    end
end

@testset "NWPR is biased at every C/N₀, a noise reference is not" begin
    # NWPR's ratio-of-sums bias settles around +0.05 dB (see `NWPRCN0Estimator`);
    # the references know the noise power instead of inferring it (see the header).
    for cn0 = 30:5:50
        @test SWEEP[cn0].nwpr.bias_db > 0.02
        @test abs(SWEEP[cn0].noncoh.bias_db) < 0.02
        @test abs(SWEEP[cn0].coh.bias_db) < 0.02
    end
end

@testset "NWPR's relative error saturates, a noise reference's keeps falling" begin
    # NWPR asymptotes to √(M/((M−1)K)) ≈ 0.49 dB (the `µ̂ → M` ceiling); the
    # references roughly halve between 45 and 50 dB-Hz.
    @test SWEEP[50].nwpr.σ_db > 0.9 * SWEEP[45].nwpr.σ_db
    @test SWEEP[50].coh.σ_db < 0.7 * SWEEP[45].coh.σ_db
    @test SWEEP[50].noncoh.σ_db < 0.7 * SWEEP[45].noncoh.σ_db
end

@testset "A coherent noise reference beats NWPR at every C/N₀" begin
    # The plan's central claim: same coherent gain, noise measured directly.
    for cn0 in CN0S_DBHZ
        @test SWEEP[cn0].coh.σ_db < SWEEP[cn0].nwpr.σ_db
    end
end

@testset "Non-coherent squaring loss is real below ~28 dB-Hz" begin
    # At 1 ms records, dropping NWPR's coherent window costs variance at low C/N₀
    # (ratio √(M−1) = 2) and pays off only above the crossover — why the plan
    # lengthens `preferred_num_code_blocks_to_integrate` instead.
    @test SWEEP[20].noncoh.σ_db > SWEEP[20].nwpr.σ_db
    @test SWEEP[40].noncoh.σ_db < SWEEP[40].nwpr.σ_db
end

# Printed so a regression shows what moved and the plan's numbers stay reproducible.
@testset "baseline table" begin
    println(
        "\nC/N₀ vs estimator, $(NUM_RECORDS) records of $(INTEGRATION_TIME), ",
        "NWPR M=$M, $(NUM_TRIALS) trials × $(NUM_SEEDS) seeds (median)",
    )
    println(
        "dB-Hz |  NWPR σ  bias  degen |  NoiseRef 1ms σ  bias |  NoiseRef $(M)ms σ  bias",
    )
    for cn0 in CN0S_DBHZ
        s = SWEEP[cn0]
        println(
            lpad(cn0, 5),
            " | ",
            lpad(round(s.nwpr.σ_db; digits = 2), 7),
            lpad(round(s.nwpr.bias_db; digits = 2), 6),
            lpad(round(100 * s.nwpr.degenerate; digits = 1), 6),
            "% | ",
            lpad(round(s.noncoh.σ_db; digits = 2), 14),
            lpad(round(s.noncoh.bias_db; digits = 2), 6),
            " | ",
            lpad(round(s.coh.σ_db; digits = 2), 14),
            lpad(round(s.coh.bias_db; digits = 2), 6),
        )
    end
    @test true
end

end
