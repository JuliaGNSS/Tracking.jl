module CN0EstimatorComparisonTest

# Baseline for docs/plans/2026-08-06-per-band-noise-estimation-for-cn0.md.
#
# Pins the four structural claims the plan's argument rests on, by Monte-Carlo over
# the post-correlation prompt model rather than by algebra:
#
#   1. NWPR reports degenerate (non-finite) estimates at low C/N₀; a noise-reference
#      estimator never does.
#   2. NWPR carries a small positive bias at every C/N₀; a noise-reference estimator
#      is unbiased.
#   3. NWPR's relative error saturates — more signal buys nothing past ~40 dB-Hz —
#      whereas a noise-reference estimator keeps improving.
#   4. A noise-reference estimator fed records as long as NWPR's narrowband window
#      beats NWPR at *every* C/N₀. This is the claim that justifies replacing it.
#
# `NWPRCN0Estimator` here is the shipped implementation, driven through its
# bare-prompt-stream `update` (free-running windows of `num_narrowband_code_blocks`,
# one block per record) — so claims 1-3 are assertions about Tracking's real
# behaviour, not about a model of it. The noise-reference columns are reference
# implementations: `NoiseRefCN0Estimator` does not exist yet, and the point of this
# file is to pin the numbers it must reproduce when it lands.
#
# NWPR is given its best case throughout: perfect carrier phase, no bit transitions,
# windows aligned to the record grid. Its real-world numbers can only be worse.

using Test: @test, @testset
using Random: Xoshiro, randn
using Statistics: mean, std
using Unitful: ms, s, dBHz, ustrip, uconvert
using Tracking: NWPRCN0Estimator, update, estimate_cn0

# Post-correlation prompt model. `λ = (C/N₀)·T` is the per-record post-correlation
# SNR; noise is CN(0,1) so the noise power is 1 and λ is also the signal power.
# Perfect phase, so the signal sits on the real axis.
_prompts(λ, num_records, rng) =
    sqrt(λ) .+ (randn(rng, num_records) .+ im .* randn(rng, num_records)) ./ sqrt(2)

# Reference non-coherent noise-reference estimator: mean prompt power minus the
# known noise power. λ̂ = ⟨|P|²⟩ − 1. No bit grid, no window, no `M`.
_noise_ref_noncoherent(prompts) = mean(abs2, prompts) - 1

# Reference coherent noise-reference estimator: sum `M` prompts coherently, then
# subtract the known noise power of the sum. Statistically identical to running
# records `M` times longer, which is the point — NWPR's narrowband window is a
# longer coherent integration, not extra information.
function _noise_ref_coherent(prompts, M)
    num_windows = div(length(prompts), M)
    acc = 0.0
    for w = 1:num_windows
        window = @view prompts[((w-1)*M+1):(w*M)]
        acc += (abs2(sum(window)) - M) / M^2
    end
    acc / num_windows
end

# The shipped estimator, folded over the same prompt stream. Returns λ̂ so all three
# are on one scale.
function _nwpr(prompts, num_records, M, integration_time)
    estimator = NWPRCN0Estimator(;
        num_records,
        num_narrowband_code_blocks = M,
        num_presync_narrowband_code_blocks = M,
    )
    for prompt in prompts
        estimator = update(estimator, prompt)
    end
    # λ = (C/N₀)·T. `estimate_cn0` returns dB-Hz, so undo the log and re-apply T.
    cn0_db = ustrip(uconvert(dBHz, estimate_cn0(estimator, integration_time)))
    10^(cn0_db / 10) * ustrip(s, integration_time)
end

# NWPR's inversion is only defined for `1 < μ̂ < M`. Outside it `estimate_cn0` returns
# 0 dB-Hz or Inf dB-Hz (src/cn0_estimation.jl:704-706), which map to λ̂ of exactly
# zero and Inf. A noise-reference estimate may legitimately be negative at low C/N₀
# — that is what makes it unbiased, not a degeneracy — so only exact zero and
# non-finite count here.
_is_degenerate(λ̂) = !isfinite(λ̂) || iszero(λ̂)

# Relative standard deviation in dB (delta method) and bias of the linear-domain
# mean, both over the non-degenerate estimates, with the degenerate fraction
# reported separately. Excluding them is the *generous* reading for NWPR — the
# discarded draws are its worst ones.
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
const NUM_TRIALS = 4_000         # keeps the file at a few seconds; σ good to ~1 %
const CN0S_DBHZ = 20:5:50

function _sweep(cn0_dbhz, seed)
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

# One sweep, reused by every testset below so the Monte Carlo runs once.
const SWEEP = Dict(cn0 => _sweep(cn0, 20260806 + cn0) for cn0 in CN0S_DBHZ)

@testset "NWPR reports degenerate estimates at low C/N₀, a noise reference never does" begin
    # NWPR's inversion is only defined for `1 < μ̂ < M`; outside it `estimate_cn0`
    # returns 0 or Inf dB-Hz (src/cn0_estimation.jl:704-706). At 20 dB-Hz that is a
    # ≈6 % of windows at 20 dB-Hz — one in sixteen, which a lock detector cannot
    # use. Threshold set well below the measured rate so the test pins the
    # structural fact, not the exact rate.
    @test SWEEP[20].nwpr.degenerate > 0.01
    for cn0 in CN0S_DBHZ
        @test SWEEP[cn0].noncoh.degenerate == 0
        @test SWEEP[cn0].coh.degenerate == 0
    end
end

@testset "NWPR is biased at every C/N₀, a noise reference is not" begin
    # The ratio-of-sums bias (src/cn0_estimation.jl:219-247) does not vanish as the
    # signal grows: it settles around +0.05 dB. The noise-reference estimators are
    # unbiased by construction, because the noise power is known rather than
    # inferred from the same prompts.
    for cn0 = 30:5:50
        @test SWEEP[cn0].nwpr.bias_db > 0.02
        @test abs(SWEEP[cn0].noncoh.bias_db) < 0.02
        @test abs(SWEEP[cn0].coh.bias_db) < 0.02
    end
end

@testset "NWPR's relative error saturates, a noise reference's keeps falling" begin
    # NWPR asymptotes to √(M/((M−1)K)) ≈ 0.49 dB regardless of signal strength —
    # the `μ̂ → M` ceiling. Between 45 and 50 dB-Hz it therefore barely moves, while
    # the noise-reference estimators roughly halve.
    @test SWEEP[50].nwpr.σ_db > 0.9 * SWEEP[45].nwpr.σ_db
    @test SWEEP[50].coh.σ_db < 0.7 * SWEEP[45].coh.σ_db
    @test SWEEP[50].noncoh.σ_db < 0.7 * SWEEP[45].noncoh.σ_db
end

@testset "A coherent noise reference beats NWPR at every C/N₀" begin
    # The plan's central claim. A coherent sum of M records *is* one M-times-longer
    # record, so this is what NWPR's narrowband window costs versus taking the same
    # coherent gain in the correlator and measuring the noise directly.
    for cn0 in CN0S_DBHZ
        @test SWEEP[cn0].coh.σ_db < SWEEP[cn0].nwpr.σ_db
    end
end

@testset "Non-coherent squaring loss is real below ~28 dB-Hz" begin
    # The honest half of the comparison: with records left at 1 ms, dropping NWPR's
    # coherent window costs variance at low C/N₀ — the low-SNR ratio is √(M−1) = 2 —
    # and only pays off above the crossover. This is why the plan recommends
    # lengthening `preferred_num_code_blocks_to_integrate` rather than presenting
    # the non-coherent form as a free win.
    @test SWEEP[20].noncoh.σ_db > SWEEP[20].nwpr.σ_db
    @test SWEEP[40].noncoh.σ_db < SWEEP[40].nwpr.σ_db
end

# Print the table so a regression shows what moved, and so the numbers quoted in
# docs/plans/2026-08-06-per-band-noise-estimation-for-cn0.md stay reproducible.
@testset "baseline table" begin
    println(
        "\nC/N₀ vs estimator, $(NUM_RECORDS) records of $(INTEGRATION_TIME), ",
        "NWPR M=$M, $(NUM_TRIALS) trials",
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
