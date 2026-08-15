module NoiseRefInTheLoopTest

# What had to be true before `NoiseRefCN0Estimator` became the default:
#
#   1. It works on the signals `NWPRCN0Estimator` cannot serve at all — GPS
#      L1C-D and Galileo E1B (one symbol per code block) and any secondary-coded
#      signal before sync. That is issue #217, and the reason for the change.
#   2. In a real loop it beats NWPR where lock decisions are made, and — unlike
#      NWPR — it is immune to the phase-noise wash at long coherent records.
#   3. The one bias it carries, the reference's self-leakage, is the size the
#      design says it is.

using Test: @test, @testset
using Random: Xoshiro, randn, rand
using Statistics: std
using Unitful: Hz, s, ms, dBHz, ustrip, uconvert
using GNSSSignals:
    GPSL1CA, GPSL1C_D, GPSL5I, GalileoE1B, gen_code, get_code_frequency, get_code_length
import Tracking
using Tracking:
    BandMeasurement,
    BitBuffer,
    CN0UpdateContext,
    CPUDownconvertAndCorrelator,
    CorrelatorNoiseEstimator,
    MomentsCN0Estimator,
    NWPRCN0Estimator,
    NoiseRefCN0Estimator,
    TrackState,
    TrackedSat,
    estimate_cn0,
    get_cn0_estimator,
    get_noise_density,
    has_bit_or_secondary_code_been_found,
    reset_loop_filters!,
    set_preferred_num_code_blocks_to_integrate!,
    track!,
    update

_db(x) = ustrip(uconvert(dBHz, x))
_median(xs) = (v = sort(collect(xs)); v[div(length(v) + 1, 2)])

# Float64, matching what the tracking loop actually builds — see the note in
# `cn0_estimators/noise_ref.jl`.
const T = 1.0ms
const N₀ = uconvert(Hz^-1, T)              # `E|P|² = N₀/T`, so a unit-noise prompt

# `λ = (C/N₀)·T` per record, CN(0,1) noise, perfect phase — the same model
# `test/cn0_estimator_comparison.jl` uses.
_prompts(λ, n, rng) = sqrt(λ) .+ (randn(rng, n) .+ im .* randn(rng, n)) ./ sqrt(2)

@testset "the signals NWPR cannot serve at all (issue #217)" begin
    # NWPR needs a coherent window of at least two records tiling a navigation
    # bit. Three configurations admit none, ever, and there it reports its
    # fallback for good — a *different* estimator with a different bias and, at
    # the default `MomentsCN0Estimator`, a ≈27.6 dB-Hz floor on pure noise. The
    # noise reference has no window, so all three are ordinary records to it.
    bit_buffer = BitBuffer{UInt64}()
    cases = (
        # (signal, blocks per bit, this record's offset in the bit)
        ("GPS L1C-D, synced", GPSL1C_D(), 1, 0),
        ("Galileo E1B, synced", GalileoE1B(), 1, 0),
        ("GPS L5I, pre-sync", GPSL5I(), 10, -1),
    )
    for (name, signal, blocks_per_bit, block_index) in cases
        @testset "$name" begin
            context =
                CN0UpdateContext(signal, 1, blocks_per_bit, block_index, bit_buffer, N₀, T)
            rng = Xoshiro(20260806)
            λ = 10^3.0 * 1e-3                          # a true 30 dB-Hz
            prompts = _prompts(λ, 400, rng)

            nwpr =
                foldl((e, p) -> update(e, p, context), prompts; init = NWPRCN0Estimator())
            reference = foldl(
                (e, p) -> update(e, p, context),
                prompts;
                init = NoiseRefCN0Estimator(; num_records = 400),
            )

            # NWPR buffered nothing at all and is reporting its fallback.
            @test Base.length(nwpr) == 0
            @test estimate_cn0(nwpr, T) ==
                  estimate_cn0(Tracking.get_fallback_cn0_estimator(nwpr), T)
            # The noise reference: inside 1 dB of the truth, with no window and
            # no fallback anywhere in sight.
            @test _db(estimate_cn0(reference, T)) ≈ 30 atol = 1.0

            # ... and on **pure noise** the difference is the floor #217 is
            # about: NWPR's fallback manufactures signal power out of the sample
            # moments, the reference divides by a floor it was told.
            noise = _prompts(0.0, 400, rng)
            nwpr_noise =
                foldl((e, p) -> update(e, p, context), noise; init = NWPRCN0Estimator())
            reference_noise = foldl(
                (e, p) -> update(e, p, context),
                noise;
                init = NoiseRefCN0Estimator(; num_records = 400),
            )
            # A single realisation, so the bound is loose — the documented
            # ≈27.6 dB-Hz is a median over seeds. The gap is the claim.
            @test _db(estimate_cn0(nwpr_noise, T)) > 18
            @test _db(estimate_cn0(reference_noise, T)) < 15
            @test _db(estimate_cn0(nwpr_noise, T)) - _db(estimate_cn0(reference_noise, T)) >
                  5
        end
    end
end

# One satellite through `track!`, data-modulated like a real L1 C/A signal.
# `blocks > 1` lengthens the coherent record once bit sync is in, which is the
# only point the cap allows it.
function _track_noisy(cn0_db, num_blocks; seed = 1, cn0_estimator, blocks = 1)
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    num_samples = 4000
    code_frequency = get_code_frequency(gpsl1)
    rng = Xoshiro(seed)
    amplitude = 10^(cn0_db / 20)
    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator)])
    dc = CPUDownconvertAndCorrelator()
    data_bits = rand(rng, (-1.0, 1.0), div(num_blocks, 20) + 2)
    promoted = false
    for block_index = 1:num_blocks
        code_phase = (block_index - 1) * num_samples * code_frequency / fs
        clean =
            amplitude .* data_bits[div(block_index-1, 20)+1] .*
            gen_code(num_samples, gpsl1, 1, fs, code_frequency, ustrip(code_phase))
        noise = randn(rng, ComplexF64, num_samples) .* sqrt(ustrip(Hz, fs))
        track!(clean .+ noise, track_state, fs; downconvert_and_correlator = dc)
        if blocks > 1 && !promoted && has_bit_or_secondary_code_been_found(track_state, 1)
            set_preferred_num_code_blocks_to_integrate!(track_state, 1, blocks)
            reset_loop_filters!(track_state, 1)
            promoted = true
        end
    end
    track_state
end

@testset "in the loop, against NWPR" begin
    # The gate on flipping the default. Medians over seeds, because every number
    # here is an estimate of a statistical quantity.
    #
    # Measured over 9 seeds and 1200 blocks (this runs 5 and 700 to stay quick):
    #
    # |true| NWPR 1-block      | NoiseRef      | NWPR 20-block  | NoiseRef      |
    # |----|-------------------|---------------|----------------|---------------|
    # | 25 | 11.2 dB, 56 % -Inf| 23.4, 0 %     | —              | —             |
    # | 30 | 29.3 ± 1.40       | 29.6 ± 0.62   | —              | —             |
    # | 40 | 39.7 ± 0.35       | 39.9 ± 0.19   | 27.1 ± 5.79    | 39.8 ± 0.37   |
    # | 45 | 44.7 ± 0.30       | 44.9 ± 0.13   | 32.9 ± 2.00    | 44.7 ± 0.30   |
    seeds = 1:5
    gpsl1 = GPSL1CA()
    run(cn0, est; blocks = 1) = [
        _db(estimate_cn0(_track_noisy(cn0, 700; seed, cn0_estimator = est(), blocks), 1)) for seed in seeds
    ]

    # 1-block records, 25 dB-Hz — the regime lock and loss decisions are made in.
    # What is asserted here is the *structural* claim: the noise reference cannot
    # produce a degenerate output, because it has no `1 < μ̂ < M` bound to fall
    # outside of, and it lands near the truth.
    #
    # NWPR's degenerate rate here is real and large — 56 % over 9 seeds × 1200
    # blocks — but five seeds of 700 blocks cannot pin it, and the RNG stream
    # differs between Julia versions, so the same seeds give 1 in 5 degenerate on
    # 1.12 and 0 in 5 on 1.10. Asserting an inequality on that count would be
    # testing the RNG. The long-record collapse below is the robust comparison.
    ref_low = run(25.0, () -> NoiseRefCN0Estimator())
    @test count(!isfinite, ref_low) == 0
    @test _median(ref_low) ≈ 25 atol = 3.0

    # 45 dB-Hz, 1-block records: both work, the reference is tighter and closer.
    nwpr_high = run(45.0, () -> NWPRCN0Estimator(gpsl1))
    ref_high = run(45.0, () -> NoiseRefCN0Estimator())
    @test _median(ref_high) ≈ 45 atol = 1.0
    @test abs(_median(ref_high) - 45) <= abs(_median(nwpr_high) - 45) + 0.1
    @test std(ref_high) < std(nwpr_high)

    # 20-block records — a whole GPS L1 C/A navigation bit. NWPR's window closes
    # on a single record, where `NBP ≡ WBP` and no window exists at all, so it
    # falls back for good. The non-coherent reference does not notice: this is
    # the phase-noise wash the plan deliberately left out of its model, and the
    # one number that had to come from a real loop.
    nwpr_long = run(45.0, () -> NWPRCN0Estimator(gpsl1); blocks = 20)
    ref_long = run(45.0, () -> NoiseRefCN0Estimator(); blocks = 20)
    @test _median(ref_long) ≈ 45 atol = 1.5
    @test _median(ref_long) - _median(nwpr_long) > 5
    # ... and the reference barely moves between the two record lengths.
    @test abs(_median(ref_long) - _median(ref_high)) < 1.0
end

@testset "the reference's self-leakage is the size the design says" begin
    # The one bias this estimator carries and NWPR does not: the reference
    # despreads with a *wrong* PRN, so it also collects the tracked satellite's
    # own power, `ε_self = C/f_chip`. Measured as a paired comparison — the same
    # noise realisation, the same PRN rotation, the same code-phase and carrier
    # draws, with and without the satellite — so the loop's own σ cancels and the
    # prediction is testable at a few per cent rather than needing hundreds of
    # seeds. The reference's `rng` is what makes the last of those hold: the
    # draws are random but seeded, so the two runs despread identically.
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    num_samples = 40_000
    f_chip = ustrip(Hz, get_code_frequency(gpsl1))
    rng = Xoshiro(20260806)
    noise = sqrt(ustrip(Hz, fs)) .* randn(rng, ComplexF64, num_samples)

    function density(signal)
        # A freshly seeded generator per run, so both borrow the same codes in
        # the same order *and* despread them at the same phases and carriers.
        ts = TrackState(
            gpsl1,
            [TrackedSat(gpsl1, 1, 0.0, 0.0Hz)];
            noise_estimators = (GPSL1CA = CorrelatorNoiseEstimator(; rng = Xoshiro(1)),),
        )
        dc = CPUDownconvertAndCorrelator()
        for _ = 1:10
            track!(signal, ts, fs; downconvert_and_correlator = dc)
        end
        ustrip(Hz^-1, get_noise_density(ts.noise_estimators.GPSL1CA))
    end

    for cn0 in (45.0, 50.0)
        carrier_power = 10^(cn0 / 10)
        satellite =
            10^(cn0 / 20) .*
            gen_code(num_samples, gpsl1, 1, fs, get_code_frequency(gpsl1), 0.0)
        ratio = density(ComplexF64.(satellite) .+ noise) / density(noise)
        # `N̂₀ → N₀ + C/f_chip`, i.e. a ratio of `1 + (C/N₀)/f_chip` against the
        # `N₀ = 1 Hz⁻¹` this noise is built to have: 1.031 at 45 dB-Hz, 1.098 at
        # 50. That is 0.13 dB and 0.40 dB of C/N₀ read low.
        @test ratio ≈ 1 + carrier_power / f_chip rtol = 0.25
        @test ratio > 1
    end
end

end
