module DiscriminatorsTest

using Test: @test, @testset, @inferred, @test_throws
using GNSSSignals:
    GPSL1CA,
    GPSL1C_P,
    GalileoE1B,
    GalileoE1B_BOC11,
    GalileoE1C,
    GalileoE5aI,
    GalileoE5aQ,
    GPSL1C_D,
    get_code,
    get_relative_power,
    get_code_frequency,
    get_code_length
using StaticArrays: SVector
using Unitful: Hz, MHz, ms, upreferred
using Tracking:
    Tracking,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    pll_disc,
    fll_disc,
    dll_disc,
    dll_disc_noise_gain,
    get_correlator_sample_shifts,
    get_early_late_sample_spacing,
    get_prompt,
    update_accumulator

@testset "PLL discriminator" begin
    correlator_minus60off = EarlyPromptLateCorrelator(
        SVector(-0.5 + sqrt(3) / 2im, -1 + sqrt(3) * 1im, -0.5 + sqrt(3) / 2im),
        0.5,
    )
    correlator_0off =
        EarlyPromptLateCorrelator(SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im), 0.5)
    correlator_plus60off = EarlyPromptLateCorrelator(
        SVector(0.5 + sqrt(3) / 2im, 1 + sqrt(3) * 1im, 0.5 + sqrt(3) / 2im),
        0.5,
    )
    gpsl1 = GPSL1CA()
    @test @inferred(pll_disc(gpsl1, correlator_minus60off)) == -π / 3  #-60°
    @test @inferred(pll_disc(gpsl1, correlator_0off)) == 0
    @test @inferred(pll_disc(gpsl1, correlator_plus60off)) == π / 3  #+60°
end

@testset "FLL discriminator" begin
    correlator_minus60off = EarlyPromptLateCorrelator(
        SVector(-0.5 + sqrt(3) / 2im, -1 + sqrt(3) * 1im, -0.5 + sqrt(3) / 2im),
        0.5,
    )
    correlator_0off =
        EarlyPromptLateCorrelator(SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im), 0.5)
    correlator_plus60off = EarlyPromptLateCorrelator(
        SVector(0.5 + sqrt(3) / 2im, 1 + sqrt(3) * 1im, 0.5 + sqrt(3) / 2im),
        0.5,
    )
    correlator_plus90off =
        EarlyPromptLateCorrelator(SVector(0.0 + 0.5im, 0.0 + 1.0im, 0.0 + 0.5im), 0.5)
    correlator_empty = EarlyPromptLateCorrelator()
    gpsl1 = GPSL1CA()
    @test @inferred(
        fll_disc(gpsl1, correlator_0off, get_prompt(correlator_minus60off), 1ms)
    ) == (166 + 2 / 3) * 1Hz
    @test @inferred(fll_disc(gpsl1, correlator_0off, get_prompt(correlator_0off), 1ms)) ==
          0.0Hz
    @test @inferred(
        fll_disc(gpsl1, correlator_0off, get_prompt(correlator_plus60off), 1ms)
    ) == -(166 + 2 / 3) * 1Hz
    @test @inferred(
        fll_disc(gpsl1, correlator_0off, get_prompt(correlator_plus90off), 1ms)
    ) == -250Hz
    @test @inferred(fll_disc(gpsl1, correlator_0off, get_prompt(correlator_empty), 1ms)) ==
          0Hz
end

@testset "DLL discriminator" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = get_code_frequency(gpsl1) * 4

    very_early_correlator =
        EarlyPromptLateCorrelator(SVector(1.0 + 0.0im, 0.5 + 0.0im, 0.0 + 0.0im), 0.5)
    early_correlator =
        EarlyPromptLateCorrelator(SVector(0.75 + 0.0im, 0.75 + 0.0im, 0.25 + 0.0im), 0.5)
    prompt_correlator =
        EarlyPromptLateCorrelator(SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im), 0.5)
    late_correlator =
        EarlyPromptLateCorrelator(SVector(0.25 + 0.0im, 0.75 + 0.0im, 0.75 + 0.0im), 0.5)
    very_late_correlator =
        EarlyPromptLateCorrelator(SVector(0.0 + 0.0im, 0.5 + 0.0im, 1.0 + 0.0im), 0.5)

    @test @inferred(
        get_early_late_sample_spacing(
            prompt_correlator,
            sampling_frequency,
            get_code_frequency(gpsl1),
        )
    ) == 4

    @test @inferred(dll_disc(gpsl1, very_early_correlator, 0.0Hz, sampling_frequency)) ==
          -0.5
    @test @inferred(dll_disc(gpsl1, early_correlator, 0.0Hz, sampling_frequency)) == -0.25
    @test @inferred(dll_disc(gpsl1, prompt_correlator, 0.0Hz, sampling_frequency)) == 0
    @test @inferred(dll_disc(gpsl1, late_correlator, 0.0Hz, sampling_frequency)) == 0.25
    @test @inferred(dll_disc(gpsl1, very_late_correlator, 0.0Hz, sampling_frequency)) == 0.5
end

# A chips-calibrated discriminator has an S-curve slope of 1 at the origin: sweep a true
# code offset τ over the linear region, build the accumulators by correlating the
# (noiseless) code against tap replicas at the correlator's actual sample-quantized
# offsets, and least-squares fit `dll_disc`'s output against τ.
function dll_disc_s_curve_slope(
    signal,
    correlator,
    sampling_frequency;
    prn = 1,
    offsets = -0.05:0.005:0.05,
)
    code_frequency = get_code_frequency(signal)
    chips_per_sample = upreferred(code_frequency / sampling_frequency)
    samples_per_code = round(Int, get_code_length(signal) / chips_per_sample)
    shifts = get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    phases = (0:(samples_per_code-1)) .* chips_per_sample
    discriminator = map(offsets) do offset
        incoming = get_code.(signal, phases .+ offset, prn)
        accumulators = map(shifts) do shift
            replica = get_code.(signal, phases .+ shift * chips_per_sample, prn)
            complex(sum(incoming .* replica) / samples_per_code, 0.0)
        end
        dll_disc(
            signal,
            update_accumulator(correlator, SVector(accumulators)),
            0.0Hz,
            sampling_frequency,
        )
    end
    sum(offsets .* discriminator) / sum(abs2, offsets)
end

@testset "DLL discriminator S-curve slope is 1 (calibrated in chips)" begin
    veml = VeryEarlyPromptLateCorrelator()
    # EPL on GPS L1 C/A: the (2 - d) / 2 normalization.
    @test dll_disc_s_curve_slope(GPSL1CA(), EarlyPromptLateCorrelator(), 20MHz) ≈ 1.0 atol =
        0.02
    # VEML on BOC(1,1) — the normalization's own correlation model — at two sampling
    # rates that quantize the taps differently (±3/±12 and ±12/±47 samples).
    @test dll_disc_s_curve_slope(GalileoE1B_BOC11(), veml, 20MHz) ≈ 1.0 atol = 0.02
    @test dll_disc_s_curve_slope(GalileoE1B_BOC11(), veml, 79.5MHz) ≈ 1.0 atol = 0.02
    # VEML on the full modulations: BOC(1,1) is the modulation model the normalization
    # assumes, so its residual gain error is the modulation mismatch — a few percent for
    # the additive CBOC⁺ (E1B) and TMBOC (L1C pilot), more for E1C's anti-phased CBOC⁻,
    # whose subtracted BOC(6,1) component flattens the correlation peak.
    @test dll_disc_s_curve_slope(GalileoE1B(), veml, 20MHz) ≈ 1.0 atol = 0.1
    @test dll_disc_s_curve_slope(GPSL1C_P(), veml, 20MHz) ≈ 1.0 atol = 0.1
    @test dll_disc_s_curve_slope(GalileoE1C(), veml, 20MHz) ≈ 1.0 atol = 0.2
end

# A correlator type with no `dll_disc_noise_gain` method: there is no generic gain
# to fall back to.
struct UnmodelledCorrelator <: Tracking.AbstractCorrelator{1} end

@testset "DLL discriminator noise gains" begin
    fs = 5e6Hz
    # A BOC VEML discriminator is several times more accurate than a 1-chip BPSK
    # early-late one at equal SNR, so its noise gain must be correspondingly
    # smaller — that is what lets a pilot outvote a legacy signal.
    veml = VeryEarlyPromptLateCorrelator()
    epl = EarlyPromptLateCorrelator()
    g_veml = dll_disc_noise_gain(GPSL1C_P(), veml, 0.0Hz, fs)
    g_epl = dll_disc_noise_gain(GPSL1CA(), epl, 0.0Hz, fs)
    @test 0 < g_veml < g_epl

    # The documented values, not just their ordering. `d / 4` is Kaplan & Hegarty's
    # tracking-jitter constant for the noncoherent early-minus-late envelope
    # discriminator with early-late spacing `d` chips, pinned at the two sampling
    # frequencies where the preferred chip shift lands on a whole number of
    # samples and `d` is therefore exactly what it was asked for — which also pins
    # the linear `d` dependence.
    @test dll_disc_noise_gain(
        GPSL1CA(),
        EarlyPromptLateCorrelator(),                  # ±0.5 chips
        0.0Hz,
        2.046e6Hz,                                    # 0.5 chips = 1 sample
    ) === 0.25
    @test dll_disc_noise_gain(
        GPSL1CA(),
        EarlyPromptLateCorrelator([complex(0.0), complex(0.0), complex(0.0)], 0.25),
        0.0Hz,
        4.092e6Hz,                                    # 0.25 chips = 1 sample
    ) === 0.125
    # The VEML gain for the default ±0.15/±0.6 taps: `1 / (2 · (3 + 1)²)`, about
    # 8× the precision of the 1-chip early-late layout.
    @test g_veml ≈ 0.03125 rtol = 1e-12
    @test dll_disc_noise_gain(GPSL1C_P(), veml, 0.0Hz, 20e6Hz) ≈ 0.03125 rtol = 1e-9

    # A degenerate tap layout — both pairs past the 1-chip correlation support —
    # carries no delay information and must come out as `Inf`, so a record gets
    # weight *zero*. Not `NaN`: the envelope sum is what the S-curve slope divides
    # by, so `-0.0 / 0.0` is one line away, and `NaN` compares equal to nothing.
    # Real accumulator values, because an all-zero correlator is a separate `0/0`.
    degenerate = VeryEarlyPromptLateCorrelator(
        [
            complex(0.60, 0.0),
            complex(0.90, 0.0),
            complex(1.00, 0.0),
            complex(0.80, 0.0),
            complex(0.50, 0.0),
        ],
        2.0,
        3.0,
    )
    @test dll_disc_noise_gain(GalileoE1C(), degenerate, 0.0Hz, 20e6Hz) === Inf
    @test Tracking._veml_discriminator_slope(2.0, 3.0) === 0.0
    # The paired guard in the discriminator itself: an uncalibrated raw value
    # rather than a division by a zero (or NaN) slope.
    @test isfinite(dll_disc(GalileoE1C(), degenerate, 0.0Hz, 20e6Hz))

    # No generic gain: a correlator type without a `dll_disc_noise_gain` method is
    # a `MethodError`, like one without `dll_disc`.
    @test_throws MethodError dll_disc_noise_gain(
        GPSL1CA(),
        UnmodelledCorrelator(),
        0.0Hz,
        fs,
    )
end

@testset "the shared tap-offset helpers are what dll_disc reads" begin
    # `dll_disc` and `dll_disc_noise_gain` must be evaluated at the same tap
    # offsets, or the gain is not the noise gain of that discriminator.
    fs = 5e6Hz
    code_doppler = 3.0Hz
    epl = EarlyPromptLateCorrelator()
    d = Tracking._early_late_spacing(GPSL1CA(), epl, code_doppler, fs)
    code_frequency = code_doppler + get_code_frequency(GPSL1CA())
    @test d ≈
          get_early_late_sample_spacing(epl, fs, code_frequency) *
          upreferred(code_frequency / fs)
    @test dll_disc_noise_gain(GPSL1CA(), epl, code_doppler, fs) === d / 4

    veml = VeryEarlyPromptLateCorrelator()
    inner, outer = Tracking._veml_tap_offsets(GalileoE1B(), veml, 0.0Hz, 20e6Hz)
    # ±3 / ±12 samples at 20 MHz and 1.023 MHz.
    @test inner ≈ 3 * 1.023e6 / 20e6
    @test outer ≈ 12 * 1.023e6 / 20e6
end

@testset "relative record SNR is the nominal power share times the sample count" begin
    # The prompt does not enter it at all: the weight is nominal.
    @test Tracking._relative_record_snr(GPSL1CA(), 100) ===
          get_relative_power(GPSL1CA()) * 100
    @test Tracking._relative_record_snr(GPSL1CA(), 200) ===
          get_relative_power(GPSL1CA()) * 200
    # GPS L1C's 75/25 pilot/data split is the one asymmetric intra-band case.
    @test get_relative_power(GPSL1C_P()) / get_relative_power(GPSL1C_D()) === 3.0
    # Galileo's intra-band pairs split evenly.
    @test get_relative_power(GalileoE1B()) === get_relative_power(GalileoE1C()) === 0.5
    @test get_relative_power(GalileoE5aI()) === get_relative_power(GalileoE5aQ()) === 0.5
end

end
