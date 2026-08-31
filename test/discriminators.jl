module DiscriminatorsTest

using Test: @test, @testset, @inferred
using GNSSSignals:
    GPSL1CA,
    GPSL1C_P,
    GalileoE1B,
    GalileoE1B_BOC11,
    GalileoE1C,
    get_code,
    get_code_frequency,
    get_code_length
using StaticArrays: SVector
using Unitful: Hz, MHz, ms, upreferred
using Tracking:
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    pll_disc,
    fll_disc,
    dll_disc,
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

# A correlator with no energy at all used to make every discriminator compute
# 0 / 0. The NaN then propagated into the loop filter and out as a NaN Doppler,
# and the next correlate iteration threw `InexactError: Int64(NaN)` while
# converting a Doppler-derived sample count -- killing the whole tracking task
# from inside, with nothing pointing back at the empty correlator. Reported from
# on-FPGA correlator bring-up (GNSSReceiver.jl#107), where a reassigned hardware
# channel can deliver a dump with `integrated_samples == 0`.
@testset "Discriminators with a zero-energy correlator" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = get_code_frequency(gpsl1) * 4
    empty_correlator = EarlyPromptLateCorrelator()
    empty_veml = VeryEarlyPromptLateCorrelator()

    @test iszero(get_prompt(empty_correlator))

    @testset "DLL returns zero rather than NaN" begin
        d = @inferred(dll_disc(gpsl1, empty_correlator, 0.0Hz, sampling_frequency))
        @test !isnan(d)
        @test d == 0
    end

    @testset "VEML DLL returns zero rather than NaN" begin
        d = @inferred(dll_disc(gpsl1, empty_veml, 0.0Hz, sampling_frequency))
        @test !isnan(d)
        @test d == 0
    end

    @testset "PLL returns zero rather than NaN" begin
        p = @inferred(pll_disc(gpsl1, empty_correlator))
        @test !isnan(p)
        @test p == 0
    end

    @testset "FLL returns zero rather than NaN for a zero current prompt" begin
        nonzero_prompt = get_prompt(
            EarlyPromptLateCorrelator(SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im), 0.5),
        )
        f = @inferred(fll_disc(gpsl1, empty_correlator, nonzero_prompt, 1ms))
        @test !isnan(f)
        @test f == 0Hz
    end

    # The guard must not swallow the legitimate quadrature cases: a purely
    # imaginary prompt is +-pi/2, not "no energy".
    @testset "A purely imaginary prompt is still +-pi/2" begin
        quadrature = EarlyPromptLateCorrelator(
            SVector(0.0 + 0.5im, 0.0 + 1.0im, 0.0 + 0.5im),
            0.5,
        )
        @test @inferred(pll_disc(gpsl1, quadrature)) == π / 2
        minus_quadrature = EarlyPromptLateCorrelator(
            SVector(0.0 - 0.5im, 0.0 - 1.0im, 0.0 - 0.5im),
            0.5,
        )
        @test @inferred(pll_disc(gpsl1, minus_quadrature)) == -π / 2
    end

    # E == L with real energy is a genuine zero code error, and must stay
    # distinguishable from the no-energy case above (both report 0, but only the
    # latter may skip the normalisation).
    @testset "A balanced correlator with energy still discriminates" begin
        balanced =
            EarlyPromptLateCorrelator(SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im), 0.5)
        @test @inferred(dll_disc(gpsl1, balanced, 0.0Hz, sampling_frequency)) == 0
        off_late =
            EarlyPromptLateCorrelator(SVector(0.25 + 0.0im, 0.75 + 0.0im, 0.75 + 0.0im), 0.5)
        @test @inferred(dll_disc(gpsl1, off_late, 0.0Hz, sampling_frequency)) == 0.25
    end
end

end
