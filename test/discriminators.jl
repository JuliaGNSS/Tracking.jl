module DiscriminatorsTest

using Test: @test, @testset, @inferred
using GNSSSignals: GPSL1, get_code_frequency
using StaticArrays: SVector
using Unitful: Hz
using Tracking: EarlyPromptLateCorrelator, pll_disc, dll_disc, get_early_late_sample_spacing

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
    gpsl1 = GPSL1()
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
    gpsl1 = GPSL1()
    @test @inferred(fll_disc(gpsl1, correlator_0off, get_prompt(correlator_minus60off), 1ms)) == (166 + 2/3) * 1Hz
    @test @inferred(fll_disc(gpsl1, correlator_0off, get_prompt(correlator_0off), 1ms)) == 0.0Hz
    @test @inferred(fll_disc(gpsl1, correlator_0off, get_prompt(correlator_plus60off), 1ms)) == - (166 + 2/3) * 1Hz
end

@testset "DLL discriminator" begin
    gpsl1 = GPSL1()
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

end
