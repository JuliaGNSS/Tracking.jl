module DiscriminatorsTest

using Test: @test, @testset, @inferred
using GNSSSignals: GPSL1
using StaticArrays: SVector
using Tracking: EarlyPromptLateCorrelator, pll_disc, dll_disc, get_early_late_sample_spacing

@testset "PLL discriminator" begin
    correlator_minus60off = EarlyPromptLateCorrelator(
        SVector(-0.5 + sqrt(3) / 2im, -1 + sqrt(3) * 1im, -0.5 + sqrt(3) / 2im),
        -1:1,
        3,
        2,
        1,
    )
    correlator_0off = EarlyPromptLateCorrelator(
        SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im),
        -1:1,
        3,
        2,
        1,
    )
    correlator_plus60off = EarlyPromptLateCorrelator(
        SVector(0.5 + sqrt(3) / 2im, 1 + sqrt(3) * 1im, 0.5 + sqrt(3) / 2im),
        -1:1,
        3,
        2,
        1,
    )
    gpsl1 = GPSL1()
    @test @inferred(pll_disc(gpsl1, correlator_minus60off)) == -π / 3  #-60°
    @test @inferred(pll_disc(gpsl1, correlator_0off)) == 0
    @test @inferred(pll_disc(gpsl1, correlator_plus60off)) == π / 3  #+60°
end

@testset "DLL discriminator" begin
    sample_shifts = SVector(-2, 0, 2)
    delta = 0.25

    gpsl1 = GPSL1()
    very_early_correlator = EarlyPromptLateCorrelator(
        SVector(1.0 + 0.0im, 0.5 + 0.0im, 0.0 + 0.0im),
        sample_shifts,
        3,
        2,
        1,
    )
    early_correlator = EarlyPromptLateCorrelator(
        SVector(0.75 + 0.0im, 0.75 + 0.0im, 0.25 + 0.0im),
        sample_shifts,
        3,
        2,
        1,
    )
    prompt_correlator = EarlyPromptLateCorrelator(
        SVector(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im),
        sample_shifts,
        3,
        2,
        1,
    )
    late_correlator = EarlyPromptLateCorrelator(
        SVector(0.25 + 0.0im, 0.75 + 0.0im, 0.75 + 0.0im),
        sample_shifts,
        3,
        2,
        1,
    )
    very_late_correlator = EarlyPromptLateCorrelator(
        SVector(0.0 + 0.0im, 0.5 + 0.0im, 1.0 + 0.0im),
        sample_shifts,
        3,
        2,
        1,
    )

    @test @inferred(get_early_late_sample_spacing(prompt_correlator)) == 4

    @test @inferred(dll_disc(gpsl1, very_early_correlator, delta)) == -0.5
    @test @inferred(dll_disc(gpsl1, early_correlator, delta)) == -0.25
    @test @inferred(dll_disc(gpsl1, prompt_correlator, delta)) == 0
    @test @inferred(dll_disc(gpsl1, late_correlator, delta)) == 0.25
    @test @inferred(dll_disc(gpsl1, very_late_correlator, delta)) == 0.5
end

end
