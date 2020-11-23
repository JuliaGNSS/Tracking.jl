@testset "PLL discriminator" begin

    correlator_minus60off = EarlyPromptLateCorrelator(
        -0.5 + sqrt(3) / 2im,
        -1 + sqrt(3) * 1im,
        -0.5 + sqrt(3) / 2im
    )
    correlator_0off = EarlyPromptLateCorrelator(
        0.5 + 0.0im,
        1 + 0.0im,
        0.5 + 0.0im
    )
    correlator_plus60off = EarlyPromptLateCorrelator(
        0.5 + sqrt(3) / 2im,
        1 + sqrt(3) * 1im,
        0.5 + sqrt(3) / 2im
    )

    @test @inferred(Tracking.pll_disc(GPSL1, correlator_minus60off)) == -π / 3  #-60°
    @test @inferred(Tracking.pll_disc(GPSL1, correlator_0off)) == 0
    @test @inferred(Tracking.pll_disc(GPSL1, correlator_plus60off)) == π / 3  #+60°
end


@testset "DLL discriminator" begin

    very_early_correlator = EarlyPromptLateCorrelator(0.0 + 0.0im, 0.5 + 0.0im, 1.0 + 0.0im)
    early_correlator = EarlyPromptLateCorrelator(0.25 + 0.0im, 0.75 + 0.0im, 0.75 + 0.0im)
    prompt_correlator = EarlyPromptLateCorrelator(0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im)
    late_correlator = EarlyPromptLateCorrelator(0.75 + 0.0im, 0.75 + 0.0im, 0.25 + 0.0im)
    very_late_correlator = EarlyPromptLateCorrelator(1.0 + 0.0im, 0.5 + 0.0im, 0.0 + 0.0im)

    sample_shift = SVector(-2, 0, 2)
    delta = 0.25

    @test @inferred(
        Tracking.dll_disc(GPSL1, very_early_correlator, sample_shift, delta)
    ) == -0.5
    @test @inferred(
        Tracking.dll_disc(GPSL1, early_correlator, sample_shift, delta)
    ) == -0.25
    @test @inferred(
        Tracking.dll_disc(GPSL1, prompt_correlator, sample_shift, delta)
    ) == 0
    @test @inferred(
        Tracking.dll_disc(GPSL1, late_correlator, sample_shift, delta)
    ) == 0.25
    @test @inferred(
        Tracking.dll_disc(GPSL1, very_late_correlator, sample_shift, delta)
    ) == 0.5
end
