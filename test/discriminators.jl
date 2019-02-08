@testset "Prompt early late" begin
    a = @SVector [i for i = 1:3]
    @test Tracking.early(a) == a[3]
    @test prompt(a) == a[2]
    @test Tracking.late(a) == a[1]

    a = @SVector [i for i = 1:5]
    @test Tracking.veryearly(a) == a[5]
    @test Tracking.early(a) == a[4]
    @test prompt(a) == a[3]
    @test Tracking.late(a) == a[2]
    @test Tracking.verylate(a) == a[1]
end

@testset "PLL discriminator" begin
    test_signal_prephase = @SVector [-0.5 + sqrt(3) / 2im, -1 + sqrt(3) * 1im, -0.5 + sqrt(3) / 2im]
    test_signal = @SVector [0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im]
    test_signal_postphase = @SVector [0.5 + sqrt(3) / 2im, 1 + sqrt(3) * 1im, 0.5 + sqrt(3) / 2im]

    @test @inferred(Tracking.pll_disc(test_signal_prephase)) == -π / 3  #-60°
    @test @inferred(Tracking.pll_disc(test_signal)) == 0
    @test @inferred(Tracking.pll_disc(test_signal_postphase)) == π / 3  #+60°
end


@testset "DLL discriminator" begin
    test_signal_very_early = @SVector [0. + 0.0im, 0.5 + 0.0im, 1.0 + 0.0im] #τ = 0.5
    test_signal_early = @SVector [0.25 + 0.0im, 0.75 + 0.0im, 0.75 + 0.0im] #τ = 0.25
    test_signal_in_time = @SVector [0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im] #τ = 0
    test_signal_late = @SVector [0.75 + 0.0im, 0.75 + 0.0im, 0.25 + 0.0im] #τ = -0.25
    test_signal_very_late = @SVector [1.0 + 0.0im, 0.5 + 0.0im, 0.0 + 0.0im] #τ = -0.5

    @test @inferred(Tracking.dll_disc(test_signal_very_early)) == 0.5
    @test @inferred(Tracking.dll_disc(test_signal_early)) == 0.25
    @test @inferred(Tracking.dll_disc(test_signal_in_time)) == 0
    @test @inferred(Tracking.dll_disc(test_signal_late)) == -0.25
    @test @inferred(Tracking.dll_disc(test_signal_very_late)) == -0.5
end
