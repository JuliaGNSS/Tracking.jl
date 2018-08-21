
@testset "PLL discriminator" begin
    test_signal_prephase = [-0.5 + sqrt(3) / 2im, -1 + sqrt(3) * 1im, -0.5 + sqrt(3) / 2im]
    test_signal = [0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im]
    test_signal_postphase = [0.5 + sqrt(3) / 2im, 1 + sqrt(3) * 1im, 0.5 + sqrt(3) / 2im]
    
    @test @inferred(Tracking.pll_disc(test_signal_prephase)) == -π / 3  #-60°
    @test @inferred(Tracking.pll_disc(test_signal)) == 0
    @test @inferred(Tracking.pll_disc(test_signal_postphase)) == π / 3  #+60°
end


@testset "DLL discriminator" begin
    test_signal_very_early = [0. + 0.0im, 0.5 + 0.0im, 1.0 + 0.0im]
    test_signal_early = [0.25 + 0.0im, 0.75 + 0.0im, 0.75 + 0.0im]
    test_signal_in_time = [0.5 + 0.0im, 1 + 0.0im, 0.5 + 0.0im]
    test_signal_late = [0.75 + 0.0im, 0.75 + 0.0im, 0.25 + 0.0im]
    test_signal_very_late = [1.0 + 0.0im, 0.5 + 0.0im, 0.0 + 0.0im]
    
    @test @inferred(Tracking.dll_disc(test_signal_very_early)) == 1
    @test @inferred(Tracking.dll_disc(test_signal_early)) == 0.5
    @test @inferred(Tracking.dll_disc(test_signal_in_time)) == 0
    @test @inferred(Tracking.dll_disc(test_signal_late)) == -0.5
    @test @inferred(Tracking.dll_disc(test_signal_very_late)) == -1
end
