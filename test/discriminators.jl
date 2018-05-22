
# use the same sampelingrate and the same timestamps and just shift the values
testSignal = [0,0.25+0.5im,0.5+1im,0.25+0,5im,0];
testSignal_early = [0,0.15+0.3im,0.4+0.8im,0.35+0.7im,0];
testSignal_late = [0,0.35+0,7im,0.4+0.8im,0.15+0.3im,0];



@Testset "PLL discriminator test" begin
    @test Tracking.pll_disc(testSignal) == 0.4636476090008061

end


@Testset "DLL discriminator test" begin
    @test Tracking.pll_disc(testSignal) == 0
    @test Tracking.pll_disc(testSignal_early) == -0.39999999999999997
    @test Tracking.pll_disc(testSignal_late) == -0.39999999999999997
end
