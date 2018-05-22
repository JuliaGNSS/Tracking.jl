@Testset CalculateOutputErrors begin
    testSignal = [];
    @Testset PLL_Test begin
        
        @test Tracking.dll_disc(testSignal) = (abs(testSignal[1]) - abs(testSignal[3])) / (abs(testSignal[1]) + abs(testSignal[3]))
    
    end
    @Testset DLL_Test begin
        
        @test Tracking.pll_disc(testSignal) = atan(imag(testSignal[2], real(testSignal[2])))
    end
end

