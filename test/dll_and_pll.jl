@testset "initialize PLL" begin
    locked_loop, sample_code, phase = Tracking.init_PLL(20.0, 50, 4000, 4e6, 18.0, 1e-3)
    @test sample_code[2000] == 0.2602417925505148 + 0.9655434787776751im
    @test phase == 1.4646033438202217
    @test locked_loop(sample_code)[3] == 1.4646033438202217

    
end
@testset "inizialize DLL" begin
    locked_loop, sample_code, phase = Tracking.init_DLL(20, 1023e3, 4000, 4e6, 1, 1e-3 1)
    @test phase == 20.0
    @test sample_code[1][90:100] == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    @test locked_loop(sample_code[2])[3] == 20.0;
end
