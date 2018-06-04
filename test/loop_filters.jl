@testset "Loop_filter first order" begin
    loop_filter = Tracking.init_1st_order_loop_filter(1) 
    loop_filter, current_y = loop_filter(1)
    @test @inferred current_y == 0.0
    loop_filter, current_y = loop_filter(2)
    @test @inferred current_y == 4.0
    loop_filter, current_y = loop_filter(3)
    @test @inferred current_y == 8.0
end

@testset "Loop_filter second order" begin
    loop_filter = Tracking.init_2nd_order_loop_filter(2 / 1.89, 2) 
    loop_filter, current_y = loop_filter(1)
    @test @inferred current_y == 2 * sqrt(2)
    loop_filter, current_y = loop_filter(2)
    @test @inferred current_y == 8.0 + 4 * sqrt(2)
    loop_filter, current_y = loop_filter(3)
    @test @inferred current_y == 24.0 + 6 * sqrt(2)
end

@testset "Loop_filter third order" begin
    loop_filter = Tracking.init_3rd_order_loop_filter(2 / 1.2, 2) 
    loop_filter, current_y = loop_filter(1.0)
    @test @inferred current_y == 4.8
    loop_filter, current_y = loop_filter(2)
    @test @inferred current_y == 8.8 + 2 * 4.8
    loop_filter, current_y = loop_filter(3)
    @test @inferred current_y == 58.4 + 3 * 4.8
end