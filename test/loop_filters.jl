@testset "First order loop filter" begin
    loop_filter = @inferred Tracking.init_1st_order_loop_filter(1) 
    loop_filter, current_y = @inferred loop_filter(1, 2)
    @test current_y == 0.0
    loop_filter, current_y = @inferred loop_filter(2, 2)
    @test current_y == 4.0
    loop_filter, current_y = @inferred loop_filter(3, 2)
    @test current_y == 8.0
end

@testset "Second order loop filter" begin
    loop_filter = @inferred Tracking.init_2nd_order_boxcar_loop_filter(2 / 1.89) 
    loop_filter, current_y = @inferred loop_filter(1, 2)
    @test current_y == 2 * sqrt(2)
    loop_filter, current_y = @inferred loop_filter(2, 2)
    @test current_y == 8.0 + 4 * sqrt(2)
    loop_filter, current_y = @inferred loop_filter(3, 2)
    @test current_y == 24.0 + 6 * sqrt(2)
end

@testset "Third order loop filter" begin
    loop_filter = @inferred Tracking.init_3rd_order_boxcar_loop_filter(2 / 1.2) 
    loop_filter, current_y = @inferred loop_filter(1.0, 2)
    @test current_y == 4.8
    loop_filter, current_y = @inferred loop_filter(2, 2)
    @test current_y == 8.8 + 2 * 4.8
    loop_filter, current_y = @inferred loop_filter(3, 2)
    @test current_y == 58.4 + 3 * 4.8
end