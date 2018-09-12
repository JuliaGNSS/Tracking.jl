@testset "First order loop filter" begin
    loop_filter = @inferred Tracking.init_1st_order_loop_filter(1Hz) 
    loop_filter, current_y = @inferred loop_filter(1, 2s)
    @test current_y == 0.0Hz
    loop_filter, current_y = @inferred loop_filter(2, 2s)
    @test current_y == 4.0Hz
    loop_filter, current_y = @inferred loop_filter(3, 2s)
    @test current_y == 8.0Hz
end

@testset "Second order loop filter" begin
    loop_filter = @inferred Tracking.init_2nd_order_boxcar_loop_filter(2Hz / 1.89) 
    loop_filter, current_y = @inferred loop_filter(1, 2s)
    @test current_y == 2Hz * sqrt(2)
    loop_filter, current_y = @inferred loop_filter(2, 2s)
    @test current_y == 8.0Hz + 4Hz * sqrt(2)
    loop_filter, current_y = @inferred loop_filter(3, 2s)
    @test current_y == 24.0Hz + 6Hz * sqrt(2)
end

@testset "Third order loop filter" begin
    loop_filter = Tracking.init_3rd_order_boxcar_loop_filter(2Hz / 1.2) # @inferred
    loop_filter, current_y = loop_filter(1.0, 2s)
    @test current_y == 4.8Hz
    loop_filter, current_y = loop_filter(2, 2s)
    @test current_y == 8.8Hz + 2Hz * 4.8
    loop_filter, current_y = loop_filter(3, 2s)
    @test current_y == 58.4Hz + 3Hz * 4.8
end