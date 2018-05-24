@testset "Loop_filter first order" begin
    new_function = Tracking.init_1st_order_loop_filter(1) 
    new_function, current_y = new_function(1)
    @test current_y[1] == 0.0
    new_function, current_y = new_function(2)
    @test current_y[1] == 4.0
    new_function, current_y = new_function(3)
    @test current_y[1] == 8.0
end

@testset "Loop_filter second order" begin
    new_function = Tracking.init_2nd_order_loop_filter(2/1.89,2) 
    new_function, current_y = new_function(1)
    @test current_y[1] == 2*sqrt(2)
    new_function, current_y = new_function(2)
    @test current_y[1] == 8.0+4*sqrt(2)
    new_function, current_y = new_function(3)
    @test current_y[1] == 24.0+6*sqrt(2)
end

@testset "Loop_filter third order" begin
    new_function = Tracking.init_3rd_order_loop_filter(2/1.2,2) 
    new_function, current_y = new_function(1.0)
    @test current_y[1] == 4.8
    new_function, current_y = new_function(2)
    @test current_y[1] == 8.8+2*4.8
    new_function, current_y = new_function(3)
    @test current_y[1] == 58.4+3*4.8
end