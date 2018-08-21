@testset "Code replica" begin
    gen_sampled_code, get_code_phase = init_gpsl1_codes()

    gen_code_replica = @inferred Tracking.init_code_replica(1023e3, 0.0, 4e6, 1, gen_sampled_code, get_code_phase)
    next_gen_code_replica, code_replicas = @inferred gen_code_replica(4000, 10)
    next1_gen_code_replica, next_code_replicas = @inferred next_gen_code_replica(4000, 10)

    @test [prompt(code_replicas); prompt(next_code_replicas)] == gen_sampled_code(1:8000, 1023e3 + 10, 0.0, 4e6, 1)
    @test Tracking.late(code_replicas) == gen_sampled_code(-1:3998, 1023e3 + 10, 0.0, 4e6, 1)
    @test Tracking.early(code_replicas) == gen_sampled_code(3:4002, 1023e3 + 10, 0.0, 4e6, 1)
end

@testset "Carrier replica" begin
    gen_carrier_replica = @inferred Tracking.init_carrier_replica(100, 1π, 4e6)
    next_gen_carrier_replica, carrier_replica = @inferred gen_carrier_replica(400, 50)
    next1_gen_carrier_replica, next_carrier_replica = @inferred next_gen_carrier_replica(400, 50)

    @test [carrier_replica; next_carrier_replica] ≈ cis.(2π * 150 / 4e6 * (1:800) + 1π)
end