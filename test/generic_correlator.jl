

@testset "Generic Correlator" begin

@testset "Empty constructor" begin
    correlator = @inferred GenericCorrelator()
    @test @inferred(get_early(correlator)) == 0.0
    @test @inferred(get_prompt(correlator)) == 0.0
    @test @inferred(get_late(correlator)) == 0.0
    @test length(get_correlators(correlator)) == 3
    @test @inferred(get_num_correlators(correlator)) == 3
    @test @inferred(get_early_index(correlator)) == 3
    @test @inferred(get_prompt_index(correlator)) == 2
    @test @inferred(get_late_index(correlator)) == 1
end

@testset "Single antenna" begin
    correlator = @inferred GenericCorrelator(num_ants = NumAnts(1))
    @test @inferred(get_early(correlator)) == 0.0
    @test @inferred(get_prompt(correlator)) == 0.0
    @test @inferred(get_late(correlator)) == 0.0
    @test length(get_correlators(correlator)) == 3
    @test @inferred(get_num_correlators(correlator)) == 3
    @test @inferred(get_early_index(correlator)) == 3
    @test @inferred(get_prompt_index(correlator)) == 2
    @test @inferred(get_late_index(correlator)) == 1
end

@testset "Multiple Antennas" begin
    correlator = @inferred GenericCorrelator(num_ants = NumAnts(2))
    @test @inferred(get_early(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    @test @inferred(get_prompt(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    @test @inferred(get_late(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    @test length(get_correlators(correlator)) == 3
    @test @inferred(get_num_correlators(correlator)) == 3
    @test @inferred(get_early_index(correlator)) == 3
    @test @inferred(get_prompt_index(correlator)) == 2
    @test @inferred(get_late_index(correlator)) == 1
end

@testset "Multiple antennas, multiple taps" begin
    correlator = @inferred GenericCorrelator(
        num_ants = NumAnts(2),
        num_correlators = NumCorrelators(5)
    )
    for k in -2:2
        @test @inferred(get_correlator(correlator, k)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    end
    @test length(get_correlators(correlator)) == 5
    @test @inferred(get_num_correlators(correlator)) == 5
    @test @inferred(get_early_index(correlator)) == 4
    @test @inferred(get_prompt_index(correlator)) == 3
    @test @inferred(get_late_index(correlator)) == 2
end

@testset "Multiple antennas, multiple taps, custom spacing" begin
    correlator = @inferred GenericCorrelator(
        num_ants = NumAnts(2),
        num_correlators = NumCorrelators(7),
        early_late_index_offset = 4
    )
    for k in -3:3
        @test @inferred(get_correlator(correlator, k)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    end
    @test length(get_correlators(correlator)) == 7
    @test @inferred(get_num_correlators(correlator)) == 7
    @test @inferred(get_early_index(correlator)) == 6
    @test @inferred(get_prompt_index(correlator)) == 4
    @test @inferred(get_late_index(correlator)) == 2
end

@testset "Default constructor" begin
    correlator = @inferred GenericCorrelator(
        [SVector(i + 1.0i * 1im, i + 2.0im) for i = 1:7],
        1, 4, 7
    )
    @test @inferred(get_early(correlator)) == SVector(1 + 1.0im, 1 + 2.0im)
    @test @inferred(get_prompt(correlator)) == SVector(4 + 4.0im, 4 + 2.0im)
    @test @inferred(get_late(correlator)) == SVector(7 + 7.0im, 7 + 2.0im)
end

@testset "Correlator zeroing" begin
    correlator = @inferred GenericCorrelator(
        [SVector(i + 1im, 2i + 3im) for i = 1:3],
        3, 2, 1
    )
    zero_correlator = @inferred Tracking.zero(correlator)
    @test all([all(t .== 0) for t in get_correlators(zero_correlator)])
end

@testset "Correlator filtering" begin
    correlator = @inferred GenericCorrelator(
        [SVector(1.0 + 0im, 1.0 + 0im) for i = 1:3],
        1, 2, 3
    )
    filtered_correlator = @inferred Tracking.filter(x -> x[1], correlator)
    @test get_correlators(filtered_correlator) == [1.0 + 0im for i = 1:3]
    filtered_correlator = @inferred Tracking.filter(x -> x, correlator)
    @test get_correlators(filtered_correlator) == [SVector(1.0 + 0im, 1.0 + 0im) for i = 1:3]
end

@testset "Sampleshift calculation" begin
    gpsl1 = GPSL1()
    correlator = @inferred GenericCorrelator(num_ants = NumAnts(2), num_correlators = NumCorrelators(7))
    correlator_sample_shifts = @inferred get_correlator_sample_shifts(
        gpsl1,
        correlator,
        4e6Hz,
        0.5
    )
    @test correlator_sample_shifts ≈ SVector(-3, -2, -1, 0, 1, 2, 3) .* round(0.5 * 4e6 / 1.023e6)
    @test typeof(correlator_sample_shifts) <: SVector
end

@testset "Correlator normalization" begin
    correlator = @inferred GenericCorrelator(
        [SVector(10i + 0im, 20i + 0im) for i = 1:3],
        1, 2, 3
    )
    normalized_correlator = @inferred Tracking.normalize(correlator, 10)
    @test get_correlators(normalized_correlator) == [SVector(1i + 0im, 2i + 0im) for i = 1:3]
end

@testset "Single antenna correlation" begin
    gpsl1 = GPSL1()
    signal = StructArray{Complex{Float32}}(
        (Float32.(get_code.(gpsl1, (1:2500) * 1023e3 / 2.5e6, 1)),
        zeros(Float32, 2500))
    )
    correlator = GenericCorrelator(num_ants = NumAnts(1), num_correlators = NumCorrelators(5))
    correlator_sample_shifts = get_correlator_sample_shifts(
        gpsl1,
        correlator,
        2.5e6Hz,
        0.5
    )
    code = get_code.(
        gpsl1,
        (1 + correlator_sample_shifts[1]:2500 + correlator_sample_shifts[end]) * 1023e3 / 2.5e6,
        1
    )
    correlator_result = Tracking.correlate(
        correlator,
        signal,
        code,
        correlator_sample_shifts,
        1,
        2500
    )
    @test get_prompt(correlator_result) == 2500
    @test get_early(correlator_result) == get_late(correlator_result)
    @test get_correlators(correlator_result) == reverse(get_correlators(correlator_result))
end

@testset "Multi antenna correlation" begin
    gpsl1 = GPSL1()
    signal = StructArray{Complex{Float32}}(
        (Float32.(get_code.(gpsl1, (1:2500) * 1023e3 / 2.5e6, 1)),
        zeros(Float32, 2500))
    )
    correlator = GenericCorrelator(num_ants = NumAnts(3), num_correlators = NumCorrelators(5))
    correlator_sample_shifts = get_correlator_sample_shifts(
        gpsl1,
        correlator,
        2.5e6Hz,
        0.5
    )
    code = get_code.(
        gpsl1,
        (1 + correlator_sample_shifts[1]:2500 + correlator_sample_shifts[end]) * 1023e3 / 2.5e6,
        1
    )
    signal_mat = repeat(signal, outer = (1,3))
    correlator_result = Tracking.correlate(
        correlator,
        signal_mat,
        code,
        correlator_sample_shifts,
        1,
        2500,
    )
    @test all(get_prompt(correlator_result) .== 2500)
    @test all(get_early(correlator_result) .== 1476)
    @test get_correlators(correlator) == reverse(get_correlators(correlator))
end

end
