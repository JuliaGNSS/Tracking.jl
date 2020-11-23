

@testset "Generic Correlator" begin

@testset "Empty constructor" begin
    correlator = @inferred GenericCorrelator()
    @test @inferred(get_early(correlator)) == 0.0
    @test @inferred(get_prompt(correlator)) == 0.0
    @test @inferred(get_late(correlator)) == 0.0
    @test length(get_taps(correlator)) == 3
    @test @inferred(get_num_taps(correlator)) == 3
    @test @inferred(get_early_index(correlator)) == 3
    @test @inferred(get_prompt_index(correlator)) == 2
    @test @inferred(get_late_index(correlator)) == 1
end

@testset "Single antenna" begin
    correlator = @inferred GenericCorrelator(NumAnts(1))
    @test @inferred(get_early(correlator)) == 0.0
    @test @inferred(get_prompt(correlator)) == 0.0
    @test @inferred(get_late(correlator)) == 0.0
    @test length(get_taps(correlator)) == 3
    @test @inferred(get_num_taps(correlator)) == 3
    @test @inferred(get_early_index(correlator)) == 3
    @test @inferred(get_prompt_index(correlator)) == 2
    @test @inferred(get_late_index(correlator)) == 1
end

@testset "Multiple Antennas" begin
    correlator = @inferred GenericCorrelator(NumAnts(2))
    @test @inferred(get_early(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    @test @inferred(get_prompt(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    @test @inferred(get_late(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    @test length(get_taps(correlator)) == 3
    @test @inferred(get_num_taps(correlator)) == 3
    @test @inferred(get_early_index(correlator)) == 3
    @test @inferred(get_prompt_index(correlator)) == 2
    @test @inferred(get_late_index(correlator)) == 1
end

@testset "Multiple antennas, multiple taps" begin
    correlator = @inferred GenericCorrelator(NumAnts(2), NumTaps(5))
    for k in -2:2
        @test @inferred(get_tap(correlator, k)) == SVector(0.0+0.0im, 0.0+0.0im)
    end
    @test length(get_taps(correlator)) == 5
    @test @inferred(get_num_taps(correlator)) == 5
    @test @inferred(get_early_index(correlator)) == 4
    @test @inferred(get_prompt_index(correlator)) == 3
    @test @inferred(get_late_index(correlator)) == 2
end

@testset "Multiple antennas, multiple taps, custom spacing" begin
    correlator = @inferred GenericCorrelator(NumAnts(2), NumTaps(7), 4)
    for k in -3:3
        @test @inferred(get_tap(correlator, k)) == SVector(0.0+0.0im, 0.0+0.0im)
    end
    @test length(get_taps(correlator)) == 7
    @test @inferred(get_num_taps(correlator)) == 7
    @test @inferred(get_early_index(correlator)) == 6
    @test @inferred(get_prompt_index(correlator)) == 4
    @test @inferred(get_late_index(correlator)) == 2
end

@testset "Default constructor" begin
    correlator = @inferred GenericCorrelator(
        [SVector(i + 1.0i*im, i + 2.0im) for i = 1:7],
        NumTaps(7),
        1, 4, 7
    )
    @test @inferred(get_early(correlator)) == SVector(1+1.0im, 1+2.0im)
    @test @inferred(get_prompt(correlator)) == SVector(4+4.0im, 4+2.0im)
    @test @inferred(get_late(correlator)) == SVector(7+7.0im, 7+2.0im)
end

@testset "Correlator zeroing" begin
    correlator = @inferred GenericCorrelator(
        [SVector(i+1im, 2i+3im) for i=1:3],
        NumTaps(3),
        3,2,1
    )
    zero_correlator = @inferred Tracking.zero(correlator)
    @test all([all(t.==0) for t in get_taps(zero_correlator)])
end

@testset "Correlator filtering" begin
    correlator = @inferred GenericCorrelator(
        [SVector(1.0+0im, 1.0+0im) for i = 1:3],
        NumTaps(3),
        1, 2, 3
    )
    filtered_correlator = @inferred Tracking.filter(x->x[1], correlator)
    @test get_taps(filtered_correlator) == [1.0+0im for i = 1:3]
    filtered_correlator = @inferred Tracking.filter(x->x, correlator)
    @test get_taps(filtered_correlator) == [SVector(1.0+0im, 1.0+0im) for i = 1:3]
end

@testset "Sampleshift calculation" begin
    correlator = @inferred GenericCorrelator(NumAnts(2), NumTaps(7))
    correlator_sample_shifts = @inferred get_correlator_sample_shifts(
        GPSL1,
        correlator,
        4e6Hz,
        0.5
    )
    @test correlator_sample_shifts ≈ SVector(-3,-2,-1,0,1,2,3).*round(0.5*4e6/1.023e6)
    @test typeof(correlator_sample_shifts) <: SVector
end

@testset "Correlator normalization" begin
    correlator = @inferred GenericCorrelator(
        [SVector(10i+0im, 20i+0im) for i = 1:3],
        NumTaps(3),
        1, 2, 3
    )
    normalized_correlator = @inferred Tracking.normalize(correlator, 10)
    @test get_taps(normalized_correlator) == [SVector(1i+0im, 2i+0im) for i = 1:3]    
end

@testset "Singe antenna correlation" begin
    signal = StructArray{Complex{Int16}}(
        (get_code.(GPSL1, (1:2500) * 1023e3 / 2.5e6, 1) * Int16(1) << (7 + 2),
        zeros(Int16, 2500))
    )
    early_late_sample_shift = 1
    code = get_code.(
        GPSL1,
        (1 - 2early_late_sample_shift:2500+2early_late_sample_shift) * 1023e3 / 2.5e6,
        1
    )
    correlator = GenericCorrelator(NumAnts(1), NumTaps(5))
    correlator_result = Tracking.correlate(
        correlator,
        signal,
        code,
        early_late_sample_shift,
        1,
        2500,
        1.0,
        2,
        Val(7)
    )
    @test get_prompt(correlator_result) == 2500
    @test get_early(correlator_result) == get_late(correlator_result)
    @test get_taps(correlator_result) == reverse(get_taps(correlator_result))
end

@testset "Multi antenna correlation" begin
    signal = StructArray{Complex{Int16}}(
        (get_code.(GPSL1, (1:2500) * 1023e3 / 2.5e6, 1) * Int16(1) << (7 + 2),
        zeros(Int16, 2500))
    )
    early_late_sample_shift = 1
    code = get_code.(
        GPSL1,
        (1 - 2early_late_sample_shift:2500+2early_late_sample_shift) * 1023e3 / 2.5e6,
        1
    )
    signal_mat = repeat(signal, outer = (1,3))
    early_late_sample_shift = 1
    correlator = GenericCorrelator(NumAnts(3), NumTaps(5))
    
    correlator_result = Tracking.correlate(
        correlator,
        signal_mat,
        code,
        early_late_sample_shift,
        1,
        2500,
        ones(Float64, 3),
        2,
        Val(7)
    )
    @test all(get_prompt(correlator_result) .== 2500)
    @test all(get_early(correlator_result) .== 1476)
    @test get_taps(correlator) == reverse(get_taps(correlator))
end

end