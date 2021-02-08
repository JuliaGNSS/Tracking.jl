@testset "Correlator" begin

    @testset "Correlator constructor" begin
        gpsl1 = GPSL1()
        correlator = @inferred EarlyPromptLateCorrelator()
        correlator_sample_shifts = SVector(-1, 0, 1)
        early_late_index_shift = 1
        @test @inferred(get_early(correlator, correlator_sample_shifts, early_late_index_shift)) == 0.0
        @test @inferred(get_prompt(correlator, correlator_sample_shifts)) == 0.0
        @test @inferred(get_late(correlator, correlator_sample_shifts, early_late_index_shift)) == 0.0

        correlator = @inferred EarlyPromptLateCorrelator(NumAnts(1))

        @test @inferred(get_early(correlator, correlator_sample_shifts, early_late_index_shift)) == 0.0
        @test @inferred(get_prompt(correlator, correlator_sample_shifts)) == 0.0
        @test @inferred(get_late(correlator, correlator_sample_shifts, early_late_index_shift)) == 0.0

        correlator = @inferred EarlyPromptLateCorrelator(NumAnts(2))

        @test @inferred(get_early(correlator, correlator_sample_shifts, early_late_index_shift)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_prompt(correlator, correlator_sample_shifts)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_late(correlator, correlator_sample_shifts, early_late_index_shift)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
    end
    @testset "Zeroing correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im)
            )
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            SVector(
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
                SVector(0.0 + 0.0im, 0.0 + 0.0im)
            )
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                1.0 + 0.0im,
                1.0 + 0.0im,
                1.0 + 0.0im
            )
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            SVector(
                0.0 + 0.0im,
                0.0 + 0.0im,
                0.0 + 0.0im
            )
        )
    end
    @testset "Filter correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im)
            )
        )
        filtered_correlator = @inferred Tracking.filter(x -> x[1], correlator)
        @test filtered_correlator == EarlyPromptLateCorrelator(
            SVector(
                1.0 + 0.0im,
                1.0 + 0.0im,
                1.0 + 0.0im
            )
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                1.0 + 0.0im,
                1.0 + 0.0im,
                1.0 + 0.0im
            )
        )
        filtered_correlator = @inferred Tracking.filter(x -> 2 * x, correlator)
        @test filtered_correlator == EarlyPromptLateCorrelator(
            SVector(
                2.0 + 0.0im,
                2.0 + 0.0im,
                2.0 + 0.0im
            )
        )
    end
    @testset "Correlator sample shifts" begin
        gpsl1 = GPSL1()
        correlator = EarlyPromptLateCorrelator()
        correlator_sample_shifts = @inferred get_correlator_sample_shifts(
            gpsl1,
            correlator,
            4e6Hz,
            0.5
        )
        @test correlator_sample_shifts ≈ round(0.5 * 4e6 / 1023e3) * SVector(-1, 0, 1)
        @test typeof(correlator_sample_shifts) <: SVector{3,<:Integer}
    end
    @testset "Early late index shift" begin
        gpsl1 = GPSL1()
        correlator = EarlyPromptLateCorrelator()
        correlator_sample_shifts = get_correlator_sample_shifts(
            gpsl1,
            correlator,
            4e6Hz,
            0.5
        )
        early_late_index_shift = @inferred get_early_late_index_shift(
           gpsl1,
           correlator_sample_shifts,
           correlator,
           4e6Hz,
           0.5
        )

        @test early_late_index_shift == 1
    end
    @testset "Early late sample spacing" begin
        @test get_early_late_sample_spacing(SVector(-2,0,2), 1) == 4
    end
    @testset "Normalize correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                1.0 + 0.0im,
                1.0 + 0.0im,
                1.0 + 0.0im
            )
        )
        normalized_correlator = @inferred Tracking.normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            SVector(
                0.1 + 0.0im,
                0.1 + 0.0im,
                0.1 + 0.0im
            )
        )

        correlator = EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im)
            )
        )
        normalized_correlator = @inferred Tracking.normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            SVector(
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
                SVector(0.1 + 0.0im, 0.1 + 0.0im)
            )
        )
    end
    @testset "Correlate" begin
        gpsl1 = GPSL1()
        signal = StructArray{Complex{Float32}}(
            (Float32.(get_code.(gpsl1, (1:2500) * 1023e3 / 2.5e6, 1)),
            zeros(Float32, 2500))
        )
        correlator_sample_shifts = SVector(-1, 0, 1)
        code = get_code.(
            gpsl1,
            (1 + correlator_sample_shifts[1]:2500 + correlator_sample_shifts[end]) * 1023e3 / 2.5e6,
            1
        )
        correlator = EarlyPromptLateCorrelator(SVector(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im))

        correlator_result = Tracking.correlate(
            correlator,
            signal,
            code,
            correlator_sample_shifts,
            1,
            2500
        )
        @test get_early(correlator_result, correlator_sample_shifts, 1) ==
            get_late(correlator_result, correlator_sample_shifts, 1)
        @test get_prompt(correlator_result, correlator_sample_shifts) == 2500


        signal_mat = repeat(signal, outer = (1,3))
        correlator_sample_shifts = SVector(-1,0,1)
        correlator = EarlyPromptLateCorrelator(NumAnts(3))

        correlator_result = Tracking.correlate(
            correlator,
            signal_mat,
            code,
            correlator_sample_shifts,
            1,
            2500,
        )
        @test get_early(correlator_result, correlator_sample_shifts, 1) == get_late(correlator_result, correlator_sample_shifts, 1)
        @test all(get_early(correlator_result, correlator_sample_shifts, 1) .== 1476)
        @test all(get_prompt(correlator_result, correlator_sample_shifts) .== 2500)
    end

end
