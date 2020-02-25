@testset "Correlator" begin

    @testset "Early prompt late correlator" begin

        correlator = @inferred EarlyPromptLateCorrelator()

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0

        correlator = @inferred EarlyPromptLateCorrelator(NumAnts(1))

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0

        correlator = @inferred EarlyPromptLateCorrelator(NumAnts(2))

        @test @inferred(get_early(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_prompt(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_late(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im),
            SVector(1.0 + 0.0im, 1.0 + 0.0im),
            SVector(1.0 + 0.0im, 1.0 + 0.0im)
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            SVector(0.0 + 0.0im, 0.0 + 0.0im),
            SVector(0.0 + 0.0im, 0.0 + 0.0im),
            SVector(0.0 + 0.0im, 0.0 + 0.0im)
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            1.0 + 0.0im,
            1.0 + 0.0im,
            1.0 + 0.0im
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            0.0 + 0.0im,
            0.0 + 0.0im,
            0.0 + 0.0im
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im),
            SVector(1.0 + 0.0im, 1.0 + 0.0im),
            SVector(1.0 + 0.0im, 1.0 + 0.0im)
        )
        filtered_correlator = @inferred Tracking.filter(x -> x[1], correlator)
        @test filtered_correlator == EarlyPromptLateCorrelator(
            1.0 + 0.0im,
            1.0 + 0.0im,
            1.0 + 0.0im
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            1.0 + 0.0im,
            1.0 + 0.0im,
            1.0 + 0.0im
        )
        filtered_correlator = @inferred Tracking.filter(x -> x, correlator)
        @test filtered_correlator == EarlyPromptLateCorrelator(
            1.0 + 0.0im,
            1.0 + 0.0im,
            1.0 + 0.0im
        )

        early_late_sample_shift = @inferred get_early_late_sample_shift(
            GPSL1,
            correlator,
            4e6Hz,
            0.5
        )
        @test early_late_sample_shift ≈ round(0.5 * 4e6 / 1023e3)
        @test typeof(early_late_sample_shift) <: Integer

        normalized_correlator = @inferred Tracking.normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            0.1 + 0.0im,
            0.1 + 0.0im,
            0.1 + 0.0im
        )

        correlator = EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im),
            SVector(1.0 + 0.0im, 1.0 + 0.0im),
            SVector(1.0 + 0.0im, 1.0 + 0.0im)
        )
        normalized_correlator = @inferred Tracking.normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            SVector(0.1 + 0.0im, 0.1 + 0.0im),
            SVector(0.1 + 0.0im, 0.1 + 0.0im),
            SVector(0.1 + 0.0im, 0.1 + 0.0im)
        )

        signal = StructArray{Complex{Int16}}(
            (get_code.(GPSL1, (1:2500) * 1023e3 / 2.5e6, 1) * Int16(1) << (7 + 2),
            zeros(Int16, 2500))
        )
        early_late_sample_shift = 1
        code = get_code.(
            GPSL1,
            (1 - early_late_sample_shift:2500 + early_late_sample_shift) * 1023e3 / 2.5e6,
            1
        )
        correlator = EarlyPromptLateCorrelator(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)

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
        @test get_early(correlator_result) == get_late(correlator_result)
        @test get_prompt(correlator_result) == 2500


        signal_mat = repeat(signal, outer = (1,3))
        early_late_sample_shift = 1
        correlator = EarlyPromptLateCorrelator(NumAnts(3))

        correlator_result = Tracking.correlate(
            correlator,
            signal,
            code,
            early_late_sample_shift,
            1,
            2500,
            [1.0, 1.0, 1.0],
            2,
            Val(7)
        )
        @test get_early(correlator_result) == get_late(correlator_result)
        @test all(get_early(correlator_result) .== 1476)
        @test all(get_prompt(correlator_result) .== 2500)
    end

end
