@testset "Correlator" begin
    @testset "Correlator constructor" begin
        gpsl1 = GPSL1()
        sampling_frequency = 5e6Hz
        correlator = @inferred EarlyPromptLateCorrelator(gpsl1, sampling_frequency)

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0
        @test correlator.shifts == -2:2:2
        @test correlator.early_shift == 3
        @test correlator.prompt_shift == 2
        @test correlator.late_shift == 1

        correlator = @inferred EarlyPromptLateCorrelator(
            gpsl1,
            sampling_frequency,
            num_ants = NumAnts(1),
            num_accumulators = 3,
        )
        @test isa(get_accumulators(correlator), Vector)
        @test correlator.shifts == -2:2:2
        @test correlator.early_shift == 3
        @test correlator.prompt_shift == 2
        @test correlator.late_shift == 1

        correlator = @inferred EarlyPromptLateCorrelator(
            gpsl1,
            sampling_frequency,
            num_ants = NumAnts(1),
        )

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0
        @test correlator.shifts == -2:2:2
        @test correlator.early_shift == 3
        @test correlator.prompt_shift == 2
        @test correlator.late_shift == 1

        correlator = @inferred EarlyPromptLateCorrelator(
            gpsl1,
            sampling_frequency,
            num_ants = NumAnts(1),
            num_accumulators = 3,
        )

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0
        @test correlator.shifts == -2:2:2
        @test correlator.early_shift == 3
        @test correlator.prompt_shift == 2
        @test correlator.late_shift == 1

        correlator = @inferred EarlyPromptLateCorrelator(
            gpsl1,
            sampling_frequency,
            num_ants = NumAnts(2),
        )

        @test @inferred(get_early(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_prompt(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_late(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test correlator.shifts == -2:2:2
        @test correlator.early_shift == 3
        @test correlator.prompt_shift == 2
        @test correlator.late_shift == 1

        correlator = @inferred EarlyPromptLateCorrelator(
            gpsl1,
            sampling_frequency,
            num_ants = NumAnts(2),
            num_accumulators = 3,
        )

        @test @inferred(get_early(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_prompt(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_late(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test correlator.shifts == -2:2:2
        @test correlator.early_shift == 3
        @test correlator.prompt_shift == 2
        @test correlator.late_shift == 1
    end
    @testset "Zeroing correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            -1:1,
            3,
            2,
            1,
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            SVector(
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
            ),
            -1:1,
            3,
            2,
            1,
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im),
            -1:1,
            3,
            2,
            1,
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            SVector(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im),
            -1:1,
            3,
            2,
            1,
        )
    end
    @testset "Filter correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 1.0im),
                SVector(1.0 + 0.0im, 2.0 + 0.0im),
                SVector(1.0 + 1.0im, 1.0 + 3.0im),
            ),
            -1:1,
            3,
            2,
            1,
        )
        default_post_corr_filter = Tracking.DefaultPostCorrFilter()
        filtered_correlator = @inferred Tracking.apply(default_post_corr_filter, correlator)
        @test filtered_correlator == EarlyPromptLateCorrelator(
            SVector(1.0 + 1.0im, 2.0 + 0.0im, 1.0 + 3.0im),
            -1:1,
            3,
            2,
            1,
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im),
            -1:1,
            3,
            2,
            1,
        )
        filtered_correlator = @inferred Tracking.apply(x -> 2 * x, correlator)
        @test filtered_correlator == EarlyPromptLateCorrelator(
            SVector(2.0 + 0.0im, 2.0 + 0.0im, 2.0 + 0.0im),
            -1:1,
            3,
            2,
            1,
        )
    end
    @testset "Early late sample spacing" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            -1:1,
            3,
            2,
            1,
        )
        @test get_early_late_sample_spacing(correlator) == 2

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            -2:2:2,
            3,
            2,
            1,
        )
        @test get_early_late_sample_spacing(correlator) == 4

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            -2:2,
            5,
            3,
            1,
        )
        @test get_early_late_sample_spacing(correlator) == 4
    end
    @testset "Normalize correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im),
            -1:1,
            3,
            2,
            1,
        )
        normalized_correlator = @inferred Tracking.normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            SVector(0.1 + 0.0im, 0.1 + 0.0im, 0.1 + 0.0im),
            -1:1,
            3,
            2,
            1,
        )

        correlator = EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            -1:1,
            3,
            2,
            1,
        )
        normalized_correlator = @inferred Tracking.normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            SVector(
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
            ),
            -1:1,
            3,
            2,
            1,
        )
    end

    @testset "Correlate" begin
        gpsl1 = GPSL1()
        sampling_frequency = 2.5e6Hz
        signal = StructArray{Complex{Float32}}((
            Float32.(get_code.(gpsl1, (1:2500) * 1023e3 / (sampling_frequency / Hz), 1)),
            zeros(Float32, 2500),
        ))
        code = get_code.(gpsl1, (1+-1:2500+1) * 1023e3 / 2.5e6, 1)

        @testset "↳ $na antennas" for na ∈ [1, 4]
            signal_mat = na == 1 ? signal : repeat(signal; outer = (1, na))

            @testset "↳ $(string(type)) accumulator" for type ∈ [:SVector, :Vector]
                if type == :SVector
                    correlator = EarlyPromptLateCorrelator(
                        gpsl1,
                        sampling_frequency;
                        num_ants = NumAnts(na),
                    )
                else
                    correlator = EarlyPromptLateCorrelator(
                        gpsl1,
                        sampling_frequency;
                        num_ants = NumAnts(na),
                        num_accumulators = 3,
                    )
                end

                correlator_result =
                    @inferred Tracking.correlate(correlator, signal_mat, code, 1, 2500)
                early = get_early(correlator_result)
                prompt = get_prompt(correlator_result)
                late = get_late(correlator_result)
                @test early == late
                @test all(early .== 1476)
                @test all(prompt .== 2500)
            end
        end
    end
end
