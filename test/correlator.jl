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

        signal = 1.0 + 0.0im
        carrier = cis(0.0)
        early_late_sample_shift = 2
        # Late Prompt Early
        early_code = 1
        code_register = 1 << 4 + 0 << 2 + 1
        correlator = EarlyPromptLateCorrelator(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im)

        next_correlator = @inferred Tracking.correlate_iteration(
            GPSL1,
            correlator,
            signal,
            carrier,
            code_register,
            early_late_sample_shift,
            early_code
        )
        @test get_early(next_correlator) ≈ 1.0 + signal * conj(carrier) * 1
        @test get_prompt(next_correlator) ≈ 1.0 + signal * conj(carrier) * -1
        @test get_late(next_correlator) ≈ 1.0 + signal * conj(carrier) * 1

        early_code = -1
        code_register = 0 << 4 + 1 << 2 + 0

        next_correlator = @inferred Tracking.correlate_iteration(
            GPSL1,
            correlator,
            signal,
            carrier,
            code_register,
            early_late_sample_shift,
            early_code
        )
        @test get_early(next_correlator) ≈ 1.0 + signal * conj(carrier) * -1
        @test get_prompt(next_correlator) ≈ 1.0 + signal * conj(carrier) * 1
        @test get_late(next_correlator) ≈ 1.0 + signal * conj(carrier) * -1

    end

end
