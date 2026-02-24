module CorrelatorTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, GalileoE1B, get_code, get_code_frequency
using StaticArrays: SVector, MVector, MMatrix
using StructArrays: StructArray
using Tracking:
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    NumAnts,
    NumAccumulators,
    add_to_previous,
    get_very_early,
    get_early,
    get_prompt,
    get_late,
    get_num_ants,
    get_very_late,
    get_accumulators,
    get_early_late_sample_spacing,
    get_early_late_code_spacing,
    get_correlator_code_shifts,
    use_fast_code_replica,
    DefaultPostCorrFilter,
    apply,
    normalize,
    correlate,
    get_correlator_sample_shifts,
    get_initial_accumulator,
    zero_accumulators

@testset "Correlator" begin
    @testset "Get initial accumulators" begin
        @test @inferred(get_initial_accumulator(NumAnts(1), NumAccumulators(3))) isa
              SVector{3,ComplexF64}
        @test @inferred(get_initial_accumulator(NumAnts(4), NumAccumulators(3))) isa
              SVector{3,SVector{4,ComplexF64}}
        @test @inferred(get_initial_accumulator(NumAnts(1), 3)) isa Vector{ComplexF64}
        @test @inferred(get_initial_accumulator(NumAnts(4), 3)) isa
              Vector{SVector{4,ComplexF64}}
    end

    @testset "Zero accumulators" begin
        accumulator =
            @inferred(zero_accumulators(ones(SVector{3,ComplexF64}), zeros(ComplexF32, 1)))
        @test accumulator isa MVector{3,Float32}
        @test accumulator == [0, 0, 0]

        accumulator = @inferred(
            zero_accumulators(
                SVector(ones(SVector{3,ComplexF64}), ones(SVector{3,ComplexF64})),
                zeros(ComplexF32, 1),
            )
        )
        @test accumulator isa MMatrix{3,2,Float32}
        @test accumulator == [0.0 0.0; 0.0 0.0; 0.0 0.0]

        accumulator =
            @inferred(zero_accumulators(ones(ComplexF64, 3), zeros(ComplexF32, 1)))
        @test accumulator isa Vector{Float32}
        @test accumulator == [0, 0, 0]

        accumulator = @inferred(
            zero_accumulators(
                [ones(SVector{3,ComplexF64}), ones(SVector{3,ComplexF64})],
                zeros(ComplexF32, 1),
            )
        )
        @test accumulator isa Matrix{Float32}
        @test accumulator == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    end

    @testset "Add to previous accumulators" begin
        accumulator = @inferred(
            add_to_previous(
                get_initial_accumulator(NumAnts(2), NumAccumulators(3)),
                ones(MMatrix{2,3,Float32}),
                ones(MMatrix{2,3,Float32}),
            )
        )
        @test accumulator isa SVector{3,SVector{2,ComplexF64}}
        @test accumulator == [
            [1.0 + 1.0im, 1.0 + 1.0im],
            [1.0 + 1.0im, 1.0 + 1.0im],
            [1.0 + 1.0im, 1.0 + 1.0im],
        ]

        accumulator = @inferred(
            add_to_previous(
                get_initial_accumulator(NumAnts(2), 3),
                ones(Float32, 2, 3),
                ones(Float32, 2, 3),
            )
        )
        @test accumulator isa Vector{SVector{2,ComplexF64}}
        @test accumulator == [
            [1.0 + 1.0im, 1.0 + 1.0im],
            [1.0 + 1.0im, 1.0 + 1.0im],
            [1.0 + 1.0im, 1.0 + 1.0im],
        ]
    end

    @testset "Correlator constructor" begin
        gpsl1 = GPSL1()
        sampling_frequency = 5e6Hz
        correlator = @inferred EarlyPromptLateCorrelator()

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0
        @test correlator.preferred_early_late_to_prompt_code_shift == 0.5
        @test get_num_ants(correlator) == 1

        correlator = @inferred EarlyPromptLateCorrelator(num_ants = NumAnts(1))

        @test @inferred(get_early(correlator)) == 0.0
        @test @inferred(get_prompt(correlator)) == 0.0
        @test @inferred(get_late(correlator)) == 0.0
        @test correlator.preferred_early_late_to_prompt_code_shift == 0.5
        @test get_num_ants(correlator) == 1

        correlator = @inferred EarlyPromptLateCorrelator(num_ants = NumAnts(2))

        @test @inferred(get_early(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_prompt(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test @inferred(get_late(correlator)) == SVector(0.0 + 0.0im, 0.0 + 0.0im)
        @test correlator.preferred_early_late_to_prompt_code_shift == 0.5
        @test get_num_ants(correlator) == 2

        correlator =
            EarlyPromptLateCorrelator(SVector(1.0 + 0.0im, 2.0 + 0.0im, 3.0 + 0.0im), 0.5)
        @test @inferred(get_early(correlator)) == 3.0
        @test @inferred(get_prompt(correlator)) == 2.0
        @test @inferred(get_late(correlator)) == 1.0
        @test get_num_ants(correlator) == 1

        correlator = VeryEarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 2.0 + 0.0im, 3.0 + 0.0im, 4.0 + 0.0im, 5.0 + 0.0im),
            0.15,
            0.6,
        )
        @test @inferred(get_very_early(correlator)) == 5.0
        @test @inferred(get_early(correlator)) == 4.0
        @test @inferred(get_prompt(correlator)) == 3.0
        @test @inferred(get_late(correlator)) == 2.0
        @test @inferred(get_very_late(correlator)) == 1.0
        @test get_num_ants(correlator) == 1
    end

    @testset "Calculate sample shift" begin
        gpsl1 = GPSL1()
        code_frequency = get_code_frequency(gpsl1)
        sampling_frequency = code_frequency * 4
        correlator = EarlyPromptLateCorrelator()
        @test @inferred(
            get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
        ) == -2:2:2

        galileoE1B = GalileoE1B()
        sampling_frequency = get_code_frequency(galileoE1B) * 4
        correlator = VeryEarlyPromptLateCorrelator()
        @test @inferred(
            get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
        ) == -2:1:2

        sampling_frequency = get_code_frequency(galileoE1B) * 8
        @test @inferred(
            get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
        ) == [-5, -1, 0, 1, 5]
    end

    @testset "Zeroing correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            0.5,
        )

        @test @inferred(zero(correlator)) == EarlyPromptLateCorrelator(
            SVector(
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
                SVector(0.0 + 0.0im, 0.0 + 0.0im),
            ),
            0.5,
        )

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im),
            0.5,
        )

        @test @inferred(zero(correlator)) ==
              EarlyPromptLateCorrelator(SVector(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im), 0.5)
    end
    @testset "Filter correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 1.0im),
                SVector(1.0 + 0.0im, 2.0 + 0.0im),
                SVector(1.0 + 1.0im, 1.0 + 3.0im),
            ),
            0.5,
        )
        default_post_corr_filter = DefaultPostCorrFilter()
        filtered_correlator = @inferred apply(default_post_corr_filter, correlator)
        @test filtered_correlator ==
              EarlyPromptLateCorrelator(SVector(1.0 + 1.0im, 2.0 + 0.0im, 1.0 + 3.0im), 0.5)

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im),
            0.5,
        )
        filtered_correlator = @inferred apply(x -> 2 * x, correlator)
        @test filtered_correlator ==
              EarlyPromptLateCorrelator(SVector(2.0 + 0.0im, 2.0 + 0.0im, 2.0 + 0.0im), 0.5)
    end
    @testset "Early late sample spacing" begin
        gpsl1 = GPSL1()
        code_frequency = get_code_frequency(gpsl1)
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            0.5,
        )
        @test get_early_late_sample_spacing(correlator, 2e6Hz, code_frequency) == 2

        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            0.5,
        )
        @test get_early_late_sample_spacing(correlator, 4e6Hz, code_frequency) == 4

        correlator = @inferred VeryEarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            0.15,
            0.6,
        )
        @test get_early_late_sample_spacing(correlator, 4e6Hz, code_frequency) == 2
    end
    @testset "Normalize correlator" begin
        correlator = @inferred EarlyPromptLateCorrelator(
            SVector(1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im),
            0.5,
        )
        normalized_correlator = @inferred normalize(correlator, 10)
        @test normalized_correlator ==
              EarlyPromptLateCorrelator(SVector(0.1 + 0.0im, 0.1 + 0.0im, 0.1 + 0.0im), 0.5)

        correlator = EarlyPromptLateCorrelator(
            SVector(
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
                SVector(1.0 + 0.0im, 1.0 + 0.0im),
            ),
            0.5,
        )
        normalized_correlator = @inferred normalize(correlator, 10)
        @test normalized_correlator == EarlyPromptLateCorrelator(
            SVector(
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
                SVector(0.1 + 0.0im, 0.1 + 0.0im),
            ),
            0.5,
        )
    end

    @testset "Correlate" begin
        gpsl1 = GPSL1()
        sampling_frequency = 2.5e6Hz
        signal = StructArray{Complex{Float32}}((
            Float32.(get_code.(gpsl1, (1:2500) * 1023e3 / (sampling_frequency / Hz), 1)),
            zeros(Float32, 2500),
        ))
        code = get_code.(gpsl1, ((1+-1):(2500+1)) * 1023e3 / 2.5e6, 1)

        @testset "↳ $na antennas" for na ∈ [1, 4]
            signal_mat = na == 1 ? signal : repeat(signal; outer = (1, na))

            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(na))

            sample_shifts = get_correlator_sample_shifts(
                correlator,
                sampling_frequency,
                get_code_frequency(gpsl1),
            )

            correlator_result =
                @inferred correlate(correlator, signal_mat, sample_shifts, code, 1, 2500)
            early = get_early(correlator_result)
            prompt = get_prompt(correlator_result)
            late = get_late(correlator_result)
            @test early == late
            @test all(early .== 1476)
            @test all(prompt .== 2500)
        end
    end

    @testset "Code shifts" begin
        correlator = EarlyPromptLateCorrelator()
        @test get_correlator_code_shifts(correlator) == SVector(-0.5, 0.0, 0.5)

        correlator = EarlyPromptLateCorrelator(;
            preferred_early_late_to_prompt_code_shift = 0.25,
        )
        @test get_correlator_code_shifts(correlator) == SVector(-0.25, 0.0, 0.25)

        correlator = VeryEarlyPromptLateCorrelator()
        @test get_correlator_code_shifts(correlator) ==
              SVector(-0.6, -0.15, 0.0, 0.15, 0.6)
    end

    @testset "Early late code spacing" begin
        correlator = EarlyPromptLateCorrelator()
        @test get_early_late_code_spacing(correlator) == 1.0

        correlator = EarlyPromptLateCorrelator(;
            preferred_early_late_to_prompt_code_shift = 0.25,
        )
        @test get_early_late_code_spacing(correlator) == 0.5

        correlator = VeryEarlyPromptLateCorrelator()
        @test get_early_late_code_spacing(correlator) == 0.3
    end

    @testset "Use fast code replica decision" begin
        gpsl1 = GPSL1()
        code_frequency = get_code_frequency(gpsl1)

        correlator = EarlyPromptLateCorrelator()

        # High oversampling (4x): sample shift = round(0.5 * 4) = 2,
        # actual = 2 / 4 = 0.5 chips, error = 0% → fast path
        @test use_fast_code_replica(correlator, code_frequency * 4, code_frequency) == true

        # Low oversampling (~1.111x): sample shift = max(1, round(0.5 * 1.111)) = 1,
        # actual = 1 / 1.111 = 0.9 chips, error = 80% → slow path
        @test use_fast_code_replica(correlator, code_frequency * 10 / 9, code_frequency) ==
              false

        # 2x oversampling: sample shift = round(0.5 * 2) = 1,
        # actual = 1 / 2 = 0.5 chips, error = 0% → fast path
        @test use_fast_code_replica(correlator, code_frequency * 2, code_frequency) == true

        # VeryEarlyPromptLateCorrelator at 4x oversampling: 0.15 shift rounds to 1 sample,
        # actual = 0.25 chips, 66% error → slow path
        vepl_correlator = VeryEarlyPromptLateCorrelator()
        @test use_fast_code_replica(vepl_correlator, code_frequency * 4, code_frequency) ==
              false

        # VeryEarlyPromptLateCorrelator at high oversampling (20x): all shifts round well
        @test use_fast_code_replica(vepl_correlator, code_frequency * 20, code_frequency) ==
              true
    end
end

end
