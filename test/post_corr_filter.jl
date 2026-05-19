module PostCorrFilterTest

using Test: @test, @testset
using StaticArrays: SVector
using Tracking:
    DefaultPostCorrFilter,
    EarlyPromptLateCorrelator,
    NumAnts,
    get_default_post_corr_filter,
    update

@testset "Post correlation filter" begin
    default_post_corr_filter = DefaultPostCorrFilter()

    @test update(default_post_corr_filter, randn(ComplexF64)) == default_post_corr_filter

    @test default_post_corr_filter(complex(1.2, 1.3)) == complex(1.2, 1.3)
    @test default_post_corr_filter(SVector(complex(1.2, 1.3), complex(2.2, 4.3))) ==
          complex(2.2, 4.3)
end

@testset "get_default_post_corr_filter fallback identity" begin
    # The generic fallback returns the identity closure, used for both
    # single-antenna (scalar accumulators) and multi-antenna (SVector
    # accumulators) correlators.
    single = EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    f_single = get_default_post_corr_filter(single)
    @test f_single(complex(1.0, 2.0)) == complex(1.0, 2.0)

    multi = EarlyPromptLateCorrelator(; num_ants = NumAnts(2))
    f_multi = get_default_post_corr_filter(multi)
    v = SVector(complex(7.0, 0.0), complex(9.0, 0.0))
    @test f_multi(v) == v
end

end
