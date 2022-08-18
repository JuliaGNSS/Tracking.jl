@testset "Post correlation filter" begin
    default_post_corr_filter = Tracking.DefaultPostCorrFilter()

    @test Tracking.update(default_post_corr_filter, randn(ComplexF64)) == default_post_corr_filter

    @test default_post_corr_filter(complex(1.2, 1.3)) == complex(1.2, 1.3)
    @test default_post_corr_filter(SVector(complex(1.2, 1.3), complex(2.2, 4.3))) == complex(2.2, 4.3)
end