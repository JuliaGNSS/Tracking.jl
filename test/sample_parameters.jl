module SampleParametersTest

using Test: @test, @testset, @inferred
using Unitful: MHz, Hz
using GNSSSignals: GPSL1CA
using Tracking: calc_num_chips_to_integrate, calc_num_samples_left_to_integrate

# `calc_num_code_blocks_to_integrate` is tested in TrackingLoops.
@testset "Calculate number of chips to integrate" begin
    gpsl1 = GPSL1CA()

    @test @inferred(calc_num_chips_to_integrate(gpsl1, 1, 10.5)) == 1023 - 10.5
    @test @inferred(calc_num_chips_to_integrate(gpsl1, 1, 1023 + 10.5)) == 1023 - 10.5
    @test @inferred(calc_num_chips_to_integrate(gpsl1, 2, 10.5)) == 2 * 1023 - 10.5
end

@testset "Calculate number of samples to integrate" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 4MHz
    code_doppler = -0.00022901219179036268Hz
    code_phase = 0.0
    preferred_num_code_blocks_to_integrate = 1

    num_samples_left_to_integrate = calc_num_samples_left_to_integrate(
        gpsl1,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        code_doppler,
        code_phase,
    )
    @test num_samples_left_to_integrate == 4001
end

end
