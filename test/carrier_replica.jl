module CarrierReplicaTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using Tracking: update_carrier_phase, get_current_carrier_frequency

@testset "Update carrier phase" begin
    carrier_phase = 0.25
    carrier_frequency = 10Hz
    sampling_frequency = 100Hz
    num_samples = 2000
    phase = @inferred update_carrier_phase(
        num_samples,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
    )
    @test phase ≈ mod(0.25 + 0.1 * 2000 + 0.5, 1) - 0.5

    carrier_phase = 0.25
    carrier_frequency = 10Hz
    sampling_frequency = 2.5e6Hz
    num_samples = 2500
    phase = @inferred update_carrier_phase(
        num_samples,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
    )
    @test phase ≈ mod(0.25 + 10 * 2500 / 2.5e6 + 0.5, 1) - 0.5
end

@testset "Current carrier frequency" begin
    @test @inferred(get_current_carrier_frequency(1.0e6Hz, 200Hz)) == 1.0002e6Hz
    @test @inferred(get_current_carrier_frequency(0Hz, -123.5Hz)) == -123.5Hz
end

end
