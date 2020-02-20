@testset "Carrier replica" begin

    carrier = StructArray(zeros(Complex{Int16}, 2500))

    Tracking.gen_carrier_replica!(
        carrier,
        1500Hz,
        2.5e6Hz,
        0.25, # π / 2
        111,
        2390,
    )

    @test sqrt(mean(abs2.(carrier.re[111:2500] ./ 1 << 7 .-
                          cos.(2π * (0:2389) * 1500Hz / 2.5e6Hz .+ π / 2)))) < 8e-3
end

@testset "Update carrier phase" begin
    carrier_phase = 0.25
    carrier_frequency = 10Hz
    sample_frequency = 100Hz
    num_samples = 2000
    phase = @inferred Tracking.update_carrier_phase(
        num_samples,
        carrier_frequency,
        sample_frequency,
        carrier_phase
    )
    @test phase ≈ mod(0.25 + 0.1 * 2000 + 0.5, 1) - 0.5 atol = 1 / 1 << 30 * 2500

    carrier_phase = 0.25
    carrier_frequency = 10Hz
    sample_frequency = 2.5e6Hz
    num_samples = 2500
    phase = @inferred Tracking.update_carrier_phase(
        num_samples,
        carrier_frequency,
        sample_frequency,
        carrier_phase
    )
    @test phase ≈ mod(0.25 + 10 * 2500 / 2.5e6 + 0.5, 1) - 0.5 atol = 1 / 1 << 30 * 2500
end
