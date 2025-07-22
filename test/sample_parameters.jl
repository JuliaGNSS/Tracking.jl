@testset "Calculate number of code blocks / chips to integrate" begin
    gpsl1 = GPSL1()

    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 1, false)) == 1
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 2, false)) == 1
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 2, true)) == 2
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 20, true)) == 20
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 21, true)) == 20

    @test @inferred(Tracking.calc_num_chips_to_integrate(gpsl1, 1, 10.5)) == 1023 - 10.5
    @test @inferred(Tracking.calc_num_chips_to_integrate(gpsl1, 1, 1023 + 10.5)) ==
          1023 - 10.5
    @test @inferred(Tracking.calc_num_chips_to_integrate(gpsl1, 2, 10.5)) == 2 * 1023 - 10.5
end

@testset "Calculate number of samples to integrate" begin
    gpsl1 = GPSL1()
    signal_start_sample = 1
    sampling_frequency = 4MHz
    code_doppler = -0.00022901219179036268Hz
    code_phase = 0.0
    preferred_num_code_blocks_to_integrate = 1
    secondary_code_or_bit_found = true
    num_samples_signal = 4000

    num_samples_left_to_integrate = Tracking.calc_num_samples_left_to_integrate(
        gpsl1,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        code_doppler,
        code_phase,
    )
    @test num_samples_left_to_integrate == 4001

    samples_to_integrate, is_integration_completed =
        Tracking.calc_signal_samples_to_integrate(
            gpsl1,
            signal_start_sample,
            sampling_frequency,
            code_doppler,
            code_phase,
            preferred_num_code_blocks_to_integrate,
            secondary_code_or_bit_found,
            num_samples_signal,
        )

    @test samples_to_integrate == num_samples_signal
    @test is_integration_completed == false
end
