@testset "Sample parameters" begin

    gpsl1 = GPSL1()

    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 1, false)) == 1
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 2, false)) == 1
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 2, true)) == 2
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 20, true)) == 20
    @test @inferred(Tracking.calc_num_code_blocks_to_integrate(gpsl1, 21, true)) == 20

    @test @inferred(Tracking.calc_num_chips_to_integrate(gpsl1, 1, 10.5)) == 1023 - 10.5
    @test @inferred(Tracking.calc_num_chips_to_integrate(gpsl1, 1, 1023 + 10.5)) == 1023 - 10.5
    @test @inferred(Tracking.calc_num_chips_to_integrate(gpsl1, 2, 10.5)) == 2 * 1023 - 10.5

    code_freq = get_code_frequency(gpsl1)
    @test @inferred(Tracking.calc_num_samples_left_to_integrate(gpsl1, 1, 5e6Hz, 0.0Hz, 10.5)) ==
        ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)

    sample_params = @inferred Tracking.SampleParams(
        gpsl1,
        1,
        false
    )

    @test sample_params.signal_samples_to_integrate == 0
    @test sample_params.signal_start_sample == 1
    @test sample_params.samples_to_integrate_until_code_end == 0
    @test sample_params.num_code_blocks_to_integrate == 1

    num_samples_signal = 5000
    sampling_frequency = 5e6Hz
    code_doppler = 0.0Hz
    code_phase = 10.5
    preferred_num_code_blocks_to_integrate = 1
    secondary_code_or_bit_found = false

    next_sample_params = @inferred Tracking.SampleParams(
        gpsl1,
        sample_params,
        num_samples_signal,
        sampling_frequency,
        code_doppler,
        code_phase,
        preferred_num_code_blocks_to_integrate,
        secondary_code_or_bit_found
    )

    @test next_sample_params.signal_samples_to_integrate == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_sample_params.signal_start_sample == 1
    @test next_sample_params.samples_to_integrate_until_code_end == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_sample_params.num_code_blocks_to_integrate == 1

    # Bit not found -> It should default to one code block to integrate
    next_sample_params = @inferred Tracking.SampleParams(
        gpsl1,
        sample_params,
        num_samples_signal,
        sampling_frequency,
        code_doppler,
        code_phase,
        2,
        secondary_code_or_bit_found
    )

    @test next_sample_params.signal_samples_to_integrate == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_sample_params.signal_start_sample == 1
    @test next_sample_params.samples_to_integrate_until_code_end == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_sample_params.num_code_blocks_to_integrate == 1

    # Bit found
    next_sample_params = @inferred Tracking.SampleParams(
        gpsl1,
        sample_params,
        num_samples_signal,
        sampling_frequency,
        code_doppler,
        code_phase,
        2,
        true
    )

    @test next_sample_params.signal_samples_to_integrate == 5000
    @test next_sample_params.signal_start_sample == 1
    @test next_sample_params.samples_to_integrate_until_code_end == ceil(Int, (2 * 1023 - 10.5) * 5e6Hz / code_freq)
    @test next_sample_params.num_code_blocks_to_integrate == 2

    system_sats_state = SystemSatsState(
        gpsl1,
        [
            SatState(gpsl1, 1, sampling_frequency, code_phase, 0.0Hz),
            SatState(gpsl1, 2, sampling_frequency, 11.5, 0.0Hz)
        ]
    )

    system_sats_sample_params = @inferred Tracking.init_sample_params(
        (system_sats_state,),
        1,
    )

    @test system_sats_sample_params[1][1].signal_samples_to_integrate == 0
    @test system_sats_sample_params[1][1].signal_start_sample == 1
    @test system_sats_sample_params[1][1].samples_to_integrate_until_code_end == 0
    @test system_sats_sample_params[1][1].num_code_blocks_to_integrate == 1

    @test system_sats_sample_params[1][2].signal_samples_to_integrate == 0
    @test system_sats_sample_params[1][2].signal_start_sample == 1
    @test system_sats_sample_params[1][2].samples_to_integrate_until_code_end == 0
    @test system_sats_sample_params[1][2].num_code_blocks_to_integrate == 1


    next_system_sats_sample_params = @inferred Tracking.calc_sample_params(
        (system_sats_state,),
        system_sats_sample_params,
        num_samples_signal,
        sampling_frequency,
        preferred_num_code_blocks_to_integrate,
    )

    @test next_system_sats_sample_params[1][1].signal_samples_to_integrate == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_system_sats_sample_params[1][1].signal_start_sample == 1
    @test next_system_sats_sample_params[1][1].samples_to_integrate_until_code_end == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_system_sats_sample_params[1][1].num_code_blocks_to_integrate == 1

    @test next_system_sats_sample_params[1][2].signal_samples_to_integrate == ceil(Int, (1023 - 11.5) * 5e6Hz / code_freq)
    @test next_system_sats_sample_params[1][2].signal_start_sample == 1
    @test next_system_sats_sample_params[1][2].samples_to_integrate_until_code_end == ceil(Int, (1023 - 11.5) * 5e6Hz / code_freq)
    @test next_system_sats_sample_params[1][2].num_code_blocks_to_integrate == 1

    next_next_system_sats_sample_params = @inferred Tracking.calc_sample_params(
        (system_sats_state,),
        next_system_sats_sample_params,
        num_samples_signal,
        sampling_frequency,
        preferred_num_code_blocks_to_integrate,
    )

    @test next_next_system_sats_sample_params[1][1].signal_samples_to_integrate == 5000 - ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_next_system_sats_sample_params[1][1].signal_start_sample == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq) + 1
    @test next_next_system_sats_sample_params[1][1].samples_to_integrate_until_code_end == ceil(Int, (1023 - 10.5) * 5e6Hz / code_freq)
    @test next_next_system_sats_sample_params[1][1].num_code_blocks_to_integrate == 1

    @test next_next_system_sats_sample_params[1][2].signal_samples_to_integrate == 5000 - ceil(Int, (1023 - 11.5) * 5e6Hz / code_freq)
    @test next_next_system_sats_sample_params[1][2].signal_start_sample == ceil(Int, (1023 - 11.5) * 5e6Hz / code_freq) + 1
    @test next_next_system_sats_sample_params[1][2].samples_to_integrate_until_code_end == ceil(Int, (1023 - 11.5) * 5e6Hz / code_freq)
    @test next_next_system_sats_sample_params[1][2].num_code_blocks_to_integrate == 1
end