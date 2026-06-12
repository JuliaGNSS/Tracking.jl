module SampleParametersTest

using Test: @test, @testset, @inferred
using Unitful: MHz, Hz
using GNSSSignals: GPSL1CA, GPSL1C_P, GPSL5I
using Tracking:
    calc_num_code_blocks_to_integrate,
    calc_num_chips_to_integrate,
    calc_num_samples_left_to_integrate,
    calc_signal_samples_to_integrate

@testset "Calculate number of code blocks / chips to integrate" begin
    gpsl1 = GPSL1CA()

    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 1, false)) == 1
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 2, false)) == 1
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 2, true)) == 2
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 20, true)) == 20
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 21, true)) == 20

    # A preferred value that doesn't divide the 20 blocks per bit is clamped
    # to the largest divisor below it — an integration that straddled a bit
    # boundary would never emit a bit again (issue #128).
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 3, true)) == 2
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 7, true)) == 5
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1, 19, true)) == 10

    # Pilot signals (data_frequency == 0) are capped by their secondary-code
    # period once the secondary code has been found — not clamped to one
    # block (issue #134; `set_preferred_num_code_blocks_to_integrate!` used
    # to be a silent no-op for GPS L1C-P). The same divisor clamp as the
    # data-bit path applies (1800 = 2³·3²·5², so 10 and 1800 are exact).
    gpsl1c_p = GPSL1C_P()
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1c_p, 10, false)) == 1
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1c_p, 10, true)) == 10
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1c_p, 1800, true)) == 1800
    @test @inferred(calc_num_code_blocks_to_integrate(gpsl1c_p, 2000, true)) == 1800
    # Data-bearing signals stay capped by the bit period.
    @test @inferred(calc_num_code_blocks_to_integrate(GPSL5I(), 11, true)) == 10

    @test @inferred(calc_num_chips_to_integrate(gpsl1, 1, 10.5)) == 1023 - 10.5
    @test @inferred(calc_num_chips_to_integrate(gpsl1, 1, 1023 + 10.5)) == 1023 - 10.5
    @test @inferred(calc_num_chips_to_integrate(gpsl1, 2, 10.5)) == 2 * 1023 - 10.5
end

@testset "Calculate number of samples to integrate" begin
    gpsl1 = GPSL1CA()
    signal_start_sample = 1
    sampling_frequency = 4MHz
    code_doppler = -0.00022901219179036268Hz
    code_phase = 0.0
    preferred_num_code_blocks_to_integrate = 1
    secondary_code_or_bit_found = true
    num_samples_signal = 4000

    num_samples_left_to_integrate = calc_num_samples_left_to_integrate(
        gpsl1,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        code_doppler,
        code_phase,
    )
    @test num_samples_left_to_integrate == 4001

    samples_to_integrate, is_integration_completed = calc_signal_samples_to_integrate(
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

end
