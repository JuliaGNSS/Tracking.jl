module ConventionalPLLAndDLLTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, get_code_center_frequency_ratio
using TrackingLoopFilters: ThirdOrderBilinearLF, SecondOrderBilinearLF
using StaticArrays: SVector
using Tracking:
    aid_dopplers,
    SatConventionalPLLAndDLL,
    EarlyPromptLateCorrelator,
    SatState,
    SystemSatsState,
    ConventionalPLLAndDLL,
    TrackState,
    estimate_dopplers_and_filter_prompt,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_filtered_prompt,
    update_accumulator,
    get_default_correlator

@testset "Doppler aiding" begin
    gpsl1 = GPSL1()
    init_carrier_doppler = 10Hz
    init_code_doppler = 1Hz
    carrier_freq_update = 2Hz
    code_freq_update = -0.5Hz

    carrier_freq, code_freq = @inferred aid_dopplers(
        gpsl1,
        init_carrier_doppler,
        init_code_doppler,
        carrier_freq_update,
        code_freq_update,
    )

    @test carrier_freq == 10Hz + 2Hz
    @test code_freq == 1Hz + 2Hz / 1540 - 0.5Hz
end

@testset "Satellite Conventional PLL and DLL" begin
    pll_and_dll = @inferred SatConventionalPLLAndDLL(
        init_carrier_doppler = 500.0Hz,
        init_code_doppler = 100.0Hz,
    )

    @test pll_and_dll.init_carrier_doppler == 500.0Hz
    @test pll_and_dll.init_code_doppler == 100.0Hz
    @test pll_and_dll.carrier_loop_filter == ThirdOrderBilinearLF()
    @test pll_and_dll.code_loop_filter == SecondOrderBilinearLF()
end

@testset "Conventional PLL and DLL" begin
    sampling_frequency = 5e6Hz

    gpsl1 = GPSL1()

    carrier_doppler = 100.0Hz
    prn = 1
    code_phase = 0.5
    preferred_num_code_blocks_to_integrate = 1

    sat_state = SatState(gpsl1, prn, code_phase, carrier_doppler)

    system_sats_state = (SystemSatsState(gpsl1, sat_state),)

    doppler_estimator = ConventionalPLLAndDLL(system_sats_state)

    # Number of samples too small to generate a new estimate for phases and dopplers
    num_samples = 2000
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    sat_state_after_small_integration =
        SatState(sat_state; integrated_samples = num_samples, correlator)

    track_state = TrackState(gpsl1, sat_state_after_small_integration; doppler_estimator)

    new_track_state = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    # Since number of samples is too small the state doesn't change
    @test get_carrier_doppler(new_track_state) == carrier_doppler
    @test get_code_doppler(new_track_state) ==
          get_code_center_frequency_ratio(gpsl1) * carrier_doppler
    @test get_last_fully_integrated_filtered_prompt(new_track_state) == 0.0

    # This time it is large enough to produce new dopplers and phases
    num_samples = 5000
    sat_state_after_full_integration = SatState(
        sat_state;
        is_integration_completed = true,
        integrated_samples = num_samples,
        correlator,
    )
    track_state = TrackState(gpsl1, sat_state_after_full_integration)

    new_track_state_after_full_integration = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    @test get_carrier_doppler(new_track_state_after_full_integration) ==
          100.52615628464486Hz
    @test get_code_doppler(new_track_state_after_full_integration) == -0.12599932765737312Hz
    @test get_last_fully_integrated_filtered_prompt(
        new_track_state_after_full_integration,
    ) == 0.4 + 0.004im
end

end
