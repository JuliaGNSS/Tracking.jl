module ConventionalPLLAndDLLTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, get_code_center_frequency_ratio
using TrackingLoopFilters: ThirdOrderBilinearLF, SecondOrderBilinearLF
using StaticArrays: SVector
using Dictionaries: Dictionary
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
    get_default_correlator,
    get_sat_state,
    merge_sats

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
    @test pll_and_dll.carrier_loop_filter_bandwidth == 18.0Hz
    @test pll_and_dll.code_loop_filter_bandwidth == 1.0Hz

    gpsl1 = GPSL1()
    sat_state = SatState(gpsl1, 1, 0.5, 100.0Hz)
    from_sat_state = @inferred SatConventionalPLLAndDLL(
        sat_state,
        ThirdOrderBilinearLF(),
        SecondOrderBilinearLF();
        carrier_loop_filter_bandwidth = 25.0Hz,
        code_loop_filter_bandwidth = 2.0Hz,
    )
    @test from_sat_state.carrier_loop_filter_bandwidth == 25.0Hz
    @test from_sat_state.code_loop_filter_bandwidth == 2.0Hz

    # Update-from-existing constructor preserves bandwidth when not overridden,
    # and overrides when provided.
    preserved = @inferred SatConventionalPLLAndDLL(from_sat_state)
    @test preserved.carrier_loop_filter_bandwidth == 25.0Hz
    @test preserved.code_loop_filter_bandwidth == 2.0Hz

    overridden = @inferred SatConventionalPLLAndDLL(
        from_sat_state;
        carrier_loop_filter_bandwidth = 30.0Hz,
    )
    @test overridden.carrier_loop_filter_bandwidth == 30.0Hz
    @test overridden.code_loop_filter_bandwidth == 2.0Hz
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
    @test get_code_doppler(new_track_state_after_full_integration) == -0.16073504885813858Hz
    @test get_last_fully_integrated_filtered_prompt(
        new_track_state_after_full_integration,
    ) == 0.4 + 0.004im
end

@testset "Per-satellite bandwidths drive the loop filters" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1()
    carrier_doppler = 100.0Hz
    code_phase = 0.5
    preferred_num_code_blocks_to_integrate = 1
    num_samples = 5000
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )

    # Two sats, both fully integrated with the same correlator: any difference
    # in the Doppler update must come from per-sat bandwidths.
    sat1_initial = SatState(gpsl1, 1, code_phase, carrier_doppler)
    sat2_initial = SatState(gpsl1, 2, code_phase, carrier_doppler)
    sat1 = SatState(
        sat1_initial;
        is_integration_completed = true,
        integrated_samples = num_samples,
        correlator,
    )
    sat2 = SatState(
        sat2_initial;
        is_integration_completed = true,
        integrated_samples = num_samples,
        correlator,
    )
    system_sats_state = (SystemSatsState(gpsl1, [sat1, sat2]),)

    # Build a doppler estimator, then bump sat2's bandwidths via the
    # update-from-existing constructor.
    base_estimator = ConventionalPLLAndDLL(system_sats_state)
    per_sat_states = only(base_estimator.states)
    sat2_state_wide = SatConventionalPLLAndDLL(
        per_sat_states[2];
        carrier_loop_filter_bandwidth = 36.0Hz,
        code_loop_filter_bandwidth = 2.0Hz,
    )
    new_sat_states = (Dictionary(
        [1, 2], [per_sat_states[1], sat2_state_wide],
    ),)
    doppler_estimator = ConventionalPLLAndDLL(base_estimator, new_sat_states)

    @test doppler_estimator.states[1][1].carrier_loop_filter_bandwidth == 18.0Hz
    @test doppler_estimator.states[1][2].carrier_loop_filter_bandwidth == 36.0Hz
    @test doppler_estimator.states[1][2].code_loop_filter_bandwidth == 2.0Hz

    track_state = TrackState(
        gpsl1, [sat1, sat2]; doppler_estimator,
    )
    new_track_state = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    # Sat 1 matches the baseline result from the previous testset.
    @test get_carrier_doppler(new_track_state, 1) == 100.52615628464486Hz
    @test get_code_doppler(new_track_state, 1) == -0.16073504885813858Hz

    # Sat 2 has different bandwidths so must produce a different update.
    @test get_carrier_doppler(new_track_state, 2) != get_carrier_doppler(new_track_state, 1)
    @test get_code_doppler(new_track_state, 2) != get_code_doppler(new_track_state, 1)
end

@testset "Bandwidths propagate through ConventionalPLLAndDLL constructor" begin
    gpsl1 = GPSL1()
    sat_state = SatState(gpsl1, 1, 0.5, 100.0Hz)
    system_sats_state = (SystemSatsState(gpsl1, sat_state),)

    estimator = ConventionalPLLAndDLL(
        system_sats_state;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test estimator.carrier_loop_filter_bandwidth == 22.0Hz
    @test estimator.code_loop_filter_bandwidth == 1.5Hz
    # Per-sat entries inherit the estimator-level bandwidths.
    @test only(estimator.states[1]).carrier_loop_filter_bandwidth == 22.0Hz
    @test only(estimator.states[1]).code_loop_filter_bandwidth == 1.5Hz
end

@testset "merge_sats propagates bandwidths to new sats" begin
    gpsl1 = GPSL1()
    sat1 = SatState(gpsl1, 1, 0.5, 100.0Hz)
    system_sats_state = (SystemSatsState(gpsl1, sat1),)
    estimator = ConventionalPLLAndDLL(
        system_sats_state;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )

    sat2 = SatState(gpsl1, 2, 0.25, 200.0Hz)
    merged = merge_sats(estimator, 1, Dictionary([2], [sat2]))

    @test merged.states[1][2].carrier_loop_filter_bandwidth == 22.0Hz
    @test merged.states[1][2].code_loop_filter_bandwidth == 1.5Hz
end

end
