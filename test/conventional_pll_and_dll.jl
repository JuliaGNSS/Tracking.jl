module ConventionalPLLAndDLLTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz
using GNSSSignals: GPSL1CA, get_code_center_frequency_ratio
using TrackingLoopFilters: ThirdOrderBilinearLF, SecondOrderBilinearLF
using StaticArrays: SVector
using Dictionaries: Dictionary
using Tracking:
    aid_dopplers,
    SatConventionalPLLAndDLL,
    EarlyPromptLateCorrelator,
    TrackedSignal,
    TrackedSat,
    init_estimator_state,
    ConventionalPLLAndDLL,
    TrackState,
    Measurement,
    estimate_dopplers_and_filter_prompt,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_sat_state,
    get_sat_states,
    update_accumulator,
    get_default_correlator,
    merge_sats

# Build a stub `(l1 = Measurement(...),)` NamedTuple to pass to the
# estimator. Samples are unused by `estimate_dopplers_and_filter_prompt`
# (it only reads `sampling_frequency` per group), so an empty buffer
# suffices.
_meas_l1(fs) = (l1 = Measurement(ComplexF64[], fs),)

@testset "Doppler aiding" begin
    gpsl1 = GPSL1CA()
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

    gpsl1 = GPSL1CA()
    sat_state = TrackedSat(gpsl1, 1, 0.5, 100.0Hz)
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

    gpsl1 = GPSL1CA()

    carrier_doppler = 100.0Hz
    prn = 1
    code_phase = 0.5
    preferred_num_code_blocks_to_integrate = 1

    doppler_estimator = ConventionalPLLAndDLL()

    sat_state = TrackedSat(gpsl1, prn, code_phase, carrier_doppler; doppler_estimator)

    # Number of samples too small to generate a new estimate for phases and dopplers
    num_samples = 2000
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    sat_state_after_small_integration = TrackedSat(
        sat_state;
        signals = (
            TrackedSignal(
                only(sat_state.signals);
                integrated_samples = num_samples,
                correlator,
            ),
        ),
    )

    track_state = TrackState(gpsl1, sat_state_after_small_integration; doppler_estimator)

    new_track_state = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        _meas_l1(sampling_frequency),
        preferred_num_code_blocks_to_integrate,
    )

    # Since number of samples is too small the state doesn't change
    @test get_carrier_doppler(new_track_state) == carrier_doppler
    @test get_code_doppler(new_track_state) ==
          get_code_center_frequency_ratio(gpsl1) * carrier_doppler
    @test get_last_fully_integrated_filtered_prompt(new_track_state) == 0.0
    # No integration completed -> no filtered prompt recorded
    @test isempty(get_filtered_prompts(get_sat_state(new_track_state, prn)))

    # This time it is large enough to produce new dopplers and phases
    num_samples = 5000
    sat_state_after_full_integration = TrackedSat(
        sat_state;
        signals = (
            TrackedSignal(
                only(sat_state.signals);
                is_integration_completed = true,
                integrated_samples = num_samples,
                correlator,
            ),
        ),
    )
    track_state = TrackState(gpsl1, sat_state_after_full_integration; doppler_estimator)

    new_track_state_after_full_integration = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        _meas_l1(sampling_frequency),
        preferred_num_code_blocks_to_integrate,
    )

    @test get_carrier_doppler(new_track_state_after_full_integration) ==
          100.52615628464486Hz
    @test get_code_doppler(new_track_state_after_full_integration) == -0.16073504885813858Hz
    @test get_last_fully_integrated_filtered_prompt(
        new_track_state_after_full_integration,
    ) == 0.4 + 0.004im
    # One integration completed -> one filtered prompt recorded, equal to the
    # last_fully_integrated_filtered_prompt value.
    let prompts = get_filtered_prompts(
            get_sat_state(new_track_state_after_full_integration, prn),
        )
        @test length(prompts) == 1
        @test prompts[1] == 0.4 + 0.004im
    end
end

@testset "Per-satellite bandwidths drive the loop filters" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
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
    doppler_estimator = ConventionalPLLAndDLL()
    sat1_initial = TrackedSat(gpsl1, 1, code_phase, carrier_doppler; doppler_estimator)
    sat2_initial = TrackedSat(gpsl1, 2, code_phase, carrier_doppler; doppler_estimator)
    sat1 = TrackedSat(
        sat1_initial;
        signals = (
            TrackedSignal(
                only(sat1_initial.signals);
                is_integration_completed = true,
                integrated_samples = num_samples,
                correlator,
            ),
        ),
    )
    sat2_pre = TrackedSat(
        sat2_initial;
        signals = (
            TrackedSignal(
                only(sat2_initial.signals);
                is_integration_completed = true,
                integrated_samples = num_samples,
                correlator,
            ),
        ),
    )
    # Bump sat2's bandwidths by replacing its doppler estimator state with a
    # custom-configured SatConventionalPLLAndDLL.
    sat1_de = sat1.doppler_estimator_state
    sat2_de = SatConventionalPLLAndDLL(
        sat2_pre.doppler_estimator_state;
        carrier_loop_filter_bandwidth = 36.0Hz,
        code_loop_filter_bandwidth = 2.0Hz,
    )
    sat2 = TrackedSat(
        sat2_pre.prn,
        sat2_pre.code_phase,
        sat2_pre.code_doppler,
        sat2_pre.carrier_phase,
        sat2_pre.carrier_doppler,
        sat2_pre.signal_start_sample,
        sat2_pre.signals,
        sat2_de,
    )
    tracked = [sat1, sat2]

    @test sat1_de.carrier_loop_filter_bandwidth == 18.0Hz
    @test sat2_de.carrier_loop_filter_bandwidth == 36.0Hz
    @test sat2_de.code_loop_filter_bandwidth == 2.0Hz

    track_state = TrackState(gpsl1, tracked; doppler_estimator)
    new_track_state = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        _meas_l1(sampling_frequency),
        preferred_num_code_blocks_to_integrate,
    )

    # Sat 1 matches the baseline result from the previous testset.
    @test get_carrier_doppler(new_track_state, 1) == 100.52615628464486Hz
    @test get_code_doppler(new_track_state, 1) == -0.16073504885813858Hz

    # Sat 2 has different bandwidths so must produce a different update.
    @test get_carrier_doppler(new_track_state, 2) != get_carrier_doppler(new_track_state, 1)
    @test get_code_doppler(new_track_state, 2) != get_code_doppler(new_track_state, 1)
end

@testset "Bandwidths propagate through ConventionalPLLAndDLL constructor" begin
    gpsl1 = GPSL1CA()
    sat_state = TrackedSat(gpsl1, 1, 0.5, 100.0Hz)
    estimator = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test estimator.carrier_loop_filter_bandwidth == 22.0Hz
    @test estimator.code_loop_filter_bandwidth == 1.5Hz
    # init_estimator_state seeds each sat with the configured bandwidths.
    de_state = init_estimator_state(estimator, sat_state)
    @test de_state.carrier_loop_filter_bandwidth == 22.0Hz
    @test de_state.code_loop_filter_bandwidth == 1.5Hz

    # Kwarg-update constructor returns a new estimator with overridden
    # bandwidths and preserves the carrier/code filter type parameters.
    bumped = ConventionalPLLAndDLL(estimator;
        carrier_loop_filter_bandwidth = 30.0Hz,
    )
    @test bumped.carrier_loop_filter_bandwidth == 30.0Hz
    @test bumped.code_loop_filter_bandwidth == 1.5Hz  # unchanged
    preserved = ConventionalPLLAndDLL(estimator)
    @test preserved.carrier_loop_filter_bandwidth == 22.0Hz
    @test preserved.code_loop_filter_bandwidth == 1.5Hz
end

@testset "merge_sats with matching estimator carries through bandwidths" begin
    gpsl1 = GPSL1CA()
    estimator = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sat1 = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator = estimator)
    track_state = TrackState(gpsl1, sat1; doppler_estimator = estimator)

    # The incoming sat must be constructed with the same estimator
    # (TrackState's slot type pins the doppler_estimator_state type).
    sat2 = TrackedSat(gpsl1, 2, 0.25, 200.0Hz; doppler_estimator = estimator)
    merged = merge_sats(track_state, sat2)

    de_state = get_sat_states(merged)[2].doppler_estimator_state
    @test de_state.carrier_loop_filter_bandwidth == 22.0Hz
    @test de_state.code_loop_filter_bandwidth == 1.5Hz
end

@testset "merge_sats errors on mismatched doppler estimator type" begin
    gpsl1 = GPSL1CA()
    estimator = ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sat1 = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator = estimator)
    track_state = TrackState(gpsl1, sat1; doppler_estimator = estimator)

    # Default estimator is ConventionalAssistedPLLAndDLL — produces a
    # different concrete state type than the TrackState's estimator.
    bad_sat = TrackedSat(gpsl1, 2, 0.25, 200.0Hz)
    @test_throws ArgumentError merge_sats(track_state, bad_sat)
end

end
