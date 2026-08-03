module VectorPLLAndDLLTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA, GalileoE1B, get_code_center_frequency_ratio
using TrackingLoopFilters: ThirdOrderAssistedBilinearLF, SecondOrderBilinearLF, filter_loop
using StaticArrays: SVector
using Dictionaries: dictionary
using Tracking:
    SatVectorPLLAndDLL,
    VectorPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    TrackedSignal,
    TrackedSat,
    TrackState,
    BandMeasurement,
    CorrelatorOutput,
    init_estimator_state,
    estimate_dopplers_and_filter_prompt,
    get_carrier_doppler,
    get_code_doppler,
    get_sat_state,
    get_sat_states,
    add_satellite!,
    get_doppler_estimator_state,
    get_default_correlator,
    update_accumulator,
    dll_disc,
    pll_disc,
    enable_vt!,
    disable_vt!,
    reset_code_discr_acc!,
    reset_carrier_discr_acc!,
    mean_code_discriminator,
    mean_carrier_discriminator,
    set_code_freq_updates!,
    set_carrier_freq_updates!,
    reset_loop_filters!

_meas_l1(fs) = (L1 = BandMeasurement(ComplexF64[], fs),)

# Rebuild `sat` with its driver signal carrying one completed integration as a
# `CorrelatorOutput` record — the estimate phase folds over `correlator_outputs`.
function _with_full_integration(sat, accumulators, num_samples)
    correlator =
        update_accumulator(get_default_correlator(GPSL1CA()), SVector(accumulators...))
    TrackedSat(
        sat;
        signals = (
            TrackedSignal(
                only(sat.signals);
                correlator_outputs = [
                    CorrelatorOutput(correlator, num_samples, num_samples),
                ],
            ),
        ),
    )
end

# Swap in a modified per-sat estimator state (same concrete type).
_with_state(sat, state) = TrackedSat(sat; doppler_estimator_state = state)

@testset "Satellite Vector PLL and DLL" begin
    state = @inferred SatVectorPLLAndDLL(
        init_carrier_doppler = 500.0Hz,
        init_code_doppler = 100.0Hz,
    )

    @test state.init_carrier_doppler == 500.0Hz
    @test state.init_code_doppler == 100.0Hz
    @test state.carrier_loop_filter == ThirdOrderAssistedBilinearLF()
    @test state.code_loop_filter == SecondOrderBilinearLF()
    @test state.carrier_loop_filter_bandwidth == 18.0Hz
    @test state.code_loop_filter_bandwidth == 1.0Hz
    @test state.code_discr_acc == (0, 0.0)
    @test state.code_freq_update == 0.0Hz
    @test state.carrier_discr_acc == (0, 0.0Hz)
    @test state.carrier_freq_update == 0.0Hz
    @test state.vt_on == false

    # Kwarg-update constructor preserves what isn't overridden.
    updated = @inferred SatVectorPLLAndDLL(
        state;
        code_discr_acc = (2, 0.5),
        carrier_freq_update = 3.0Hz,
        vt_on = true,
    )
    @test updated.code_discr_acc == (2, 0.5)
    @test updated.carrier_freq_update == 3.0Hz
    @test updated.vt_on == true
    @test updated.init_carrier_doppler == 500.0Hz
    @test updated.carrier_loop_filter_bandwidth == 18.0Hz
end

@testset "Vector PLL and DLL estimator seeding" begin
    gpsl1 = GPSL1CA()

    # Auto bandwidths resolve from the driver signal at seeding time, sized
    # exactly like the conventional estimator (the scalar fallback loop).
    estimator = VectorPLLAndDLL()
    sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator = estimator)
    state = get_doppler_estimator_state(sat)
    @test state isa SatVectorPLLAndDLL
    @test state.carrier_loop_filter_bandwidth == 18.0Hz
    @test state.code_loop_filter_bandwidth == 1.0Hz
    @test state.init_carrier_doppler == 100.0Hz
    @test state.vt_on == false

    # Explicit bandwidths override the auto-sizing; the kwarg-update
    # constructor tweaks them afterwards.
    explicit = VectorPLLAndDLL(;
        carrier_loop_filter_bandwidth = 12.0Hz,
        code_loop_filter_bandwidth = 1.0Hz,
    )
    explicit_state = @inferred init_estimator_state(explicit, sat)
    @test explicit_state.carrier_loop_filter_bandwidth == 12.0Hz
    @test explicit_state.code_loop_filter_bandwidth == 1.0Hz
    bumped = VectorPLLAndDLL(explicit; carrier_loop_filter_bandwidth = 15.0Hz)
    @test bumped.carrier_loop_filter_bandwidth == 15.0Hz
    @test bumped.code_loop_filter_bandwidth == 1.0Hz
end

@testset "Scalar fallback matches the conventional assisted PLL and DLL" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    prn = 1
    code_phase = 0.5
    num_samples = 5000
    accumulators = (1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im)

    # With vt_on = false (the initial state), the vector estimator must
    # reproduce the conventional assisted estimator's Doppler update
    # exactly — both auto-size their bandwidths identically.
    vector_estimator = VectorPLLAndDLL()
    conventional_estimator = ConventionalAssistedPLLAndDLL()

    results = map((vector_estimator, conventional_estimator)) do doppler_estimator
        sat = TrackedSat(gpsl1, prn, code_phase, carrier_doppler; doppler_estimator)
        sat = _with_full_integration(sat, accumulators, num_samples)
        track_state = TrackState(gpsl1, sat; doppler_estimator)
        new_track_state = @inferred estimate_dopplers_and_filter_prompt(
            track_state,
            _meas_l1(sampling_frequency),
        )
        (get_carrier_doppler(new_track_state, prn), get_code_doppler(new_track_state, prn))
    end
    @test results[1] == results[2]

    # The scalar fallback accumulates nothing, and neither NCO-correction
    # field is written by the estimator (both are navigation-filter-owned).
    vector_sat = TrackedSat(
        gpsl1,
        prn,
        code_phase,
        carrier_doppler;
        doppler_estimator = vector_estimator,
    )
    vector_sat = _with_full_integration(vector_sat, accumulators, num_samples)
    track_state = TrackState(gpsl1, vector_sat; doppler_estimator = vector_estimator)
    new_track_state =
        estimate_dopplers_and_filter_prompt(track_state, _meas_l1(sampling_frequency))
    state = get_doppler_estimator_state(get_sat_state(new_track_state, prn))
    @test state.code_discr_acc == (0, 0.0)
    @test state.carrier_discr_acc == (0, 0.0Hz)
    @test state.code_freq_update == 0.0Hz
    @test state.carrier_freq_update == 0.0Hz
    # But the scalar loop still updates the Doppler (the DLL/PLL ran).
    @test get_carrier_doppler(new_track_state, prn) != carrier_doppler
end

@testset "Vector loop closure applies the NCO corrections" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    prn = 1
    code_phase = 0.5
    num_samples = 5000
    accumulators = (1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im)
    nav_carrier_freq_update = 5.0Hz
    nav_code_freq_update = -0.25Hz

    doppler_estimator = VectorPLLAndDLL()
    sat = TrackedSat(gpsl1, prn, code_phase, carrier_doppler; doppler_estimator)
    init_code_doppler = get_code_doppler(sat)
    sat = _with_full_integration(sat, accumulators, num_samples)
    track_state = TrackState(gpsl1, sat; doppler_estimator)
    enable_vt!(track_state, (prn,))
    set_carrier_freq_updates!(track_state, dictionary((prn => nav_carrier_freq_update,)))
    set_code_freq_updates!(track_state, dictionary((prn => nav_code_freq_update,)))

    new_track_state = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        _meas_l1(sampling_frequency),
    )
    state = get_doppler_estimator_state(get_sat_state(new_track_state, prn))

    # Expected values, computed on the normalized correlator (the post-corr
    # filter is an identity for the first prompt).
    normalized_correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(accumulators...) ./ num_samples,
    )
    integration_time = num_samples / sampling_frequency
    pll_discriminator = pll_disc(gpsl1, normalized_correlator)
    dll_discriminator =
        dll_disc(gpsl1, normalized_correlator, init_code_doppler, sampling_frequency)
    # The FLL branch is driven by the navigation filter's carrier update.
    expected_carrier_freq_update, _ = filter_loop(
        ThirdOrderAssistedBilinearLF(),
        (pll_discriminator, nav_carrier_freq_update),
        integration_time,
        18.0Hz,
    )

    @test get_carrier_doppler(new_track_state, prn) ==
          carrier_doppler + expected_carrier_freq_update
    # The code Doppler follows the navigation filter's update directly
    # (plus carrier aiding) — the code loop filter is bypassed.
    @test get_code_doppler(new_track_state, prn) ==
          init_code_doppler +
          nav_code_freq_update +
          expected_carrier_freq_update * get_code_center_frequency_ratio(gpsl1)

    # Discriminator outputs are accumulated for the navigation filter; the
    # FLL discriminator is zero here since there is no previous prompt.
    @test state.code_discr_acc == (1, dll_discriminator)
    @test state.carrier_discr_acc == (1, 0.0Hz)
    # The mean accessors divide sum by count (one sample here).
    @test mean_code_discriminator(state) == dll_discriminator
    @test mean_carrier_discriminator(state) == 0.0Hz
    # The navigation filter's corrections survive the update untouched.
    @test state.code_freq_update == nav_code_freq_update
    @test state.carrier_freq_update == nav_carrier_freq_update

    # The accumulator-reset functions bring the accumulators back to zero.
    reset_code_discr_acc!(new_track_state)
    reset_carrier_discr_acc!(new_track_state)
    state = get_doppler_estimator_state(get_sat_state(new_track_state, prn))
    @test state.code_discr_acc == (0, 0.0)
    @test state.carrier_discr_acc == (0, 0.0Hz)
    # With count == 0 the mean accessors return `nothing`.
    @test mean_code_discriminator(state) === nothing
    @test mean_carrier_discriminator(state) === nothing
end

@testset "Vector tracking state management" begin
    gpsl1 = GPSL1CA()
    doppler_estimator = VectorPLLAndDLL()
    sats = [
        TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator),
        TrackedSat(gpsl1, 2, 0.25, 200.0Hz; doppler_estimator),
    ]
    track_state = TrackState(gpsl1, sats; doppler_estimator)
    _state(prn) = get_doppler_estimator_state(get_sat_state(track_state, prn))

    # PRN 1 in lock joins the vector loop; PRN 2 keeps the scalar fallback.
    enable_vt!(track_state, (1,))
    @test _state(1).vt_on == true
    @test _state(2).vt_on == false

    # Satellites outside the addressed set are left as they are — an empty set
    # drops nobody, so a satellite in an outage keeps being steered until it is
    # explicitly disabled.
    enable_vt!(track_state, ())
    @test _state(1).vt_on == true
    @test _state(2).vt_on == false

    # PRN 2 joins once in lock; re-enabling PRN 1 is a no-op.
    enable_vt!(track_state, (1, 2))
    @test _state(1).vt_on == true
    @test _state(2).vt_on == true

    # `disable_vt!` is the inverse and equally addressed.
    disable_vt!(track_state, (1,))
    @test _state(1).vt_on == false
    @test _state(2).vt_on == true

    # Disabling a satellite that is already out is a no-op too.
    disable_vt!(track_state, (1,))
    @test _state(1).vt_on == false
    @test _state(2).vt_on == true

    # Per-PRN NCO corrections only reach satellites in the vector loop.
    set_code_freq_updates!(track_state, dictionary((1 => 1.0Hz, 2 => 2.0Hz)))
    set_carrier_freq_updates!(track_state, dictionary((1 => 10.0Hz, 2 => 20.0Hz)))
    @test _state(1).code_freq_update == 0.0Hz
    @test _state(1).carrier_freq_update == 0.0Hz
    @test _state(2).code_freq_update == 2.0Hz
    @test _state(2).carrier_freq_update == 20.0Hz
end

@testset "Group-scoped vector tracking state management" begin
    # Two constellations sharing a PRN number: the group-scoped managers must
    # address exactly one group — a PRN alone is ambiguous across groups.
    gpsl1 = GPSL1CA()
    galileo_e1b = GalileoE1B()
    doppler_estimator = VectorPLLAndDLL()
    track_state = TrackState(;
        signals = (gps = (gpsl1,), galileo = (galileo_e1b,)),
        doppler_estimator,
    )
    add_satellite!(
        track_state;
        prn = 5,
        group = :gps,
        code_phase = 0.5,
        carrier_doppler = 100.0Hz,
    )
    add_satellite!(
        track_state;
        prn = 5,
        group = :galileo,
        code_phase = 0.25,
        carrier_doppler = 200.0Hz,
    )
    _state(group) = get_doppler_estimator_state(get_sat_state(track_state, group, 5))

    # Enabling only reaches the addressed group.
    enable_vt!(track_state, :gps, (5,))
    @test _state(:gps).vt_on == true
    @test _state(:galileo).vt_on == false

    # Group-scoped NCO corrections: the other group's dictionary is never
    # consulted (no KeyError for its vt-on sats) and its state is untouched.
    enable_vt!(track_state, :galileo, (5,))
    set_code_freq_updates!(track_state, :gps, dictionary((5 => 1.0Hz,)))
    set_carrier_freq_updates!(track_state, :gps, dictionary((5 => 10.0Hz,)))
    @test _state(:gps).code_freq_update == 1.0Hz
    @test _state(:gps).carrier_freq_update == 10.0Hz
    @test _state(:galileo).code_freq_update == 0.0Hz
    @test _state(:galileo).carrier_freq_update == 0.0Hz

    # Membership changes are group-scoped too.
    disable_vt!(track_state, :gps, (5,))
    @test _state(:gps).vt_on == false
    @test _state(:galileo).vt_on == true

    # Enabling the same PRN in one group must not touch the other group's
    # membership.
    disable_vt!(track_state, :galileo, (5,))
    enable_vt!(track_state, :gps, (5,))
    @test _state(:gps).vt_on == true
    @test _state(:galileo).vt_on == false

    # Group-scoped accumulator resets: only the addressed group is cleared.
    enable_vt!(track_state, :galileo, (5,))
    for group in (:gps, :galileo)
        sats = get_sat_states(track_state, group)
        sats[5] = _with_state(
            sats[5],
            SatVectorPLLAndDLL(
                get_doppler_estimator_state(sats[5]);
                code_discr_acc = (2, 0.5),
                carrier_discr_acc = (2, 4.0Hz),
            ),
        )
    end
    reset_code_discr_acc!(track_state, :gps)
    reset_carrier_discr_acc!(track_state, :gps)
    @test _state(:gps).code_discr_acc == (0, 0.0)
    @test _state(:gps).carrier_discr_acc == (0, 0.0Hz)
    @test _state(:galileo).code_discr_acc == (2, 0.5)
    @test _state(:galileo).carrier_discr_acc == (2, 4.0Hz)
end

@testset "reset_loop_filters! preserves VT flags and zeroes corrections" begin
    gpsl1 = GPSL1CA()
    doppler_estimator = VectorPLLAndDLL()
    sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator)
    dirty_state = SatVectorPLLAndDLL(
        get_doppler_estimator_state(sat);
        code_discr_acc = (3, 0.7),
        carrier_discr_acc = (3, 2.0Hz),
        code_freq_update = 0.5Hz,
        carrier_freq_update = 4.0Hz,
        carrier_loop_filter_bandwidth = 12.0Hz,
        vt_on = true,
    )
    track_state = TrackState(gpsl1, _with_state(sat, dirty_state); doppler_estimator)
    reset_loop_filters!(track_state, 1)
    state = get_doppler_estimator_state(get_sat_state(track_state, 1))
    @test state.code_discr_acc == (0, 0.0)
    @test state.carrier_discr_acc == (0, 0.0Hz)
    @test state.code_freq_update == 0.0Hz
    @test state.carrier_freq_update == 0.0Hz
    # Per-sat bandwidth override and the vt_on flag survive the reset.
    @test state.carrier_loop_filter_bandwidth == 12.0Hz
    @test state.vt_on == true
    # Init Dopplers are re-seeded from the sat's current Dopplers.
    @test state.init_carrier_doppler == 100.0Hz
end

end
