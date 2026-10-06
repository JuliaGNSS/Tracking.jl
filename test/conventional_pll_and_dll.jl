module ConventionalPLLAndDLLTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz, ms
import Tracking
using GNSSSignals: GPSL1CA, get_code_center_frequency_ratio
using TrackingLoopFilters:
    ThirdOrderBilinearLF, ThirdOrderAssistedBilinearLF, SecondOrderBilinearLF, filter_loop
using Random: MersenneTwister
using StaticArrays: SVector
using Dictionaries: Dictionary
using Tracking:
    aid_dopplers,
    _close_loops,
    pll_disc,
    fll_disc,
    SatConventionalPLLAndDLL,
    EarlyPromptLateCorrelator,
    TrackedSignal,
    TrackedSat,
    init_estimator_state,
    ConventionalPLLAndDLL,
    TrackState,
    BandMeasurement,
    estimate_dopplers_and_filter_prompt,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_sat_state,
    get_sat_states,
    update_accumulator,
    get_default_correlator,
    CorrelatorOutput,
    merge_sats

# Stub band measurement: the estimator only reads `sampling_frequency`.
_meas_l1(fs) = (L1 = BandMeasurement(ComplexF64[], fs),)

# A signal carrying one completed integration as a `CorrelatorOutput` record.
_completed_signal(sig, correlator, num_samples) = TrackedSignal(
    sig;
    correlator_outputs = [CorrelatorOutput(correlator, num_samples, num_samples)],
)

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
        signals = (_completed_signal(only(sat_state.signals), correlator, num_samples),),
    )
    track_state = TrackState(gpsl1, sat_state_after_full_integration; doppler_estimator)

    new_track_state_after_full_integration = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        _meas_l1(sampling_frequency),
    )

    # The carrier update is the radian-fed one of before #244 divided by 2π.
    @test get_carrier_doppler(new_track_state_after_full_integration) == 100.0837403735401Hz
    @test get_code_doppler(new_track_state_after_full_integration) == -0.1610223319172975Hz
    @test get_last_fully_integrated_filtered_prompt(
        new_track_state_after_full_integration,
    ) == 0.4 + 0.004im
    # One integration completed -> one filtered prompt recorded, equal to the
    # last_fully_integrated_filtered_prompt value.
    let prompts =
            get_filtered_prompts(get_sat_state(new_track_state_after_full_integration, prn))
        @test length(prompts) == 1
        @test prompts[1] == 0.4 + 0.004im
    end
end

@testset "loop bandwidth follows the record's actual integration length" begin
    # The stability caps must pair with the time a record ACTUALLY
    # integrated (recovered from its sample count), not the intended integration
    # length: a single-block record folded when the bit buffer already reports
    # sync — a mid-fold sync detection with an enlarged
    # `doppler_update_interval`, or the truncated first post-sync integration —
    # must be filtered at the full single-period bandwidth. The Doppler update
    # therefore depends only on the record itself, not on the preferred
    # integration length or the sync state.
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    code_phase = 0.5
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    doppler_estimator = ConventionalPLLAndDLL()

    # Same-typed bit buffer with the bit/secondary sync already found
    # (secondary phase 0, polarity +1).
    synced(::Tracking.BitBuffer{B}) where {B} = Tracking.BitBuffer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(+1),
        complex(0.0, 0.0),
        0,
        Float32[],
        Tracking.PhaseAccumulators(),
    )

    doppler_after(preferred, found, integrated_samples) = begin
        sat = TrackedSat(gpsl1, 1, code_phase, carrier_doppler; doppler_estimator)
        base = only(sat.signals)
        sig = TrackedSignal(
            base;
            bit_buffer = found ? synced(base.bit_buffer) : base.bit_buffer,
            preferred_num_code_blocks_to_integrate = preferred,
            correlator_outputs = [
                CorrelatorOutput(correlator, integrated_samples, integrated_samples),
            ],
        )
        ts = TrackState(gpsl1, TrackedSat(sat; signals = (sig,)); doppler_estimator)
        get_carrier_doppler(
            estimate_dopplers_and_filter_prompt(ts, _meas_l1(sampling_frequency)),
        )
    end

    one_block = 5000                      # one 1 ms L1 C/A code period at 5 MHz
    # Synced with a 20-block preferred length, but the record covered a single
    # block: full bandwidth — exactly the plain single-block baseline. (Capping
    # against the intended length would cut the bandwidth to 4.5 Hz here.)
    @test doppler_after(20, true, one_block) == doppler_after(1, false, one_block)

    # A record that actually covered 20 blocks is capped regardless of the
    # preferred integration length.
    @test doppler_after(20, true, 20 * one_block) == doppler_after(1, true, 20 * one_block)
end

# The DLL bandwidth is capped by its stability product, not scaled by the block
# count (see `effective_code_loop_filter_bandwidth`).
@testset "code loop bandwidth is not narrowed by the integration length" begin
    bw = 1.0Hz
    l1ca_period = 1ms  # 1023 chips at 1.023 Mcps

    # Single 1 ms block, and 10 ms / 20 ms coherent integrations: the cap starts
    # binding only past 18 ms, so the first two run at the full reference.
    @test Tracking.effective_code_loop_filter_bandwidth(bw, l1ca_period) == bw
    @test Tracking.effective_code_loop_filter_bandwidth(bw, 10 * l1ca_period) == bw
    @test Tracking.effective_code_loop_filter_bandwidth(bw, 20 * l1ca_period) ≈ 0.9Hz

    # Where it binds it holds the stability product, not a block ratio.
    for num_blocks in (20, 100, 1500)
        integration_time = num_blocks * l1ca_period
        @test Tracking.effective_code_loop_filter_bandwidth(bw, integration_time) *
              integration_time ≈ Tracking.MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT
    end

    # An explicit bandwidth below the cap is used verbatim, at any length.
    @test Tracking.effective_code_loop_filter_bandwidth(0.25Hz, 20 * l1ca_period) == 0.25Hz

    # End to end through the estimator, on a 4 ms (4-block) integration whose cap
    # sits at 4.5 Hz: two configured bandwidths above the cap must produce the
    # same code Doppler, because both saturate to it. Under a `1/N` scaling they
    # would come out at 5/4 and 10/4 Hz and differ — so this pins "capped" rather
    # than "scaled" through the real fold.
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    code_doppler_after(integrated_samples, code_bandwidth) = begin
        doppler_estimator =
            ConventionalPLLAndDLL(; code_loop_filter_bandwidth = code_bandwidth)
        sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator)
        sig = TrackedSignal(
            only(sat.signals);
            correlator_outputs = [
                CorrelatorOutput(correlator, integrated_samples, integrated_samples),
            ],
        )
        ts = TrackState(gpsl1, TrackedSat(sat; signals = (sig,)); doppler_estimator)
        get_code_doppler(
            estimate_dopplers_and_filter_prompt(ts, _meas_l1(sampling_frequency)),
        )
    end
    one_block = 5000
    @test code_doppler_after(4 * one_block, 5.0Hz) ==
          code_doppler_after(4 * one_block, 10.0Hz)
    # Below the cap the configured value still drives the loop, so it is a cap
    # and not a ceiling on everything.
    @test code_doppler_after(4 * one_block, 1.0Hz) !=
          code_doppler_after(4 * one_block, 4.0Hz)
end

@testset "Per-satellite bandwidths drive the loop filters" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    code_phase = 0.5
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
        signals = (_completed_signal(only(sat1_initial.signals), correlator, num_samples),),
    )
    sat2_pre = TrackedSat(
        sat2_initial;
        signals = (_completed_signal(only(sat2_initial.signals), correlator, num_samples),),
    )
    # Give sat2 different bandwidths via its estimator state.
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
    )

    # Sat 1 matches the baseline result from the previous testset.
    @test get_carrier_doppler(new_track_state, 1) == 100.0837403735401Hz
    @test get_code_doppler(new_track_state, 1) == -0.1610223319172975Hz

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
    bumped = ConventionalPLLAndDLL(estimator; carrier_loop_filter_bandwidth = 30.0Hz)
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

    # Same estimator: the slot type pins the estimator-state type.
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

    # The default (assisted) estimator yields a different state type.
    bad_sat = TrackedSat(gpsl1, 2, 0.25, 200.0Hz)
    @test_throws ArgumentError merge_sats(track_state, bad_sat)
end

# Likewise the carrier bandwidth, at 0.09 (issue #245).
@testset "carrier loop bandwidth is capped, not scaled, by the integration length" begin
    bw = 18.0Hz
    @test Tracking.effective_carrier_loop_filter_bandwidth(bw, 1ms) == bw
    @test Tracking.effective_carrier_loop_filter_bandwidth(bw, 4ms) == bw
    @test Tracking.effective_carrier_loop_filter_bandwidth(bw, 5ms) ≈ bw
    @test Tracking.effective_carrier_loop_filter_bandwidth(bw, 10ms) ≈ 9.0Hz
    @test Tracking.effective_carrier_loop_filter_bandwidth(bw, 20ms) ≈ 4.5Hz
    for integration_time in (10ms, 20ms, 100ms, 1500ms)
        @test Tracking.effective_carrier_loop_filter_bandwidth(bw, integration_time) *
              integration_time ≈ Tracking.MAX_CARRIER_LOOP_BANDWIDTH_TIME_PRODUCT
    end
    @test Tracking.effective_carrier_loop_filter_bandwidth(2.0Hz, 20ms) == 2.0Hz

    # Through the estimator: a 4 ms integration runs at the full bandwidth (a
    # `1/N` scaling would give 4.5 Hz); at 20 ms bandwidths above the cap saturate.
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    carrier_doppler_after(integrated_samples, carrier_bandwidth) = begin
        doppler_estimator =
            ConventionalPLLAndDLL(; carrier_loop_filter_bandwidth = carrier_bandwidth)
        sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator)
        sig = TrackedSignal(
            only(sat.signals);
            correlator_outputs = [
                CorrelatorOutput(correlator, integrated_samples, integrated_samples),
            ],
        )
        ts = TrackState(gpsl1, TrackedSat(sat; signals = (sig,)); doppler_estimator)
        get_carrier_doppler(
            estimate_dopplers_and_filter_prompt(ts, _meas_l1(sampling_frequency)),
        )
    end
    one_block = 5000
    # `pll_disc` is scale invariant, so the raw correlator will do.
    freq_update, _ =
        filter_loop(ThirdOrderBilinearLF(), pll_disc(gpsl1, correlator), 4ms, 18.0Hz)
    @test carrier_doppler_after(4 * one_block, 18.0Hz) ≈ 100.0Hz + freq_update
    @test carrier_doppler_after(20 * one_block, 10.0Hz) ==
          carrier_doppler_after(20 * one_block, 20.0Hz)
end

@testset "Carrier loop has the configured noise bandwidth (issue #244)" begin
    # Closed loop with white phase noise σ_n: the NCO jitter must be
    # σ² = σ_n² · 2 · BL · T. Fed radians, it was ≈5.6× wider.
    gpsl1 = GPSL1CA()
    T = 1ms
    bandwidth = 18.0Hz
    σ_n = 0.05
    rng = MersenneTwister(1)
    # The carrier loop closes on the phase error in cycles: the closure's
    # carrier update, which the code loop does not touch.
    close_carrier_loop(carrier_loop_filter, correlator, previous_prompt) = _close_loops(
        SatConventionalPLLAndDLL(;
            init_carrier_doppler = 0.0Hz,
            init_code_doppler = 0.0Hz,
            carrier_loop_filter,
        ),
        gpsl1,
        correlator,
        previous_prompt,
        0.0Hz,
        5e6Hz,
        T,
        bandwidth,
        1.0Hz,
    )
    lf = ThirdOrderBilinearLF()
    φ = 0.0          # signal phase minus NCO phase, rad
    acc = 0.0
    n = 200_000
    settle = 1_000
    for i = 1:n
        prompt = cis(φ + σ_n * randn(rng))
        correlator = EarlyPromptLateCorrelator(SVector(prompt, prompt, prompt), 0.5)
        freq_update, _, state = close_carrier_loop(lf, correlator, prompt)
        lf = state.carrier_loop_filter
        φ -= 2π * Float64(freq_update * T)
        i > settle && (acc += φ^2)
    end
    effective_bandwidth = acc / (n - settle) / σ_n^2 / (2 * Float64(T * Hz))
    @test 0.85 < effective_bandwidth / Float64(bandwidth / Hz) < 1.15

    # The FLL-assisted filter gets the phase error in cycles too.
    prompt = cis(0.3)
    previous_prompt = cis(0.1)
    correlator = EarlyPromptLateCorrelator(SVector(prompt, prompt, prompt), 0.5)
    @test first(
        close_carrier_loop(ThirdOrderAssistedBilinearLF(), correlator, previous_prompt),
    ) == first(
        filter_loop(
            ThirdOrderAssistedBilinearLF(),
            (pll_disc(gpsl1, correlator), fll_disc(gpsl1, correlator, previous_prompt, T)),
            T,
            bandwidth,
        ),
    )
end

end
