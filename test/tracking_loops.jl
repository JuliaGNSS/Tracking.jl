module TrackingLoopsTest

using Test: @test, @testset, @inferred, @test_throws
using GNSSSignals: CBOC
using Logging: NullLogger, with_logger
import Random
using Unitful: Hz, ms
import Tracking
import GNSSSignals
using GNSSSignals:
    GPSL1CA,
    GalileoE1C,
    get_chip_rate_to_center_frequency_ratio,
    get_chip_rate,
    get_code_length,
    get_secondary_code
using TrackingLoopFilters:
    ThirdOrderAssistedBilinearLF, ThirdOrderBilinearLF, SecondOrderBilinearLF
using StaticArrays: SVector
using Tracking:
    get_bit_synchronizer,
    TrackedSignal,
    TrackedSat,
    TrackState,
    BandMeasurement,
    step_tracking_loops,
    step_tracking_loops!,
    is_bit_or_secondary_code_synced,
    add_satellite!,
    append_correlator_bank_output!,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_beamformed_prompt,
    get_beamformed_prompts,
    get_tracked_sat,
    get_tracked_sats,
    merge_sats
import TrackingLoops
using TrackingLoops:
    ScalarLoopState,
    init_loop_state,
    ScalarTrackingLoop,
    LoopBandwidths,
    with_correlator_outputs,
    get_default_correlator_bank,
    CorrelatorBankOutput

# The estimate step only reads `sampling_frequency`.
_meas_l1(fs) = (L1 = BandMeasurement(ComplexF64[], fs),)

# A signal carrying one completed integration as a `CorrelatorBankOutput`.
_completed_signal(sig, correlator_bank, num_samples) = TrackedSignal(
    sig;
    correlator_bank_outputs = [
        CorrelatorBankOutput(correlator_bank, num_samples, num_samples),
    ],
)

@testset "Satellite Conventional PLL and DLL" begin
    pll_and_dll = @inferred ScalarLoopState(
        reference_carrier_doppler = 500.0Hz,
        reference_code_doppler = 100.0Hz,
    )

    @test pll_and_dll.reference_carrier_doppler == 500.0Hz
    @test pll_and_dll.reference_code_doppler == 100.0Hz
    @test pll_and_dll.carrier_loop_filter == ThirdOrderBilinearLF()
    @test pll_and_dll.code_loop_filter == SecondOrderBilinearLF()
    @test pll_and_dll.bandwidths == LoopBandwidths()

    gpsl1 = GPSL1CA()
    tracked_sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz)
    from_tracked_sat = @inferred init_loop_state(
        ScalarTrackingLoop(
            ThirdOrderBilinearLF,
            SecondOrderBilinearLF;
            wide_carrier_loop_filter_bandwidth = 25.0Hz,
            code_loop_filter_bandwidth = 2.0Hz,
        ),
        tracked_sat,
    )
    @test from_tracked_sat.bandwidths.wide_carrier == 25.0Hz
    @test from_tracked_sat.bandwidths.code == 2.0Hz

    preserved = @inferred ScalarLoopState(from_tracked_sat)
    @test preserved.bandwidths == from_tracked_sat.bandwidths

    overridden = @inferred ScalarLoopState(
        from_tracked_sat;
        bandwidths = LoopBandwidths(; wide_carrier = 30.0Hz, code = 2.0Hz),
    )
    @test overridden.bandwidths.wide_carrier == 30.0Hz
    @test overridden.bandwidths.code == 2.0Hz
end

@testset "Conventional PLL and DLL" begin
    sampling_frequency = 5e6Hz

    gpsl1 = GPSL1CA()

    carrier_doppler = 100.0Hz
    prn = 1
    code_phase = 0.5

    tracking_loop = ScalarTrackingLoop()

    tracked_sat = TrackedSat(gpsl1, prn, code_phase, carrier_doppler; tracking_loop)

    # No completed integration: nothing changes.
    num_samples = 2000
    correlator_bank = with_correlator_outputs(
        get_default_correlator_bank(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    sat_state_after_small_integration = TrackedSat(
        tracked_sat;
        signals = (
            TrackedSignal(
                only(tracked_sat.signals);
                num_integrated_samples = num_samples,
                correlator_bank,
            ),
        ),
    )

    track_state = TrackState(gpsl1, sat_state_after_small_integration; tracking_loop)

    new_track_state =
        @inferred step_tracking_loops(track_state, _meas_l1(sampling_frequency))

    @test get_carrier_doppler(new_track_state) == carrier_doppler
    @test get_code_doppler(new_track_state) ==
          get_chip_rate_to_center_frequency_ratio(gpsl1) * carrier_doppler
    @test get_last_fully_integrated_beamformed_prompt(new_track_state) == 0.0
    @test isempty(get_beamformed_prompts(get_tracked_sat(new_track_state, prn)))

    # One completed integration.
    num_samples = 5000
    sat_state_after_full_integration = TrackedSat(
        tracked_sat;
        signals = (
            _completed_signal(only(tracked_sat.signals), correlator_bank, num_samples),
        ),
    )
    track_state = TrackState(gpsl1, sat_state_after_full_integration; tracking_loop)

    new_track_state_after_full_integration =
        @inferred step_tracking_loops(track_state, _meas_l1(sampling_frequency))

    @test get_carrier_doppler(new_track_state_after_full_integration) ==
          100.23888678999913Hz
    @test get_code_doppler(new_track_state_after_full_integration) == -0.1609215874910254Hz
    @test get_last_fully_integrated_beamformed_prompt(
        new_track_state_after_full_integration,
    ) == 0.4 + 0.004im
    let prompts = get_beamformed_prompts(
            get_tracked_sat(new_track_state_after_full_integration, prn),
        )
        @test length(prompts) == 1
        @test prompts[1] == 0.4 + 0.004im
    end
end

@testset "loop bandwidth follows the actual integration length" begin
    # The carrier cap pairs with the time the integration covered (from its sample
    # count), not the intended length: a single-code-period integration after sync (a
    # mid-update sync, or the truncated first post-sync integration) gets the full
    # single-period bandwidth.
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    code_phase = 0.5
    correlator_bank = with_correlator_outputs(
        get_default_correlator_bank(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    tracking_loop = ScalarTrackingLoop()

    # Synced at secondary phase 0, polarity +1.
    synced(::TrackingLoops.BitSynchronizer{B}) where {B} = TrackingLoops.BitSynchronizer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(+1),
        complex(0.0, 0.0),
        0,
        Float32[],
        TrackingLoops.SyncHypothesisAccumulators(),
    )

    doppler_after(preferred, found, num_integrated_samples) = begin
        sat = TrackedSat(gpsl1, 1, code_phase, carrier_doppler; tracking_loop)
        base = only(sat.signals)
        sig = TrackedSignal(
            base;
            bit_synchronizer = found ? synced(get_bit_synchronizer(base)) :
                               get_bit_synchronizer(base),
            preferred_num_code_periods_to_integrate = preferred,
            correlator_bank_outputs = [
                CorrelatorBankOutput(
                    correlator_bank,
                    num_integrated_samples,
                    num_integrated_samples,
                ),
            ],
        )
        ts = TrackState(gpsl1, TrackedSat(sat; signals = (sig,)); tracking_loop)
        get_carrier_doppler(step_tracking_loops(ts, _meas_l1(sampling_frequency)))
    end

    one_code_period = 5000                      # one 1 ms L1 C/A code period at 5 MHz
    # Capping against the intended 20 code periods would give 4.5 Hz.
    @test doppler_after(20, true, one_code_period) ==
          doppler_after(1, false, one_code_period)

    # 20 covered code periods are capped at 0.09 / 20 ms either way.
    @test doppler_after(20, true, 20 * one_code_period) ==
          doppler_after(1, true, 20 * one_code_period)
end

# The DLL bandwidth is absolute: nothing that sets it scales with the
# integration length, so only its `BL·Δt` stability cap may narrow it. A `1/N`
# scaling would take a 20 ms L1 C/A integration to 0.05 Hz where 0.9 Hz is
# stable.
@testset "code loop bandwidth is not narrowed by the integration length" begin
    bw = 1.0Hz
    l1ca_period = 1ms  # 1023 chips at 1.023 Mcps

    # The cap binds only past 18 ms.
    @test TrackingLoops.effective_code_loop_filter_bandwidth(bw, l1ca_period) == bw
    @test TrackingLoops.effective_code_loop_filter_bandwidth(bw, 10 * l1ca_period) == bw
    @test TrackingLoops.effective_code_loop_filter_bandwidth(bw, 20 * l1ca_period) ≈ 0.9Hz

    for num_code_periods in (20, 100, 1500)
        integration_time = num_code_periods * l1ca_period
        @test TrackingLoops.effective_code_loop_filter_bandwidth(bw, integration_time) *
              integration_time ≈ TrackingLoops.MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT
    end

    @test TrackingLoops.effective_code_loop_filter_bandwidth(0.25Hz, 20 * l1ca_period) ==
          0.25Hz

    # Through the tracking loop at 4 ms, where the cap is 4.5 Hz: 5 and 10 Hz both
    # saturate (a `1/N` scaling would make them differ).
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    correlator_bank = with_correlator_outputs(
        get_default_correlator_bank(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    code_doppler_after(num_integrated_samples, code_bandwidth) = begin
        tracking_loop = ScalarTrackingLoop(; code_loop_filter_bandwidth = code_bandwidth)
        sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; tracking_loop)
        sig = TrackedSignal(
            only(sat.signals);
            correlator_bank_outputs = [
                CorrelatorBankOutput(
                    correlator_bank,
                    num_integrated_samples,
                    num_integrated_samples,
                ),
            ],
        )
        ts = TrackState(gpsl1, TrackedSat(sat; signals = (sig,)); tracking_loop)
        get_code_doppler(step_tracking_loops(ts, _meas_l1(sampling_frequency)))
    end
    one_code_period = 5000
    @test code_doppler_after(4 * one_code_period, 5.0Hz) ==
          code_doppler_after(4 * one_code_period, 10.0Hz)
    @test code_doppler_after(4 * one_code_period, 1.0Hz) !=
          code_doppler_after(4 * one_code_period, 4.0Hz)
end

@testset "Per-satellite bandwidths drive the loop filters" begin
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    code_phase = 0.5
    num_samples = 5000
    correlator_bank = with_correlator_outputs(
        get_default_correlator_bank(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )

    # Same correlator for both: only the bandwidths differ.
    tracking_loop = ScalarTrackingLoop()
    sat1_initial = TrackedSat(gpsl1, 1, code_phase, carrier_doppler; tracking_loop)
    sat2_initial = TrackedSat(gpsl1, 2, code_phase, carrier_doppler; tracking_loop)
    sat1 = TrackedSat(
        sat1_initial;
        signals = (
            _completed_signal(only(sat1_initial.signals), correlator_bank, num_samples),
        ),
    )
    sat2_pre = TrackedSat(
        sat2_initial;
        signals = (
            _completed_signal(only(sat2_initial.signals), correlator_bank, num_samples),
        ),
    )
    sat1_de = sat1.loop_state
    sat2_de = ScalarLoopState(
        sat2_pre.loop_state;
        bandwidths = LoopBandwidths(; wide_carrier = 36.0Hz, code = 2.0Hz),
    )
    sat2 = TrackedSat(
        sat2_pre.prn,
        sat2_pre.code_phase,
        sat2_pre.code_doppler,
        sat2_pre.carrier_phase,
        sat2_pre.carrier_doppler,
        sat2_pre.next_sample_index,
        sat2_pre.signals,
        sat2_de,
    )
    tracked = [sat1, sat2]

    @test sat1_de.bandwidths.wide_carrier == 50.0Hz
    @test sat2_de.bandwidths.wide_carrier == 36.0Hz
    @test sat2_de.bandwidths.code == 2.0Hz

    track_state = TrackState(gpsl1, tracked; tracking_loop)
    new_track_state =
        @inferred step_tracking_loops(track_state, _meas_l1(sampling_frequency))

    # Sat 1 matches the testset above.
    @test get_carrier_doppler(new_track_state, 1) == 100.23888678999913Hz
    @test get_code_doppler(new_track_state, 1) == -0.1609215874910254Hz

    @test get_carrier_doppler(new_track_state, 2) != get_carrier_doppler(new_track_state, 1)
    @test get_code_doppler(new_track_state, 2) != get_code_doppler(new_track_state, 1)
end

@testset "Bandwidths propagate through ScalarTrackingLoop constructor" begin
    gpsl1 = GPSL1CA()
    tracked_sat = TrackedSat(gpsl1, 1, 0.5, 100.0Hz)
    tracking_loop = ScalarTrackingLoop(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    de_state = init_loop_state(tracking_loop, tracked_sat)
    @test de_state.bandwidths.wide_carrier == 22.0Hz
    @test de_state.bandwidths.code == 1.5Hz

    bumped = ScalarTrackingLoop(tracking_loop; wide_carrier_loop_filter_bandwidth = 30.0Hz)
    @test init_loop_state(bumped, tracked_sat).bandwidths.wide_carrier == 30.0Hz
    @test init_loop_state(bumped, tracked_sat).bandwidths.code == 1.5Hz  # unchanged
    preserved = ScalarTrackingLoop(tracking_loop)
    @test init_loop_state(preserved, tracked_sat).bandwidths == de_state.bandwidths
end

@testset "merge_sats with matching tracking loop carries through bandwidths" begin
    gpsl1 = GPSL1CA()
    tracking_loop = ScalarTrackingLoop(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sat1 = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; tracking_loop)
    track_state = TrackState(gpsl1, sat1; tracking_loop)

    # The slot type pins the loop-state type.
    sat2 = TrackedSat(gpsl1, 2, 0.25, 200.0Hz; tracking_loop)
    merged = merge_sats(track_state, sat2)

    de_state = get_tracked_sats(merged)[2].loop_state
    @test de_state.bandwidths.wide_carrier == 22.0Hz
    @test de_state.bandwidths.code == 1.5Hz
end

@testset "merge_sats errors on mismatched tracking loop type" begin
    gpsl1 = GPSL1CA()
    tracking_loop = ScalarTrackingLoop(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sat1 = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; tracking_loop)
    track_state = TrackState(gpsl1, sat1; tracking_loop)

    # The default `ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)` has another state type.
    bad_sat = TrackedSat(gpsl1, 2, 0.25, 200.0Hz)
    @test_throws ArgumentError merge_sats(track_state, bad_sat)
end

# A tracking loop that keeps what each integration hands the carrier discriminators,
# around the loop it wraps.
struct PromptRecordingTrackingLoop{E<:TrackingLoops.AbstractTrackingLoop} <:
       TrackingLoops.AbstractTrackingLoop
    inner::E
    integrations::Vector{Tuple{ComplexF64,Int8}}
end
TrackingLoops.init_loop_state(
    e::PromptRecordingTrackingLoop,
    driver_signal,
    carrier_doppler,
    code_doppler,
) = TrackingLoops.init_loop_state(e.inner, driver_signal, carrier_doppler, code_doppler)
function TrackingLoops.step_loop!(
    e::PromptRecordingTrackingLoop,
    state,
    integration::TrackingLoops.PredetectionIntegration,
    nco_history,
    application_sample_index::Union{Nothing,Integer},
)
    push!(e.integrations, (integration.previous_prompt, integration.carrier_phase_polarity))
    TrackingLoops.step_loop!(
        e.inner,
        state,
        integration,
        nco_history,
        application_sample_index,
    )
end

@testset "An integration hands the FLL its previous prompt only at an unchanged polarity" begin
    sampling_frequency = 5e6Hz
    # One integration per call; returns the integrations' (previous prompt, polarity)
    # and the prompts left on the signal.
    function step_integrations(signal, synced, code_periods_per_integration)
        tracking_loop = PromptRecordingTrackingLoop(
            ScalarTrackingLoop(ThirdOrderAssistedBilinearLF),
            Tuple{ComplexF64,Int8}[],
        )
        track_state = TrackState(; signal, tracking_loop)
        track_state = add_satellite!(track_state; prn = 1, carrier_doppler = 100.0Hz)
        sats = get_tracked_sats(track_state)
        sat = sats[1]
        sats[1] = TrackedSat(
            sat;
            signals = map(
                s -> TrackedSignal(s; bit_synchronizer = synced(get_bit_synchronizer(s))),
                sat.signals,
            ),
        )
        correlator_bank = get_default_correlator_bank(signal)
        correlator_bank = with_correlator_outputs(
            correlator_bank,
            map(_ -> 1000.0 + 10im, TrackingLoops.get_correlator_outputs(correlator_bank)),
        )
        one_code_period =
            Int(sampling_frequency * get_code_length(signal) / get_chip_rate(signal))
        prompts = ComplexF64[]
        for num_code_periods in code_periods_per_integration
            sat = get_tracked_sats(track_state)[1]
            samples = num_code_periods * one_code_period
            append_correlator_bank_output!(
                only(sat.signals),
                CorrelatorBankOutput(correlator_bank, samples, samples),
            )
            step_tracking_loops!(track_state, _meas_l1(sampling_frequency))
            push!(
                prompts,
                get_last_fully_integrated_beamformed_prompt(
                    get_tracked_sats(track_state)[1],
                ),
            )
        end
        tracking_loop.integrations, prompts
    end

    # A bit buffer synced at polarity -1, on the first code period of a bit.
    synced(::TrackingLoops.BitSynchronizer{B}) where {B} = TrackingLoops.BitSynchronizer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(-1),
        complex(0.0, 0.0),
        0,
        Float32[],
        TrackingLoops.SyncHypothesisAccumulators(),
    )

    # A data signal's prompt sign stays unknown, so every integration but the first
    # chains from the one before, across a change of integration length too.
    integrations, prompts = step_integrations(GPSL1CA(), synced, (1, 1, 2, 2))
    @test integrations[1][1] == 0
    @test integrations[2][1] == prompts[1]
    @test integrations[3][1] == prompts[2]
    @test integrations[4][1] == prompts[3]
    @test all(r -> r[2] == 0, integrations)

    # A synced pilot knows its prompt's sign, the sync polarity. The FLL
    # compares no prompts across the switch from the unknown sign.
    integrations, prompts = step_integrations(GalileoE1C(CBOC), synced, (1, 1))
    @test integrations[1] == (0, -1)
    @test integrations[2] == (prompts[1], -1)
end

# A tracking loop that keeps each integration's carrier phase polarity and prompt.
struct PolarityRecordingTrackingLoop{E<:TrackingLoops.AbstractTrackingLoop} <:
       TrackingLoops.AbstractTrackingLoop
    inner::E
    integrations::Vector{Tuple{Int8,ComplexF64}}
end
TrackingLoops.init_loop_state(
    e::PolarityRecordingTrackingLoop,
    driver_signal,
    carrier_doppler,
    code_doppler,
) = TrackingLoops.init_loop_state(e.inner, driver_signal, carrier_doppler, code_doppler)
function TrackingLoops.step_loop!(
    e::PolarityRecordingTrackingLoop,
    state,
    integration::TrackingLoops.PredetectionIntegration,
    nco_history,
    application_sample_index::Union{Nothing,Integer},
)
    push!(
        e.integrations,
        (
            integration.carrier_phase_polarity,
            TrackingLoops.get_prompt(integration.beamformed_correlator_bank),
        ),
    )
    TrackingLoops.step_loop!(
        e.inner,
        state,
        integration,
        nco_history,
        application_sample_index,
    )
end

@testset "A synced pilot's carrier phase polarity is its prompt's sign" begin
    # Pre-sync the replica carries chip 0 on every code period, so the prompt is
    # sign · chip(k) · chip(0); synced, it is the sign. E1C's CS25 starts with a
    # negative chip, so reading the pre-sync prompts as they come flips the
    # sign, a half-cycle PLL error.
    e1c = GalileoE1C(CBOC)
    prn = 1
    chip(k) = GNSSSignals.secondary_value(get_secondary_code(e1c), prn, k)
    @test chip(0) == -1
    sampling_frequency = 16.368e6Hz
    num_samples = 65472 # 4 ms, one code period
    correlator_output(prompt, sample_index) = CorrelatorBankOutput(
        with_correlator_outputs(
            get_default_correlator_bank(e1c),
            ComplexF64(prompt * num_samples) .* SVector(0.2, 0.7, 1.0, 0.7, 0.2),
        ),
        num_samples,
        sample_index,
    )
    rng = Random.Xoshiro(1)
    for sign in (1, -1)
        tracking_loop = PolarityRecordingTrackingLoop(
            ScalarTrackingLoop(ThirdOrderAssistedBilinearLF),
            Tuple{Int8,ComplexF64}[],
        )
        track_state = TrackState(; signal = e1c, tracking_loop)
        track_state = add_satellite!(track_state; prn, carrier_doppler = 0.0Hz)
        noisy(prompt) = prompt + 0.05 * randn(rng, ComplexF64)
        step!(prompt, k) = begin
            append_correlator_bank_output!(
                get_tracked_sat(track_state, prn),
                correlator_output(prompt, k),
            )
            with_logger(NullLogger()) do
                step_tracking_loops!(track_state, _meas_l1(sampling_frequency))
            end
        end
        k = 0
        while !is_bit_or_secondary_code_synced(get_tracked_sat(track_state, prn)) && k < 200
            step!(noisy(sign * chip(k) * chip(0)), (k += 1) * num_samples)
        end
        @test is_bit_or_secondary_code_synced(get_tracked_sat(track_state, prn))
        num_pre_sync = length(tracking_loop.integrations)
        @test all(r -> r[1] == 0, tracking_loop.integrations)
        for _ = 1:10
            step!(noisy(sign), (k += 1) * num_samples)
        end
        post_sync = tracking_loop.integrations[(num_pre_sync+1):end]
        @test all(r -> r[1] == sign, post_sync)
        @test all(r -> Base.sign(real(r[2])) == sign, post_sync)
    end
end

end
