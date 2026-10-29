module ConventionalPLLAndDLLTest

using Test: @test, @testset, @inferred, @test_throws
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
using TrackingLoopFilters: ThirdOrderBilinearLF, SecondOrderBilinearLF
using StaticArrays: SVector
using Tracking:
    get_bit_buffer,
    TrackedSignal,
    TrackedSat,
    TrackState,
    BandMeasurement,
    estimate_dopplers_and_filter_prompt,
    estimate_dopplers_and_filter_prompt!,
    has_bit_or_secondary_code_been_found,
    add_satellite!,
    append_correlator_output!,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_sat_state,
    get_sat_states,
    merge_sats
import TrackingLoops
using TrackingLoops:
    SatConventionalPLLAndDLL,
    init_estimator_state,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    LoopBandwidths,
    update_accumulator,
    get_default_correlator,
    CorrelatorOutput

# The estimate step only reads `sampling_frequency`.
_meas_l1(fs) = (L1 = BandMeasurement(ComplexF64[], fs),)

# A signal carrying one completed integration as a `CorrelatorOutput`.
_completed_signal(sig, correlator, num_samples) = TrackedSignal(
    sig;
    correlator_outputs = [CorrelatorOutput(correlator, num_samples, num_samples)],
)

@testset "Satellite Conventional PLL and DLL" begin
    pll_and_dll = @inferred SatConventionalPLLAndDLL(
        init_carrier_doppler = 500.0Hz,
        init_code_doppler = 100.0Hz,
    )

    @test pll_and_dll.init_carrier_doppler == 500.0Hz
    @test pll_and_dll.init_code_doppler == 100.0Hz
    @test pll_and_dll.carrier_loop_filter == ThirdOrderBilinearLF()
    @test pll_and_dll.code_loop_filter == SecondOrderBilinearLF()
    @test pll_and_dll.bandwidths == LoopBandwidths()

    gpsl1 = GPSL1CA()
    sat_state = TrackedSat(gpsl1, 1, 0.5, 100.0Hz)
    from_sat_state = @inferred init_estimator_state(
        ConventionalPLLAndDLL(
            ThirdOrderBilinearLF,
            SecondOrderBilinearLF;
            wide_carrier_loop_filter_bandwidth = 25.0Hz,
            code_loop_filter_bandwidth = 2.0Hz,
        ),
        sat_state,
    )
    @test from_sat_state.bandwidths.wide_carrier == 25.0Hz
    @test from_sat_state.bandwidths.code == 2.0Hz

    preserved = @inferred SatConventionalPLLAndDLL(from_sat_state)
    @test preserved.bandwidths == from_sat_state.bandwidths

    overridden = @inferred SatConventionalPLLAndDLL(
        from_sat_state;
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

    doppler_estimator = ConventionalPLLAndDLL()

    sat_state = TrackedSat(gpsl1, prn, code_phase, carrier_doppler; doppler_estimator)

    # No completed integration: nothing changes.
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

    @test get_carrier_doppler(new_track_state) == carrier_doppler
    @test get_code_doppler(new_track_state) ==
          get_chip_rate_to_center_frequency_ratio(gpsl1) * carrier_doppler
    @test get_last_fully_integrated_filtered_prompt(new_track_state) == 0.0
    @test isempty(get_filtered_prompts(get_sat_state(new_track_state, prn)))

    # One completed integration.
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

    @test get_carrier_doppler(new_track_state_after_full_integration) ==
          100.23888678999913Hz
    @test get_code_doppler(new_track_state_after_full_integration) == -0.1609215874910254Hz
    @test get_last_fully_integrated_filtered_prompt(
        new_track_state_after_full_integration,
    ) == 0.4 + 0.004im
    let prompts =
            get_filtered_prompts(get_sat_state(new_track_state_after_full_integration, prn))
        @test length(prompts) == 1
        @test prompts[1] == 0.4 + 0.004im
    end
end

@testset "loop bandwidth follows the record's actual integration length" begin
    # The carrier cap pairs with the time the record covered (from its sample
    # count), not the intended length: a single-block record after sync (a
    # mid-fold sync, or the truncated first post-sync integration) gets the full
    # single-period bandwidth.
    sampling_frequency = 5e6Hz
    gpsl1 = GPSL1CA()
    carrier_doppler = 100.0Hz
    code_phase = 0.5
    correlator = update_accumulator(
        get_default_correlator(gpsl1),
        SVector(1000.0 + 10im, 2000.0 + 20im, 750.0 + 10im),
    )
    doppler_estimator = ConventionalPLLAndDLL()

    # Synced at secondary phase 0, polarity +1.
    synced(::TrackingLoops.BitBuffer{B}) where {B} = TrackingLoops.BitBuffer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(+1),
        complex(0.0, 0.0),
        0,
        Float32[],
        TrackingLoops.PhaseAccumulators(),
    )

    doppler_after(preferred, found, integrated_samples) = begin
        sat = TrackedSat(gpsl1, 1, code_phase, carrier_doppler; doppler_estimator)
        base = only(sat.signals)
        sig = TrackedSignal(
            base;
            bit_buffer = found ? synced(get_bit_buffer(base)) : get_bit_buffer(base),
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
    # Capping against the intended 20 blocks would give 4.5 Hz.
    @test doppler_after(20, true, one_block) == doppler_after(1, false, one_block)

    # 20 covered blocks are capped at 0.09 / 20 ms either way.
    @test doppler_after(20, true, 20 * one_block) == doppler_after(1, true, 20 * one_block)
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

    for num_blocks in (20, 100, 1500)
        integration_time = num_blocks * l1ca_period
        @test TrackingLoops.effective_code_loop_filter_bandwidth(bw, integration_time) *
              integration_time ≈ TrackingLoops.MAX_CODE_LOOP_BANDWIDTH_TIME_PRODUCT
    end

    @test TrackingLoops.effective_code_loop_filter_bandwidth(0.25Hz, 20 * l1ca_period) ==
          0.25Hz

    # Through the estimator at 4 ms, where the cap is 4.5 Hz: 5 and 10 Hz both
    # saturate (a `1/N` scaling would make them differ).
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

    # Same correlator for both: only the bandwidths differ.
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
    sat1_de = sat1.doppler_estimator_state
    sat2_de = SatConventionalPLLAndDLL(
        sat2_pre.doppler_estimator_state;
        bandwidths = LoopBandwidths(; wide_carrier = 36.0Hz, code = 2.0Hz),
    )
    sat2 = TrackedSat(
        sat2_pre.prn,
        sat2_pre.code_phase,
        sat2_pre.code_doppler,
        sat2_pre.carrier_phase,
        sat2_pre.carrier_doppler,
        sat2_pre.resume_index,
        sat2_pre.signals,
        sat2_de,
    )
    tracked = [sat1, sat2]

    @test sat1_de.bandwidths.wide_carrier == 50.0Hz
    @test sat2_de.bandwidths.wide_carrier == 36.0Hz
    @test sat2_de.bandwidths.code == 2.0Hz

    track_state = TrackState(gpsl1, tracked; doppler_estimator)
    new_track_state = @inferred estimate_dopplers_and_filter_prompt(
        track_state,
        _meas_l1(sampling_frequency),
    )

    # Sat 1 matches the testset above.
    @test get_carrier_doppler(new_track_state, 1) == 100.23888678999913Hz
    @test get_code_doppler(new_track_state, 1) == -0.1609215874910254Hz

    @test get_carrier_doppler(new_track_state, 2) != get_carrier_doppler(new_track_state, 1)
    @test get_code_doppler(new_track_state, 2) != get_code_doppler(new_track_state, 1)
end

@testset "Bandwidths propagate through ConventionalPLLAndDLL constructor" begin
    gpsl1 = GPSL1CA()
    sat_state = TrackedSat(gpsl1, 1, 0.5, 100.0Hz)
    estimator = ConventionalPLLAndDLL(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    @test estimator.wide_carrier_loop_filter_bandwidth == 22.0Hz
    @test estimator.code_loop_filter_bandwidth == 1.5Hz
    de_state = init_estimator_state(estimator, sat_state)
    @test de_state.bandwidths.wide_carrier == 22.0Hz
    @test de_state.bandwidths.code == 1.5Hz

    bumped = ConventionalPLLAndDLL(estimator; wide_carrier_loop_filter_bandwidth = 30.0Hz)
    @test bumped.wide_carrier_loop_filter_bandwidth == 30.0Hz
    @test bumped.code_loop_filter_bandwidth == 1.5Hz  # unchanged
    preserved = ConventionalPLLAndDLL(estimator)
    @test preserved.wide_carrier_loop_filter_bandwidth == 22.0Hz
    @test preserved.code_loop_filter_bandwidth == 1.5Hz
end

@testset "merge_sats with matching estimator carries through bandwidths" begin
    gpsl1 = GPSL1CA()
    estimator = ConventionalPLLAndDLL(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sat1 = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator = estimator)
    track_state = TrackState(gpsl1, sat1; doppler_estimator = estimator)

    # The slot type pins the estimator-state type.
    sat2 = TrackedSat(gpsl1, 2, 0.25, 200.0Hz; doppler_estimator = estimator)
    merged = merge_sats(track_state, sat2)

    de_state = get_sat_states(merged)[2].doppler_estimator_state
    @test de_state.bandwidths.wide_carrier == 22.0Hz
    @test de_state.bandwidths.code == 1.5Hz
end

@testset "merge_sats errors on mismatched doppler estimator type" begin
    gpsl1 = GPSL1CA()
    estimator = ConventionalPLLAndDLL(;
        wide_carrier_loop_filter_bandwidth = 22.0Hz,
        code_loop_filter_bandwidth = 1.5Hz,
    )
    sat1 = TrackedSat(gpsl1, 1, 0.5, 100.0Hz; doppler_estimator = estimator)
    track_state = TrackState(gpsl1, sat1; doppler_estimator = estimator)

    # The default `ConventionalAssistedPLLAndDLL` has another state type.
    bad_sat = TrackedSat(gpsl1, 2, 0.25, 200.0Hz)
    @test_throws ArgumentError merge_sats(track_state, bad_sat)
end

# An estimator that keeps what each record hands the carrier discriminators,
# around the loop it wraps.
struct PromptRecordingEstimator{E<:TrackingLoops.AbstractDopplerEstimator} <:
       TrackingLoops.AbstractDopplerEstimator
    inner::E
    records::Vector{Tuple{ComplexF64,Int8}}
end
TrackingLoops.init_estimator_state(
    e::PromptRecordingEstimator,
    driver_signal,
    carrier_doppler,
    code_doppler,
) = TrackingLoops.init_estimator_state(
    e.inner,
    driver_signal,
    carrier_doppler,
    code_doppler,
)
function TrackingLoops.step_loop(
    e::PromptRecordingEstimator,
    state,
    record::TrackingLoops.LoopRecord,
    words,
    landing_sample::Int64,
)
    push!(e.records, (record.previous_prompt, record.carrier_phase_polarity))
    TrackingLoops.step_loop(e.inner, state, record, words, landing_sample)
end

@testset "A record hands the FLL its previous prompt only at an unchanged polarity" begin
    sampling_frequency = 5e6Hz
    # One record per call; returns the records' (previous prompt, polarity)
    # and the prompts left on the signal.
    function step_records(signal, synced, num_blocks_per_record)
        estimator = PromptRecordingEstimator(
            ConventionalAssistedPLLAndDLL(),
            Tuple{ComplexF64,Int8}[],
        )
        track_state = TrackState(; signal, doppler_estimator = estimator)
        track_state = add_satellite!(track_state; prn = 1, carrier_doppler = 100.0Hz)
        sats = get_sat_states(track_state)
        sat = sats[1]
        sats[1] = TrackedSat(
            sat;
            signals = map(
                s -> TrackedSignal(s; bit_buffer = synced(get_bit_buffer(s))),
                sat.signals,
            ),
        )
        correlator = get_default_correlator(signal)
        correlator = update_accumulator(
            correlator,
            map(_ -> 1000.0 + 10im, TrackingLoops.get_accumulators(correlator)),
        )
        one_block =
            Int(sampling_frequency * get_code_length(signal) / get_chip_rate(signal))
        prompts = ComplexF64[]
        for num_blocks in num_blocks_per_record
            sat = get_sat_states(track_state)[1]
            samples = num_blocks * one_block
            append_correlator_output!(
                only(sat.signals),
                CorrelatorOutput(correlator, samples, samples),
            )
            estimate_dopplers_and_filter_prompt!(track_state, _meas_l1(sampling_frequency))
            push!(
                prompts,
                get_last_fully_integrated_filtered_prompt(get_sat_states(track_state)[1]),
            )
        end
        estimator.records, prompts
    end

    # A bit buffer synced at polarity -1, on the first block of a bit.
    synced(::TrackingLoops.BitBuffer{B}) where {B} = TrackingLoops.BitBuffer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(-1),
        complex(0.0, 0.0),
        0,
        Float32[],
        TrackingLoops.PhaseAccumulators(),
    )

    # A data signal's prompt sign stays unknown, so every record but the first
    # chains from the one before, across a change of record length too.
    records, prompts = step_records(GPSL1CA(), synced, (1, 1, 2, 2))
    @test records[1][1] == 0
    @test records[2][1] == prompts[1]
    @test records[3][1] == prompts[2]
    @test records[4][1] == prompts[3]
    @test all(r -> r[2] == 0, records)

    # A synced pilot knows its prompt's sign, the sync polarity. The FLL
    # compares no prompts across the switch from the unknown sign.
    records, prompts = step_records(GalileoE1C(), synced, (1, 1))
    @test records[1] == (0, -1)
    @test records[2] == (prompts[1], -1)
end

# An estimator that keeps each record's carrier phase polarity and prompt.
struct PolarityRecordingEstimator{E<:TrackingLoops.AbstractDopplerEstimator} <:
       TrackingLoops.AbstractDopplerEstimator
    inner::E
    records::Vector{Tuple{Int8,ComplexF64}}
end
TrackingLoops.init_estimator_state(
    e::PolarityRecordingEstimator,
    driver_signal,
    carrier_doppler,
    code_doppler,
) = TrackingLoops.init_estimator_state(
    e.inner,
    driver_signal,
    carrier_doppler,
    code_doppler,
)
function TrackingLoops.step_loop(
    e::PolarityRecordingEstimator,
    state,
    record::TrackingLoops.LoopRecord,
    words,
    landing_sample::Int64,
)
    push!(
        e.records,
        (
            record.carrier_phase_polarity,
            TrackingLoops.get_prompt(record.filtered_correlator),
        ),
    )
    TrackingLoops.step_loop(e.inner, state, record, words, landing_sample)
end

@testset "A synced pilot's carrier phase polarity is its prompt's sign" begin
    # Pre-sync the replica carries chip 0 on every block, so the prompt is
    # sign · chip(k) · chip(0); synced, it is the sign. E1C's CS25 starts with a
    # negative chip, so reading the pre-sync prompts as they come flips the
    # sign, a half-cycle PLL error.
    e1c = GalileoE1C()
    prn = 1
    chip(k) = GNSSSignals.secondary_value(get_secondary_code(e1c), prn, k)
    @test chip(0) == -1
    sampling_frequency = 16.368e6Hz
    num_samples = 65472 # 4 ms, one code period
    record(prompt, sample_index) = CorrelatorOutput(
        update_accumulator(
            get_default_correlator(e1c),
            ComplexF64(prompt * num_samples) .* SVector(0.2, 0.7, 1.0, 0.7, 0.2),
        ),
        num_samples,
        sample_index,
    )
    rng = Random.Xoshiro(1)
    for sign in (1, -1)
        estimator = PolarityRecordingEstimator(
            ConventionalAssistedPLLAndDLL(),
            Tuple{Int8,ComplexF64}[],
        )
        track_state = TrackState(; signal = e1c, doppler_estimator = estimator)
        track_state = add_satellite!(track_state; prn, carrier_doppler = 0.0Hz)
        noisy(prompt) = prompt + 0.05 * randn(rng, ComplexF64)
        step!(prompt, k) = begin
            append_correlator_output!(get_sat_state(track_state, prn), record(prompt, k))
            with_logger(NullLogger()) do
                estimate_dopplers_and_filter_prompt!(track_state, _meas_l1(sampling_frequency))
            end
        end
        k = 0
        while !has_bit_or_secondary_code_been_found(get_sat_state(track_state, prn)) &&
              k < 200
            step!(noisy(sign * chip(k) * chip(0)), (k += 1) * num_samples)
        end
        @test has_bit_or_secondary_code_been_found(get_sat_state(track_state, prn))
        num_pre_sync = length(estimator.records)
        @test all(r -> r[1] == 0, estimator.records)
        for _ = 1:10
            step!(noisy(sign), (k += 1) * num_samples)
        end
        post_sync = estimator.records[(num_pre_sync+1):end]
        @test all(r -> r[1] == sign, post_sync)
        @test all(r -> Base.sign(real(r[2])) == sign, post_sync)
    end
end

end
