module ExternalCorrelatorProducerTest

using Test: @test, @testset
using Unitful: Hz, s
using StaticArrays: SVector
import Tracking
using GNSSSignals: GPSL1CA, gen_code, get_code_frequency
using Tracking:
    TrackedSat,
    TrackState,
    BandMeasurement,
    CPUDownconvertAndCorrelator,
    downconvert_and_correlate!,
    estimate_dopplers_and_filter_prompt!,
    estimate_dopplers_and_filter_prompt,
    get_correlator_outputs,
    append_correlator_output!,
    get_sat_state,
    get_signals,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_filtered_prompt,
    get_default_correlator,
    NumAnts,
    EarlyPromptLateCorrelator,
    CorrelatorOutput

# A track state with one L1 C/A sat, an offset initial Doppler, and its
# `correlator_outputs` buffer filled by the software correlate phase (so the
# outputs are realistic). Returns the state plus the measurements used.
function correlated_state(; fs = 4e6Hz, prn = 1, num_samples = 12_000)
    gpsl1 = GPSL1CA()
    code_frequency = get_code_frequency(gpsl1)
    signal = ComplexF32.(gen_code(num_samples, gpsl1, prn, fs, code_frequency, 0.0))
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, 0.0, 100.0Hz)])
    dc = CPUDownconvertAndCorrelator()
    measurements = (L1 = BandMeasurement(signal, fs),)
    # One big chunk spanning the whole buffer, estimator NOT run, so the
    # collected outputs survive for the estimator calls below.
    downconvert_and_correlate!(
        dc,
        measurements,
        ts;
        chunk_index = 0,
        chunk_duration = (num_samples / (fs / Hz)) * s,
    )
    return ts, measurements
end

@testset "fs source matches the BandMeasurements path exactly (in-place)" begin
    fs = 4e6Hz
    ts_meas, measurements = correlated_state(; fs)
    ts_fs, _ = correlated_state(; fs)
    # Same inputs went into both, so the collected outputs are identical.
    @test length(get_correlator_outputs(only(get_signals(get_sat_state(ts_meas, 1))))) ==
          length(get_correlator_outputs(only(get_signals(get_sat_state(ts_fs, 1)))))

    estimate_dopplers_and_filter_prompt!(ts_meas, measurements)          # BandMeasurements
    estimate_dopplers_and_filter_prompt!(ts_fs, (L1 = fs,))              # per-band NamedTuple

    @test get_carrier_doppler(ts_fs, 1) == get_carrier_doppler(ts_meas, 1)
    @test get_code_doppler(ts_fs, 1) == get_code_doppler(ts_meas, 1)
    @test get_last_fully_integrated_filtered_prompt(ts_fs, 1) ==
          get_last_fully_integrated_filtered_prompt(ts_meas, 1)
    # Both consumed and cleared the buffer.
    @test isempty(get_correlator_outputs(only(get_signals(get_sat_state(ts_fs, 1)))))
    @test isempty(get_correlator_outputs(only(get_signals(get_sat_state(ts_meas, 1)))))
end

@testset "Dict fs source matches the NamedTuple fs source" begin
    fs = 4e6Hz
    ts_nt, _ = correlated_state(; fs)
    ts_dict, _ = correlated_state(; fs)
    estimate_dopplers_and_filter_prompt!(ts_nt, (L1 = fs,))
    estimate_dopplers_and_filter_prompt!(ts_dict, Dict(:L1 => fs))
    @test get_carrier_doppler(ts_dict, 1) == get_carrier_doppler(ts_nt, 1)
    @test get_code_doppler(ts_dict, 1) == get_code_doppler(ts_nt, 1)
end

@testset "immutable fs form matches the in-place fs form" begin
    fs = 4e6Hz
    ts_inplace, _ = correlated_state(; fs)
    ts_src, _ = correlated_state(; fs)
    ts_imm = estimate_dopplers_and_filter_prompt(ts_src, (L1 = fs,))
    estimate_dopplers_and_filter_prompt!(ts_inplace, (L1 = fs,))
    @test get_carrier_doppler(ts_imm, 1) == get_carrier_doppler(ts_inplace, 1)
    @test get_code_doppler(ts_imm, 1) == get_code_doppler(ts_inplace, 1)
end

@testset "append_correlator_output! at every addressing level" begin
    corr = get_default_correlator(GPSL1CA(), NumAnts(1))
    shift = corr.preferred_early_late_to_prompt_code_shift
    raw = EarlyPromptLateCorrelator(
        SVector(complex(2000.0), complex(4000.0), complex(2000.0)),
        shift,
    )
    out = CorrelatorOutput(raw, 4000, 3999)

    ts = TrackState(; signal = GPSL1CA())
    ts = Tracking.add_satellite!(ts; prn = 1, code_phase = 0.0, carrier_doppler = 0.0Hz)

    # TrackState-level, single group/sat/signal.
    append_correlator_output!(ts, out)
    @test length(get_correlator_outputs(ts)) == 1
    # TrackState-level with a prn selector.
    append_correlator_output!(ts, out, 1)
    @test length(get_correlator_outputs(ts, 1)) == 2
    # TrackState-level per-signal form (group, prn, signal selector).
    append_correlator_output!(ts, out, :default, 1, 1)
    @test length(get_correlator_outputs(ts, :default, 1, 1)) == 3
    # TrackedSat-level.
    sat = get_sat_state(ts, 1)
    append_correlator_output!(sat, out)
    @test length(get_correlator_outputs(sat)) == 4
    # TrackedSignal-level (returns the signal).
    sig = only(get_signals(sat))
    @test append_correlator_output!(sig, out) === sig
    @test length(get_correlator_outputs(sig)) == 5

    # The estimator consumes and clears whatever an external producer appended.
    estimate_dopplers_and_filter_prompt!(ts, (L1 = 4e6Hz,))
    @test isempty(get_correlator_outputs(ts))
end

@testset "external-producer standalone estimate (no downconvert_and_correlate!)" begin
    # Emulate an FPGA: build outputs by hand, append per signal in sample_index
    # order, then fold them through the estimator with only a per-band rate.
    fs = 4e6Hz
    corr = get_default_correlator(GPSL1CA(), NumAnts(1))
    shift = corr.preferred_early_late_to_prompt_code_shift
    ts = TrackState(; signal = GPSL1CA())
    ts = Tracking.add_satellite!(ts; prn = 3, code_phase = 0.0, carrier_doppler = 50.0Hz)

    period = 4000                                   # 1 ms at 4 MHz
    for k = 1:3
        raw = EarlyPromptLateCorrelator(
            SVector(complex(1500.0), complex(3000.0), complex(1500.0)),
            shift,
        )
        # sample_index is chunk-relative (end sample of the k-th integration).
        append_correlator_output!(ts, CorrelatorOutput(raw, period, k * period - 1), 3)
    end
    @test length(get_correlator_outputs(ts, 3)) == 3
    ts = estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
    @test isempty(get_correlator_outputs(ts, 3))
    @test isfinite(get_carrier_doppler(ts, 3) / Hz)
    @test isfinite(get_code_doppler(ts, 3) / Hz)
end

end
