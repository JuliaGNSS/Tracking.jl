module ExternalCorrelatorProducerTest

using Test: @test, @testset, @test_logs
using Unitful: Hz, s, dBHz, ustrip, uconvert
using StaticArrays: SVector
import Tracking
using GNSSSignals: GPSL1CA, gen_code, get_chip_rate
using Tracking:
    TrackedSat,
    TrackedSignal,
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
    get_last_fully_integrated_filtered_prompt
import TrackingLoops
using TrackingLoops:
    BitBuffer,
    has_bit_or_secondary_code_been_found,
    get_default_correlator,
    NumAnts,
    EarlyPromptLateCorrelator,
    CorrelatorOutput,
    append_noise_observation!,
    noise_observation_from_samples

# One L1 C/A satellite, its `correlator_outputs` filled by the software
# correlate step. Returns the state and the measurements.
function correlated_state(; fs = 4e6Hz, prn = 1, num_samples = 12_000)
    gpsl1 = GPSL1CA()
    code_frequency = get_chip_rate(gpsl1)
    signal = ComplexF32.(gen_code(num_samples, gpsl1, prn, fs, code_frequency, 0.0))
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, 0.0, 100.0Hz)])
    dc = CPUDownconvertAndCorrelator()
    measurements = (L1 = BandMeasurement(signal, fs),)
    # Correlate only, so the outputs survive for the estimate calls below.
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
    @test length(get_correlator_outputs(only(get_signals(get_sat_state(ts_meas, 1))))) ==
          length(get_correlator_outputs(only(get_signals(get_sat_state(ts_fs, 1)))))

    estimate_dopplers_and_filter_prompt!(ts_meas, measurements)          # BandMeasurements
    estimate_dopplers_and_filter_prompt!(ts_fs, (L1 = fs,))              # per-band NamedTuple

    @test get_carrier_doppler(ts_fs, 1) == get_carrier_doppler(ts_meas, 1)
    @test get_code_doppler(ts_fs, 1) == get_code_doppler(ts_meas, 1)
    @test get_last_fully_integrated_filtered_prompt(ts_fs, 1) ==
          get_last_fully_integrated_filtered_prompt(ts_meas, 1)
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

    append_correlator_output!(ts, out)
    @test length(get_correlator_outputs(ts)) == 1
    append_correlator_output!(ts, out, 1)
    @test length(get_correlator_outputs(ts, 1)) == 2
    append_correlator_output!(ts, out, :default, 1, 1)
    @test length(get_correlator_outputs(ts, :default, 1, 1)) == 3
    sat = get_sat_state(ts, 1)
    append_correlator_output!(sat, out)
    @test length(get_correlator_outputs(sat)) == 4
    sig = only(get_signals(sat))
    @test append_correlator_output!(sig, out) === sig
    @test length(get_correlator_outputs(sig)) == 5

    estimate_dopplers_and_filter_prompt!(ts, (L1 = 4e6Hz,))
    @test isempty(get_correlator_outputs(ts))
end

@testset "external-producer standalone estimate (no downconvert_and_correlate!)" begin
    # As an FPGA would: hand-built outputs in `sample_index` order, folded with
    # only a per-band rate.
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

# One producer and one buffer of hand-built outputs, folded once.
function _ingest_and_fold!(ts; fs = 4e6Hz, prn = 3, num_records = 3, period = 4000)
    corr = get_default_correlator(GPSL1CA(), NumAnts(1))
    shift = corr.preferred_early_late_to_prompt_code_shift
    for k = 1:num_records
        raw = EarlyPromptLateCorrelator(
            SVector(complex(1500.0), complex(3000.0), complex(1500.0)),
            shift,
        )
        append_correlator_output!(ts, CorrelatorOutput(raw, period, k * period - 1), prn)
    end
    estimate_dopplers_and_filter_prompt!(ts, (L1 = fs,))
end

@testset "the ingest path's C/N₀ needs a noise observation under the default" begin
    # The default C/N₀ estimator is noise-referenced, so a correlator-ingest
    # producer must also feed noise: without it the satellite reports
    # `-Inf dB-Hz` for good, which reads as loss of lock. Nothing throws and the
    # one warning is `maxlog = 1`, hence the assertion.
    fs = 4e6Hz
    starved = Tracking.add_satellite!(
        TrackState(; signal = GPSL1CA());
        prn = 3,
        code_phase = 0.0,
        carrier_doppler = 50.0Hz,
    )
    @test keys(starved.noise_estimators) == (:GPSL1CA,)
    _ingest_and_fold!(starved; fs)
    @test TrackingLoops.estimate_cn0(starved, 3) == -Inf * dBHz
    @test Base.length(Tracking.get_cn0_estimator(starved, 3)) == 0

    # With noise fed: `σ² = 1e-6·f_s` per sample gives `N₀ = 1e-6 Hz⁻¹`,
    # against a normalised prompt power of `(3000/4000)² = 0.5625`.
    fed = Tracking.add_satellite!(
        TrackState(; signal = GPSL1CA());
        prn = 3,
        code_phase = 0.0,
        carrier_doppler = 50.0Hz,
    )
    n = 4000
    append_noise_observation!(
        fed,
        noise_observation_from_samples(n * 1e-6 * ustrip(Hz, fs), n, fs),
        :GPSL1CA,
    )
    _ingest_and_fold!(fed; fs)
    @test Base.length(Tracking.get_cn0_estimator(fed, 3)) == 3
    @test ustrip(uconvert(dBHz, TrackingLoops.estimate_cn0(fed, 3))) ≈
          10log10((3000 / 4000)^2 / 1e-6 - 1 / 1e-3) atol = 1e-6
end

@testset "a producer record straddling a bit boundary warns and drops sync" begin
    # The bit buffer reports the overshoot and this package logs it: 18 blocks
    # into a bit, a three-block record from a producer off the bit grid.
    fs = 4e6Hz
    ts = Tracking.add_satellite!(
        TrackState(; signal = GPSL1CA());
        prn = 1,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    satellites = first(ts.groups).satellites.values
    sat = satellites[1]
    synced = BitBuffer(UInt64(0xff), 8, true, complex(9.0, 1.0), 18, Float32[1.0])
    satellites[1] = TrackedSat(
        sat;
        signals = (TrackedSignal(only(get_signals(sat)); bit_buffer = synced),),
    )
    shift =
        get_default_correlator(GPSL1CA(), NumAnts(1)).preferred_early_late_to_prompt_code_shift
    prompt = complex(12000.0)
    raw = EarlyPromptLateCorrelator(SVector(prompt / 6, prompt, prompt / 6), shift)
    append_correlator_output!(ts, CorrelatorOutput(raw, 12_000, 11_999))
    @test_logs (:warn, r"past the navigation-bit boundary") match_mode = :any estimate_dopplers_and_filter_prompt!(
        ts,
        (L1 = fs,),
    )
    signal = only(get_signals(first(ts.groups).satellites.values[1]))
    @test !has_bit_or_secondary_code_been_found(signal)
end

end
