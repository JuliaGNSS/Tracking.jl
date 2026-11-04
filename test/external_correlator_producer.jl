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
    step_tracking_loops!,
    step_tracking_loops,
    get_correlator_bank_outputs,
    append_correlator_bank_output!,
    get_sat_state,
    get_signals,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_beamformed_prompt
import TrackingLoops
using TrackingLoops:
    BitSynchronizer,
    is_bit_or_secondary_code_synced,
    get_default_correlator_bank,
    NumAntennas,
    EarlyPromptLateCorrelatorBank,
    CorrelatorBankOutput,
    append_noise_observation!,
    calc_noise_observation_from_samples

# One L1 C/A satellite, its `correlator_bank_outputs` filled by the software
# correlate step. Returns the state and the measurements.
function correlated_state(; fs = 4e6Hz, prn = 1, num_samples = 12_000)
    gpsl1 = GPSL1CA()
    code_frequency = get_chip_rate(gpsl1)
    samples = ComplexF32.(gen_code(num_samples, gpsl1, prn, fs, code_frequency, 0.0))
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, 0.0, 100.0Hz)])
    dc = CPUDownconvertAndCorrelator()
    measurements = (L1 = BandMeasurement(samples, fs),)
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
    @test length(
        get_correlator_bank_outputs(only(get_signals(get_sat_state(ts_meas, 1)))),
    ) == length(get_correlator_bank_outputs(only(get_signals(get_sat_state(ts_fs, 1)))))

    step_tracking_loops!(ts_meas, measurements)          # BandMeasurements
    step_tracking_loops!(ts_fs, (L1 = fs,))              # per-band NamedTuple

    @test get_carrier_doppler(ts_fs, 1) == get_carrier_doppler(ts_meas, 1)
    @test get_code_doppler(ts_fs, 1) == get_code_doppler(ts_meas, 1)
    @test get_last_fully_integrated_beamformed_prompt(ts_fs, 1) ==
          get_last_fully_integrated_beamformed_prompt(ts_meas, 1)
    @test isempty(get_correlator_bank_outputs(only(get_signals(get_sat_state(ts_fs, 1)))))
    @test isempty(get_correlator_bank_outputs(only(get_signals(get_sat_state(ts_meas, 1)))))
end

@testset "Dict fs source matches the NamedTuple fs source" begin
    fs = 4e6Hz
    ts_nt, _ = correlated_state(; fs)
    ts_dict, _ = correlated_state(; fs)
    step_tracking_loops!(ts_nt, (L1 = fs,))
    step_tracking_loops!(ts_dict, Dict(:L1 => fs))
    @test get_carrier_doppler(ts_dict, 1) == get_carrier_doppler(ts_nt, 1)
    @test get_code_doppler(ts_dict, 1) == get_code_doppler(ts_nt, 1)
end

@testset "immutable fs form matches the in-place fs form" begin
    fs = 4e6Hz
    ts_inplace, _ = correlated_state(; fs)
    ts_src, _ = correlated_state(; fs)
    ts_imm = step_tracking_loops(ts_src, (L1 = fs,))
    step_tracking_loops!(ts_inplace, (L1 = fs,))
    @test get_carrier_doppler(ts_imm, 1) == get_carrier_doppler(ts_inplace, 1)
    @test get_code_doppler(ts_imm, 1) == get_code_doppler(ts_inplace, 1)
end

@testset "append_correlator_bank_output! at every addressing level" begin
    corr = get_default_correlator_bank(GPSL1CA(), NumAntennas(1))
    shift = corr.preferred_early_late_to_prompt_spacing
    raw = EarlyPromptLateCorrelatorBank(
        SVector(complex(2000.0), complex(4000.0), complex(2000.0)),
        shift,
    )
    out = CorrelatorBankOutput(raw, 4000, 3999)

    ts = TrackState(; signal = GPSL1CA())
    ts = Tracking.add_satellite!(ts; prn = 1, code_phase = 0.0, carrier_doppler = 0.0Hz)

    append_correlator_bank_output!(ts, out)
    @test length(get_correlator_bank_outputs(ts)) == 1
    append_correlator_bank_output!(ts, out, 1)
    @test length(get_correlator_bank_outputs(ts, 1)) == 2
    append_correlator_bank_output!(ts, out, :default, 1, 1)
    @test length(get_correlator_bank_outputs(ts, :default, 1, 1)) == 3
    sat = get_sat_state(ts, 1)
    append_correlator_bank_output!(sat, out)
    @test length(get_correlator_bank_outputs(sat)) == 4
    sig = only(get_signals(sat))
    @test append_correlator_bank_output!(sig, out) === sig
    @test length(get_correlator_bank_outputs(sig)) == 5

    step_tracking_loops!(ts, (L1 = 4e6Hz,))
    @test isempty(get_correlator_bank_outputs(ts))
end

@testset "external-producer standalone estimate (no downconvert_and_correlate!)" begin
    # As an FPGA would: hand-built outputs in `sample_index` order, stepped with
    # only a per-band rate.
    fs = 4e6Hz
    corr = get_default_correlator_bank(GPSL1CA(), NumAntennas(1))
    shift = corr.preferred_early_late_to_prompt_spacing
    ts = TrackState(; signal = GPSL1CA())
    ts = Tracking.add_satellite!(ts; prn = 3, code_phase = 0.0, carrier_doppler = 50.0Hz)

    period = 4000                                   # 1 ms at 4 MHz
    for k = 1:3
        raw = EarlyPromptLateCorrelatorBank(
            SVector(complex(1500.0), complex(3000.0), complex(1500.0)),
            shift,
        )
        # sample_index is chunk-relative (end sample of the k-th integration).
        append_correlator_bank_output!(
            ts,
            CorrelatorBankOutput(raw, period, k * period - 1),
            3,
        )
    end
    @test length(get_correlator_bank_outputs(ts, 3)) == 3
    ts = step_tracking_loops!(ts, (L1 = fs,))
    @test isempty(get_correlator_bank_outputs(ts, 3))
    @test isfinite(get_carrier_doppler(ts, 3) / Hz)
    @test isfinite(get_code_doppler(ts, 3) / Hz)
end

# One producer and one buffer of hand-built outputs, stepped once.
function _ingest_and_step!(ts; fs = 4e6Hz, prn = 3, num_integrations = 3, period = 4000)
    corr = get_default_correlator_bank(GPSL1CA(), NumAntennas(1))
    shift = corr.preferred_early_late_to_prompt_spacing
    for k = 1:num_integrations
        raw = EarlyPromptLateCorrelatorBank(
            SVector(complex(1500.0), complex(3000.0), complex(1500.0)),
            shift,
        )
        append_correlator_bank_output!(
            ts,
            CorrelatorBankOutput(raw, period, k * period - 1),
            prn,
        )
    end
    step_tracking_loops!(ts, (L1 = fs,))
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
    _ingest_and_step!(starved; fs)
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
        calc_noise_observation_from_samples(n * 1e-6 * ustrip(Hz, fs), n, fs),
        :GPSL1CA,
    )
    _ingest_and_step!(fed; fs)
    @test Base.length(Tracking.get_cn0_estimator(fed, 3)) == 3
    @test ustrip(uconvert(dBHz, TrackingLoops.estimate_cn0(fed, 3))) ≈
          10log10((3000 / 4000)^2 / 1e-6 - 1 / 1e-3) atol = 1e-6
end

@testset "a producer's integration straddling a bit boundary warns and drops sync" begin
    # The bit buffer reports the overshoot and this package logs it: 18 code periods
    # into a bit, a three-code-period integration from a producer off the bit grid.
    fs = 4e6Hz
    ts = Tracking.add_satellite!(
        TrackState(; signal = GPSL1CA());
        prn = 1,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    satellites = first(ts.groups).satellites.values
    sat = satellites[1]
    synced = BitSynchronizer(UInt64(0xff), 8, true, complex(9.0, 1.0), 18, Float32[1.0])
    satellites[1] = TrackedSat(
        sat;
        signals = (TrackedSignal(only(get_signals(sat)); bit_synchronizer = synced),),
    )
    shift =
        get_default_correlator_bank(GPSL1CA(), NumAntennas(1)).preferred_early_late_to_prompt_spacing
    prompt = complex(12000.0)
    raw = EarlyPromptLateCorrelatorBank(SVector(prompt / 6, prompt, prompt / 6), shift)
    append_correlator_bank_output!(ts, CorrelatorBankOutput(raw, 12_000, 11_999))
    @test_logs (:warn, r"past the navigation-bit boundary") match_mode = :any step_tracking_loops!(
        ts,
        (L1 = fs,),
    )
    signal = only(get_signals(first(ts.groups).satellites.values[1]))
    @test !is_bit_or_secondary_code_synced(signal)
end

end
