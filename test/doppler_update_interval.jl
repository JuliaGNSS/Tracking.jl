module UpdateIntervalTest

using Test: @test, @testset, @test_throws
using Unitful: Hz, s, uconvert, NoUnits
import Tracking
using GNSSSignals:
    GPSL1CA,
    GalileoE1B,
    gen_code,
    get_chip_rate,
    get_code_length,
    get_chip_rate_to_center_frequency_ratio
using Tracking:
    TrackedSat,
    TrackState,
    BandMeasurement,
    track,
    track!,
    CPUDownconvertAndCorrelator,
    downconvert_and_correlate!,
    get_correlator_outputs,
    get_sat_state,
    get_carrier_doppler
using TrackingLoops: CorrelatorOutput

@testset "smallest-code-period default resolution" begin
    # `nothing`: the shortest primary-code period of all signals.
    ts = TrackState(; signal = GPSL1CA())
    dt = Tracking._smallest_code_period(ts)
    @test uconvert(NoUnits, dt * get_chip_rate(GPSL1CA())) ≈ get_code_length(GPSL1CA())
    @test dt ≈ 1e-3 * s

    ts_mixed = TrackState(; signals = (l1ca = (GPSL1CA(),), e1b = (GalileoE1B(),)))
    @test Tracking._smallest_code_period(ts_mixed) ≈ 1e-3 * s
end

@testset "chunk-grid arithmetic (`_chunk_last_sample`)" begin
    # Boundaries are computed from the absolute chunk index, so bands stay
    # aligned without drift.
    @test Tracking._chunk_last_sample(1e-3 * s, 0, 5e6Hz, 100_000) == 5000
    @test Tracking._chunk_last_sample(1e-3 * s, 1, 5e6Hz, 100_000) == 10_000
    @test Tracking._chunk_last_sample(1e-3 * s, 9, 5e6Hz, 100_000) == 50_000
    # A second band at a different sampling rate, same 1 ms time grid.
    @test Tracking._chunk_last_sample(1e-3 * s, 0, 13e6Hz, 200_000) == 13_000
    @test Tracking._chunk_last_sample(1e-3 * s, 4, 13e6Hz, 200_000) == 65_000
    @test Tracking._chunk_last_sample(1e-3 * s, 20, 5e6Hz, 12_345) == 12_345
    @test Tracking._chunk_last_sample(nothing, 0, 5e6Hz, 12_345) == 12_345
end

@testset "doppler_update_interval shorter than a sample throws" begin
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])
    samples = zeros(ComplexF32, 4000)
    @test_throws ArgumentError track(samples, ts, 4e6Hz; doppler_update_interval = 1e-9 * s)
end

@testset "doppler_update_interval without time units throws a clear error" begin
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])
    samples = zeros(ComplexF32, 4000)
    @test_throws ArgumentError track(samples, ts, 4e6Hz; doppler_update_interval = 1e-3)
    err = try
        track(samples, ts, 4e6Hz; doppler_update_interval = 1e-3)
    catch e
        e
    end
    @test occursin("time quantity", err.msg)
end

@testset "a chunk collects every completed correlator output, sample-indexed" begin
    gpsl1 = GPSL1CA()
    fs = 5e6Hz
    prn = 1
    start_code_phase = 0.0
    n_periods = 3
    period = 5000                       # 1 ms at 5 MHz = one L1 C/A period
    num_samples = n_periods * period
    samples = ComplexF32.(
        gen_code(num_samples, gpsl1, prn, fs, get_chip_rate(gpsl1), start_code_phase),
    )
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, 0.0Hz)])
    dc = CPUDownconvertAndCorrelator()
    measurements = (L1 = BandMeasurement(samples, fs),)

    # Correlate only, so the outputs survive.
    downconvert_and_correlate!(
        dc,
        measurements,
        ts;
        chunk_index = 0,
        chunk_duration = n_periods * 1e-3 * s,
    )
    outputs = get_correlator_outputs(only(get_sat_state(ts, prn).signals))

    @test length(outputs) == n_periods
    @test all(o -> o isa CorrelatorOutput, outputs)
    for (i, o) in enumerate(outputs)
        @test isapprox(o.sample_index, i * period; atol = 2)
        @test isapprox(o.integrated_samples, period; atol = 2)
    end
end

@testset "estimator empties correlator_outputs after each chunk" begin
    gpsl1 = GPSL1CA()
    fs = 5e6Hz
    prn = 1
    num_samples = 15_000
    samples = ComplexF32.(gen_code(num_samples, gpsl1, prn, fs, get_chip_rate(gpsl1), 0.0))
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, 0.0, 0.0Hz)])
    track!(samples, ts, fs)
    @test isempty(get_correlator_outputs(only(get_sat_state(ts, prn).signals)))
end

@testset "enlarged doppler_update_interval still converges" begin
    # Both the per-period default and a 4 ms interval pull the Doppler in.
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    prn = 1
    true_doppler = 300.0Hz
    start_code_phase = 0.0
    code_frequency =
        true_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1) + get_chip_rate(gpsl1)
    num_samples = 400_000                   # 100 ms
    range = 0:(num_samples-1)
    samples = ComplexF32.(
        cis.(2π .* true_doppler .* range ./ fs) .*
        gen_code(num_samples, gpsl1, prn, fs, code_frequency, start_code_phase),
    )

    converged(doppler_update_interval) = begin
        ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, true_doppler - 30Hz)])
        ts = track(samples, ts, fs; doppler_update_interval)
        abs(get_carrier_doppler(ts) / Hz - true_doppler / Hz)
    end

    @test converged(nothing) < 5      # default = 1 ms
    @test converged(4e-3 * s) < 5     # 4 ms => 4 correlator outputs per chunk
end

end
