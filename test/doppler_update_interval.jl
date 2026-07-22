module UpdateIntervalTest

using Test: @test, @testset, @test_throws
using Unitful: Hz, s, uconvert, NoUnits
import Tracking
using GNSSSignals:
    GPSL1CA,
    GalileoE1B,
    gen_code,
    get_code_frequency,
    get_code_length,
    get_code_center_frequency_ratio
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
    get_carrier_doppler,
    get_code_doppler,
    CorrelatorOutput

@testset "smallest-code-period default resolution" begin
    # `nothing` doppler_update_interval => smallest primary-code period across all
    # signals in all groups. GPS L1 C/A: 1023 chips / 1.023 MHz = 1 ms.
    ts = TrackState(; signal = GPSL1CA())
    dt = Tracking._smallest_code_period(ts)
    @test uconvert(NoUnits, dt * get_code_frequency(GPSL1CA())) ≈ get_code_length(GPSL1CA())
    @test dt ≈ 1e-3 * s

    # Mixed L1 C/A (1 ms) + Galileo E1B (4 ms): the smaller wins.
    ts_mixed = TrackState(; signals = (l1ca = (GPSL1CA(),), e1b = (GalileoE1B(),)))
    @test Tracking._smallest_code_period(ts_mixed) ≈ 1e-3 * s
end

@testset "chunk-grid arithmetic (`_chunk_last_sample`)" begin
    # Boundaries lie on a shared time grid, re-anchored to the absolute chunk
    # index so different bands stay time-aligned without cumulative drift.
    @test Tracking._chunk_last_sample(1e-3 * s, 0, 5e6Hz, 100_000) == 5000
    @test Tracking._chunk_last_sample(1e-3 * s, 1, 5e6Hz, 100_000) == 10_000
    @test Tracking._chunk_last_sample(1e-3 * s, 9, 5e6Hz, 100_000) == 50_000
    # A second band at a different sampling rate, same 1 ms time grid.
    @test Tracking._chunk_last_sample(1e-3 * s, 0, 13e6Hz, 200_000) == 13_000
    @test Tracking._chunk_last_sample(1e-3 * s, 4, 13e6Hz, 200_000) == 65_000
    # Clamps to the buffer end on the final chunk.
    @test Tracking._chunk_last_sample(1e-3 * s, 20, 5e6Hz, 12_345) == 12_345
    # `nothing` => no chunking, whole buffer at once.
    @test Tracking._chunk_last_sample(nothing, 0, 5e6Hz, 12_345) == 12_345
end

@testset "doppler_update_interval shorter than a sample throws" begin
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])
    signal = zeros(ComplexF32, 4000)
    @test_throws ArgumentError track(signal, ts, 4e6Hz; doppler_update_interval = 1e-9 * s)
end

@testset "doppler_update_interval without time units throws a clear error" begin
    gpsl1 = GPSL1CA()
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, 1, 0.0, 100.0Hz)])
    signal = zeros(ComplexF32, 4000)
    @test_throws ArgumentError track(signal, ts, 4e6Hz; doppler_update_interval = 1e-3)
    err = try
        track(signal, ts, 4e6Hz; doppler_update_interval = 1e-3)
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
    signal = ComplexF32.(
        gen_code(num_samples, gpsl1, prn, fs, get_code_frequency(gpsl1), start_code_phase),
    )
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, 0.0Hz)])
    dc = CPUDownconvertAndCorrelator()
    measurements = (L1 = BandMeasurement(signal, fs),)

    # One chunk spanning all three periods, NCO held fixed, estimator not run,
    # so the collected outputs survive for inspection.
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
        # Each integration ends on its code-period boundary; the recorded end
        # sample index and integrated-sample count reflect that.
        @test isapprox(o.sample_index, i * period; atol = 2)
        @test isapprox(o.integrated_samples, period; atol = 2)
    end
end

@testset "estimator empties correlator_outputs after each chunk" begin
    gpsl1 = GPSL1CA()
    fs = 5e6Hz
    prn = 1
    num_samples = 15_000
    signal =
        ComplexF32.(gen_code(num_samples, gpsl1, prn, fs, get_code_frequency(gpsl1), 0.0))
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, 0.0, 0.0Hz)])
    track!(signal, ts, fs)
    @test isempty(get_correlator_outputs(only(get_sat_state(ts, prn).signals)))
end

@testset "enlarged doppler_update_interval still converges" begin
    # A long continuous buffer with a constant true Doppler; the initial
    # estimate is offset. Both the default (per-code-period) and a 4 ms
    # (4-periods-per-chunk) update interval must pull the carrier Doppler in.
    gpsl1 = GPSL1CA()
    fs = 4e6Hz
    prn = 1
    true_doppler = 300.0Hz
    start_code_phase = 0.0
    code_frequency =
        true_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    num_samples = 400_000                   # 100 ms
    range = 0:(num_samples-1)
    signal = ComplexF32.(
        cis.(2π .* true_doppler .* range ./ fs) .*
        gen_code(num_samples, gpsl1, prn, fs, code_frequency, start_code_phase),
    )

    converged(doppler_update_interval) = begin
        ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, true_doppler - 30Hz)])
        ts = track(signal, ts, fs; doppler_update_interval)
        abs(get_carrier_doppler(ts) / Hz - true_doppler / Hz)
    end

    @test converged(nothing) < 5      # default = 1 ms
    @test converged(4e-3 * s) < 5     # 4 ms => 4 correlator outputs per chunk
end

end
