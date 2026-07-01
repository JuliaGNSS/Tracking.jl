module DownconvertAndCorrelateOneBitTest

using Test: @test, @testset, @test_throws
using Unitful: Hz
import GNSSSignals
using GNSSSignals:
    GPSL1CA, GPSL5I, GalileoE1B, gen_code, get_code_center_frequency_ratio, get_code_frequency
using Tracking:
    TrackedSat,
    TrackState,
    track,
    downconvert_and_correlate,
    BandMeasurement,
    get_sat_state,
    get_carrier_doppler,
    get_prompt,
    get_early,
    get_late,
    get_accumulators,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    NumAnts,
    TrackedSignal,
    DefaultPostCorrFilter,
    ConventionalAssistedPLLAndDLL,
    CPUThreadedDownconvertAndCorrelator,
    OneBitDownconvertAndCorrelator,
    OneBitThreadedDownconvertAndCorrelator

# 12-bit-ADC-style Complex{Int16} capture (carrier × unit-normalised code, scaled to `peak`).
function make_capture(sig, prn, fs, nsamp, cdopp, cphase; peak = 2000)
    fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
    code = gen_code(nsamp, sig, prn, fs, fc, cphase)
    code = code ./ maximum(abs, code)
    s = cis.(2π .* cdopp .* (0:(nsamp-1)) ./ fs .+ 0.6) .* code
    complex.(round.(Int16, real.(s) .* peak), round.(Int16, imag.(s) .* peak))
end

band_key_for(sig) = sig isa GPSL5I ? :l5 : :l1
make_capture_mat(sig, fs, nsamp, cdopp, cphase, M; peak = 2000) =
    repeat(make_capture(sig, 1, fs, nsamp, cdopp, cphase; peak); outer = (1, M))

function correlate_once(dc, sig, fs, nsamp, cdopp, cphase; correlator = nothing, mat = false, M = 1)
    cap =
        mat ? make_capture_mat(sig, fs, nsamp, cdopp, cphase, M) :
        make_capture(sig, 1, fs, nsamp, cdopp, cphase)
    meas = NamedTuple{(band_key_for(sig),)}((BandMeasurement(cap, fs, 0.0Hz),))
    sat =
        correlator === nothing ? TrackedSat(sig, 1, cphase, cdopp) :
        TrackedSat(sig, 1, cphase, cdopp; correlator)
    ts = TrackState(sig, [sat])
    ts2 = downconvert_and_correlate(dc, meas, ts)
    first(get_sat_state(ts2, 1).signals).correlator
end

@testset "One-bit downconvert and correlate" begin
    # The code is ±1 in both the float and 1-bit pipelines, so the E/P·L/P magnitude
    # ratios (the correlation-triangle shape) survive 1-bit quantisation at high SNR;
    # only the amplitude and (square-wave) carrier phase change. Compare ratios to the
    # Float32 backend; the prompt phase gets a loose bound (the 1-bit carrier adds a
    # systematic phase bias the PLL absorbs — see the convergence test).
    @testset "ratios match Float32: $(nameof(typeof(sig))) @ $(fs/1e6Hz) MHz" for (sig, fs) in
                ((GPSL1CA(), 5e6Hz), (GPSL1CA(), 2e6Hz), (GPSL5I(), 12e6Hz))
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cf = correlate_once(CPUThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        cb = correlate_once(OneBitThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        @test abs(get_early(cb)) / abs(get_prompt(cb)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 3e-2
        @test abs(get_late(cb)) / abs(get_prompt(cb)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 3e-2
        # prompt is the correlation peak; early ≈ late (symmetric)
        @test abs(get_prompt(cb)) > abs(get_early(cb))
        @test abs(get_prompt(cb)) > abs(get_late(cb))
        @test abs(get_early(cb)) / abs(get_late(cb)) ≈ 1 atol = 5e-2
        # prompt phase within the 1-bit carrier bias
        @test abs(mod2pi(angle(get_prompt(cb)) - angle(get_prompt(cf)) + π) - π) < 0.6
    end

    @testset "VeryEarlyPromptLateCorrelator (NC=5)" begin
        sig, fs = GPSL5I(), 12e6Hz          # BPSK (Int8 code); the 1-bit backend is BPSK-only
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        c = correlate_once(
            OneBitThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0;
            correlator = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
        )
        a = get_accumulators(c)
        @test length(a) == 5
        # tap 3 is prompt (the peak); it dominates the outer very-early/very-late taps
        @test abs(a[3]) > abs(a[1]) && abs(a[3]) > abs(a[5])
        @test abs(a[3]) ≥ abs(a[2]) && abs(a[3]) ≥ abs(a[4])
    end

    @testset "multiple antennas (M=$M)" for M in (2, 4)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cm = correlate_once(OneBitThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(M)), mat = true, M)
        c1 = correlate_once(OneBitThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        pm = get_prompt(cm)
        @test length(pm) == M
        # every antenna sees the identical capture → identical to the M=1 result
        for j = 1:M
            @test get_prompt(cm)[j] == get_prompt(c1)
            @test get_early(cm)[j] == get_early(c1)
            @test get_late(cm)[j] == get_late(c1)
        end
    end

    @testset "single-threaded and threaded backends agree" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        c1 = correlate_once(OneBitDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        ct = correlate_once(OneBitThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        @test get_prompt(c1) == get_prompt(ct)
        @test get_early(c1) == get_early(ct)
        @test get_late(c1) == get_late(ct)
    end

    @testset "multi-signal-per-sat (N=$N) matches single signal" for N in (2, 3)
        # A sat carrying N GPS L1 signals; the 1-bit backend correlates each in turn,
        # so each signal's correlator is identical to correlating it alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = OneBitThreadedDownconvertAndCorrelator()
        est = ConventionalAssistedPLLAndDLL()
        mksig() = TrackedSignal(
            sig;
            num_ants = NumAnts(1),
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            post_corr_filter = DefaultPostCorrFilter(),
        )
        cs = first(
            get_sat_state(
                downconvert_and_correlate(
                    dc, meas,
                    TrackState(
                        sig,
                        TrackedSat((mksig(),), 1, 100.0, 200Hz; doppler_estimator = est);
                        doppler_estimator = est,
                    ),
                ),
                1,
            ).signals,
        ).correlator
        satN = TrackedSat(ntuple(_ -> mksig(), N), 1, 100.0, 200Hz; doppler_estimator = est)
        tsN = downconvert_and_correlate(dc, meas, TrackState(sig, satN; doppler_estimator = est))
        for s in get_sat_state(tsN, 1).signals
            @test get_prompt(s.correlator) == get_prompt(cs)
            @test get_early(s.correlator) == get_early(cs)
            @test get_late(s.correlator) == get_late(cs)
        end
    end

    @testset "errors on non-Complex{Int16} (Float) measurement" begin
        sig, fs = GPSL1CA(), 5e6Hz
        capf = Complex{Float32}.(make_capture(sig, 1, fs, 5000, 200Hz, 100.0))
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        meas = (l1 = BandMeasurement(capf, fs, 0.0Hz),)
        @test_throws ArgumentError downconvert_and_correlate(
            OneBitThreadedDownconvertAndCorrelator(), meas, ts,
        )
    end

    @testset "full track converges (GPS L1 C/A)" begin
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase = 300Hz, 230.0
        nsamp = round(Int, (fs / 1Hz) * 1e-3) * 5
        cap = make_capture(sig, 1, fs, nsamp, cdopp, cphase)
        ts = TrackState(sig, [TrackedSat(sig, 1, cphase, cdopp - 20Hz)])
        ts = track(cap, ts, fs; downconvert_and_correlator = OneBitThreadedDownconvertAndCorrelator())
        # 1-bit tracking is noisier than the integer/float paths; still pulls in from
        # the 20 Hz offset toward the true Doppler.
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 15Hz
    end
end

end
