module DownconvertAndCorrelateOneBitTest

using Test: @test, @testset, @test_throws
using Random: MersenneTwister
using Unitful: Hz, ustrip
import GNSSSignals
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    GalileoE1B,
    gen_code,
    get_code_center_frequency_ratio,
    get_code_frequency
using Tracking:
    TrackedSat,
    TrackState,
    track,
    downconvert_and_correlate,
    BandMeasurement,
    get_sat_state,
    get_carrier_doppler,
    estimate_cn0,
    get_prompt,
    get_early,
    get_late,
    get_accumulators,
    get_correlator_sample_shifts,
    AbstractCorrelator,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    NumAnts,
    TrackedSignal,
    DefaultPostCorrFilter,
    ConventionalAssistedPLLAndDLL,
    CPUThreadedDownconvertAndCorrelator,
    OneBitDownconvertAndCorrelator,
    OneBitThreadedDownconvertAndCorrelator
import Tracking
using StaticArrays: SVector

# Dynamic-tap-count correlator: sample shifts are a runtime Vector and the
# accumulators are a Vector (issue #126 (b) extension point), to exercise the
# one-bit backend's AbstractVector-shifts fallback. Parametric over the antenna
# count M so both the M=1 and M>1 fallback branches can be tested.
struct DynShiftsCorrelator{M} <: AbstractCorrelator{M}
    accumulators::Vector
    shifts::Vector{Int}
end
Tracking.get_accumulators(c::DynShiftsCorrelator) = c.accumulators
Tracking.update_accumulator(c::DynShiftsCorrelator{M}, acc) where {M} =
    DynShiftsCorrelator{M}(collect(acc), c.shifts)
Tracking.get_correlator_sample_shifts(
    c::DynShiftsCorrelator,
    sampling_frequency,
    code_frequency,
) = c.shifts

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

function correlate_once(
    dc,
    sig,
    fs,
    nsamp,
    cdopp,
    cphase;
    correlator = nothing,
    mat = false,
    M = 1,
)
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

# Track a noisy GPS L1CA capture at a known C/N0 with backend `dc`, returning the
# post-settle carrier-Doppler samples and the final C/N0 estimate. Carrier and code
# phase run continuously across 1 ms epochs; the capture is Complex{Int16} so the
# Float32 and one-bit backends see the identical bits. `amp`/`noise_std` follow the
# C/N0 recipe of the CN0-estimation tests (signal ×10^(C/N0/20), noise ×√fs).
function _track_noisy(
    dc,
    cn0_dbhz,
    seed;
    fs = 5e6Hz,
    cdopp = 300Hz,
    nepoch = 110,
    settle = 10,
)
    gpsl1 = GPSL1CA()
    nsamp = round(Int, (fs / 1Hz) * 1e-3)
    prn = 1
    cfreq = cdopp * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    range = 0:(nsamp-1)
    start_code_phase = 100.0
    start_carrier_phase = π / 2
    amp = 10^(cn0_dbhz / 20)
    noise_std = sqrt(fs / 1Hz)
    rng = MersenneTwister(seed)
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, cdopp)])
    dopplers = Float64[]
    for i = 0:(nepoch-1)
        carrier_phase = 2π * (cdopp / 1Hz) * nsamp * i / (fs / 1Hz) + start_carrier_phase
        code_phase = mod((cfreq / 1Hz) * nsamp * i / (fs / 1Hz) + start_code_phase, 1023)
        clean =
            cis.(2π .* (cdopp / 1Hz) .* range ./ (fs / 1Hz) .+ carrier_phase) .*
            gen_code(nsamp, gpsl1, prn, fs, cfreq, code_phase) .* amp
        noisy = clean .+ randn(rng, ComplexF64, nsamp) .* noise_std
        cap = complex.(
            round.(Int16, clamp.(real.(noisy), -32000, 32000)),
            round.(Int16, clamp.(imag.(noisy), -32000, 32000)),
        )
        ts = track(cap, ts, fs; downconvert_and_correlator = dc)
        i >= settle && push!(dopplers, get_carrier_doppler(ts) / Hz)
    end
    (; dopplers, cn0 = ustrip(estimate_cn0(ts)))
end

_mean(x) = sum(x) / length(x)
_std(x) = (m = _mean(x); sqrt(sum(v -> abs2(v - m), x) / (length(x) - 1)))

@testset "One-bit downconvert and correlate" begin
    # The code is ±1 in both the float and 1-bit pipelines, so the E/P·L/P magnitude
    # ratios (the correlation-triangle shape) survive 1-bit quantisation at high SNR;
    # only the amplitude and (square-wave) carrier phase change. Compare ratios to the
    # Float32 backend; the prompt phase gets a loose bound (the 1-bit carrier adds a
    # systematic phase bias the PLL absorbs — see the convergence test).
    @testset "ratios match Float32: $(nameof(typeof(sig))) @ $(fs/1e6Hz) MHz" for (
        sig,
        fs,
    ) in (
        (GPSL1CA(), 5e6Hz),
        (GPSL1CA(), 2e6Hz),
        (GPSL5I(), 12e6Hz),
    )
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cf = correlate_once(
            CPUThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        cb = correlate_once(
            OneBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
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
            OneBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
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
        cm = correlate_once(
            OneBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(M)),
            mat = true,
            M,
        )
        c1 = correlate_once(
            OneBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
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
        ct = correlate_once(
            OneBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
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
                    dc,
                    meas,
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
        tsN = downconvert_and_correlate(
            dc,
            meas,
            TrackState(sig, satN; doppler_estimator = est),
        )
        for s in get_sat_state(tsN, 1).signals
            @test get_prompt(s.correlator) == get_prompt(cs)
            @test get_early(s.correlator) == get_early(cs)
            @test get_late(s.correlator) == get_late(cs)
        end
    end

    @testset "multi-signal-per-sat × multi-antenna (N=$N, M=$M) matches single" for N in
                                                                                    (2, 3),
        M in (2, 4)
        # Exercises the tile-share kernel's M>1 path: N signals sharing one carrier +
        # measurement, each with M antennas. The shared carrier/measurement downconvert
        # is bit-identical, so every signal's per-antenna correlator (an SVector{M}) must
        # equal correlating that signal alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        dc = OneBitThreadedDownconvertAndCorrelator()
        cs = correlate_once(
            dc,
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(M)),
            mat = true,
            M,
        )
        cap = make_capture_mat(sig, fs, nsamp, 200Hz, 100.0, M)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        est = ConventionalAssistedPLLAndDLL()
        mksig() = TrackedSignal(
            sig;
            num_ants = NumAnts(M),
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(M)),
            post_corr_filter = DefaultPostCorrFilter(),
        )
        satN = TrackedSat(ntuple(_ -> mksig(), N), 1, 100.0, 200Hz; doppler_estimator = est)
        tsN = downconvert_and_correlate(
            dc,
            meas,
            TrackState(sig, satN; doppler_estimator = est),
        )
        for s in get_sat_state(tsN, 1).signals
            @test length(get_prompt(s.correlator)) == M
            @test get_prompt(s.correlator) == get_prompt(cs)
            @test get_early(s.correlator) == get_early(cs)
            @test get_late(s.correlator) == get_late(cs)
        end
    end

    @testset "dynamic Vector-shifts correlator matches static EPL" begin
        # A DynShiftsCorrelator with the same shifts as EPL must produce identical
        # accumulators through the one-bit backend's AbstractVector-shifts fallback as
        # the static @generated EPL kernel — the fallback sums each popcount chunk into
        # Int64 totals, which is bit-exact with the @generated kernel's Vec-then-sum.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        dc = OneBitThreadedDownconvertAndCorrelator()

        # M=1: dynamic fallback == static @generated kernel, exactly.
        cs = correlate_once(dc, sig, fs, nsamp, 200Hz, 100.0)
        cd = correlate_once(
            dc,
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = DynShiftsCorrelator{1}(zeros(ComplexF64, 3), shifts),
        )
        @test length(get_accumulators(cd)) == 3
        @test get_accumulators(cd) == get_accumulators(cs)

        # M=2: multi-antenna dynamic fallback matches the static EPL M=2 kernel.
        csm = correlate_once(
            dc,
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(2)),
            mat = true,
            M = 2,
        )
        cdm = correlate_once(
            dc,
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = DynShiftsCorrelator{2}(
                fill(zero(SVector{2,ComplexF64}), 3),
                shifts,
            ),
            mat = true,
            M = 2,
        )
        @test get_accumulators(cdm) == get_accumulators(csm)
    end

    @testset "dynamic Vector-shifts, band-shared (>1 sat) matches per-sat" begin
        # ≥2 sats with dynamic correlators trip the band-shared measurement path in the
        # AbstractVector fallback (`_ob_realign_meas!`); the shared pack is bit-identical,
        # so each sat must equal correlating that PRN alone (1 sat → direct pack).
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = OneBitThreadedDownconvertAndCorrelator()
        mkcorr() = DynShiftsCorrelator{1}(zeros(ComplexF64, 3), shifts)

        prns = (1, 5, 12)
        tsN = downconvert_and_correlate(
            dc,
            meas,
            TrackState(
                sig,
                [TrackedSat(sig, prn, 100.0, 200Hz; correlator = mkcorr()) for prn in prns],
            ),
        )
        for prn in prns
            alone = first(
                get_sat_state(
                    downconvert_and_correlate(
                        dc,
                        meas,
                        TrackState(
                            sig,
                            [TrackedSat(sig, prn, 100.0, 200Hz; correlator = mkcorr())],
                        ),
                    ),
                    prn,
                ).signals,
            ).correlator
            shared = first(get_sat_state(tsN, prn).signals).correlator
            @test get_accumulators(shared) == get_accumulators(alone)
        end
    end

    @testset "band-shared measurement (>1 sat) matches per-sat packing" begin
        # ≥2 sats on one band trip the pack-measurement-once-per-band path
        # (`_ob_pack_band!` + per-sat `_ob_realign_meas!`); a single sat uses the
        # per-sat pack. The shared pack is bit-identical, so every sat's correlator
        # must exactly equal correlating that PRN alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = OneBitThreadedDownconvertAndCorrelator()

        prns = (1, 5, 12)
        tsN = downconvert_and_correlate(
            dc,
            meas,
            TrackState(sig, [TrackedSat(sig, prn, 100.0, 200Hz) for prn in prns]),
        )
        for prn in prns
            alone = first(
                get_sat_state(
                    downconvert_and_correlate(
                        dc,
                        meas,
                        TrackState(sig, [TrackedSat(sig, prn, 100.0, 200Hz)]),
                    ),
                    prn,
                ).signals,
            ).correlator
            shared = first(get_sat_state(tsN, prn).signals).correlator
            @test get_prompt(shared) == get_prompt(alone)
            @test get_early(shared) == get_early(alone)
            @test get_late(shared) == get_late(alone)
        end
    end

    @testset "errors on non-Complex{Int16} (Float) measurement" begin
        sig, fs = GPSL1CA(), 5e6Hz
        capf = Complex{Float32}.(make_capture(sig, 1, fs, 5000, 200Hz, 100.0))
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        meas = (l1 = BandMeasurement(capf, fs, 0.0Hz),)
        @test_throws ArgumentError downconvert_and_correlate(
            OneBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "errors on CBOC (non-binary) code" begin
        # CBOC (e.g. Galileo E1B) carries an amplitude — a multi-level weighted sum of
        # two BOCs — so keeping only the code sign loses information; the modulation gate
        # rejects it (regardless of whether the replica is Float or a quantised integer).
        # Binary ±1 codes (BPSK, BOC, TMBOC) are accepted; the converging GPS L1 C/A test
        # below covers the BPSK case. E1B's BOC(6,1) needs fs ≥ code_freq·12 = 12.276 MHz.
        sig, fs = GalileoE1B(), 15e6Hz
        cap = make_capture(sig, 1, fs, 5000, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        @test_throws ArgumentError downconvert_and_correlate(
            OneBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "1-bit SNR loss vs Float32: tracking jitter and C/N0" begin
        # At a fixed 45 dB-Hz C/N0, the one-bit backend locks the same true Doppler as the
        # Float32 backend but with more jitter and a lower C/N0 estimate — the ≈2–3 dB cost
        # of 1-bit quantisation (≈2 dB measurement + ≈1 dB square-wave carrier). Pin the
        # *bounded* degradation (a fixed seed keeps it deterministic), not exact numbers.
        cn0_in = 45.0
        f = _track_noisy(CPUThreadedDownconvertAndCorrelator(), cn0_in, 1234)
        b = _track_noisy(OneBitThreadedDownconvertAndCorrelator(), cn0_in, 1234)

        # Both lock to the true 300 Hz Doppler — the 1-bit loss is jitter, not a bias.
        @test abs(_mean(f.dopplers) - 300) < 3
        @test abs(_mean(b.dopplers) - 300) < 3

        # One-bit carrier-Doppler jitter is worse but bounded (measured ≈1.55×; assert <3×).
        jitter_ratio = _std(b.dopplers) / _std(f.dopplers)
        @test 1.0 < jitter_ratio < 3.0

        # Float32 recovers the input C/N0; the one-bit estimate is biased low by the
        # quantisation loss (measured ≈2.5–2.9 dB), not wildly off.
        @test abs(f.cn0 - cn0_in) < 2.0
        @test b.cn0 < f.cn0
        @test 1.0 < (f.cn0 - b.cn0) < 5.0
    end

    @testset "full track converges (GPS L1 C/A)" begin
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase = 300Hz, 230.0
        nsamp = round(Int, (fs / 1Hz) * 1e-3) * 5
        cap = make_capture(sig, 1, fs, nsamp, cdopp, cphase)
        ts = TrackState(sig, [TrackedSat(sig, 1, cphase, cdopp - 20Hz)])
        ts = track(
            cap,
            ts,
            fs;
            downconvert_and_correlator = OneBitThreadedDownconvertAndCorrelator(),
        )
        # 1-bit tracking is noisier than the integer/float paths; still pulls in from
        # the 20 Hz offset toward the true Doppler.
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 15Hz
    end
end

end
