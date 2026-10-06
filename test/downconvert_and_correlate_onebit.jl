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
    get_band_id,
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
    MomentsCN0Estimator,
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
include("repeated_signal_test_helpers.jl")  # `_repeated_signal_track_state`
using StaticArrays: SVector

# Runtime-`Vector` shifts and accumulators (issue #126 (b)), to exercise the one-bit
# backend's AbstractVector-shifts fallback; parametric in M for both its branches.
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

band_key_for(sig) = get_band_id(sig)
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
    _completed_or_partial_correlator(first(get_sat_state(ts2, 1).signals))
end

# First completed integration, or the live partial correlator if none completed
# (see the same helper in downconvert_and_correlate_int16.jl).
function _completed_or_partial_correlator(sig)
    outs = sig.correlator_outputs
    isempty(outs) ? sig.correlator : first(outs).correlator
end

# Track a noisy GPS L1CA capture at a known C/N0 with backend `dc`; returns the
# post-settle carrier Dopplers and the final C/N0 estimate. The capture is
# Complex{Int16}, so Float32 and one-bit backends see identical bits; the C/N0
# scaling follows the CN0-estimation tests.
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
    # Pin the moment estimator: NWPR also reacts to the quantisation's residual
    # phase noise, which would fold a second effect into the measured C/N0 loss.
    ts = TrackState(
        gpsl1,
        [
            TrackedSat(
                gpsl1,
                prn,
                start_code_phase,
                cdopp;
                cn0_estimator = MomentsCN0Estimator(100),
            ),
        ],
    )
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
    # With a ±1 code the E/P·L/P ratios survive 1-bit quantisation at high SNR; the
    # prompt phase gets a loose bound (the square-wave carrier adds a phase bias).
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

    @testset "high sampling rate: tap span ≥ 64 words uncorrupted (regression)" begin
        # At 100 MHz the tap offsets exceed one UInt64 word, exercising the whole-word
        # part of `_ob_shift_plane!`.
        sig, fs = GPSL1CA(), 100e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = get_correlator_sample_shifts(
            EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            fs,
            fc,
        )
        @test maximum(shifts) - minimum(shifts) ≥ 64        # the regime the bug lived in
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
        # early and late are symmetric about the peak → equal magnitude
        @test abs(get_early(cb)) / abs(get_late(cb)) ≈ 1 atol = 5e-2
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
        # Each signal's correlator must equal correlating it alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        tsN = downconvert_and_correlate(dc, meas, _repeated_signal_track_state(satN, est))
        for s in get_sat_state(tsN, 1).signals
            @test get_prompt(s.correlator) == get_prompt(cs)
            @test get_early(s.correlator) == get_early(cs)
            @test get_late(s.correlator) == get_late(cs)
        end
    end

    @testset "multi-signal-per-sat × multi-antenna (N=$N, M=$M) matches single" for N in
                                                                                    (2, 3),
        M in (2, 4)
        # Tile-share kernel's M>1 path: the shared downconvert is bit-identical, so
        # each signal's per-antenna correlator must equal correlating it alone.
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
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        est = ConventionalAssistedPLLAndDLL()
        mksig() = TrackedSignal(
            sig;
            num_ants = NumAnts(M),
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(M)),
            post_corr_filter = DefaultPostCorrFilter(),
        )
        satN = TrackedSat(ntuple(_ -> mksig(), N), 1, 100.0, 200Hz; doppler_estimator = est)
        tsN = downconvert_and_correlate(dc, meas, _repeated_signal_track_state(satN, est))
        for s in get_sat_state(tsN, 1).signals
            c = _completed_or_partial_correlator(s)
            @test length(get_prompt(c)) == M
            @test get_prompt(c) == get_prompt(cs)
            @test get_early(c) == get_early(cs)
            @test get_late(c) == get_late(cs)
        end
    end

    @testset "dynamic Vector-shifts correlator matches static EPL" begin
        # The AbstractVector-shifts fallback must be bit-exact with the static
        # @generated EPL kernel (both sum integer popcounts).
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
        # ≥2 sats trip the fallback's band-shared measurement path
        # (`_ob_realign_meas!`); each sat must equal correlating that PRN alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        # ≥2 sats trip the pack-once-per-band path (`_ob_pack_band!` +
        # `_ob_realign_meas!`); each sat must equal correlating that PRN alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        meas = (L1 = BandMeasurement(capf, fs, 0.0Hz),)
        @test_throws ArgumentError downconvert_and_correlate(
            OneBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "errors on CBOC (non-binary) code" begin
        # The CBOC modulation gate (see `_onebit_hybrid_blocked!`). E1B's BOC(6,1)
        # needs fs ≥ 12.276 MHz.
        sig, fs = GalileoE1B(), 15e6Hz
        cap = make_capture(sig, 1, fs, 5000, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        @test_throws ArgumentError downconvert_and_correlate(
            OneBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
        # Named via `get_signal_name`, not the parametrised type.
        @test_throws "Galileo E1B" downconvert_and_correlate(
            OneBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "1-bit SNR loss vs Float32: tracking jitter and C/N0" begin
        # Pin the bounded ≈2–3 dB 1-bit loss (see the header of
        # src/downconvert_and_correlate_onebit.jl), not exact numbers; fixed seed.
        cn0_in = 45.0
        f = _track_noisy(CPUThreadedDownconvertAndCorrelator(), cn0_in, 1234)
        b = _track_noisy(OneBitThreadedDownconvertAndCorrelator(), cn0_in, 1234)

        # The 1-bit loss is jitter, not a bias.
        @test abs(_mean(f.dopplers) - 300) < 3
        @test abs(_mean(b.dopplers) - 300) < 3

        # One-bit carrier-Doppler jitter is worse but bounded (measured ≈1.55×; assert <3×).
        jitter_ratio = _std(b.dopplers) / _std(f.dopplers)
        @test 1.0 < jitter_ratio < 3.0

        # One-bit C/N0 is biased low by the loss (measured ≈2.5–2.9 dB).
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
        # Noisier than the integer/float paths, hence the wider tolerance.
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 15Hz
    end
end

end
