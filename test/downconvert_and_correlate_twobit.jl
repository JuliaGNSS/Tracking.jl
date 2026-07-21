module DownconvertAndCorrelateTwoBitTest

using Test: @test, @testset, @test_throws
using Random: MersenneTwister, Random
using Statistics: var
using Unitful: Hz, ustrip
import GNSSSignals
import SinCosLUT
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    GalileoE1B,
    gen_code,
    gen_code!,
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
    OneBitThreadedDownconvertAndCorrelator,
    TwoBitDownconvertAndCorrelator,
    TwoBitThreadedDownconvertAndCorrelator
import Tracking
using StaticArrays: SVector

# Dynamic-tap-count correlator (see the one-bit test file): runtime Vector shifts and
# accumulators, exercising the AbstractVector-shifts fallback; parametric over M.
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

# A bare `downconvert_and_correlate` treats the whole buffer as one chunk. If a
# code period completed, its (raw) correlator was snapshotted into
# `correlator_outputs` and the live correlator holds only the residue — return
# the first completed integration. If the buffer was shorter than one code
# period (e.g. Galileo E1B's 4 ms period in a 1 ms buffer) nothing completed, so
# the live correlator holds the whole partial integration — return that. Either
# way this matches the value the old single-step call left in the live correlator.
function _completed_or_partial_correlator(sig)
    outs = sig.correlator_outputs
    isempty(outs) ? sig.correlator : first(outs).correlator
end

# Direct (per-sample) reference of the two-bit quantised correlation: quantise the
# measurement to `{±1,±3}` at `thr`, reconstruct the `{±1,±3}` carrier from the SAME
# SinCosLUT bit planes the backend reads, take the sign of the same code replica, and sum
# the products sample by sample — following the kernel's `blk`-sample blocking so the
# per-block carrier-NCO restarts and code re-fills see identical inputs. The kernel's
# masked-popcount expansion must reproduce this exactly (all values are small integers,
# exact in Float64).
function ref_twobit(
    cap,
    sig,
    prn,
    shifts,
    code_phase,
    carrier_phase,
    code_freq,
    carrier_freq,
    fs,
    start_sample,
    nsamp,
    thr;
    blk = 8192,
)
    NC = length(shifts)
    min_shift = Int(minimum(shifts))
    span = Int(maximum(shifts)) - min_shift
    cps_car = Float64(ustrip(carrier_freq / Hz)) / Float64(ustrip(fs / Hz))
    cps_code = Float64(ustrip(code_freq / Hz)) / Float64(ustrip(fs / Hz))
    acc = zeros(ComplexF64, NC)
    bit(v, n) = (v[(n>>6)+1] >> (n & 63)) & 1
    blk_off = 0
    while blk_off < nsamp
        len = min(blk, nsamp - blk_off)
        extb = zeros(Int8, len + span + 64)
        gen_code!(
            view(extb, 1:(len+span)),
            sig,
            prn,
            fs,
            code_freq,
            code_phase + cps_code * blk_off,
            min_shift,
        )
        nwb = cld(len, 64)
        ss = zeros(UInt64, nwb)
        cs = zeros(UInt64, nwb)
        sm = zeros(UInt64, nwb)
        cm = zeros(UInt64, nwb)
        SinCosLUT.generate_carrier_signs_mags!(
            ss,
            cs,
            sm,
            cm,
            len,
            cps_car;
            phase = carrier_phase + cps_car * blk_off,
        )
        for n = 0:(len-1)
            sv = (1 - 2 * Int(bit(ss, n))) * (1 + 2 * Int(bit(sm, n)))
            cv = (1 - 2 * Int(bit(cs, n))) * (1 + 2 * Int(bit(cm, n)))
            z = cap[start_sample+blk_off+n]
            re = real(z)
            im = imag(z)
            mr = (re < 0 ? -1 : 1) * (1 + 2 * (abs(Int(re)) >= thr))
            mi = (im < 0 ? -1 : 1) * (1 + 2 * (abs(Int(im)) >= thr))
            DI = mr * cv + mi * sv
            DQ = mi * cv - mr * sv
            for k = 1:NC
                code = extb[n+(Int(shifts[k])-min_shift)+1] < 0 ? -1 : 1
                acc[k] += complex(code * DI, code * DQ)
            end
        end
        blk_off += len
    end
    acc
end

# Track a noisy GPS L1CA capture at a known C/N0 with backend `dc` (see the one-bit test
# file for the recipe), returning post-settle carrier-Doppler samples and the C/N0 estimate.
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
    amp = 10^(cn0_dbhz / 20)
    noise_std = sqrt(fs / 1Hz)
    rng = MersenneTwister(seed)
    ts = TrackState(gpsl1, [TrackedSat(gpsl1, prn, start_code_phase, cdopp)])
    dopplers = Float64[]
    for i = 0:(nepoch-1)
        carrier_phase = 2π * (cdopp / 1Hz) * nsamp * i / (fs / 1Hz) + π / 2
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

@testset "Two-bit downconvert and correlate" begin
    @testset "bit-exact vs direct quantised sum: $(nameof(typeof(sig))) @ $(fs/1e6Hz) MHz, NC=$NC, thr=$thr" for (
            sig,
            fs,
        ) in (
            (GPSL1CA(), 5e6Hz),        # single block, partial last word
            (GPSL5I(), 12e6Hz),        # two blocks
            (GPSL1CA(), 100e6Hz),      # tap offsets ≥ 64 (whole-word funnel shifts)
        ),
        NC in (3, 5),
        thr in (1, 700, 2001)

        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cdopp, cphase, prn = 250Hz, 100.0, 1
        cap = make_capture(sig, prn, fs, nsamp, cdopp, cphase)
        fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        corr =
            NC == 3 ? EarlyPromptLateCorrelator(; num_ants = NumAnts(1)) :
            VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
        shifts = get_correlator_sample_shifts(corr, fs, fc)
        dc = TwoBitDownconvertAndCorrelator(; threshold = thr)
        got = Tracking._twobit_hybrid_blocked!(
            dc,
            cap,
            NumAnts(1),
            sig,
            prn,
            shifts,
            cphase,
            0.35,
            fc,
            cdopp,
            fs,
            1,
            nsamp,
        )
        want = ref_twobit(cap, sig, prn, shifts, cphase, 0.35, fc, cdopp, fs, 1, nsamp, thr)
        @test got == want

        # the dynamic (Vector shifts) fallback is bit-identical too
        got_dyn = Tracking._twobit_hybrid_blocked!(
            dc,
            cap,
            NumAnts(1),
            sig,
            prn,
            collect(shifts),
            cphase,
            0.35,
            fc,
            cdopp,
            fs,
            1,
            nsamp,
        )
        @test got_dyn == want
    end

    @testset "bit-exact on full-range Int16 noise (incl. ±32767 / −32768)" begin
        # The magnitude plane uses a wrap-exact unsigned-range compare, so even the
        # Int16 extremes quantise correctly (an |v|-based test wraps at −32768).
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = 5000
        rng = MersenneTwister(42)
        cap = complex.(rand(rng, Int16, nsamp), rand(rng, Int16, nsamp))
        cap[1] = complex(typemin(Int16), typemax(Int16))
        cap[2] = complex(typemax(Int16), typemin(Int16))
        cdopp, cphase, prn = 250Hz, 100.0, 1
        fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = get_correlator_sample_shifts(
            EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            fs,
            fc,
        )
        for thr in (1, 512, 32767)
            dc = TwoBitDownconvertAndCorrelator(; threshold = thr)
            got = Tracking._twobit_hybrid_blocked!(
                dc,
                cap,
                NumAnts(1),
                sig,
                prn,
                shifts,
                cphase,
                0.35,
                fc,
                cdopp,
                fs,
                1,
                nsamp,
            )
            want = ref_twobit(
                cap,
                sig,
                prn,
                shifts,
                cphase,
                0.35,
                fc,
                cdopp,
                fs,
                1,
                nsamp,
                thr,
            )
            @test got == want
        end
    end

    @testset "rejects invalid threshold" begin
        @test_throws ArgumentError TwoBitDownconvertAndCorrelator(; threshold = 0)
        @test_throws ArgumentError TwoBitDownconvertAndCorrelator(; threshold = 32768)
        @test_throws ArgumentError TwoBitThreadedDownconvertAndCorrelator(; threshold = -5)
    end

    # The two-bit correlation preserves the correlation-triangle shape and its ratios
    # (E/P, L/P) at high SNR; a small tolerance absorbs the quantisation coarseness.
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
            TwoBitThreadedDownconvertAndCorrelator(; threshold = 700),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        @test abs(get_early(cb)) / abs(get_prompt(cb)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 5e-2
        @test abs(get_late(cb)) / abs(get_prompt(cb)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 5e-2
        # prompt is the correlation peak; early ≈ late (symmetric)
        @test abs(get_prompt(cb)) > abs(get_early(cb))
        @test abs(get_prompt(cb)) > abs(get_late(cb))
        @test abs(get_early(cb)) / abs(get_late(cb)) ≈ 1 atol = 5e-2
        # prompt phase within the 2-bit carrier bias (smaller than the 1-bit one)
        @test abs(mod2pi(angle(get_prompt(cb)) - angle(get_prompt(cf)) + π) - π) < 0.6
    end

    @testset "high sampling rate: tap span ≥ 64 words uncorrupted" begin
        # Tap sample-shift offsets > 64 need the whole-word + sub-word funnel shift (the
        # one-bit regression); the bit-exact reference above covers correctness, this
        # pins the correlation shape at 100 MHz through the full plumbing.
        sig, fs = GPSL1CA(), 100e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = get_correlator_sample_shifts(
            EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            fs,
            fc,
        )
        @test maximum(shifts) - minimum(shifts) ≥ 64
        cf = correlate_once(
            CPUThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        cb = correlate_once(
            TwoBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        @test abs(get_early(cb)) / abs(get_prompt(cb)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 5e-2
        @test abs(get_late(cb)) / abs(get_prompt(cb)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 5e-2
        @test abs(get_early(cb)) / abs(get_late(cb)) ≈ 1 atol = 5e-2
    end

    @testset "VeryEarlyPromptLateCorrelator (NC=5)" begin
        sig, fs = GPSL5I(), 12e6Hz          # BPSK; the bit-wise backend is BPSK/BOC-only
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        c = correlate_once(
            TwoBitThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
        )
        a = get_accumulators(c)
        @test length(a) == 5
        @test abs(a[3]) > abs(a[1]) && abs(a[3]) > abs(a[5])
        @test abs(a[3]) ≥ abs(a[2]) && abs(a[3]) ≥ abs(a[4])
    end

    @testset "multiple antennas (M=$M)" for M in (2, 4)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        dc = TwoBitThreadedDownconvertAndCorrelator()
        cm = correlate_once(
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
        c1 = correlate_once(dc, sig, fs, nsamp, 200Hz, 100.0)
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
        c1 = correlate_once(TwoBitDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        ct = correlate_once(
            TwoBitThreadedDownconvertAndCorrelator(),
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
        # A sat carrying N GPS L1 signals shares one carrier + measurement downconvert
        # (and the mask totals) via the tile-share kernel; each signal's correlator must
        # be identical to correlating it alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = TwoBitThreadedDownconvertAndCorrelator()
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

        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        dc = TwoBitThreadedDownconvertAndCorrelator()
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
        tsN = downconvert_and_correlate(
            dc,
            meas,
            TrackState(sig, satN; doppler_estimator = est),
        )
        for s in get_sat_state(tsN, 1).signals
            c = _completed_or_partial_correlator(s)
            @test length(get_prompt(c)) == M
            @test get_prompt(c) == get_prompt(cs)
            @test get_early(c) == get_early(cs)
            @test get_late(c) == get_late(cs)
        end
    end

    @testset "dynamic Vector-shifts correlator matches static EPL" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        dc = TwoBitThreadedDownconvertAndCorrelator()

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
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = TwoBitThreadedDownconvertAndCorrelator()
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

    @testset "band-shared measurement (>1 sat) matches per-sat packing @ $(fs/1e6Hz) MHz" for fs in
                                                                                              (
        5e6Hz,      # small capture → serial band pack
        40e6Hz,     # ≥ _TB_BAND_PAR_MIN full words → the band pack itself runs threaded
    )
        # ≥2 sats on one band trip the pack-measurement-once-per-band path
        # (`_tb_pack_band!` + per-sat `_tb_realign_meas!`, sign AND magnitude planes);
        # each sat's correlator must exactly equal correlating that PRN alone (which
        # packs directly, always serially) — for the threaded backend this also pins
        # the chunked threaded band pack against the serial one.
        sig = GPSL1CA()
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = TwoBitThreadedDownconvertAndCorrelator()

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
            TwoBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "errors on CBOC (non-binary) code" begin
        # CBOC (Galileo E1B) carries an amplitude, so keeping only the code sign loses
        # information; the modulation gate rejects it. E1B's BOC(6,1) needs fs ≥ 12.276 MHz.
        sig, fs = GalileoE1B(), 15e6Hz
        cap = make_capture(sig, 1, fs, 5000, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        @test_throws ArgumentError downconvert_and_correlate(
            TwoBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "two-bit tracks Float32 better than one-bit in noise" begin
        # Weak BPSK signal in Gaussian noise: the two-bit measurement + carrier recover
        # quantisation information the one-bit backend discards, so the recovered prompt
        # phase has lower variance about the noiseless-reference angle.
        Random.seed!(3)
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase, prn = 1000.0Hz, 0.0, 1
        fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        N = 20000
        A, σ = 40.0, 200.0
        d1 = OneBitThreadedDownconvertAndCorrelator()
        d2 = TwoBitThreadedDownconvertAndCorrelator(; threshold = round(Int, 0.98σ))
        ref = angle(
            get_prompt(
                correlate_once(
                    CPUThreadedDownconvertAndCorrelator(),
                    sig,
                    fs,
                    N,
                    cdopp,
                    cphase,
                ),
            ),
        )
        err1 = Float64[]
        err2 = Float64[]
        for _ = 1:40
            code = gen_code(N, sig, prn, fs, fc, cphase)
            code = code ./ maximum(abs, code)
            s = cis.(2π .* cdopp .* (0:(N-1)) ./ fs .+ 0.6) .* code
            cap = complex.(
                round.(Int16, A .* real.(s) .+ σ .* randn(N)),
                round.(Int16, A .* imag.(s) .+ σ .* randn(N)),
            )
            meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
            # Separate track states per backend: `downconvert_and_correlate`
            # shares each signal's (reused) `correlator_outputs` buffer with its
            # input, so running both backends on one `ts` would append d2's
            # records after d1's and make `first(...)` return d1's for both.
            ts1 = TrackState(sig, [TrackedSat(sig, prn, cphase, cdopp)])
            ts2 = TrackState(sig, [TrackedSat(sig, prn, cphase, cdopp)])
            c1 = _completed_or_partial_correlator(
                first(get_sat_state(downconvert_and_correlate(d1, meas, ts1), 1).signals),
            )
            c2 = _completed_or_partial_correlator(
                first(get_sat_state(downconvert_and_correlate(d2, meas, ts2), 1).signals),
            )
            push!(err1, angle(get_prompt(c1)) - ref)
            push!(err2, angle(get_prompt(c2)) - ref)
        end
        @test var(err2) < var(err1)
    end

    @testset "2-bit SNR vs Float32 and one-bit: tracking jitter and C/N0" begin
        # At a fixed 45 dB-Hz C/N0, all backends lock the same true Doppler; the losses
        # show up as C/N0-estimate bias and Doppler jitter. Measured over 8 seeds (at 45
        # and 40 dB-Hz alike): one-bit loses ≈2.8–2.9 dB of C/N0, two-bit only ≈0.85 dB
        # (theory: 0.55 dB from the 2-bit measurement + 0.25 dB from the 2-bit carrier)
        # — a ≈2 dB recovery; Doppler jitter ≈1.55× Float32 for one-bit vs ≈1.16× for
        # two-bit. Pin the *bounded* degradation and the ≥1 dB recovery over one-bit
        # with a fixed seed (this seed measures: two-bit loss 0.74 dB, one-bit loss
        # 2.76 dB, recovery 2.02 dB, jitter ratio 1.15×), not exact numbers.
        cn0_in = 45.0
        f = _track_noisy(CPUThreadedDownconvertAndCorrelator(), cn0_in, 1234)
        b1 = _track_noisy(OneBitThreadedDownconvertAndCorrelator(), cn0_in, 1234)
        b2 = _track_noisy(
            TwoBitThreadedDownconvertAndCorrelator(; threshold = 2214),  # ≈1σ of √fs noise
            cn0_in,
            1234,
        )

        @test abs(_mean(f.dopplers) - 300) < 3
        @test abs(_mean(b2.dopplers) - 300) < 3

        # two-bit jitter: worse than Float32 but bounded, and no worse than one-bit
        jitter_ratio = _std(b2.dopplers) / _std(f.dopplers)
        @test 1.0 < jitter_ratio < 2.0
        @test _std(b2.dopplers) <= _std(b1.dopplers) * 1.1

        # C/N0 estimate: biased low by ≲1.5 dB (vs ≈2.5–3 dB for one-bit), recovering
        # at least 1 dB of the one-bit backend's quantisation loss
        @test b2.cn0 < f.cn0
        @test (f.cn0 - b2.cn0) < 1.5
        @test (b2.cn0 - b1.cn0) > 1.0
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
            downconvert_and_correlator = TwoBitThreadedDownconvertAndCorrelator(),
        )
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 15Hz
    end
end

end
