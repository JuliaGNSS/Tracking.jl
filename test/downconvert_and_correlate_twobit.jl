module DownconvertAndCorrelateTwoBitTest

using Test: @test, @testset, @test_throws
using Random: Random
using Statistics: var
using Unitful: Hz
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
    OneBitThreadedDownconvertAndCorrelator,
    TwoBitDownconvertAndCorrelator,
    TwoBitThreadedDownconvertAndCorrelator

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

@testset "Two-bit downconvert and correlate" begin
    # The two-bit sign+magnitude correlation preserves the correlation-triangle shape and
    # its ratios (E/P, L/P) at high SNR, for both a 1-bit carrier (CB=1) and a 2-bit
    # carrier (CB=2); a small tolerance absorbs the quantisation coarseness.
    @testset "ratios match Float32: $(nameof(typeof(sig))) @ $(fs/1e6Hz) MHz, CB=$CB" for (
            sig,
            fs,
        ) in (
            (GPSL1CA(), 5e6Hz),
            (GPSL1CA(), 20e6Hz),
            (GPSL5I(), 12e6Hz),
        ),
        CB in (1, 2)

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
            TwoBitDownconvertAndCorrelator(; carrier_bits = CB, threshold = 700),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        @test abs(get_early(cb)) / abs(get_prompt(cb)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 8e-2
        @test abs(get_late(cb)) / abs(get_prompt(cb)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 8e-2
        @test abs(get_prompt(cb)) > abs(get_early(cb))
        @test abs(get_prompt(cb)) > abs(get_late(cb))
        @test abs(get_early(cb)) / abs(get_late(cb)) ≈ 1 atol = 8e-2
        # prompt phase agrees with Float32 (bit-wise loses amplitude, not phase)
        @test abs(mod2pi(angle(get_prompt(cb)) - angle(get_prompt(cf)) + π) - π) < 0.6
    end

    @testset "VeryEarlyPromptLateCorrelator (NC=5), CB=$CB" for CB in (1, 2)
        sig, fs = GPSL5I(), 12e6Hz          # BPSK; the bit-wise backend is BPSK-only
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        c = correlate_once(
            TwoBitThreadedDownconvertAndCorrelator(; carrier_bits = CB),
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

    @testset "multiple antennas (M=$M), CB=$CB" for M in (2, 4), CB in (1, 2)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        dc = TwoBitThreadedDownconvertAndCorrelator(; carrier_bits = CB)
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

    @testset "single-threaded and threaded backends agree (exact), CB=$CB" for CB in (1, 2)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        c1 = correlate_once(
            TwoBitDownconvertAndCorrelator(; carrier_bits = CB),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        ct = correlate_once(
            TwoBitThreadedDownconvertAndCorrelator(; carrier_bits = CB),
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
        # A sat carrying N GPS L1 signals; the backend correlates each in turn, so each
        # signal's correlator is identical to correlating it alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
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

    @testset "band-shared measurement (>1 sat) matches per-sat packing" begin
        # ≥2 sats on one band trip the pack-measurement-once-per-band path
        # (`_tb_pack_band!` + per-sat `_tb_realign_meas!`); a single sat uses the per-sat
        # pack. The shared pack is bit-identical, so every sat's correlator must exactly
        # equal correlating that PRN alone.
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        meas = (l1 = BandMeasurement(capf, fs, 0.0Hz),)
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
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        @test_throws ArgumentError downconvert_and_correlate(
            TwoBitThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "rejects invalid carrier_bits" begin
        @test_throws ArgumentError TwoBitDownconvertAndCorrelator(; carrier_bits = 0)
        @test_throws ArgumentError TwoBitDownconvertAndCorrelator(; carrier_bits = 3)
        @test_throws ArgumentError TwoBitThreadedDownconvertAndCorrelator(;
            carrier_bits = 3,
        )
    end

    @testset "two-bit tracks Float32 better than one-bit in noise" begin
        # Weak BPSK signal in Gaussian noise: the two-bit measurement recovers input
        # quantisation the one-bit backend discards, so its recovered prompt phase has
        # lower variance about the noiseless-reference angle (higher effective SNR).
        Random.seed!(3)
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase, prn = 1000.0Hz, 0.0, 1
        fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        N = 20000
        A, σ = 40.0, 200.0
        d1 = OneBitThreadedDownconvertAndCorrelator()
        d2 = TwoBitThreadedDownconvertAndCorrelator(;
            carrier_bits = 2,
            threshold = round(Int, 0.98σ),
        )
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
            meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
            ts = TrackState(sig, [TrackedSat(sig, prn, cphase, cdopp)])
            c1 =
                first(get_sat_state(downconvert_and_correlate(d1, meas, ts), 1).signals).correlator
            c2 =
                first(get_sat_state(downconvert_and_correlate(d2, meas, ts), 1).signals).correlator
            push!(err1, angle(get_prompt(c1)) - ref)
            push!(err2, angle(get_prompt(c2)) - ref)
        end
        @test var(err2) < var(err1)
    end

    @testset "full track converges (GPS L1 C/A), CB=$CB" for CB in (1, 2)
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase = 300Hz, 230.0
        nsamp = round(Int, (fs / 1Hz) * 1e-3) * 5
        cap = make_capture(sig, 1, fs, nsamp, cdopp, cphase)
        ts = TrackState(sig, [TrackedSat(sig, 1, cphase, cdopp - 20Hz)])
        ts = track(
            cap,
            ts,
            fs;
            downconvert_and_correlator = TwoBitThreadedDownconvertAndCorrelator(;
                carrier_bits = CB,
            ),
        )
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 15Hz
    end
end

end
