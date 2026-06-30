module DownconvertAndCorrelateInt16Test

using Test: @test, @testset, @test_throws
using Unitful: Hz
import GNSSSignals
using GNSSSignals:
    GPSL1CA, GPSL5I, GalileoE1B, gen_code, get_code_center_frequency_ratio, get_code_frequency
using Tracking:
    TrackedSat, TrackState, track, downconvert_and_correlate, BandMeasurement,
    get_sat_state, get_carrier_doppler, get_prompt, get_early, get_late,
    get_accumulators, EarlyPromptLateCorrelator, VeryEarlyPromptLateCorrelator, NumAnts,
    CPUThreadedDownconvertAndCorrelator, Int16DownconvertAndCorrelator,
    Int16ThreadedDownconvertAndCorrelator

# 12-bit-ADC-style Complex{Int16} capture: carrier × unit-normalised code, scaled
# to a `peak` magnitude (≤ 2^11), rounded. Normalising the code keeps CBOC's
# ±13/±25 sub-carrier amplitudes within the 12-bit range the Int16 wipe assumes.
function make_capture(sig, prn, fs, nsamp, cdopp, cphase; peak = 2000)
    fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
    code = gen_code(nsamp, sig, prn, fs, fc, cphase)
    code = code ./ maximum(abs, code)
    s = cis.(2π .* cdopp .* (0:(nsamp-1)) ./ fs .+ 0.6) .* code
    complex.(round.(Int16, real.(s) .* peak), round.(Int16, imag.(s) .* peak))
end

band_key_for(sig) = sig isa GPSL5I ? :l5 : :l1

# M-antenna capture: a dense samples×M `Matrix` (same signal across antennas).
make_capture_mat(sig, fs, nsamp, cdopp, cphase, M; peak = 2000) =
    repeat(make_capture(sig, 1, fs, nsamp, cdopp, cphase; peak); outer = (1, M))

# One downconvert+correlate step with the given backend; returns the correlator.
function correlate_once(dc, sig, fs, nsamp, cdopp, cphase; correlator = nothing)
    cap = make_capture(sig, 1, fs, nsamp, cdopp, cphase)
    meas = NamedTuple{(band_key_for(sig),)}((BandMeasurement(cap, fs, 0.0Hz),))
    sat =
        correlator === nothing ? TrackedSat(sig, 1, cphase, cdopp) :
        TrackedSat(sig, 1, cphase, cdopp; correlator)
    ts = TrackState(sig, [sat])
    ts2 = downconvert_and_correlate(dc, meas, ts)
    first(get_sat_state(ts2, 1).signals).correlator
end

@testset "Int16 downconvert and correlate" begin
    # Correlator agreement with the Float32 backend. The integer pipeline is
    # bit-faithful on the code correlation (E/P/L magnitude ratios) and differs
    # only by the Int8 carrier's phase quantisation (a small constant prompt-phase
    # offset, ~π/64). So compare amplitude-invariant ratios tightly and the prompt
    # phase loosely.
    @testset "matches Float32 backend: $(nameof(typeof(sig))) @ $(fs/1e6Hz) MHz" for (sig, fs) in (
        (GPSL1CA(), 5e6Hz),
        (GPSL1CA(), 2e6Hz),
        (GPSL5I(), 12e6Hz),
        (GalileoE1B(), 15e6Hz),
    )
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cf = correlate_once(CPUThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        ci = correlate_once(Int16ThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        # E/P and L/P magnitude ratios match closely.
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
        # Prompt phase agrees up to the Int8 carrier quantisation.
        @test mod2pi(angle(get_prompt(ci)) - angle(get_prompt(cf)) + π) - π ≈ 0 atol = 0.1
    end

    # Parameterized (static) tap count: the kernel is @generated over
    # NC = length(sample_shifts), so the 5-tap VeryEarlyPromptLateCorrelator works
    # with no kernel change. Compare every accumulator's magnitude (normalised to
    # the prompt) against the Float32 backend.
    @testset "VeryEarlyPromptLateCorrelator (NC=5) matches Float32" begin
        sig, fs = GalileoE1B(), 15e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        mk() = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
        cf = correlate_once(CPUThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0; correlator = mk())
        ci = correlate_once(Int16ThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0; correlator = mk())
        af = get_accumulators(cf)
        ai = get_accumulators(ci)
        @test length(ai) == 5
        pf = abs(get_prompt(cf))
        pii = abs(get_prompt(ci))
        for k = 1:5
            @test abs(ai[k]) / pii ≈ abs(af[k]) / pf atol = 1e-2
        end
        @test mod2pi(angle(get_prompt(ci)) - angle(get_prompt(cf)) + π) - π ≈ 0 atol = 0.1
    end

    # Multiple antennas: the kernel runs M antenna-outer passes over the shared
    # code+carrier block and returns SVector{M,ComplexF64} accumulators. Compare
    # each antenna's E/P/L magnitude ratios against the Float32 backend.
    @testset "multiple antennas (M=$M) matches Float32" for M in (2, 4)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture_mat(sig, fs, nsamp, 200Hz, 100.0, M)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        mk() = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz; num_ants = NumAnts(M))])
        run(dc) = first(get_sat_state(downconvert_and_correlate(dc, meas, mk()), 1).signals).correlator
        cf = run(CPUThreadedDownconvertAndCorrelator())
        ci = run(Int16ThreadedDownconvertAndCorrelator())
        ef, pf, lf = get_early(cf), get_prompt(cf), get_late(cf)   # SVector{M,Complex}
        ei, pii, li = get_early(ci), get_prompt(ci), get_late(ci)
        @test length(pf) == M && length(pii) == M
        for j = 1:M
            @test abs(ei[j]) / abs(pii[j]) ≈ abs(ef[j]) / abs(pf[j]) atol = 1e-2
            @test abs(li[j]) / abs(pii[j]) ≈ abs(lf[j]) / abs(pf[j]) atol = 1e-2
            @test mod2pi(angle(pii[j]) - angle(pf[j]) + π) - π ≈ 0 atol = 0.1
        end
    end

    @testset "single-threaded and threaded backends agree" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        c1 = correlate_once(Int16DownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        ct = correlate_once(Int16ThreadedDownconvertAndCorrelator(), sig, fs, nsamp, 200Hz, 100.0)
        @test get_prompt(c1) ≈ get_prompt(ct)
        @test get_early(c1) ≈ get_early(ct)
        @test get_late(c1) ≈ get_late(ct)
    end

    @testset "errors on non-Complex{Int16} (Float) measurement" begin
        sig, fs = GPSL1CA(), 5e6Hz
        capf = Complex{Float32}.(make_capture(sig, 1, fs, 5000, 200Hz, 100.0))
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        meas = (l1 = BandMeasurement(capf, fs, 0.0Hz),)
        @test_throws ArgumentError downconvert_and_correlate(
            Int16ThreadedDownconvertAndCorrelator(), meas, ts,
        )
    end

    @testset "full track converges (GPS L1 C/A)" begin
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase = 300Hz, 230.0
        nsamp = round(Int, (fs / 1Hz) * 1e-3) * 5
        cap = make_capture(sig, 1, fs, nsamp, cdopp, cphase)
        ts = TrackState(sig, [TrackedSat(sig, 1, cphase, cdopp - 20Hz)])
        ts = track(cap, ts, fs; downconvert_and_correlator = Int16ThreadedDownconvertAndCorrelator())
        # Converges back toward the true Doppler from the 20 Hz initial offset.
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 10Hz
    end
end

end
