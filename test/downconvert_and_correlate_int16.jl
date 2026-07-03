module DownconvertAndCorrelateInt16Test

using Test: @test, @testset, @test_throws
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
    get_correlator_sample_shifts,
    update_accumulator,
    AbstractCorrelator,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    NumAnts,
    TrackedSignal,
    DefaultPostCorrFilter,
    ConventionalAssistedPLLAndDLL,
    CPUThreadedDownconvertAndCorrelator,
    Int16DownconvertAndCorrelator,
    Int16ThreadedDownconvertAndCorrelator
import Tracking

# Dynamic-tap-count correlator: sample shifts are a runtime Vector and the
# accumulators are a Vector (issue #126 (b) extension point), to exercise the
# Int16 backend's AbstractVector-shifts fallback.
struct DynShiftsCorrelator <: AbstractCorrelator{1}
    accumulators::Vector{ComplexF64}
    shifts::Vector{Int}
end
Tracking.get_accumulators(c::DynShiftsCorrelator) = c.accumulators
Tracking.update_accumulator(c::DynShiftsCorrelator, acc) =
    DynShiftsCorrelator(collect(acc), c.shifts)
Tracking.get_correlator_sample_shifts(
    c::DynShiftsCorrelator,
    sampling_frequency,
    code_frequency,
) = c.shifts

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
    @testset "matches Float32 backend: $(nameof(typeof(sig))) @ $(fs/1e6Hz) MHz" for (
        sig,
        fs,
    ) in (
        (GPSL1CA(), 5e6Hz),
        (GPSL1CA(), 2e6Hz),
        (GPSL5I(), 12e6Hz),
        (GalileoE1B(), 15e6Hz),
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
        ci = correlate_once(
            Int16ThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
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
        cf = correlate_once(
            CPUThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = mk(),
        )
        ci = correlate_once(
            Int16ThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = mk(),
        )
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
        run(dc) =
            first(get_sat_state(downconvert_and_correlate(dc, meas, mk()), 1).signals).correlator
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
        ct = correlate_once(
            Int16ThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        @test get_prompt(c1) ≈ get_prompt(ct)
        @test get_early(c1) ≈ get_early(ct)
        @test get_late(c1) ≈ get_late(ct)
    end

    # Configurable measurement amplitude: the carrier amplitude (and wipe type)
    # are derived from `max_meas`. For each declared amplitude, a capture scaled to
    # that range must still correlate correctly (E/P/L ratios are amplitude-
    # invariant, so they match the Float32 backend regardless of `max_meas`).
    @testset "configurable max_meas = $mm" for mm in (512, 2048, 8192)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0; peak = mm)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        corr(dc) = first(
            get_sat_state(
                downconvert_and_correlate(
                    dc,
                    meas,
                    TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)]),
                ),
                1,
            ).signals,
        ).correlator
        cf = corr(CPUThreadedDownconvertAndCorrelator())
        ci = corr(Int16ThreadedDownconvertAndCorrelator(; max_meas = mm))
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
    end

    # Exact-Int32 wipe fallback (used on backends without `vpmaddwd`: NEON /
    # scalar). A `max_meas` too large for any Int16-safe amplitude forces the
    # Int32 path even on x86, so this exercises and validates that path here.
    @testset "Int32 wipe fallback (forced via large max_meas)" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        dc32 = Int16ThreadedDownconvertAndCorrelator(; max_meas = 2^15)
        @test dc32 isa Int16ThreadedDownconvertAndCorrelator{<:Any,Int32}   # Int32 wipe selected
        cf = correlate_once(
            CPUThreadedDownconvertAndCorrelator(),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        ci = correlate_once(dc32, sig, fs, nsamp, 200Hz, 100.0)
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
    end

    # #167: a full-scale 16-bit `max_meas` must not silently overflow the
    # correlation accumulators. The large-`max_meas` fallback picks the exact
    # Int32 wipe with amp = 127, so the `|DI| ≤ typemax(Int16)` invariant no
    # longer holds; with near-full-scale samples and multi-level CBOC code (±25)
    # the per-block Σ codeₓ·DI runs ~25× past typemax(Int32) unless the block
    # length is shrunk to keep it bounded. Feed near-full-scale Galileo E1B and
    # check the E/P/L ratios still match the Float32 backend (they are garbage if
    # the Int32 accumulators wrapped).
    @testset "full-scale max_meas = 2^15 does not overflow (#167)" begin
        sig, fs = GalileoE1B(), 15e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0; peak = 2^15 - 1)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        corr(dc) = first(
            get_sat_state(
                downconvert_and_correlate(
                    dc,
                    meas,
                    TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)]),
                ),
                1,
            ).signals,
        ).correlator
        cf = corr(CPUThreadedDownconvertAndCorrelator())
        ci = corr(Int16ThreadedDownconvertAndCorrelator(; max_meas = 2^15))
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
    end

    # Dynamic (runtime Vector) tap count: a DynShiftsCorrelator with the same
    # shifts as EPL must produce the same accumulators through the Int16 backend's
    # AbstractVector-shifts fallback as the static @generated EPL kernel.
    @testset "dynamic Vector-shifts correlator matches static EPL" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        dc = Int16ThreadedDownconvertAndCorrelator()
        cs = correlate_once(dc, sig, fs, nsamp, 200Hz, 100.0)                       # static EPL
        cd = correlate_once(
            dc,
            sig,
            fs,
            nsamp,
            200Hz,
            100.0;
            correlator = DynShiftsCorrelator(zeros(ComplexF64, 3), shifts),
        )         # dynamic
        ad = get_accumulators(cd)
        as = get_accumulators(cs)
        @test length(ad) == 3
        # Same shift order ([-d, 0, +d]) and the same exact integer sums (same code,
        # carrier, DI), so they agree element-wise to floating-point round-off.
        for k = 1:3
            @test ad[k] ≈ as[k] rtol = 1e-9
        end
    end

    # Multi-signal-per-sat tile-share: a sat carrying N identical GPS L1 C/A
    # signals shares one carrier downconvert. Each signal's correlator must equal
    # the single-signal result (same code/carrier; only the downconvert is shared).
    @testset "multi-signal-per-sat tile-share (N=$N) matches single signal" for N in (2, 3)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = Int16ThreadedDownconvertAndCorrelator()
        est = ConventionalAssistedPLLAndDLL()
        mksig() = TrackedSignal(
            sig;
            num_ants = NumAnts(1),
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            post_corr_filter = DefaultPostCorrFilter(),
        )
        # single-signal reference
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
        # N-signal sat
        satN = TrackedSat(ntuple(_ -> mksig(), N), 1, 100.0, 200Hz; doppler_estimator = est)
        tsN = downconvert_and_correlate(
            dc,
            meas,
            TrackState(sig, satN; doppler_estimator = est),
        )
        for s in get_sat_state(tsN, 1).signals
            @test get_prompt(s.correlator) ≈ get_prompt(cs)
            @test get_early(s.correlator) ≈ get_early(cs)
            @test get_late(s.correlator) ≈ get_late(cs)
        end
    end

    @testset "errors on non-Complex{Int16} (Float) measurement" begin
        sig, fs = GPSL1CA(), 5e6Hz
        capf = Complex{Float32}.(make_capture(sig, 1, fs, 5000, 200Hz, 100.0))
        ts = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz)])
        meas = (l1 = BandMeasurement(capf, fs, 0.0Hz),)
        # Both integer backends run the sample-type check through the shared
        # `_check_sample_type` hook, so both must reject Float samples.
        @test_throws ArgumentError downconvert_and_correlate(
            Int16DownconvertAndCorrelator(),
            meas,
            ts,
        )
        @test_throws ArgumentError downconvert_and_correlate(
            Int16ThreadedDownconvertAndCorrelator(),
            meas,
            ts,
        )
    end

    @testset "rejects a `steps` that changes the engine width (#168)" begin
        # The kernels stride by the compile-time `_INT16_W`; a `steps` whose
        # SinCosLUT backend has a different SIMD width would silently corrupt the
        # carrier. `steps = 48` is never a supported permute size → Portable
        # (width 1), so on any non-portable host it mismatches and must be
        # rejected at construction rather than producing garbage correlations.
        if Tracking._INT16_W != 1
            @test_throws ArgumentError Int16DownconvertAndCorrelator(steps = 48)
            @test_throws ArgumentError Int16ThreadedDownconvertAndCorrelator(steps = 48)
        end
        # The default `steps = 64` always matches the kernel stride.
        @test Int16DownconvertAndCorrelator(steps = 64) isa Int16DownconvertAndCorrelator
        @test Int16ThreadedDownconvertAndCorrelator(steps = 64) isa
              Int16ThreadedDownconvertAndCorrelator
    end

    # Block-flush horizontal sum must not wrap in Int32 (#165). The per-lane Int32
    # bound holds on the fast path, but the block-flush reduction ACROSS lanes is the
    # full block total, which for a strong, code-aligned, full-scale CBOC capture at
    # DEFAULT settings exceeds typemax(Int32). If that reduction happens in Int32 the
    # prompt wraps and the E/P·L/P triangle is destroyed. Worst case from the issue:
    # full-scale (|m| = max_meas on every sample) real signal carrying the E1B code
    # sign, zero Doppler → Σ code·DI is same-sign across the whole 8192-sample block.
    @testset "block-flush accumulator does not wrap in Int32 (#165)" begin
        sig, fs = GalileoE1B(), 15e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)   # 15000 > 8192-sample block
        cphase = 0.0
        code = gen_code(nsamp, sig, 1, fs, get_code_frequency(sig), cphase)
        cap = complex.(round.(Int16, 2048 .* sign.(real.(code))), zero(Int16))
        meas = (l1 = BandMeasurement(cap, fs, 0.0Hz),)
        corr(dc) = first(
            get_sat_state(
                downconvert_and_correlate(
                    dc,
                    meas,
                    TrackState(sig, [TrackedSat(sig, 1, cphase, 0.0Hz)]),
                ),
                1,
            ).signals,
        ).correlator
        cf = corr(CPUThreadedDownconvertAndCorrelator())
        ci = corr(Int16ThreadedDownconvertAndCorrelator())
        # The capture is strong enough that the integer prompt total exceeds
        # typemax(Int32): this is exactly what wraps an Int32 block-flush reduction.
        @test abs(get_prompt(ci)) > typemax(Int32)
        # Amplitude-invariant E/P and L/P ratios still match the (overflow-free)
        # Float32 backend — they would be garbage if the prompt had wrapped.
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
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
            downconvert_and_correlator = Int16ThreadedDownconvertAndCorrelator(),
        )
        # Converges back toward the true Doppler from the 20 Hz initial offset.
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 10Hz
    end
end

end
