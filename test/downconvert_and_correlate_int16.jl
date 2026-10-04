module DownconvertAndCorrelateInt16Test

using Test: @test, @testset, @test_throws
using Unitful: Hz
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

# Runtime-`Vector` shifts and accumulators (issue #126 (b)), to exercise the Int16
# backend's AbstractVector-shifts fallback.
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

# 12-bit-ADC-style Complex{Int16} capture: carrier × unit-normalised code (keeps
# CBOC's ±13/±25 levels within range), scaled to `peak` and rounded.
function make_capture(sig, prn, fs, nsamp, cdopp, cphase; peak = 2000)
    fc = cdopp * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
    code = gen_code(nsamp, sig, prn, fs, fc, cphase)
    code = code ./ maximum(abs, code)
    s = cis.(2π .* cdopp .* (0:(nsamp-1)) ./ fs .+ 0.6) .* code
    complex.(round.(Int16, real.(s) .* peak), round.(Int16, imag.(s) .* peak))
end

# Tightest `max_meas` covering the default `peak = 2000` captures.
const TEST_MAX_MEAS = 2^11

band_key_for(sig) = get_band_id(sig)

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
    _completed_or_partial_correlator(first(get_sat_state(ts2, 1).signals))
end

# First completed integration from `correlator_outputs`, or the live (partial)
# correlator if the buffer was shorter than one code period (e.g. E1B's 4 ms in a
# 1 ms buffer). See `update` in src/downconvert_and_correlate.jl.
function _completed_or_partial_correlator(sig)
    outs = sig.correlator_outputs
    isempty(outs) ? sig.correlator : first(outs).correlator
end

@testset "Int16 downconvert and correlate" begin
    # The integer pipeline differs from Float32 only by the Int8 carrier's
    # quantisation (~π/64 phase offset): compare E/P/L ratios tightly, the prompt
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
            Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS),
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
        # Absolute magnitude agrees too (finalize divides out the carrier amplitude);
        # the residual is the Int8 amplitude quantisation.
        @test abs(get_prompt(ci)) ≈ abs(get_prompt(cf)) rtol = 0.05
        # Prompt phase agrees up to the Int8 carrier quantisation.
        @test mod2pi(angle(get_prompt(ci)) - angle(get_prompt(cf)) + π) - π ≈ 0 atol = 0.1
    end

    # The kernel is @generated over the tap count NC: a 5-tap correlator must match
    # Float32 on every prompt-normalised accumulator magnitude.
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
            Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS),
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

    # Multi-antenna path (shared code+carrier block, SVector{M} accumulators):
    # per-antenna E/P/L ratios must match Float32.
    @testset "multiple antennas (M=$M) matches Float32" for M in (2, 4)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture_mat(sig, fs, nsamp, 200Hz, 100.0, M)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        mk() = TrackState(sig, [TrackedSat(sig, 1, 100.0, 200Hz; num_ants = NumAnts(M))])
        run(dc) = _completed_or_partial_correlator(
            first(get_sat_state(downconvert_and_correlate(dc, meas, mk()), 1).signals),
        )
        cf = run(CPUThreadedDownconvertAndCorrelator())
        ci = run(Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS))
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
        c1 = correlate_once(
            Int16DownconvertAndCorrelator(TEST_MAX_MEAS),
            sig,
            fs,
            nsamp,
            200Hz,
            100.0,
        )
        ct = correlate_once(
            Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS),
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

    # Carrier amplitude and wipe type derive from `max_meas` (see
    # `_int16_choose_carrier`); a capture at each declared range must still match.
    @testset "configurable max_meas = $mm" for mm in (512, 2048, 8192)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0; peak = mm)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        ci = corr(Int16ThreadedDownconvertAndCorrelator(mm))
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
    end

    # A `max_meas` too large for an Int16-safe amplitude forces the exact-Int32 wipe
    # on any host, so this covers that path even on x86.
    @testset "Int32 wipe fallback (forced via large max_meas)" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        dc32 = Int16ThreadedDownconvertAndCorrelator(2^15)
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

    # #167: near-full-scale E1B with `max_meas = 2^15` (amp = 127 fallback) must not
    # wrap the Int32 accumulators (see `_int16_safe_blk`); wrapping would garble the
    # E/P/L ratios.
    @testset "full-scale max_meas = 2^15 does not overflow (#167)" begin
        sig, fs = GalileoE1B(), 15e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0; peak = 2^15 - 1)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        ci = corr(Int16ThreadedDownconvertAndCorrelator(2^15))
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
    end

    # The AbstractVector-shifts fallback must match the static @generated EPL kernel.
    @testset "dynamic Vector-shifts correlator matches static EPL" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts = collect(get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc))
        dc = Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS)
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
        # Same shifts and exact integer sums, so equal up to round-off.
        for k = 1:3
            @test ad[k] ≈ as[k] rtol = 1e-9
        end
    end

    # The AbstractVector-shifts fallback should allocate only its length-NC result
    # (~48 B). Measured as the gap to the static path, which cancels the shared
    # per-block `gen_code!`/`fill_carrier!` cost (allocating on some Julia
    # versions); `minimum` strips GC-sampling noise.
    @testset "dynamic Vector-shifts fallback adds no per-call scratch vs static path" begin
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        fc = 200Hz * get_code_center_frequency_ratio(sig) + get_code_frequency(sig)
        shifts_s = get_correlator_sample_shifts(EarlyPromptLateCorrelator(), fs, fc)  # SVector → @generated
        shifts_v = collect(shifts_s)                                                  # Vector → fallback
        dc = Int16DownconvertAndCorrelator(TEST_MAX_MEAS)
        runN(shifts, cap, ns) = Tracking._int16_hybrid_blocked!(
            dc,
            cap,
            NumAnts(1),
            sig,
            1,
            shifts,
            100.0,
            0.6,
            fc,
            200Hz + 0.0Hz,
            fs,
            1,
            ns,
        )
        function minalloc(shifts, mult)
            ns = mult * nsamp
            cap = make_capture(sig, 1, fs, ns, 200Hz, 100.0)
            runN(shifts, cap, ns)
            runN(shifts, cap, ns)               # warm up (grow scratch once)
            minimum(@allocated(runN(shifts, cap, ns)) for _ = 1:8)
        end
        # 1 vs ~13 blocks: the gap must be small and not grow with block count.
        for mult in (1, 13)
            gap = minalloc(shifts_v, mult) - minalloc(shifts_s, mult)
            @test gap <= 128
        end
    end

    # N identical signals on one sat share the carrier downconvert; each must equal
    # the single-signal result.
    @testset "multi-signal-per-sat tile-share (N=$N) matches single signal" for N in (2, 3)
        sig, fs = GPSL1CA(), 5e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)
        cap = make_capture(sig, 1, fs, nsamp, 200Hz, 100.0)
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
        dc = Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS)
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
        meas = (L1 = BandMeasurement(capf, fs, 0.0Hz),)
        # Both backends share the `_check_sample_type` hook.
        @test_throws ArgumentError downconvert_and_correlate(
            Int16DownconvertAndCorrelator(TEST_MAX_MEAS),
            meas,
            ts,
        )
        @test_throws ArgumentError downconvert_and_correlate(
            Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS),
            meas,
            ts,
        )
    end

    @testset "max_meas is a required positional argument (no default)" begin
        # No default, so under-declaring can't happen silently (see
        # `_int16_choose_carrier`).
        @test_throws MethodError Int16DownconvertAndCorrelator()
        @test_throws MethodError Int16ThreadedDownconvertAndCorrelator()
    end

    @testset "rejects a `steps` that changes the engine width (#168)" begin
        # `steps = 48` always selects the Portable (width 1) table, which mismatches
        # `_INT16_W` on any non-portable host (see `_int16_assert_engine_width`).
        if Tracking._INT16_W != 1
            @test_throws ArgumentError Int16DownconvertAndCorrelator(
                TEST_MAX_MEAS;
                steps = 48,
            )
            @test_throws ArgumentError Int16ThreadedDownconvertAndCorrelator(
                TEST_MAX_MEAS;
                steps = 48,
            )
        end
        # The default `steps = 64` always matches the kernel stride.
        @test Int16DownconvertAndCorrelator(TEST_MAX_MEAS; steps = 64) isa
              Int16DownconvertAndCorrelator
        @test Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS; steps = 64) isa
              Int16ThreadedDownconvertAndCorrelator
    end

    # #165: the cross-lane block-flush sum must not wrap Int32 (see `_wide64`).
    # Worst case from the issue: full-scale, code-aligned E1B at zero Doppler, so
    # Σ code·DI is same-sign across the whole default block.
    @testset "block-flush accumulator does not wrap in Int32 (#165)" begin
        sig, fs = GalileoE1B(), 15e6Hz
        nsamp = round(Int, (fs / 1Hz) * 1e-3)   # 15000 > 8192-sample block
        cphase = 0.0
        code = gen_code(nsamp, sig, 1, fs, get_code_frequency(sig), cphase)
        cap = complex.(round.(Int16, 2048 .* sign.(real.(code))), zero(Int16))
        meas = (L1 = BandMeasurement(cap, fs, 0.0Hz),)
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
        ci_dc = Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS)
        cf = corr(CPUThreadedDownconvertAndCorrelator())
        ci = corr(ci_dc)
        # Precondition: the raw block total (finalize divided out the carrier
        # amplitude) really exceeds typemax(Int32).
        @test abs(get_prompt(ci)) * ci_dc.carrier_amplitude > typemax(Int32)
        # Ratios would be garbage if the prompt had wrapped.
        @test abs(get_early(ci)) / abs(get_prompt(ci)) ≈
              abs(get_early(cf)) / abs(get_prompt(cf)) atol = 1e-2
        @test abs(get_late(ci)) / abs(get_prompt(ci)) ≈
              abs(get_late(cf)) / abs(get_prompt(cf)) atol = 1e-2
    end

    # Per-lane Int32 overflow guard `_int16_flush_len` (issue #166), checked at every
    # host SIMD width.
    @testset "flush cadence bounds the Int32 accumulator (issue #166)" begin
        max_code = 25                       # Galileo E1B CBOC sub-carrier peak
        max_wipe = Int(typemax(Int16))      # Int16-safe wipe bound (|DI| ≤ typemax(Int16))
        max_product = max_code * max_wipe
        blk = Tracking._INT16_BLK
        for W in (1, 16, 32, 64)
            L = Tracking._int16_flush_len(W, max_wipe, blk)
            @test 0 < L <= blk
            # 2× covers the vpmaddwd path.
            @test 2 * cld(L, W) * max_product <= typemax(Int32)
        end
        # The Portable (W == 1) path must flush before a full default block …
        @test Tracking._int16_flush_len(1, max_wipe, blk) < blk
        # … while wide-SIMD hosts are unaffected (no benchmark impact).
        @test Tracking._int16_flush_len(16, max_wipe, blk) == blk
        @test Tracking._int16_flush_len(32, max_wipe, blk) == blk
    end

    @testset "blk ≤ 0 is rejected at construction (issue #169)" begin
        # blk ≤ 0 would hang track! (see `_int16_validate_blk`).
        for bad in (0, -1, -8192)
            @test_throws ArgumentError Int16DownconvertAndCorrelator(
                TEST_MAX_MEAS;
                blk = bad,
            )
            @test_throws ArgumentError Int16ThreadedDownconvertAndCorrelator(
                TEST_MAX_MEAS;
                blk = bad,
            )
        end

        # A positive blk is accepted and stored verbatim (default and custom).
        @test Int16DownconvertAndCorrelator(TEST_MAX_MEAS).blk == 8192
        @test Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS).blk == 8192
        @test Int16DownconvertAndCorrelator(TEST_MAX_MEAS; blk = 4096).blk == 4096

        # A large blk is accepted and clamped at run time (see `_int16_validate_blk`);
        # the correlation must match the default-blk result.
        @test Int16DownconvertAndCorrelator(TEST_MAX_MEAS; blk = 200_000).blk == 200_000
        @test Int16ThreadedDownconvertAndCorrelator(TEST_MAX_MEAS; blk = 200_000).blk ==
              200_000
        let sig = GalileoE1B(), fs = 15e6Hz
            nsamp = round(Int, (fs / 1Hz) * 1e-3)
            cbig = correlate_once(
                Int16DownconvertAndCorrelator(TEST_MAX_MEAS; blk = 200_000),
                sig,
                fs,
                nsamp,
                200Hz,
                100.0,
            )
            cdef = correlate_once(
                Int16DownconvertAndCorrelator(TEST_MAX_MEAS),
                sig,
                fs,
                nsamp,
                200Hz,
                100.0,
            )
            @test abs(get_early(cbig)) / abs(get_prompt(cbig)) ≈
                  abs(get_early(cdef)) / abs(get_prompt(cdef)) atol = 1e-2
            @test abs(get_late(cbig)) / abs(get_prompt(cbig)) ≈
                  abs(get_late(cdef)) / abs(get_prompt(cdef)) atol = 1e-2
        end
    end

    @testset "full track converges (GPS L1 C/A)" begin
        sig, fs = GPSL1CA(), 5e6Hz
        cdopp, cphase = 300Hz, 230.0
        # 100 ms: long enough for the 18 Hz carrier loop to settle a 20 Hz step.
        nsamp = round(Int, (fs / 1Hz) * 1e-3) * 100
        cap = make_capture(sig, 1, fs, nsamp, cdopp, cphase)
        ts = TrackState(sig, [TrackedSat(sig, 1, cphase, cdopp - 20Hz)])
        ts = track(
            cap,
            ts,
            fs;
            downconvert_and_correlator = Int16ThreadedDownconvertAndCorrelator(
                TEST_MAX_MEAS,
            ),
        )
        # Converges back toward the true Doppler from the 20 Hz initial offset.
        @test get_carrier_doppler(get_sat_state(ts, 1)) ≈ cdopp atol = 10Hz
    end
end

end
