# Integer (`Complex{Int16}`) hybrid-blocked downconvert + correlate backend.
#
# This is the integer counterpart of the Float32 fused path
# (`downconvert_and_correlate_fused.jl`). It targets `Complex{Int16}` (12-bit
# ADC) sample buffers and gets its speed from an all-integer pipeline, ported
# from the GNSSSignals integration benchmark's `correlate_epl_hybrid_blocked!`
# (the fastest variant there) and generalised to Tracking's arbitrary correlator
# tap count, antenna count, and `ComplexF64` accumulators:
#
#   * code   — Int8 ±1 (CBOC ±13/±25) replica from GNSSSignals' embedded SIMD
#              LUT, generated per block by the value-threaded `CodeFillEngine`
#              (`code_engine(signal, prn, fs, fc; …)` + `code_state` +
#              `gen_code!(out, eng, st)`): allocation-free, thread-safe (the
#              caller threads the state), exact across block seams, and applies
#              non-baked secondaries (e.g. GPS L5I NH10).
#   * carrier — Int8 sin/cos from SinCosLUT's register-resident table lookup,
#              filled per block 4-way-unrolled (`_int16_fill_carrier!`).
#   * wipe + correlate — the shared-carrier wipe DI = mᵣ·cos + mᵢ·sin,
#              DQ = mᵢ·cos − mᵣ·sin runs in Int16 when `2·max|meas|·amp ≤
#              typemax(Int16)` (`choose_carrier`), and the code accumulate
#              Σ codeₓ·DI uses the fused widening `vpmaddwd` (Int16×Int16→Int32)
#              on x86; otherwise the exact Int32 path. Per-tap lane accumulators
#              are flushed into Int64 totals every block so no integration length
#              or code amplitude can overflow.
#
# Multiple antennas (M>1): the sample buffer is a dense `Matrix` (rows = samples,
# columns = antennas). One code + carrier block is filled per strip-mine block
# and shared across `M` antenna-outer correlate passes (each pass keeps only its
# own NC tap accumulators live, mirroring the Float32 tile-share kernel), so the
# carrier wipe-off is computed once-per-antenna against the shared code/carrier.
#
# Strategy is "hybrid-blocked": strip-mine the integration into `blk`-sample
# blocks, regenerating BOTH code and carrier into small, L1-resident scratch
# reused across blocks. See docs/plans/2026-06-30-int16-hybrid-blocked-...
#
# Accumulators are converted to `ComplexF64` (M=1) / `SVector{M,ComplexF64}`
# (M>1), scaled by the constant carrier amplitude, at finalize — so every
# downstream consumer (discriminators, C/N0, bit buffer — all ratio/normalised)
# is unaffected.
#
# Scope: one signal per sat, static correlator tap counts (EPL NC=3, VEPL NC=5);
# the kernel is `@generated` over (NC, M).

import SinCosLUT
using SinCosLUT: SinCosTable, carrier_engine, carrier_state, carrier_lookup, carrier_advance

# SIMD width of the Int8 carrier/code LUT backend on this host (compile-time
# const), mirroring the benchmark's `_CORR_W`. SinCosLUT's Int8 table and the
# GNSSSignals Int8 code LUT use the same width per backend.
const _INT16_W = let be = SinCosLUT.default_backend(Int8, 64)
    be isa SinCosLUT.AVX512 ? 64 :
    be isa SinCosLUT.AVX2 ? 32 :
    be isa SinCosLUT.Neon ? 16 : 1
end

# The fused widening `vpmaddwd` accumulate (Int16×Int16→Int32 pairwise) is x86
# (AVX2/AVX-512) only. Other backends use the exact Int32 multiply path.
const _INT16_HAS_MADDWD = Sys.ARCH in (:x86_64, :i686) && _INT16_W in (32, 64)

# Maximum |measurement component| we size the Int16 wipe against: a 12-bit ADC
# (|m| ≤ 2^11). Pick the LARGEST Int16-safe carrier amplitude (carrier rounding
# error is ±0.5 regardless of amplitude, so a bigger amplitude is a finer
# carrier), capped at the Int8 storage limit; if none is safe (or no SIMD wipe),
# fall back to the exact Int32 wipe at full Int8 amplitude. Mirrors the
# benchmark's `choose_carrier`.
const _INT16_MAX_MEAS = 1 << 11
function _int16_choose_carrier(max_meas::Integer)
    _INT16_HAS_MADDWD || return (Int(typemax(Int8)), Int32)
    a = Int(typemax(Int16)) ÷ (2 * Int(max_meas))
    a >= 1 ? (min(a, Int(typemax(Int8))), Int16) : (Int(typemax(Int8)), Int32)
end
const _INT16_AMP, _INT16_WIPE_TI = _int16_choose_carrier(_INT16_MAX_MEAS)

# Default strip-mine block length (samples); multiple of every backend `W`,
# sized so the per-block L1 scratch (code + sin + cos) stays in L1.
const _INT16_BLK = 8192

# ── Widening / vpmaddwd helpers (ported from the GNSSSignals benchmark) ───────
@inline _wide32(v::SIMD.Vec{W,Int8}) where {W} = convert(SIMD.Vec{W,Int32}, v)
@inline _wide32(v::SIMD.Vec{W,Int16}) where {W} = convert(SIMD.Vec{W,Int32}, v)
@inline _wide16(v::SIMD.Vec{W,Int8}) where {W} = convert(SIMD.Vec{W,Int16}, v)

# vpmaddwd: Int16×Int16 → Int32 pairwise-add — the dot-product primitive for the
# code accumulate Σ codeₓ·DI. Native 512-/256-bit intrinsics tiled to width W;
# x86-only (gated on `_INT16_HAS_MADDWD`). Output is Vec{W÷2,Int32} (adjacent
# samples pre-summed); the final `sum` is bit-exact with a per-lane Int32 reduce.
@static if Sys.ARCH in (:x86_64, :i686)
    @inline _madd_tile(a::SIMD.Vec{M,Int16}, ::Val{o}, ::Val{t}) where {M,o,t} =
        shufflevector(a, Val(ntuple(i -> i - 1 + o, Val(t))))
    @inline function _madd512(a::SIMD.Vec{32,Int16}, b::SIMD.Vec{32,Int16})
        SIMD.Vec(Base.llvmcall(("""
            declare <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16>, <32 x i16>)
            define <16 x i32> @entry(<32 x i16> %a, <32 x i16> %b) #0 {
              %r = call <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16> %a, <32 x i16> %b)
              ret <16 x i32> %r }
            attributes #0 = { alwaysinline }""", "entry"), NTuple{16,Base.VecElement{Int32}},
            Tuple{NTuple{32,Base.VecElement{Int16}},NTuple{32,Base.VecElement{Int16}}}, a.data, b.data))
    end
    @inline function _madd256(a::SIMD.Vec{16,Int16}, b::SIMD.Vec{16,Int16})
        SIMD.Vec(Base.llvmcall(("""
            declare <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16>, <16 x i16>)
            define <8 x i32> @entry(<16 x i16> %a, <16 x i16> %b) #0 {
              %r = call <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16> %a, <16 x i16> %b)
              ret <8 x i32> %r }
            attributes #0 = { alwaysinline }""", "entry"), NTuple{8,Base.VecElement{Int32}},
            Tuple{NTuple{16,Base.VecElement{Int16}},NTuple{16,Base.VecElement{Int16}}}, a.data, b.data))
    end
    @inline function _maddacc(a::SIMD.Vec{64,Int16}, b::SIMD.Vec{64,Int16})   # AVX-512: 2×512-bit tiles
        lo = _madd512(_madd_tile(a, Val(0), Val(32)), _madd_tile(b, Val(0), Val(32)))
        hi = _madd512(_madd_tile(a, Val(32), Val(32)), _madd_tile(b, Val(32), Val(32)))
        shufflevector(lo, hi, Val(ntuple(i -> i - 1, Val(32))))
    end
    @inline function _maddacc(a::SIMD.Vec{32,Int16}, b::SIMD.Vec{32,Int16})   # AVX2: 2×256-bit tiles
        lo = _madd256(_madd_tile(a, Val(0), Val(16)), _madd_tile(b, Val(0), Val(16)))
        hi = _madd256(_madd_tile(a, Val(16), Val(16)), _madd_tile(b, Val(16), Val(16)))
        shufflevector(lo, hi, Val(ntuple(i -> i - 1, Val(16))))
    end
end

# Per-(thread) scratch: the strip-mine block code buffer (`extb`, Int8, sized
# `blk + tap-span`) plus the carrier sin/cos blocks (`csb`/`ccb`, Int8). Grown
# lazily and reused, so a hoisted backend is allocation-free in steady state.
mutable struct Int16ScratchBuffers
    extb::Vector{Int8}
    csb::Vector{Int8}
    ccb::Vector{Int8}
    # Shared DI/DQ tile (the carrier-wiped measurement) for the multi-signal-per-
    # sat tile-share path: filled once per block and reused across the sat's
    # signals (one downconvert per sat, not per signal). Wipe element type.
    dib::Vector{_INT16_WIPE_TI}
    dqb::Vector{_INT16_WIPE_TI}
end
Int16ScratchBuffers() =
    Int16ScratchBuffers(Int8[], Int8[], Int8[], _INT16_WIPE_TI[], _INT16_WIPE_TI[])

"""
$(SIGNATURES)

Integer (`Complex{Int16}`) hybrid-blocked CPU downconvert + correlate backend
(single-threaded). Opt-in alternative to [`CPUDownconvertAndCorrelator`] for
`Complex{Int16}` (12-bit ADC) sample buffers; errors on any other sample
element type. Construct **once outside** the `track!` loop and pass it via the
`downconvert_and_correlator` keyword for an allocation-free steady state.
"""
struct Int16DownconvertAndCorrelator{TBL<:SinCosTable} <: AbstractDownconvertAndCorrelator
    buffers::Int16ScratchBuffers
    table::TBL
    blk::Int
end

"""
$(SIGNATURES)

Multi-threaded integer (`Complex{Int16}`) hybrid-blocked backend. One
[`Int16ScratchBuffers`](@ref) per thread (indexed by `Threads.threadid()` inside
`@batch`); the carrier `table` is immutable and shared.
"""
struct Int16ThreadedDownconvertAndCorrelator{TBL<:SinCosTable} <:
       AbstractDownconvertAndCorrelator
    buffers::Vector{Int16ScratchBuffers}
    table::TBL
    blk::Int
end

function Int16DownconvertAndCorrelator(;
    steps::Integer = 64,
    amplitude::Integer = _INT16_AMP,
    blk::Integer = _INT16_BLK,
)
    Int16DownconvertAndCorrelator(
        Int16ScratchBuffers(),
        SinCosTable(Int8; steps, amplitude),
        Int(blk),
    )
end

function Int16ThreadedDownconvertAndCorrelator(;
    steps::Integer = 64,
    amplitude::Integer = _INT16_AMP,
    blk::Integer = _INT16_BLK,
)
    Int16ThreadedDownconvertAndCorrelator(
        [Int16ScratchBuffers() for _ = 1:Threads.maxthreadid()],
        SinCosTable(Int8; steps, amplitude),
        Int(blk),
    )
end

const _Int16DC =
    Union{Int16DownconvertAndCorrelator,Int16ThreadedDownconvertAndCorrelator}

@inline _scratch_buffers(dc::Int16DownconvertAndCorrelator) = dc.buffers
@inline _scratch_buffers(dc::Int16ThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

# Type-stable `NumAnts{M}` from a correlator (M from its type parameter).
@inline _num_ants_val(::AbstractCorrelator{M}) where {M} = NumAnts{M}()

# Fill `len` carrier samples (sin→csb, cos→ccb) starting at absolute sample
# `start`, 4-way unrolled so the permute lookups pipeline (the value engine is
# latency-bound single-stream). Ported from the benchmark's `_epl_fill_carrier!`.
@inline function _int16_fill_carrier!(
    csb::Vector{Int8},
    ccb::Vector{Int8},
    reng,
    start::Int,
    len::Int,
    phase,
    ::Val{W},
) where {W}
    s0 = carrier_state(reng, start; phase)
    s1 = carrier_state(reng, start + W; phase)
    s2 = carrier_state(reng, start + 2W; phase)
    s3 = carrier_state(reng, start + 3W; phase)
    n = 0
    @inbounds while n + 4W <= len
        a0, b0 = carrier_lookup(reng, s0)
        a1, b1 = carrier_lookup(reng, s1)
        a2, b2 = carrier_lookup(reng, s2)
        a3, b3 = carrier_lookup(reng, s3)
        vstore(a0, csb, n + 1);      vstore(b0, ccb, n + 1)
        vstore(a1, csb, n + W + 1);  vstore(b1, ccb, n + W + 1)
        vstore(a2, csb, n + 2W + 1); vstore(b2, ccb, n + 2W + 1)
        vstore(a3, csb, n + 3W + 1); vstore(b3, ccb, n + 3W + 1)
        s0 = carrier_advance(reng, s0, 4)
        s1 = carrier_advance(reng, s1, 4)
        s2 = carrier_advance(reng, s2, 4)
        s3 = carrier_advance(reng, s3, 4)
        n += 4W
    end
    @inbounds while n < len                     # < 4W tail, one chunk at a time
        a0, b0 = carrier_lookup(reng, s0)
        vstore(a0, csb, n + 1)
        vstore(b0, ccb, n + 1)
        s0 = carrier_advance(reng, s0, 1)
        n += W
    end
    nothing
end

# ── The integer hybrid-blocked kernel ────────────────────────────────────────
# Returns this integration's correlation contribution: `SVector{NC,ComplexF64}`
# (M=1) or `SVector{NC,SVector{M,ComplexF64}}` (M>1), one (multi-antenna) complex
# sum per tap, to be added to the correlator's running accumulators. `@generated`
# over (NC = length(sample_shifts), M = antenna count): NC tap accumulators per
# antenna live in named locals, and the M antenna-outer correlate passes are
# emitted explicitly (each keeps only its own NC accumulators live). Path (Int16
# vpmaddwd vs exact Int32) is chosen at generation time from host-derived consts.
@generated function _int16_hybrid_blocked!(
    dc::_Int16DC,
    signal::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    signal_type,
    prn::Integer,
    sample_shifts::SVector{NC},
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample::Integer,
    num_samples::Integer,
) where {M,NC}
    W = _INT16_W
    TI = _INT16_WIPE_TI
    use_madd = _INT16_HAS_MADDWD && TI === Int16
    AW = use_madd ? W ÷ 2 : W          # vpmaddwd pre-sums adjacent pairs → W÷2
    MT = TI === Int16 ? Int16 : Int32  # meas/wipe SIMD element type
    cc_load = TI === Int16 ? :(_wide16(vload(SIMD.Vec{$W,Int8}, ccb, n))) :
              :(_wide32(vload(SIMD.Vec{$W,Int8}, ccb, n)))
    sn_load = TI === Int16 ? :(_wide16(vload(SIMD.Vec{$W,Int8}, csb, n))) :
              :(_wide32(vload(SIMD.Vec{$W,Int8}, csb, n)))

    init = Expr(:block)
    for k = 1:NC
        push!(init.args, :($(Symbol("off_$k")) = sample_shifts[$k] - min_shift))
    end
    for j = 1:M, k = 1:NC
        push!(init.args, :($(Symbol("tI_$(j)_$k")) = zero(Int64)))
        push!(init.args, :($(Symbol("tQ_$(j)_$k")) = zero(Int64)))
    end

    # One antenna-outer correlate pass over the current block (emitted per j).
    function antenna_pass(j)
        blk_init = Expr(:block); flush = Expr(:block); chunk = Expr(:block); scalar = Expr(:block)
        for k = 1:NC
            aI = Symbol("aI_$(j)_$k"); aQ = Symbol("aQ_$(j)_$k")
            tI = Symbol("tI_$(j)_$k"); tQ = Symbol("tQ_$(j)_$k")
            offk = Symbol("off_$k")
            push!(blk_init.args, :($aI = zero(SIMD.Vec{$AW,Int32})))
            push!(blk_init.args, :($aQ = zero(SIMD.Vec{$AW,Int32})))
            push!(flush.args, :($tI += Int64(sum($aI))))
            push!(flush.args, :($tQ += Int64(sum($aQ))))
            if use_madd
                push!(chunk.args, quote
                    codew = _wide16(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                    $aI += _maddacc(codew, DI)
                    $aQ += _maddacc(codew, DQ)
                end)
            else
                push!(chunk.args, quote
                    codew = _wide32(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                    $aI += codew * DI
                    $aQ += codew * DQ
                end)
            end
            push!(scalar.args, quote
                c = Int32(extb[n+$offk])
                $tI += Int64(c * di)
                $tQ += Int64(c * dq)
            end)
        end
        quote
            $blk_init
            colbase = $(j - 1) * num_rows         # 0-based column (antenna) offset, in samples
            n = 1
            while n + W - 1 <= len
                byte_off = (colbase + base + n - 2) * 2 * sizeof(Int16)
                mr, mi = _deinterleave_load(SIMD.Vec{W,$MT}, p_sig, byte_off)
                cc = $cc_load
                sn = $sn_load
                DI = mr * cc + mi * sn
                DQ = mi * cc - mr * sn
                $chunk
                n += W
            end
            $flush
            while n <= len
                sig = signal[base+n-1, $j]
                mr_s = Int32(real(sig)); mi_s = Int32(imag(sig))
                cc_s = Int32(ccb[n]); sn_s = Int32(csb[n])
                di = mr_s * cc_s + mi_s * sn_s
                dq = mi_s * cc_s - mr_s * sn_s
                $scalar
                n += 1
            end
        end
    end
    correlate_passes = Expr(:block)
    for j = 1:M
        push!(correlate_passes.args, antenna_pass(j))
    end

    # Finalize: one accumulator per tap — ComplexF64 (M=1) or SVector{M} (M>1).
    function tap_expr(k)
        if M == 1
            :(complex(Float64($(Symbol("tI_1_$k"))), Float64($(Symbol("tQ_1_$k")))))
        else
            ant = [:(complex(Float64($(Symbol("tI_$(j)_$k"))), Float64($(Symbol("tQ_$(j)_$k"))))) for j = 1:M]
            :(SVector{$M,ComplexF64}(tuple($(ant...))))
        end
    end
    taps = [tap_expr(k) for k = 1:NC]
    result = :(SVector{$NC}(tuple($(taps...))))

    quote
        W = $W
        min_shift = minimum(sample_shifts)
        max_shift = maximum(sample_shifts)
        span = max_shift - min_shift
        num_rows = size(signal, 1)

        bufs = _scratch_buffers(dc)
        blk = dc.blk
        ncar = (cld(blk, W) + 4) * W            # room for the 4-way carrier fill tail
        length(bufs.extb) < blk + span && resize!(bufs.extb, blk + span)
        length(bufs.csb) < ncar && resize!(bufs.csb, ncar)
        length(bufs.ccb) < ncar && resize!(bufs.ccb, ncar)
        extb = bufs.extb
        csb = bufs.csb
        ccb = bufs.ccb

        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        reng = carrier_engine(dc.table, carrier_freq / sampling_freq)
        phase0 = Float64(carrier_phase)

        # Continuing code-fill engine: first emitted sample is output sample
        # `min_shift`, so the earliest tap reads real code (no zero edge), and
        # tap k at output n reads `extb[n + shift_k - min_shift]`.
        ceng = code_engine(
            signal_type, prn, sampling_frequency, code_frequency;
            start_phase = Float64(code_phase), start_index_shift = min_shift,
        )
        cst = code_state(ceng)

        $init
        p_sig = Ptr{Int16}(pointer(signal))

        blk_off = 0
        prevlen = 0
        @inbounds while blk_off < num_samples
            len = min(blk, num_samples - blk_off)

            # 1. Fill this block's code into extb[1 .. len+span].
            if prevlen == 0
                cst = gen_code!(view(extb, 1:(len+span)), ceng, cst)
            else
                for i = 1:span
                    extb[i] = extb[prevlen+i]
                end
                cst = gen_code!(view(extb, (span+1):(span+len)), ceng, cst)
            end

            # 2. Fill this block's carrier sin/cos (shared across antennas).
            _int16_fill_carrier!(csb, ccb, reng, blk_off, len, phase0, Val(W))

            # 3. M antenna-outer correlate passes (wipe + per-tap accumulate),
            #    flushed per block. `base` = 1-based signal row of local n=0.
            base = signal_start_sample + blk_off
            $correlate_passes

            blk_off += len
            prevlen = len
        end

        $result
    end
end

# Dynamic (runtime tap count) fallback: correlators whose sample shifts are a
# runtime-sized `AbstractVector` (e.g. a `Vector`-accumulator correlator —
# issue #126 (b)). Mirrors the Float32 backend's `AbstractVector`-shifts path.
# Loops over taps/antennas at runtime and uses the exact Int32 wipe + per-chunk
# horizontal sum, so it is not the hot path but stays correct for any tap count;
# returns `Vector{ComplexF64}` (M=1) or `Vector{SVector{M,ComplexF64}}` (M>1).
# `SVector` shifts are more specific and dispatch to the `@generated` method
# above, so EPL/VEPL keep the fast unrolled kernel.
function _int16_hybrid_blocked!(
    dc::_Int16DC,
    signal::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    signal_type,
    prn::Integer,
    sample_shifts::AbstractVector,
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample::Integer,
    num_samples::Integer,
) where {M}
    W = _INT16_W
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    max_shift = maximum(sample_shifts)
    span = max_shift - min_shift
    offs = [Int(sample_shifts[k]) - min_shift for k = 1:NC]   # extb read offsets
    num_rows = size(signal, 1)

    bufs = _scratch_buffers(dc)
    blk = dc.blk
    ncar = (cld(blk, W) + 4) * W
    length(bufs.extb) < blk + span && resize!(bufs.extb, blk + span)
    length(bufs.csb) < ncar && resize!(bufs.csb, ncar)
    length(bufs.ccb) < ncar && resize!(bufs.ccb, ncar)
    extb = bufs.extb
    csb = bufs.csb
    ccb = bufs.ccb

    carrier_freq = Float64(upreferred(carrier_frequency / Hz))
    sampling_freq = Float64(upreferred(sampling_frequency / Hz))
    reng = carrier_engine(dc.table, carrier_freq / sampling_freq)
    phase0 = Float64(carrier_phase)

    ceng = code_engine(
        signal_type, prn, sampling_frequency, code_frequency;
        start_phase = Float64(code_phase), start_index_shift = min_shift,
    )
    cst = code_state(ceng)

    tI = zeros(Int64, M, NC)
    tQ = zeros(Int64, M, NC)
    p_sig = Ptr{Int16}(pointer(signal))

    blk_off = 0
    prevlen = 0
    @inbounds while blk_off < num_samples
        len = min(blk, num_samples - blk_off)
        if prevlen == 0
            cst = gen_code!(view(extb, 1:(len+span)), ceng, cst)
        else
            for i = 1:span
                extb[i] = extb[prevlen+i]
            end
            cst = gen_code!(view(extb, (span+1):(span+len)), ceng, cst)
        end
        _int16_fill_carrier!(csb, ccb, reng, blk_off, len, phase0, Val(W))
        base = signal_start_sample + blk_off
        for j = 1:M
            colbase = (j - 1) * num_rows
            n = 1
            while n + W - 1 <= len
                byte_off = (colbase + base + n - 2) * 2 * sizeof(Int16)
                mr, mi = _deinterleave_load(SIMD.Vec{W,Int32}, p_sig, byte_off)
                cc = _wide32(vload(SIMD.Vec{W,Int8}, ccb, n))
                sn = _wide32(vload(SIMD.Vec{W,Int8}, csb, n))
                DI = mr * cc + mi * sn
                DQ = mi * cc - mr * sn
                for k = 1:NC
                    codev = _wide32(vload(SIMD.Vec{W,Int8}, extb, n + offs[k]))
                    tI[j, k] += Int64(sum(codev * DI))
                    tQ[j, k] += Int64(sum(codev * DQ))
                end
                n += W
            end
            while n <= len
                sig = signal[base+n-1, j]
                mr_s = Int32(real(sig)); mi_s = Int32(imag(sig))
                cc_s = Int32(ccb[n]); sn_s = Int32(csb[n])
                di = mr_s * cc_s + mi_s * sn_s
                dq = mi_s * cc_s - mr_s * sn_s
                for k = 1:NC
                    c = Int32(extb[n+offs[k]])
                    tI[j, k] += Int64(c * di)
                    tQ[j, k] += Int64(c * dq)
                end
                n += 1
            end
        end
        blk_off += len
        prevlen = len
    end

    if M == 1
        return [complex(Float64(tI[1, k]), Float64(tQ[1, k])) for k = 1:NC]
    else
        return [
            SVector{M,ComplexF64}(ntuple(j -> complex(Float64(tI[j, k]), Float64(tQ[j, k])), M))
            for k = 1:NC
        ]
    end
end

# ── Multi-signal-per-sat tile-share ───────────────────────────────────────────
# A satellite carrying several signals on one carrier (e.g. GPS L1 C/A + L1C-D +
# L1C-P) shares the carrier and therefore the per-sample carrier wipe-off. So
# per strip-mine block we fill the carrier once and materialise the shared DI/DQ
# tile (the carrier-wiped measurement) once — ONE downconvert per sat, not per
# signal — then correlate each signal's own code against that tile. Mirrors the
# Float32 `downconvert_and_correlate_fused_tuple!`, adapted to the hybrid-blocked
# strip-mine + integer pipeline.

# Fill the shared DI/DQ tile for the current block: per antenna `j`, the slice
# `dib/dqb[(j-1)*len + n]` holds the carrier-wiped measurement at sample n.
# `@generated` so the M antenna slices unroll.
@generated function _int16_fill_ditile!(
    dib, dqb, signal, ::NumAnts{M}, csb, ccb, base, len, num_rows, p_sig,
) where {M}
    W = _INT16_W
    TI = _INT16_WIPE_TI
    widen = TI === Int16 ? :_wide16 : :_wide32
    passes = Expr(:block)
    for j = 1:M
        push!(passes.args, quote
            colbase = $(j - 1) * num_rows
            toff = $(j - 1) * len
            n = 1
            while n + $W - 1 <= len
                byte_off = (colbase + base + n - 2) * 2 * sizeof(Int16)
                mr, mi = _deinterleave_load(SIMD.Vec{$W,$TI}, p_sig, byte_off)
                cc = $widen(vload(SIMD.Vec{$W,Int8}, ccb, n))
                sn = $widen(vload(SIMD.Vec{$W,Int8}, csb, n))
                vstore(mr * cc + mi * sn, dib, toff + n)
                vstore(mi * cc - mr * sn, dqb, toff + n)
                n += $W
            end
            while n <= len
                sig = signal[base+n-1, $j]
                mr_s = $TI(real(sig)); mi_s = $TI(imag(sig))
                cc_s = $TI(ccb[n]); sn_s = $TI(csb[n])
                dib[toff+n] = mr_s * cc_s + mi_s * sn_s
                dqb[toff+n] = mi_s * cc_s - mr_s * sn_s
                n += 1
            end
        end)
    end
    quote
        @inbounds begin
            $passes
        end
        nothing
    end
end

# Tile-share kernel for N signals sharing a carrier. Returns a tuple of N
# per-signal accumulators (`SVector{NCᵢ,ComplexF64}` or `…{SVector{M,…}}`).
# `@generated` over (M, the signals tuple) so each signal's tap count NCᵢ and the
# M antenna passes unroll. Per block: fill carrier + DI/DQ tile once, then for
# each signal fill its code (one-shot `gen_code!` at the block's analytically
# advanced phase, into the reused `extb`) and accumulate its taps from the tile.
@generated function _int16_hybrid_blocked_multi!(
    dc::_Int16DC,
    signal::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    signal_types::Tuple,
    prn::Integer,
    all_shifts::Tuple{Vararg{SVector}},
    code_phases::Tuple,
    code_freqs::Tuple,
    carrier_phase,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample::Integer,
    num_samples::Integer,
) where {M}
    W = _INT16_W
    TI = _INT16_WIPE_TI
    use_madd = _INT16_HAS_MADDWD && TI === Int16
    AW = use_madd ? W ÷ 2 : W
    N = length(all_shifts.parameters)
    NCs = [length(all_shifts.parameters[i]) for i = 1:N]

    setup = Expr(:block)
    maxspan = Expr(:call, :max)
    for i = 1:N
        push!(setup.args, :($(Symbol("sh_$i")) = all_shifts[$i]))
        push!(setup.args, :($(Symbol("mins_$i")) = minimum($(Symbol("sh_$i")))))
        push!(setup.args,
            :($(Symbol("span_$i")) = maximum($(Symbol("sh_$i"))) - $(Symbol("mins_$i"))))
        push!(setup.args,
            :($(Symbol("cps_$i")) = Float64(upreferred(code_freqs[$i] / Hz)) / sampling_freq))
        for k = 1:NCs[i]
            push!(setup.args,
                :($(Symbol("off_$(i)_$k")) = $(Symbol("sh_$i"))[$k] - $(Symbol("mins_$i"))))
        end
        push!(maxspan.args, Symbol("span_$i"))
        for j = 1:M, k = 1:NCs[i]
            push!(setup.args, :($(Symbol("tI_$(i)_$(j)_$k")) = zero(Int64)))
            push!(setup.args, :($(Symbol("tQ_$(i)_$(j)_$k")) = zero(Int64)))
        end
    end

    # Per-signal: one-shot code fill into `extb`, then M antenna passes that
    # accumulate the signal's taps from the shared DI/DQ tile.
    function signal_corr(i)
        b = Expr(:block)
        push!(b.args, :(blk_phase = code_phases[$i] + $(Symbol("cps_$i")) * blk_off))
        push!(b.args, :(gen_code!(
            view(extb, 1:(len+$(Symbol("span_$i")))),
            signal_types[$i], prn, sampling_frequency, code_freqs[$i],
            blk_phase, $(Symbol("mins_$i")),
        )))
        for j = 1:M
            blk_init = Expr(:block); flush = Expr(:block); chunk = Expr(:block); scalar = Expr(:block)
            for k = 1:NCs[i]
                aI = Symbol("aI_$k"); aQ = Symbol("aQ_$k")
                tI = Symbol("tI_$(i)_$(j)_$k"); tQ = Symbol("tQ_$(i)_$(j)_$k")
                offk = Symbol("off_$(i)_$k")
                push!(blk_init.args, :($aI = zero(SIMD.Vec{$AW,Int32})))
                push!(blk_init.args, :($aQ = zero(SIMD.Vec{$AW,Int32})))
                push!(flush.args, :($tI += Int64(sum($aI))))
                push!(flush.args, :($tQ += Int64(sum($aQ))))
                if use_madd
                    push!(chunk.args, quote
                        codew = _wide16(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                        $aI += _maddacc(codew, DI)
                        $aQ += _maddacc(codew, DQ)
                    end)
                else
                    push!(chunk.args, quote
                        codew = _wide32(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                        $aI += codew * DI
                        $aQ += codew * DQ
                    end)
                end
                push!(scalar.args, quote
                    c = Int32(extb[n+$offk])
                    $tI += Int64(c * Int32(dib[toff+n]))
                    $tQ += Int64(c * Int32(dqb[toff+n]))
                end)
            end
            push!(b.args, quote
                $blk_init
                toff = $(j - 1) * len
                n = 1
                while n + $W - 1 <= len
                    DI = vload(SIMD.Vec{$W,$TI}, dib, toff + n)
                    DQ = vload(SIMD.Vec{$W,$TI}, dqb, toff + n)
                    $chunk
                    n += $W
                end
                $flush
                while n <= len
                    $scalar
                    n += 1
                end
            end)
        end
        b
    end
    sigs = Expr(:block)
    for i = 1:N
        push!(sigs.args, signal_corr(i))
    end

    function tap(i, k)
        if M == 1
            :(complex(Float64($(Symbol("tI_$(i)_1_$k"))), Float64($(Symbol("tQ_$(i)_1_$k")))))
        else
            ant = [:(complex(Float64($(Symbol("tI_$(i)_$(j)_$k"))), Float64($(Symbol("tQ_$(i)_$(j)_$k"))))) for j = 1:M]
            :(SVector{$M,ComplexF64}(tuple($(ant...))))
        end
    end
    finals = Expr(:tuple)
    for i = 1:N
        taps = [tap(i, k) for k = 1:NCs[i]]
        push!(finals.args, :(SVector{$(NCs[i])}(tuple($(taps...)))))
    end

    quote
        W = $W
        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        reng = carrier_engine(dc.table, carrier_freq / sampling_freq)
        phase0 = Float64(carrier_phase)
        num_rows = size(signal, 1)
        $setup
        maxspan_v = $maxspan

        bufs = _scratch_buffers(dc)
        blk = dc.blk
        ncar = (cld(blk, W) + 4) * W
        length(bufs.extb) < blk + maxspan_v && resize!(bufs.extb, blk + maxspan_v)
        length(bufs.csb) < ncar && resize!(bufs.csb, ncar)
        length(bufs.ccb) < ncar && resize!(bufs.ccb, ncar)
        ntile = blk * $M
        length(bufs.dib) < ntile && resize!(bufs.dib, ntile)
        length(bufs.dqb) < ntile && resize!(bufs.dqb, ntile)
        extb = bufs.extb; csb = bufs.csb; ccb = bufs.ccb; dib = bufs.dib; dqb = bufs.dqb

        p_sig = Ptr{Int16}(pointer(signal))
        blk_off = 0
        @inbounds while blk_off < num_samples
            len = min(blk, num_samples - blk_off)
            _int16_fill_carrier!(csb, ccb, reng, blk_off, len, phase0, Val(W))
            base = signal_start_sample + blk_off
            _int16_fill_ditile!(dib, dqb, signal, NumAnts{$M}(), csb, ccb, base, len, num_rows, p_sig)
            $sigs
            blk_off += len
        end
        $finals
    end
end

# Multi-signal-per-sat correlate: shares one carrier downconvert across the sat's
# signals via the tile-share kernel. Returns the per-signal
# `(new_correlator, is_integration_completed)` tuples.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal,TrackedSignal,Vararg{TrackedSignal}},
    per_signal_completed::Tuple,
    dc::_Int16DC,
    signal,
    code_doppler,
    code_phase,
    carrier_frequency,
    carrier_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_samples_signal,
)
    params = map(signals) do head
        _signal_replica_params(head, code_doppler, code_phase, sampling_frequency, num_samples_signal)
    end
    correlators = map(s -> s.correlator, signals)
    signal_types = map(s -> s.signal, signals)
    all_shifts = map(p -> p.sample_shifts, params)
    code_phases = map(p -> p.signal_code_phase, params)
    code_freqs = map(p -> p.code_frequency, params)
    new_accs = _int16_hybrid_blocked_multi!(
        dc,
        signal,
        _num_ants_val(correlators[1]),
        signal_types,
        prn,
        all_shifts,
        code_phases,
        code_freqs,
        carrier_phase,
        carrier_frequency,
        sampling_frequency,
        signal_start_sample,
        samples_to_integrate,
    )
    new_corrs = map(correlators, new_accs) do c, a
        update_accumulator(c, get_accumulators(c) .+ a)
    end
    map(tuple, new_corrs, per_signal_completed)
end

# Per-signal correlate for the integer backend: single signal, static tap count,
# any antenna count M. Returns the one-tuple
# `((new_correlator, is_integration_completed),)` matching the Float32 backend's
# `_correlate_signals` contract.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal},
    per_signal_completed::Tuple{Bool},
    dc::_Int16DC,
    signal,
    code_doppler,
    code_phase,
    carrier_frequency,
    carrier_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_samples_signal,
)
    head = signals[1]
    s = head.signal
    correlator = head.correlator
    p = _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    )
    new_acc = _int16_hybrid_blocked!(
        dc,
        signal,
        _num_ants_val(correlator),
        s,
        prn,
        p.sample_shifts,
        p.signal_code_phase,
        carrier_phase,
        p.code_frequency,
        carrier_frequency,
        sampling_frequency,
        signal_start_sample,
        samples_to_integrate,
    )
    prev = get_accumulators(correlator)
    new_corr = update_accumulator(correlator, prev .+ new_acc)
    ((new_corr, per_signal_completed[1]),)
end

# Per-sat downconvert+correlate for the integer backend. Mirrors the Float32
# backend, reusing the shared boundary / completion / update helpers; only
# `_correlate_signals` (the kernel) differs.
function _update_tracked_sat_correlator(
    sat::TrackedSat,
    dc::_Int16DC,
    signal,
    num_samples_signal,
    sampling_frequency,
    intermediate_frequency,
)
    samples_to_integrate, per_signal_completed = _calc_min_samples_and_completed(
        sat.signals,
        sat.signal_start_sample,
        sampling_frequency,
        sat.code_doppler,
        sat.code_phase,
        num_samples_signal,
    )
    if samples_to_integrate == 0
        return sat
    end
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    new_signals_data = _correlate_signals(
        sat.signals,
        per_signal_completed,
        dc,
        signal,
        sat.code_doppler,
        sat.code_phase,
        carrier_frequency,
        sat.carrier_phase,
        sampling_frequency,
        sat.signal_start_sample,
        samples_to_integrate,
        sat.prn,
        num_samples_signal,
    )
    update(
        sat,
        samples_to_integrate,
        intermediate_frequency,
        sampling_frequency,
        new_signals_data,
    )
end

# ── Group/measurement plumbing (parallels the Float32 backend) ────────────────
@inline function _dc_one_group!(
    g::SignalGroup,
    dc::_Int16DC,
    measurements::BandMeasurements,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    m = measurements[band_key(g.band)]
    eltype(m.samples) === Complex{Int16} || throw(
        ArgumentError(
            string(
                "Int16DownconvertAndCorrelator requires `Complex{Int16}` measurement ",
                "samples (12-bit ADC); got element type ",
                eltype(m.samples),
                ". Use a CPU(Threaded)DownconvertAndCorrelator for floating-point samples.",
            ),
        ),
    )
    _dc_group_loop!(
        dc,
        vals,
        m.samples,
        get_num_samples(m),
        m.sampling_frequency,
        m.intermediate_frequency,
    )
end

@inline function _dc_group_loop!(dc::Int16DownconvertAndCorrelator, vals, args::Vararg{Any,4})
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

@inline function _dc_group_loop!(
    dc::Int16ThreadedDownconvertAndCorrelator,
    vals,
    args::Vararg{Any,4},
)
    @batch for i = 1:length(vals)
        @inbounds vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

"""
$(SIGNATURES)

Downconvert and correlate all satellites with the integer (`Complex{Int16}`)
backend. Detaches the per-sat slot values from the input (shares the key set),
then writes new values in place — see the Float32 backend's `downconvert_and_correlate`.
"""
function downconvert_and_correlate(
    dc::_Int16DC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    downconvert_and_correlate!(dc, measurements, new_track_state)
end

"""
$(SIGNATURES)

In-place integer (`Complex{Int16}`) downconvert and correlate. Returns the same
`track_state`.
"""
function downconvert_and_correlate!(
    dc::_Int16DC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    _foreach_group!(_dc_one_group!, track_state.groups, dc, measurements)
    return track_state
end
