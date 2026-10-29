# Integer (`Complex{Int16}`, 12-bit ADC) hybrid-blocked downconvert + correlate backend,
# the integer counterpart of the Float32 fused path. The integration is strip-mined into
# `blk`-sample blocks whose code and carrier are regenerated into L1-resident scratch:
#
#   * code — Int8 ±1 (CBOC ±13/±25) replica, one-shot `gen_code!` per block (exact
#     across block seams, non-baked secondaries such as GPS L5I NH10 applied).
#   * carrier — Int8 sin/cos from SinCosLUT's table lookup (`_int16_fill_carrier!`).
#   * wipe + correlate — DI = mᵣ·cos + mᵢ·sin, DQ = mᵢ·cos − mᵣ·sin, in Int16 where it
#     cannot overflow (`_int16_choose_carrier`) and exact Int32 otherwise; Σ codeₓ·DI
#     accumulates in Int32 lanes flushed into Int64 totals every block.
#
# With M > 1 antennas (a dense `Matrix`, one column per antenna) one code + carrier block
# serves M antenna-outer correlate passes; a satellite's signals on one carrier also share
# the wipe (tile-share, below). Results are `ComplexF64` (M=1) / `SVector{M,ComplexF64}`
# (M>1) with the carrier amplitude divided out. The kernels are `@generated` over the tap
# count NC and M, with a runtime-tap-count fallback.

import SinCosLUT
using SinCosLUT: SinCosTable, carrier_engine, carrier_state, carrier_lookup, carrier_advance

# SIMD width of the Int8 carrier and code LUTs on this host.
const _INT16_W = let be = SinCosLUT.default_backend(Int8, 64)
    be isa SinCosLUT.AVX512 ? 64 : be isa SinCosLUT.AVX2 ? 32 : be isa SinCosLUT.Neon ? 16 : 1
end

# The fused widening `vpmaddwd` accumulate (Int16×Int16→Int32 pairwise) is x86-only.
const _INT16_HAS_MADDWD = Sys.ARCH in (:x86_64, :i686) && _INT16_W in (32, 64)

# Whether to run the carrier wipe in Int16 rather than exact Int32: half the SIMD ops and
# DI/DQ traffic, and the code accumulate becomes a widening MAC — `vpmaddwd` on x86, NEON
# `SMLAL`/`SMLAL2` on aarch64 (LLVM selects it from a `sext(i16)·sext(i16)` accumulate).
# Elsewhere there is no widening MAC to gain, so the wipe stays Int32.
const _INT16_WIDE_WIPE = _INT16_HAS_MADDWD || Sys.ARCH === :aarch64

# Choose the Int8 carrier amplitude and wipe type from `max_meas`, the declared largest
# `|real|`/`|imag|` of a sample (e.g. 2^11 for a 12-bit ADC): the largest `amp` with
# `2·max_meas·amp ≤ typemax(Int16)`, capped at `typemax(Int8)`. The factor 2 is because
# the complex wipe sums two products, `|DI| ≤ |mᵣ·cos| + |mᵢ·sin|`. The amplitude is the
# same on every backend, only the wipe type differs (see `_INT16_WIDE_WIPE`): a larger
# Int32-only amplitude overflowed the CBOC accumulate on the scalar path.
#
# Three guards keep the per-block Σ codeₓ·DI from wrapping:
#   * the cross-lane flush reduce is widened to Int64 (`_wide64`, #165);
#   * for `max_meas ≥ 2^14` no `amp ≥ 1` is Int16-safe, so the wipe falls back to amp 127
#     in Int32 and `_int16_safe_blk` shrinks the block (#167);
#   * on the SinCosLUT `Portable` backend (`_INT16_W == 1`) one lane sums a whole block,
#     so `_int16_flush_len` shrinks it (#166).
function _int16_choose_carrier(max_meas::Integer)
    a = Int(typemax(Int16)) ÷ (2 * Int(max_meas))
    a >= 1 || return (Int(typemax(Int8)), Int32)
    (min(a, Int(typemax(Int8))), _INT16_WIDE_WIPE ? Int16 : Int32)
end

# Largest code magnitude (Galileo E1B CBOC ±25), bounding the Int32 lane sums.
const _INT16_MAX_CODE = 25

# Block length keeping a lane's per-block Σ codeₓ·DI within Int32 on the `max_meas ≥ 2^14`
# fallback (amp 127, Int32 wipe, #167); `blk` unchanged on the Int16-safe path.
function _int16_safe_blk(blk::Integer, max_meas::Integer, amp::Integer)
    di_max = 2 * Int(max_meas) * Int(amp)                     # |DI| = |mᵣ·cos + mᵢ·sin|
    di_max <= Int(typemax(Int16)) && return Int(blk)          # Int16-safe path: unchanged
    safe = Int(typemax(Int32)) ÷ (_INT16_MAX_CODE * di_max)   # samples/block keeping Σ within Int32
    min(Int(blk), max(1, safe))
end

# Default strip-mine block length (samples): a multiple of every backend `W`, small enough
# for the code + sin + cos scratch to stay in L1.
const _INT16_BLK = 8192

# Block length keeping an Int32 lane from wrapping before its per-block Int64 flush. A lane
# sums `fld(L, W)` products (twice that with `vpmaddwd`'s pair pre-sum) of at most
# `_INT16_MAX_CODE·max_wipe`. For `W ≥ 16` the cap exceeds the default block; on the
# SinCosLUT `Portable` fallback (`W == 1`) it shrinks it (#166). Composed by `min` with
# `_int16_safe_blk` (#167).
function _int16_flush_len(W::Integer, max_wipe::Integer, blk::Integer)
    max_product = _INT16_MAX_CODE * Int(max_wipe)
    products_per_lane = Int(typemax(Int32)) ÷ (2 * max_product)   # 2× covers the vpmaddwd path
    min(Int(blk), max(Int(W), products_per_lane * Int(W)))
end

# `blk ≤ 0` would never advance the strip-mine loop and hang `track!` (#169), so reject it;
# an oversized `blk` is safe, as `_int16_flush_len` and `_int16_safe_blk` clamp it.
function _int16_validate_blk(blk::Integer, max_meas::Integer, amp::Integer)
    blk >= 1 || throw(
        ArgumentError(
            "Int16 backend: blk must be ≥ 1 (got $blk); blk ≤ 0 makes the " *
            "strip-mine loop never advance and hangs track! forever.",
        ),
    )
    _int16_safe_blk(blk, max_meas, amp)   # clamp on the amp = 127 fallback (#167)
end

# ── Widening / vpmaddwd helpers ───────────────────────────────────────────────
@inline _wide32(v::SIMD.Vec{W,Int8}) where {W} = convert(SIMD.Vec{W,Int32}, v)
@inline _wide32(v::SIMD.Vec{W,Int16}) where {W} = convert(SIMD.Vec{W,Int32}, v)
@inline _wide32(v::SIMD.Vec{W,Int32}) where {W} = v   # already widened (Int32 wipe tile)
@inline _wide16(v::SIMD.Vec{W,Int8}) where {W} = convert(SIMD.Vec{W,Int16}, v)
# Widen the Int32 lanes to Int64 before the flush's horizontal `sum`: each lane is safe,
# but the block total can exceed typemax(Int32) (~2.7× for a strong, code-aligned,
# full-scale CBOC capture at default settings, #165).
@inline _wide64(v::SIMD.Vec{W,Int32}) where {W} = convert(SIMD.Vec{W,Int64}, v)

# vpmaddwd (Int16×Int16 → Int32, adjacent pairs summed) for Σ codeₓ·DI, from the native
# 512-/256-bit intrinsics tiled to width W; returns Vec{W÷2,Int32}.
@static if Sys.ARCH in (:x86_64, :i686)
    @inline _madd_tile(a::SIMD.Vec{M,Int16}, ::Val{o}, ::Val{t}) where {M,o,t} =
        shufflevector(a, Val(ntuple(i -> i - 1 + o, Val(t))))
    @inline function _madd512(a::SIMD.Vec{32,Int16}, b::SIMD.Vec{32,Int16})
        SIMD.Vec(
            Base.llvmcall(
                (
                    """
declare <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16>, <32 x i16>)
define <16 x i32> @entry(<32 x i16> %a, <32 x i16> %b) #0 {
  %r = call <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16> %a, <32 x i16> %b)
  ret <16 x i32> %r }
attributes #0 = { alwaysinline }""",
                    "entry",
                ),
                NTuple{16,Base.VecElement{Int32}},
                Tuple{NTuple{32,Base.VecElement{Int16}},NTuple{32,Base.VecElement{Int16}}},
                a.data,
                b.data,
            ),
        )
    end
    @inline function _madd256(a::SIMD.Vec{16,Int16}, b::SIMD.Vec{16,Int16})
        SIMD.Vec(
            Base.llvmcall(
                (
                    """
declare <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16>, <16 x i16>)
define <8 x i32> @entry(<16 x i16> %a, <16 x i16> %b) #0 {
  %r = call <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16> %a, <16 x i16> %b)
  ret <8 x i32> %r }
attributes #0 = { alwaysinline }""",
                    "entry",
                ),
                NTuple{8,Base.VecElement{Int32}},
                Tuple{NTuple{16,Base.VecElement{Int16}},NTuple{16,Base.VecElement{Int16}}},
                a.data,
                b.data,
            ),
        )
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

# Per-thread scratch, grown lazily and reused (allocation-free in steady state): the block
# code `extb` (`blk + tap-span`) and carrier sin/cos `csb`/`ccb`, all Int8. `TI` is the
# wipe element type chosen by `_int16_choose_carrier`.
struct Int16ScratchBuffers{TI}
    extb::Vector{Int8}
    csb::Vector{Int8}
    ccb::Vector{Int8}
    # Shared DI/DQ tile (the carrier-wiped measurement) of the tile-share and dynamic paths.
    dib::Vector{TI}
    dqb::Vector{TI}
    # Int64 tap totals (`M·NC`) of the runtime-tap-count fallback.
    tsumI::Vector{Int64}
    tsumQ::Vector{Int64}
end
Int16ScratchBuffers{TI}() where {TI} =
    Int16ScratchBuffers{TI}(Int8[], Int8[], Int8[], TI[], TI[], Int64[], Int64[])

"""
$(SIGNATURES)

Integer (`Complex{Int16}`) hybrid-blocked CPU downconvert + correlate backend
(single-threaded). Opt-in alternative to [`CPUDownconvertAndCorrelator`] for
`Complex{Int16}` sample buffers; errors on any other sample element type.
Construct it **once outside** the `track!` loop and pass it via the
`downconvert_and_correlator` keyword: allocation-free in steady state for `SVector`
sample shifts (EPL/VEPL); runtime `AbstractVector` shifts allocate only the returned
`Vector`.

# Arguments

`max_meas` (required, no default) is the largest `|real|`/`|imag|` any sample will
take, i.e. the front end's full-scale (e.g. `2^11` for a 12-bit ADC). The carrier
amplitude is the largest with `2·max_meas·amplitude ≤ typemax(Int16)`, since the complex
wipe `DI = mᵣ·cos + mᵢ·sin` sums two products. **Under-declaring `max_meas` silently
overflows the `Int16` wipe and corrupts the correlation**; over-declaring is safe and
only coarsens the carrier quantisation.

!!! note "Performance: keep `max_meas < 2^14`"

    For `max_meas ≥ 2^14` no `Int16`-safe amplitude `≥ 1` exists: the backend falls
    back to an exact `Int32` wipe (off the `vpmaddwd` fast path) and a shorter
    strip-mine block. Correct, but measurably slower.

`blk` is the strip-mine block length in samples. It must be `≥ 1` (`blk ≤ 0` would
hang `track!`, #169); a larger-than-safe `blk` is clamped so the accumulators never wrap.
"""
struct Int16DownconvertAndCorrelator{TBL<:SinCosTable,TI} <:
       AbstractDownconvertAndCorrelator
    buffers::Int16ScratchBuffers{TI}
    table::TBL
    blk::Int
    # Peak of the Int8 carrier replica, divided out at finalize so the magnitudes match
    # the unit-carrier Float32 backend.
    carrier_amplitude::Int
end

"""
$(SIGNATURES)

Multi-threaded [`Int16DownconvertAndCorrelator`](@ref): one `Int16ScratchBuffers` per
thread, the immutable carrier `table` shared.
"""
struct Int16ThreadedDownconvertAndCorrelator{TBL<:SinCosTable,TI} <:
       AbstractDownconvertAndCorrelator
    buffers::Vector{Int16ScratchBuffers{TI}}
    table::TBL
    blk::Int
    # Peak of the Int8 carrier replica, divided out at finalize.
    carrier_amplitude::Int
end

# The kernels stride by the compile-time `_INT16_W`, but the carrier engine's width follows
# the backend SinCosLUT picks from `steps`; a mismatch silently corrupts the carrier (#168).
# `steps` values keeping the same width are allowed.
function _int16_assert_engine_width(table::SinCosTable)
    W = SinCosLUT.carrier_width(carrier_engine(table, 0))
    W == _INT16_W || throw(
        ArgumentError(
            string(
                "Int16DownconvertAndCorrelator: the carrier table's SIMD engine ",
                "width (",
                W,
                ") does not match the kernel stride `_INT16_W` (",
                _INT16_W,
                "). This happens when `steps` selects a different ",
                "SinCosLUT backend than the default `steps = 64` on this host; a ",
                "mismatched width would silently corrupt the carrier replica. ",
                "Use `steps = 64`.",
            ),
        ),
    )
    return nothing
end

function Int16DownconvertAndCorrelator(
    max_meas::Integer;
    steps::Integer = 64,
    blk::Integer = _INT16_BLK,
)
    amplitude, TI = _int16_choose_carrier(max_meas)
    table = SinCosTable(Int8; steps, amplitude)
    _int16_assert_engine_width(table)
    Int16DownconvertAndCorrelator{typeof(table),TI}(
        Int16ScratchBuffers{TI}(),
        table,
        _int16_validate_blk(blk, max_meas, amplitude),
        amplitude,
    )
end

function Int16ThreadedDownconvertAndCorrelator(
    max_meas::Integer;
    steps::Integer = 64,
    blk::Integer = _INT16_BLK,
)
    amplitude, TI = _int16_choose_carrier(max_meas)
    table = SinCosTable(Int8; steps, amplitude)
    _int16_assert_engine_width(table)
    Int16ThreadedDownconvertAndCorrelator{typeof(table),TI}(
        [Int16ScratchBuffers{TI}() for _ = 1:Threads.maxthreadid()],
        table,
        _int16_validate_blk(blk, max_meas, amplitude),
        amplitude,
    )
end

const _Int16DC = Union{Int16DownconvertAndCorrelator,Int16ThreadedDownconvertAndCorrelator}

@inline _scratch_buffers(dc::Int16DownconvertAndCorrelator) = dc.buffers
@inline _scratch_buffers(dc::Int16ThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

# Fill `len` carrier samples (sin→csb, cos→ccb) from absolute sample `start`, 4-way
# unrolled so the latency-bound lookups pipeline.
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
        vstore(a0, csb, n + 1);
        vstore(b0, ccb, n + 1)
        vstore(a1, csb, n + W + 1);
        vstore(b1, ccb, n + W + 1)
        vstore(a2, csb, n + 2W + 1);
        vstore(b2, ccb, n + 2W + 1)
        vstore(a3, csb, n + 3W + 1);
        vstore(b3, ccb, n + 3W + 1)
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
# Returns this integration's contribution per tap, `SVector{NC}` of `ComplexF64` (M=1) or
# `SVector{M,ComplexF64}` (M>1). `@generated` over (NC, M): each of the M antenna passes
# keeps only its own NC accumulators live; the accumulate (vpmaddwd / SMLAL / exact Int32)
# is chosen at generation time.
@generated function _int16_hybrid_blocked!(
    dc::_Int16DC,
    samples::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    signal,
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
    TI = dc.parameters[2]              # carrier-wipe element type, per backend instance
    use_madd = _INT16_HAS_MADDWD && TI === Int16
    # aarch64 Int16 wipe: widen to Int32 and MAC (→ SMLAL); no pair pre-sum, so `AW = W`.
    use_smlal = TI === Int16 && !use_madd
    AW = use_madd ? W ÷ 2 : W          # vpmaddwd pre-sums adjacent pairs → W÷2
    MT = TI === Int16 ? Int16 : Int32  # meas/wipe SIMD element type
    cc_load =
        TI === Int16 ? :(_wide16(vload(SIMD.Vec{$W,Int8}, ccb, n))) :
        :(_wide32(vload(SIMD.Vec{$W,Int8}, ccb, n)))
    sn_load =
        TI === Int16 ? :(_wide16(vload(SIMD.Vec{$W,Int8}, csb, n))) :
        :(_wide32(vload(SIMD.Vec{$W,Int8}, csb, n)))

    init = Expr(:block)
    for k = 1:NC
        push!(init.args, :($(Symbol("off_$k")) = sample_shifts[$k] - min_shift))
    end
    for j = 1:M, k = 1:NC
        push!(init.args, :($(Symbol("tI_$(j)_$k")) = zero(Int64)))
        push!(init.args, :($(Symbol("tQ_$(j)_$k")) = zero(Int64)))
    end

    # One antenna pass over the current block.
    function antenna_pass(j)
        blk_init = Expr(:block);
        flush = Expr(:block);
        chunk = Expr(:block);
        scalar = Expr(:block)
        for k = 1:NC
            aI = Symbol("aI_$(j)_$k");
            aQ = Symbol("aQ_$(j)_$k")
            tI = Symbol("tI_$(j)_$k");
            tQ = Symbol("tQ_$(j)_$k")
            offk = Symbol("off_$k")
            push!(blk_init.args, :($aI = zero(SIMD.Vec{$AW,Int32})))
            push!(blk_init.args, :($aQ = zero(SIMD.Vec{$AW,Int32})))
            push!(flush.args, :($tI += sum(_wide64($aI))))
            push!(flush.args, :($tQ += sum(_wide64($aQ))))
            if use_madd
                push!(
                    chunk.args,
                    quote
                        codew = _wide16(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                        $aI += _maddacc(codew, DI)
                        $aQ += _maddacc(codew, DQ)
                    end,
                )
            elseif use_smlal
                push!(
                    chunk.args,
                    quote
                        codew = _wide16(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                        $aI += _wide32(codew) * _wide32(DI)
                        $aQ += _wide32(codew) * _wide32(DQ)
                    end,
                )
            else
                push!(
                    chunk.args,
                    quote
                        codew = _wide32(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                        $aI += codew * DI
                        $aQ += codew * DQ
                    end,
                )
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
                sig = samples[base+n-1, $j]
                mr_s = Int32(real(sig));
                mi_s = Int32(imag(sig))
                cc_s = Int32(ccb[n]);
                sn_s = Int32(csb[n])
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

    # Finalize, dividing out the carrier amplitude.
    function tap_expr(k)
        if M == 1
            :(complex(
                Float64($(Symbol("tI_1_$k"))) / carrier_amp,
                Float64($(Symbol("tQ_1_$k"))) / carrier_amp,
            ))
        else
            ant = [
                :(complex(
                    Float64($(Symbol("tI_$(j)_$k"))) / carrier_amp,
                    Float64($(Symbol("tQ_$(j)_$k"))) / carrier_amp,
                )) for j = 1:M
            ]
            :(SVector{$M,ComplexF64}(tuple($(ant...))))
        end
    end
    taps = [tap_expr(k) for k = 1:NC]
    result = :(SVector{$NC}(tuple($(taps...))))

    quote
        W = $W
        carrier_amp = Float64(dc.carrier_amplitude)
        min_shift = minimum(sample_shifts)
        max_shift = maximum(sample_shifts)
        span = max_shift - min_shift
        num_rows = size(samples, 1)

        bufs = _scratch_buffers(dc)
        # Lane-overflow cap for the `W == 1` Portable path (#166).
        blk = _int16_flush_len(W, typemax(Int16), dc.blk)
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

        # Code is generated one-shot per block from its analytically advanced phase (a
        # continuing engine would allocate per sat per integration). The block starts at
        # output sample `min_shift`, so tap k at output n reads `extb[n + shift_k - min_shift]`.
        cps = Float64(upreferred(code_frequency / Hz)) / sampling_freq
        code_phase0 = Float64(code_phase)

        $init
        p_sig = Ptr{Int16}(pointer(samples))

        blk_off = 0
        @inbounds while blk_off < num_samples
            len = min(blk, num_samples - blk_off)

            # 1. This block's code into extb[1 .. len+span].
            gen_code!(
                view(extb, 1:(len+span)),
                signal,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase0 + cps * blk_off,
                min_shift,
            )

            # 2. This block's carrier sin/cos, shared across antennas.
            _int16_fill_carrier!(csb, ccb, reng, blk_off, len, phase0, Val(W))

            # 3. M antenna passes, flushed per block; `base` = signal row of local n=0.
            base = signal_start_sample + blk_off
            $correlate_passes

            blk_off += len
        end

        $result
    end
end

# One block of the dynamic fallback: accumulate every antenna/tap from the DI/DQ tile into
# the Int64 totals. A `Val{W}` function barrier so `W` and `TIw` are constants in the SIMD
# types (a bare local `W` would box the accumulators and allocate per block).
@inline function _int16_dyn_accumulate_block!(
    tI::Vector{Int64},
    tQ::Vector{Int64},
    dib::Vector{TIw},
    dqb::Vector{TIw},
    extb::Vector{Int8},
    sample_shifts::AbstractVector,
    min_shift::Integer,
    ::NumAnts{M},
    len::Int,
    ::Val{W},
) where {TIw,M,W}
    NC = length(sample_shifts)
    @inbounds for j = 1:M
        toff = (j - 1) * len
        for k = 1:NC
            offk = Int(sample_shifts[k]) - min_shift
            idx = (j - 1) * NC + k
            aI = zero(SIMD.Vec{W,Int32})
            aQ = zero(SIMD.Vec{W,Int32})
            n = 1
            while n + W - 1 <= len
                # Int16 tile: widen code through Int16 so the MAC lowers to SMLAL.
                if TIw === Int16
                    codev = _wide32(_wide16(vload(SIMD.Vec{W,Int8}, extb, n + offk)))
                else
                    codev = _wide32(vload(SIMD.Vec{W,Int8}, extb, n + offk))
                end
                DI = _wide32(vload(SIMD.Vec{W,TIw}, dib, toff + n))
                DQ = _wide32(vload(SIMD.Vec{W,TIw}, dqb, toff + n))
                aI += codev * DI
                aQ += codev * DQ
                n += W
            end
            tI[idx] += Int64(sum(aI))
            tQ[idx] += Int64(sum(aQ))
            while n <= len
                c = Int32(extb[n+offk])
                tI[idx] += Int64(c * Int32(dib[toff+n]))
                tQ[idx] += Int64(c * Int32(dqb[toff+n]))
                n += 1
            end
        end
    end
    nothing
end

# Dynamic fallback for runtime-sized `AbstractVector` sample shifts (#126); `SVector`
# shifts dispatch to the `@generated` kernel. Per block it fills the DI/DQ tile once
# (`_int16_fill_ditile!`) and accumulates each antenna/tap with the exact Int32 accumulate.
# Returns a `Vector` of `ComplexF64` (M=1) or `SVector{M,ComplexF64}` (M>1).
function _int16_hybrid_blocked!(
    dc::_Int16DC,
    samples::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    signal,
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
    num_rows = size(samples, 1)

    bufs = _scratch_buffers(dc)
    blk = dc.blk
    ncar = (cld(blk, W) + 4) * W
    ntile = blk * M
    length(bufs.extb) < blk + span && resize!(bufs.extb, blk + span)
    length(bufs.csb) < ncar && resize!(bufs.csb, ncar)
    length(bufs.ccb) < ncar && resize!(bufs.ccb, ncar)
    length(bufs.dib) < ntile && resize!(bufs.dib, ntile)
    length(bufs.dqb) < ntile && resize!(bufs.dqb, ntile)
    length(bufs.tsumI) < M * NC && resize!(bufs.tsumI, M * NC)
    length(bufs.tsumQ) < M * NC && resize!(bufs.tsumQ, M * NC)
    extb = bufs.extb
    csb = bufs.csb
    ccb = bufs.ccb
    dib = bufs.dib
    dqb = bufs.dqb
    tI = bufs.tsumI
    tQ = bufs.tsumQ
    @inbounds for i = 1:(M*NC)
        tI[i] = zero(Int64)
        tQ[i] = zero(Int64)
    end

    carrier_freq = Float64(upreferred(carrier_frequency / Hz))
    sampling_freq = Float64(upreferred(sampling_frequency / Hz))
    reng = carrier_engine(dc.table, carrier_freq / sampling_freq)
    phase0 = Float64(carrier_phase)

    cps = Float64(upreferred(code_frequency / Hz)) / sampling_freq
    code_phase0 = Float64(code_phase)

    p_sig = Ptr{Int16}(pointer(samples))

    blk_off = 0
    @inbounds while blk_off < num_samples
        len = min(blk, num_samples - blk_off)
        gen_code!(
            view(extb, 1:(len+span)),
            signal,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase0 + cps * blk_off,
            min_shift,
        )
        _int16_fill_carrier!(csb, ccb, reng, blk_off, len, phase0, Val(W))
        base = signal_start_sample + blk_off
        _int16_fill_ditile!(
            dib,
            dqb,
            samples,
            NumAnts{M}(),
            csb,
            ccb,
            base,
            len,
            num_rows,
            p_sig,
        )
        _int16_dyn_accumulate_block!(
            tI,
            tQ,
            dib,
            dqb,
            extb,
            sample_shifts,
            min_shift,
            NumAnts{M}(),
            len,
            Val(_INT16_W),
        )
        blk_off += len
    end

    # Divide out the carrier amplitude.
    carrier_amp = Float64(dc.carrier_amplitude)
    if M == 1
        return [
            complex(Float64(tI[k]) / carrier_amp, Float64(tQ[k]) / carrier_amp) for k = 1:NC
        ]
    else
        return [
            SVector{M,ComplexF64}(
                ntuple(
                    j -> complex(
                        Float64(tI[(j-1)*NC+k]) / carrier_amp,
                        Float64(tQ[(j-1)*NC+k]) / carrier_amp,
                    ),
                    M,
                ),
            ) for k = 1:NC
        ]
    end
end

# ── Multi-signal-per-sat tile-share ───────────────────────────────────────────
# Signals on one carrier (e.g. GPS L1 C/A + L1C-D + L1C-P) share the carrier wipe, so per
# block the carrier and DI/DQ tile are filled once and each signal correlates its own code
# against the tile: one downconvert per sat, as in the Float32
# `downconvert_and_correlate_fused_tuple!`.

# Fill the DI/DQ tile for the current block: antenna `j`'s slice is `dib/dqb[(j-1)*len + n]`.
@generated function _int16_fill_ditile!(
    dib,
    dqb,
    samples,
    ::NumAnts{M},
    csb,
    ccb,
    base,
    len,
    num_rows,
    p_sig,
) where {M}
    W = _INT16_W
    TI = eltype(dib)                   # DI/DQ tile element type = the wipe type
    widen = TI === Int16 ? :_wide16 : :_wide32
    passes = Expr(:block)
    for j = 1:M
        push!(
            passes.args,
            quote
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
                    sig = samples[base+n-1, $j]
                    mr_s = $TI(real(sig));
                    mi_s = $TI(imag(sig))
                    cc_s = $TI(ccb[n]);
                    sn_s = $TI(csb[n])
                    dib[toff+n] = mr_s * cc_s + mi_s * sn_s
                    dqb[toff+n] = mi_s * cc_s - mr_s * sn_s
                    n += 1
                end
            end,
        )
    end
    quote
        @inbounds begin
            $passes
        end
        nothing
    end
end

# Tile-share kernel for the signals of one satellite. Returns a tuple of per-signal
# `SVector{NCᵢ}`s, `@generated` over (M, the signals tuple).
@generated function _int16_hybrid_blocked_multi!(
    dc::_Int16DC,
    samples::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    gnss_signals::Tuple,
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
    TI = dc.parameters[2]              # carrier-wipe element type, per backend instance
    use_madd = _INT16_HAS_MADDWD && TI === Int16
    use_smlal = TI === Int16 && !use_madd   # see the single-signal kernel
    AW = use_madd ? W ÷ 2 : W
    N = length(all_shifts.parameters)
    NCs = [length(all_shifts.parameters[i]) for i = 1:N]

    setup = Expr(:block)
    maxspan = Expr(:call, :max)
    for i = 1:N
        push!(setup.args, :($(Symbol("sh_$i")) = all_shifts[$i]))
        push!(setup.args, :($(Symbol("mins_$i")) = minimum($(Symbol("sh_$i")))))
        push!(
            setup.args,
            :($(Symbol("span_$i")) = maximum($(Symbol("sh_$i"))) - $(Symbol("mins_$i"))),
        )
        push!(
            setup.args,
            :(
                $(Symbol("cps_$i")) =
                    Float64(upreferred(code_freqs[$i] / Hz)) / sampling_freq
            ),
        )
        for k = 1:NCs[i]
            push!(
                setup.args,
                :($(Symbol("off_$(i)_$k")) = $(Symbol("sh_$i"))[$k] - $(Symbol("mins_$i"))),
            )
        end
        push!(maxspan.args, Symbol("span_$i"))
        for j = 1:M, k = 1:NCs[i]
            push!(setup.args, :($(Symbol("tI_$(i)_$(j)_$k")) = zero(Int64)))
            push!(setup.args, :($(Symbol("tQ_$(i)_$(j)_$k")) = zero(Int64)))
        end
    end

    # Per signal: fill its code, then M antenna passes over the shared DI/DQ tile.
    function signal_corr(i)
        b = Expr(:block)
        push!(b.args, :(blk_phase = code_phases[$i] + $(Symbol("cps_$i")) * blk_off))
        push!(
            b.args,
            :(gen_code!(
                view(extb, 1:(len+$(Symbol("span_$i")))),
                gnss_signals[$i],
                prn,
                sampling_frequency,
                code_freqs[$i],
                blk_phase,
                $(Symbol("mins_$i")),
            )),
        )
        for j = 1:M
            blk_init = Expr(:block);
            flush = Expr(:block);
            chunk = Expr(:block);
            scalar = Expr(:block)
            for k = 1:NCs[i]
                aI = Symbol("aI_$k");
                aQ = Symbol("aQ_$k")
                tI = Symbol("tI_$(i)_$(j)_$k");
                tQ = Symbol("tQ_$(i)_$(j)_$k")
                offk = Symbol("off_$(i)_$k")
                push!(blk_init.args, :($aI = zero(SIMD.Vec{$AW,Int32})))
                push!(blk_init.args, :($aQ = zero(SIMD.Vec{$AW,Int32})))
                push!(flush.args, :($tI += sum(_wide64($aI))))
                push!(flush.args, :($tQ += sum(_wide64($aQ))))
                if use_madd
                    push!(
                        chunk.args,
                        quote
                            codew = _wide16(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                            $aI += _maddacc(codew, DI)
                            $aQ += _maddacc(codew, DQ)
                        end,
                    )
                elseif use_smlal
                    push!(
                        chunk.args,
                        quote
                            codew = _wide16(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                            $aI += _wide32(codew) * _wide32(DI)
                            $aQ += _wide32(codew) * _wide32(DQ)
                        end,
                    )
                else
                    push!(
                        chunk.args,
                        quote
                            codew = _wide32(vload(SIMD.Vec{$W,Int8}, extb, n + $offk))
                            $aI += codew * DI
                            $aQ += codew * DQ
                        end,
                    )
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

    # Divide out the carrier amplitude.
    function tap(i, k)
        if M == 1
            :(complex(
                Float64($(Symbol("tI_$(i)_1_$k"))) / carrier_amp,
                Float64($(Symbol("tQ_$(i)_1_$k"))) / carrier_amp,
            ))
        else
            ant = [
                :(complex(
                    Float64($(Symbol("tI_$(i)_$(j)_$k"))) / carrier_amp,
                    Float64($(Symbol("tQ_$(i)_$(j)_$k"))) / carrier_amp,
                )) for j = 1:M
            ]
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
        carrier_amp = Float64(dc.carrier_amplitude)
        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        reng = carrier_engine(dc.table, carrier_freq / sampling_freq)
        phase0 = Float64(carrier_phase)
        num_rows = size(samples, 1)
        $setup
        maxspan_v = $maxspan

        bufs = _scratch_buffers(dc)
        # Lane-overflow cap for the `W == 1` Portable path (#166).
        blk = _int16_flush_len(W, typemax(Int16), dc.blk)
        ncar = (cld(blk, W) + 4) * W
        length(bufs.extb) < blk + maxspan_v && resize!(bufs.extb, blk + maxspan_v)
        length(bufs.csb) < ncar && resize!(bufs.csb, ncar)
        length(bufs.ccb) < ncar && resize!(bufs.ccb, ncar)
        ntile = blk * $M
        length(bufs.dib) < ntile && resize!(bufs.dib, ntile)
        length(bufs.dqb) < ntile && resize!(bufs.dqb, ntile)
        extb = bufs.extb;
        csb = bufs.csb;
        ccb = bufs.ccb;
        dib = bufs.dib;
        dqb = bufs.dqb

        p_sig = Ptr{Int16}(pointer(samples))
        blk_off = 0
        @inbounds while blk_off < num_samples
            len = min(blk, num_samples - blk_off)
            _int16_fill_carrier!(csb, ccb, reng, blk_off, len, phase0, Val(W))
            base = signal_start_sample + blk_off
            _int16_fill_ditile!(
                dib,
                dqb,
                samples,
                NumAnts{$M}(),
                csb,
                ccb,
                base,
                len,
                num_rows,
                p_sig,
            )
            $sigs
            blk_off += len
        end
        $finals
    end
end

# Several signals per sat: one carrier downconvert shared through the tile-share kernel.
# Returns the per-signal `(new_correlator, is_integration_completed)` tuples.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal,TrackedSignal,Vararg{TrackedSignal}},
    per_signal_completed::Tuple,
    dc::_Int16DC,
    samples,
    code_doppler,
    code_phase,
    carrier_frequency,
    carrier_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_rows,
)
    params = map(signals) do head
        _signal_replica_params(head, code_doppler, code_phase, sampling_frequency, num_rows)
    end
    correlators = map(s -> s.correlator, signals)
    gnss_signals = map(s -> s.signal, signals)
    all_shifts = map(p -> p.sample_shifts, params)
    code_phases = map(p -> p.signal_code_phase, params)
    code_freqs = map(p -> p.code_frequency, params)
    new_accs = _int16_hybrid_blocked_multi!(
        dc,
        samples,
        NumAnts(correlators[1]),
        gnss_signals,
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

# The despread primitive on this backend's kernel; see `_despread_one_signal!` in
# downconvert_and_correlate_cpu.jl. `code_replica_size` and `use_band_cache` are ignored:
# the kernel generates its code block in its own scratch and has no band cache.
@inline _despread_one_signal!(
    dc::_Int16DC,
    correlator,
    samples,
    signal,
    prn,
    sample_shifts,
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    start_sample,
    num_samples,
    ::Any,
    ::Bool = true,
) = update_accumulator(
    correlator,
    get_accumulators(correlator) .+ _int16_hybrid_blocked!(
        dc,
        samples,
        NumAnts(correlator),
        signal,
        prn,
        sample_shifts,
        code_phase,
        carrier_phase,
        code_frequency,
        carrier_frequency,
        sampling_frequency,
        start_sample,
        num_samples,
    ),
)

# ── Group/measurement plumbing ────────────────────────────────────────────────
# The per-sat loop, group body and `downconvert_and_correlate(!)` entry points are the
# backend-agnostic ones of downconvert_and_correlate_cpu.jl; beyond the kernel hooks above
# this backend adds only the sample-type check and the threading trait.

# Reject non-`Complex{Int16}` sample buffers up front (12-bit ADC contract).
@inline _check_sample_type(::_Int16DC, m) =
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

@inline _threading(::Int16ThreadedDownconvertAndCorrelator) = _BatchLoop()
