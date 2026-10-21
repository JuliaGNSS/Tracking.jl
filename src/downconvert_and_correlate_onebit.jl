# One-bit (hard-limited) bit-wise downconvert + correlate backend for `Complex{Int16}`
# sample buffers. Measurement, carrier and code keep only their SIGN (+1 → bit 0, −1 →
# bit 1), so a ±1 product is an XOR and Σ over N samples is N − 2·popcount:
#
#   DI = mᵣ·cos + mᵢ·sin,  DQ = mᵢ·cos − mᵣ·sin      (all ±1) ⇒
#   Iₓ = Σ codeₓ·DI = (N − 2·pc(codeₓ⊻mᵣ⊻cos)) + (N − 2·pc(codeₓ⊻mᵢ⊻sin))
#   Qₓ = Σ codeₓ·DQ = (N − 2·pc(codeₓ⊻mᵢ⊻cos)) − (N − 2·pc(codeₓ⊻mᵣ⊻sin))
#
# Per `blk`-sample block (L1-resident scratch, as in the Int16 backend): pack the
# measurement sign planes per antenna, take the carrier sign planes off SinCosLUT's NCO,
# pack the code sign plane once and funnel-shift it per tap, then XOR + `count_ones`.
# This file also holds the rationale for the machinery the two-bit backend shares.
#
# 1-bit quantisation costs ≈2 dB of SNR (≈1 dB more for the 1-bit carrier). Results are
# `ComplexF64` (M=1) / `SVector{M,ComplexF64}` (M>1) on the popcount scale; downstream
# consumers are ratio-normalised, so the scale is immaterial.

import SinCosLUT
using SinCosLUT: generate_carrier_signs!

# Default strip-mine block length (samples), a multiple of 64 (whole UInt64 words).
const _ONEBIT_BLK = 8192

# UInt64 lanes per SIMD popcount step. Block planes are padded to a multiple of VW words
# with zeros, so the correlate loop runs whole VW chunks with no scalar tail.
const _OB_VW = 8

# ── sign-mask helpers: pack the sign bit of 64 lanes into a UInt64 (bit j ⇔ lane j < 0) ──
# The portable `icmp slt` + `bitcast <64 x i1> to i64` is optimal on x86 (`vpmovmskb` /
# `vpmov{b,w}2m`), but on aarch64 it becomes serial `addv`→GPR moves. There each lane's
# sign byte is isolated with the `{1,2,…,128}` mask and reduced by a 4-deep `addp` tree,
# one GPR move in all (≈2.4× Int8, ≈1.9× Int16 faster on Cortex-A78); the Int16 form
# first narrows the predicate to bytes. Both lowerings are bit-identical.
@static if Sys.ARCH === :aarch64
    @inline _ob_mm(v::SIMD.Vec{64,Int8}) = Base.llvmcall(
        (
            """
declare <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8>, <16 x i8>)
define i64 @entry(<64 x i8> %v) #0 {
  %c = icmp slt <64 x i8> %v, zeroinitializer
  %s = sext <64 x i1> %c to <64 x i8>
  %v0 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 0,i32 1,i32 2,i32 3,i32 4,i32 5,i32 6,i32 7,i32 8,i32 9,i32 10,i32 11,i32 12,i32 13,i32 14,i32 15>
  %v1 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 16,i32 17,i32 18,i32 19,i32 20,i32 21,i32 22,i32 23,i32 24,i32 25,i32 26,i32 27,i32 28,i32 29,i32 30,i32 31>
  %v2 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 32,i32 33,i32 34,i32 35,i32 36,i32 37,i32 38,i32 39,i32 40,i32 41,i32 42,i32 43,i32 44,i32 45,i32 46,i32 47>
  %v3 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 48,i32 49,i32 50,i32 51,i32 52,i32 53,i32 54,i32 55,i32 56,i32 57,i32 58,i32 59,i32 60,i32 61,i32 62,i32 63>
  %bm0 = and <16 x i8> %v0, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %bm1 = and <16 x i8> %v1, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %bm2 = and <16 x i8> %v2, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %bm3 = and <16 x i8> %v3, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %p0 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %bm0, <16 x i8> %bm1)
  %p1 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %bm2, <16 x i8> %bm3)
  %p2 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %p0, <16 x i8> %p1)
  %p3 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %p2, <16 x i8> %p2)
  %r = bitcast <16 x i8> %p3 to <2 x i64>
  %lo = extractelement <2 x i64> %r, i32 0
  ret i64 %lo }
attributes #0 = { alwaysinline }""",
            "entry",
        ),
        UInt64,
        Tuple{NTuple{64,Base.VecElement{Int8}}},
        v.data,
    )
    @inline _ob_mm(v::SIMD.Vec{64,Int16}) = Base.llvmcall(
        (
            """
declare <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8>, <16 x i8>)
define i64 @entry(<64 x i16> %v) #0 {
  %c = icmp slt <64 x i16> %v, zeroinitializer
  %s = sext <64 x i1> %c to <64 x i8>
  %v0 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 0,i32 1,i32 2,i32 3,i32 4,i32 5,i32 6,i32 7,i32 8,i32 9,i32 10,i32 11,i32 12,i32 13,i32 14,i32 15>
  %v1 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 16,i32 17,i32 18,i32 19,i32 20,i32 21,i32 22,i32 23,i32 24,i32 25,i32 26,i32 27,i32 28,i32 29,i32 30,i32 31>
  %v2 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 32,i32 33,i32 34,i32 35,i32 36,i32 37,i32 38,i32 39,i32 40,i32 41,i32 42,i32 43,i32 44,i32 45,i32 46,i32 47>
  %v3 = shufflevector <64 x i8> %s, <64 x i8> undef, <16 x i32> <i32 48,i32 49,i32 50,i32 51,i32 52,i32 53,i32 54,i32 55,i32 56,i32 57,i32 58,i32 59,i32 60,i32 61,i32 62,i32 63>
  %bm0 = and <16 x i8> %v0, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %bm1 = and <16 x i8> %v1, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %bm2 = and <16 x i8> %v2, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %bm3 = and <16 x i8> %v3, <i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128,i8 1,i8 2,i8 4,i8 8,i8 16,i8 32,i8 64,i8 -128>
  %p0 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %bm0, <16 x i8> %bm1)
  %p1 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %bm2, <16 x i8> %bm3)
  %p2 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %p0, <16 x i8> %p1)
  %p3 = call <16 x i8> @llvm.aarch64.neon.addp.v16i8(<16 x i8> %p2, <16 x i8> %p2)
  %r = bitcast <16 x i8> %p3 to <2 x i64>
  %lo = extractelement <2 x i64> %r, i32 0
  ret i64 %lo }
attributes #0 = { alwaysinline }""",
            "entry",
        ),
        UInt64,
        Tuple{NTuple{64,Base.VecElement{Int16}}},
        v.data,
    )
else
    @inline _ob_mm(v::SIMD.Vec{64,Int8}) = Base.llvmcall(
        (
            """
define i64 @entry(<64 x i8> %v) #0 { %c = icmp slt <64 x i8> %v, zeroinitializer
  %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""",
            "entry",
        ),
        UInt64,
        Tuple{NTuple{64,Base.VecElement{Int8}}},
        v.data,
    )
    @inline _ob_mm(v::SIMD.Vec{64,Int16}) = Base.llvmcall(
        (
            """
define i64 @entry(<64 x i16> %v) #0 { %c = icmp slt <64 x i16> %v, zeroinitializer
  %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""",
            "entry",
        ),
        UInt64,
        Tuple{NTuple{64,Base.VecElement{Int16}}},
        v.data,
    )
end
@inline _ob_masklast(w::UInt64, r::Int) = r == 0 ? w : w & ((UInt64(1) << r) - UInt64(1))

# Pack a tap's code sign plane: bit n = sign(extb[byteoff+n]). `extb` is padded so the last
# (masked) word can be read whole.
@inline function _ob_pack_code!(
    codeb::Vector{UInt64},
    koff::Int,
    extb::Vector{Int8},
    byteoff::Int,
    nwb::Int,
    r::Int,
)
    @inbounds for w = 0:(nwb-1)
        wd = _ob_mm(SIMD.vload(SIMD.Vec{64,Int8}, extb, byteoff + (w << 6) + 1))
        codeb[koff+w+1] = w == nwb - 1 ? _ob_masklast(wd, r) : wd
    end
    nothing
end

# Derive a tap plane from the prompt-extended plane by a funnel bit-shift, so one movemask
# of the code block feeds all taps: `dst[w] = pe[w+off]` in bit terms, for any `off ≥ 0`
# (`off ÷ 64` whole words + `off mod 64` bits). The caller zero-pads `pe` so
# `pe[nwb+wsh+1]` is valid (`pepad` in the kernels).
@inline function _ob_shift_plane!(
    dst::Vector{UInt64},
    doff::Int,
    pe::Vector{UInt64},
    off::Int,
    nwb::Int,
)
    wsh = off >> 6
    bit = off & 63
    if bit == 0
        @inbounds for w = 1:nwb
            dst[doff+w] = pe[w+wsh]
        end
    else
        lo = bit;
        hi = 64 - bit
        @inbounds for w = 1:nwb
            dst[doff+w] = (pe[w+wsh] >> lo) | (pe[w+wsh+1] << hi)
        end
    end
    nothing
end

# Mask the last word to `r` bits (if partial) and zero the VW pad words [nwb+1, nwv], which
# then contribute nothing to the correlate.
@inline function _ob_finish!(buf::Vector{UInt64}, off::Int, nwb::Int, nwv::Int, r::Int)
    @inbounds (r != 0) && (buf[off+nwb] &= (UInt64(1) << r) - UInt64(1))
    @inbounds for w = (nwb+1):nwv
        buf[off+w] = 0
    end
    nothing
end
# The packer already masked the last word; just zero the VW pad.
@inline function _ob_zeropad!(buf::Vector{UInt64}, off::Int, nwb::Int, nwv::Int)
    @inbounds for w = (nwb+1):nwv
        buf[off+w] = 0
    end
    nothing
end

# Direct measurement sign pack for a single-sat group (cheaper than the band pack plus
# realign): sign(real)→mrb, sign(imag)→mib.
@inline function _ob_pack_meas!(
    mrb::Vector{UInt64},
    mib::Vector{UInt64},
    joff::Int,
    signal,
    j::Int,
    p_sig::Ptr{Int16},
    colbase::Int,
    base::Int,
    len::Int,
    nwb::Int,
    r::Int,
)
    full = fld(len, 64)
    @inbounds for w = 0:(full-1)
        byte_off = (colbase + base + (w << 6) - 1) * 2 * sizeof(Int16)
        re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
        mrb[joff+w+1] = _ob_mm(re)
        mib[joff+w+1] = _ob_mm(im)
    end
    if r != 0
        wr = zero(UInt64);
        wi = zero(UInt64)
        @inbounds for i = 0:(r-1)
            sig = signal[base+(full<<6)+i, j]
            (real(sig) < 0) && (wr |= UInt64(1) << i)
            (imag(sig) < 0) && (wi |= UInt64(1) << i)
        end
        mrb[joff+full+1] = wr;
        mib[joff+full+1] = wi
    end
    nothing
end

# Per-thread scratch, grown lazily and reused (allocation-free in steady state): code block
# `extb`, prompt-extended code plane `peb`, carrier sign planes, per-tap code planes
# (`codeb`, NC·wpbv) and per-antenna measurement planes (`mrb`/`mib`, M·wpbv). Plane
# strides `wpbv` are padded to a multiple of `_OB_VW`.
struct OneBitScratchBuffers
    extb::Vector{Int8}
    peb::Vector{UInt64}
    sinw::Vector{UInt64}
    cosw::Vector{UInt64}
    codeb::Vector{UInt64}
    mrb::Vector{UInt64}
    mib::Vector{UInt64}
end
OneBitScratchBuffers() =
    OneBitScratchBuffers(Int8[], UInt64[], UInt64[], UInt64[], UInt64[], UInt64[], UInt64[])

# The measurement sign planes are the same for every satellite on a band, so with more than
# one satellite `_dc_one_group!` packs them once per group (aligned to sample 1, with funnel
# pad words), before the per-sat loop; each sat funnel-realigns its block from them.
struct OneBitBandBuffers
    mrband::Vector{UInt64}
    miband::Vector{UInt64}
end
OneBitBandBuffers() = OneBitBandBuffers(UInt64[], UInt64[])

# Pack the band's `M` measurement sign planes; plane `j` occupies `stride` words at
# `(j-1)*stride`, ending in zero funnel pad words.
function _ob_pack_band!(
    band::OneBitBandBuffers,
    signal,
    num_samples::Int,
    M::Int,
    p_sig::Ptr{Int16},
    num_rows::Int,
)
    nwb = cld(num_samples, 64)
    stride = nwb + 2                                # 2 funnel-pad words (realign reads past nwb)
    length(band.mrband) < M * stride && resize!(band.mrband, M * stride)
    length(band.miband) < M * stride && resize!(band.miband, M * stride)
    r = num_samples & 63
    full = fld(num_samples, 64)
    @inbounds for j = 1:M
        off = (j - 1) * stride
        colbase = (j - 1) * num_rows
        for w = 0:(full-1)
            byte_off = (colbase + w * 64) * 2 * sizeof(Int16)
            re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
            band.mrband[off+w+1] = _ob_mm(re)
            band.miband[off+w+1] = _ob_mm(im)
        end
        if r != 0
            wr = zero(UInt64);
            wi = zero(UInt64)
            for i = 0:(r-1)
                sig = signal[full*64+i+1, j]
                (real(sig) < 0) && (wr |= UInt64(1) << i)
                (imag(sig) < 0) && (wi |= UInt64(1) << i)
            end
            band.mrband[off+full+1] = wr;
            band.miband[off+full+1] = wi
        end
        band.mrband[off+nwb+1] = 0;
        band.miband[off+nwb+1] = 0     # funnel pad
        band.mrband[off+nwb+2] = 0;
        band.miband[off+nwb+2] = 0
    end
    nothing
end

# Funnel-realign antenna `j`'s block (`len` samples from 1-based sample `base`) from the
# band planes into the thread-local `mrb`/`mib` at word offset `moff`.
@inline function _ob_realign_meas!(
    mrb::Vector{UInt64},
    mib::Vector{UInt64},
    moff::Int,
    band::OneBitBandBuffers,
    jbase::Int,
    base::Int,
    nwb::Int,
)
    b0 = base - 1                                  # 0-based absolute start bit
    sw = b0 >> 6
    bit = b0 & 63
    @inbounds if bit == 0
        for w = 1:nwb
            mrb[moff+w] = band.mrband[jbase+sw+w]
            mib[moff+w] = band.miband[jbase+sw+w]
        end
    else
        hi = 64 - bit
        for w = 1:nwb
            mrb[moff+w] =
                (band.mrband[jbase+sw+w] >> bit) | (band.mrband[jbase+sw+w+1] << hi)
            mib[moff+w] =
                (band.miband[jbase+sw+w] >> bit) | (band.miband[jbase+sw+w+1] << hi)
        end
    end
    nothing
end

"""
$(SIGNATURES)

One-bit (hard-limited) bit-wise CPU downconvert + correlate backend
(single-threaded). Opt-in alternative to [`CPUDownconvertAndCorrelator`] for
`Complex{Int16}` sample buffers: it 1-bit-quantises measurement, carrier and code and
correlates with XOR + popcount, trading ≈2–3 dB of SNR for bit-wise speed. Construct it
**once outside** the `track!` loop and pass it via the `downconvert_and_correlator`
keyword for an allocation-free steady state.
"""
struct OneBitDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::OneBitScratchBuffers
    band::OneBitBandBuffers
    blk::Int
end

"""
$(SIGNATURES)

Multi-threaded [`OneBitDownconvertAndCorrelator`](@ref), one `OneBitScratchBuffers` per
thread.
"""
struct OneBitThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{OneBitScratchBuffers}
    band::OneBitBandBuffers
    blk::Int
end

OneBitDownconvertAndCorrelator(; blk::Integer = _ONEBIT_BLK) =
    OneBitDownconvertAndCorrelator(OneBitScratchBuffers(), OneBitBandBuffers(), Int(blk))

OneBitThreadedDownconvertAndCorrelator(; blk::Integer = _ONEBIT_BLK) =
    OneBitThreadedDownconvertAndCorrelator(
        [OneBitScratchBuffers() for _ = 1:Threads.maxthreadid()],
        OneBitBandBuffers(),
        Int(blk),
    )

const _OneBitDC =
    Union{OneBitDownconvertAndCorrelator,OneBitThreadedDownconvertAndCorrelator}

@inline _ob_scratch(dc::OneBitDownconvertAndCorrelator) = dc.buffers
@inline _ob_scratch(dc::OneBitThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

@inline _ob_num_ants_val(::AbstractCorrelator{M}) where {M} = NumAnts{M}()

# ── The one-bit hybrid-blocked kernel ─────────────────────────────────────────
# Returns this integration's contribution per tap, `SVector{NC}` of `ComplexF64` (M=1) or
# `SVector{M,ComplexF64}` (M>1). `@generated` over (NC, M), the per-(antenna, tap) popcount
# accumulators in named locals. An `SVector`, not a `Vector`: the compiler does not elide a
# `Vector` at the noise reference's call site (112 B per despread).
@generated function _onebit_hybrid_blocked!(
    dc::_OneBitDC,
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
    use_band_cache::Bool = true,
) where {M,NC}
    # Accumulator declarations: A/B (→ I) and C/E (→ Q) popcounts per (antenna j, tap k).
    init = Expr(:block)
    for k = 1:NC
        push!(init.args, :($(Symbol("off_$k")) = Int(sample_shifts[$k]) - min_shift))
    end
    for j = 1:M, k = 1:NC
        for s in ("A", "B", "C", "E")
            push!(init.args, :($(Symbol("$(s)_$(j)_$k")) = zero(SIMD.Vec{_OB_VW,UInt64})))
        end
    end

    # Per-VW-chunk unrolled accumulate: shared sample⊻carrier products per antenna, then one
    # XOR + `count_ones` (VPOPCNTQ) per (tap, product), accumulated into Vec lane counters.
    body = Expr(:block)
    for j = 1:M
        push!(body.args, :(MR = vload(SIMD.Vec{_OB_VW,UInt64}, mrb, $(j - 1) * wpbv + i)))
        push!(body.args, :(MI = vload(SIMD.Vec{_OB_VW,UInt64}, mib, $(j - 1) * wpbv + i)))
        push!(body.args, :(prc = MR ⊻ CCv))
        push!(body.args, :(pis = MI ⊻ CSv))
        push!(body.args, :(pic = MI ⊻ CCv))
        push!(body.args, :(prs = MR ⊻ CSv))
        for k = 1:NC
            push!(
                body.args,
                :(CW = vload(SIMD.Vec{_OB_VW,UInt64}, codeb, $(k - 1) * wpbv + i)),
            )
            push!(body.args, :($(Symbol("A_$(j)_$k")) += count_ones(CW ⊻ prc)))
            push!(body.args, :($(Symbol("B_$(j)_$k")) += count_ones(CW ⊻ pis)))
            push!(body.args, :($(Symbol("C_$(j)_$k")) += count_ones(CW ⊻ pic)))
            push!(body.args, :($(Symbol("E_$(j)_$k")) += count_ones(CW ⊻ prs)))
        end
    end

    # Finalize: reduce lane counters, then Iₓ = 2N − 2(A+B); Qₓ = 2(E − C).
    s(sym) = :(Int64(sum($sym)))
    function tapval(k)
        if M == 1
            :(complex(
                Float64(2 * N - 2 * ($(s(Symbol("A_1_$k"))) + $(s(Symbol("B_1_$k"))))),
                Float64(2 * ($(s(Symbol("E_1_$k"))) - $(s(Symbol("C_1_$k"))))),
            ))
        else
            ant = [
                :(complex(
                    Float64(
                        2 * N - 2 * ($(s(Symbol("A_$(j)_$k"))) + $(s(Symbol("B_$(j)_$k")))),
                    ),
                    Float64(2 * ($(s(Symbol("E_$(j)_$k"))) - $(s(Symbol("C_$(j)_$k"))))),
                )) for j = 1:M
            ]
            :(SVector{$M,ComplexF64}(tuple($(ant...))))
        end
    end
    ret = :(SVector{$NC}(tuple($([tapval(k) for k = 1:NC]...))))

    quote
        # The code plane keeps only the sign, exact for binary codes (BPSK, BOC, TMBOC) but
        # not for CBOC's multi-level amplitude. Gate on the modulation: GNSSSignals'
        # CBOC replica is integer-typed too.
        get_modulation(signal_type) isa GNSSSignals.CBOC && throw(
            ArgumentError(
                string(
                    "OneBitDownconvertAndCorrelator supports binary (±1) codes only ",
                    "(BPSK, BOC, TMBOC); got ",
                    get_signal_name(signal_type),
                    " with ",
                    typeof(get_modulation(signal_type)),
                    " modulation. Bit-wise correlation keeps only the code sign, so it ",
                    "cannot represent CBOC — a multi-level, amplitude-carrying code.",
                ),
            ),
        )
        N = Int(num_samples)
        min_shift = minimum(sample_shifts)
        span = maximum(sample_shifts) - min_shift
        $init

        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        cps_car = carrier_freq / sampling_freq
        carphase = Float64(carrier_phase)
        cps_code = Float64(upreferred(code_frequency / Hz)) / sampling_freq
        code_phase0 = Float64(code_phase)
        num_rows = size(signal, 1)
        p_sig = Ptr{Int16}(pointer(signal))
        band = dc.band                              # band-shared measurement sign planes (>1 sat)
        band_shared = use_band_cache && !isempty(band.mrband)  # pre-packed per group when it pays off
        bandstride = cld(num_rows, 64) + 2          # per-antenna plane stride in `band` (see _ob_pack_band!)

        bufs = _ob_scratch(dc)
        blk = dc.blk
        wpb = cld(blk, 64)                          # words per block plane
        wpbv = cld(wpb, _OB_VW) * _OB_VW            # …padded to a whole number of VW chunks
        pepad = (span >> 6) + 1                     # zero pad words the funnel shift reads past nwe
        pewords = cld(blk + span, 64) + pepad + _OB_VW   # extended prompt plane + funnel pad
        length(bufs.extb) < blk + span + 64 && resize!(bufs.extb, blk + span + 64)
        length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
        length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
        length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
        length(bufs.codeb) < $NC * wpbv && resize!(bufs.codeb, $NC * wpbv)
        length(bufs.mrb) < $M * wpbv && resize!(bufs.mrb, $M * wpbv)
        length(bufs.mib) < $M * wpbv && resize!(bufs.mib, $M * wpbv)
        extb = bufs.extb;
        peb = bufs.peb;
        sinw = bufs.sinw;
        cosw = bufs.cosw
        codeb = bufs.codeb;
        mrb = bufs.mrb;
        mib = bufs.mib

        blk_off = 0
        @inbounds while blk_off < N
            len = min(blk, N - blk_off)
            nwb = cld(len, 64)
            nwv = cld(nwb, _OB_VW) * _OB_VW
            r = len & 63
            # code (Int8 ±1) for samples [min_shift, len+span): pack once, funnel-shift per tap
            gen_code!(
                view(extb, 1:(len+span)),
                signal_type,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase0 + cps_code * blk_off,
                min_shift,
            )
            nwe = cld(len + span, 64)
            _ob_pack_code!(peb, 0, extb, 0, nwe, (len + span) & 63)
            for _pw = 1:pepad                           # zero the words the funnel shift reads past nwe
                peb[nwe+_pw] = 0
            end
            $(Expr(
                :block,
                (
                    quote
                        _ob_shift_plane!(
                            codeb,
                            $(k - 1) * wpbv,
                            peb,
                            $(Symbol("off_$k")),
                            nwb,
                        )
                        _ob_finish!(codeb, $(k - 1) * wpbv, nwb, nwv, r)
                    end for k = 1:NC
                )...,
            ))
            # carrier sign planes, shared across antennas and taps
            generate_carrier_signs!(
                sinw,
                cosw,
                len,
                cps_car;
                phase = carphase + cps_car * blk_off,
            )
            _ob_zeropad!(sinw, 0, nwb, nwv)
            _ob_zeropad!(cosw, 0, nwb, nwv)
            # measurement sign planes per antenna: realigned from the band planes (>1 sat)
            # or packed directly (1 sat)
            base = signal_start_sample + blk_off
            if band_shared
                $(Expr(
                    :block,
                    (
                        quote
                            _ob_realign_meas!(
                                mrb,
                                mib,
                                $(j - 1) * wpbv,
                                band,
                                $(j - 1) * bandstride,
                                base,
                                nwb,
                            )
                            _ob_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                            _ob_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
                        end for j = 1:M
                    )...,
                ))
            else
                $(Expr(
                    :block,
                    (
                        quote
                            _ob_pack_meas!(
                                mrb,
                                mib,
                                $(j - 1) * wpbv,
                                signal,
                                $j,
                                p_sig,
                                $(j - 1) * num_rows,
                                base,
                                len,
                                nwb,
                                r,
                            )
                            _ob_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                            _ob_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
                        end for j = 1:M
                    )...,
                ))
            end
            # XOR + popcount accumulate over whole VW chunks
            i = 1
            while i <= nwv
                CCv = vload(SIMD.Vec{_OB_VW,UInt64}, cosw, i)
                CSv = vload(SIMD.Vec{_OB_VW,UInt64}, sinw, i)
                $body
                i += _OB_VW
            end
            blk_off += len
        end
        $ret
    end
end

# Dynamic fallback for runtime-sized `AbstractVector` sample shifts (#126); `SVector`
# shifts dispatch to the `@generated` kernel. Same block pipeline, with taps and antennas
# looped at runtime into Int64 totals; bit-identical to the `@generated` kernel. Returns a
# `Vector` of `ComplexF64` (M=1) or `SVector{M,ComplexF64}` (M>1).
function _onebit_hybrid_blocked!(
    dc::_OneBitDC,
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
    use_band_cache::Bool = true,
) where {M}
    get_modulation(signal_type) isa GNSSSignals.CBOC && throw(
        ArgumentError(
            string(
                "OneBitDownconvertAndCorrelator supports binary (±1) codes only ",
                "(BPSK, BOC, TMBOC); got ",
                get_signal_name(signal_type),
                " with ",
                typeof(get_modulation(signal_type)),
                " modulation. Bit-wise correlation keeps only the code sign, so it ",
                "cannot represent CBOC — a multi-level, amplitude-carrying code.",
            ),
        ),
    )
    NC = length(sample_shifts)
    N = Int(num_samples)
    min_shift = Int(minimum(sample_shifts))
    span = Int(maximum(sample_shifts)) - min_shift
    offs = [Int(sample_shifts[k]) - min_shift for k = 1:NC]

    sampling_freq = Float64(upreferred(sampling_frequency / Hz))
    carrier_freq = Float64(upreferred(carrier_frequency / Hz))
    cps_car = carrier_freq / sampling_freq
    carphase = Float64(carrier_phase)
    cps_code = Float64(upreferred(code_frequency / Hz)) / sampling_freq
    code_phase0 = Float64(code_phase)
    num_rows = size(signal, 1)
    p_sig = Ptr{Int16}(pointer(signal))
    band = dc.band                              # band-shared measurement sign planes (>1 sat)
    band_shared = use_band_cache && !isempty(band.mrband)
    bandstride = cld(num_rows, 64) + 2

    bufs = _ob_scratch(dc)
    blk = dc.blk
    wpb = cld(blk, 64)
    wpbv = cld(wpb, _OB_VW) * _OB_VW
    pepad = (span >> 6) + 1                     # zero pad words the funnel shift reads past nwe
    pewords = cld(blk + span, 64) + pepad + _OB_VW
    length(bufs.extb) < blk + span + 64 && resize!(bufs.extb, blk + span + 64)
    length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
    length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
    length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
    length(bufs.codeb) < NC * wpbv && resize!(bufs.codeb, NC * wpbv)
    length(bufs.mrb) < M * wpbv && resize!(bufs.mrb, M * wpbv)
    length(bufs.mib) < M * wpbv && resize!(bufs.mib, M * wpbv)
    extb = bufs.extb
    peb = bufs.peb
    sinw = bufs.sinw
    cosw = bufs.cosw
    codeb = bufs.codeb
    mrb = bufs.mrb
    mib = bufs.mib

    A = zeros(Int64, M, NC)
    B = zeros(Int64, M, NC)
    C = zeros(Int64, M, NC)
    E = zeros(Int64, M, NC)

    blk_off = 0
    @inbounds while blk_off < N
        len = min(blk, N - blk_off)
        nwb = cld(len, 64)
        nwv = cld(nwb, _OB_VW) * _OB_VW
        r = len & 63
        gen_code!(
            view(extb, 1:(len+span)),
            signal_type,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase0 + cps_code * blk_off,
            min_shift,
        )
        nwe = cld(len + span, 64)
        _ob_pack_code!(peb, 0, extb, 0, nwe, (len + span) & 63)
        for _pw = 1:pepad
            peb[nwe+_pw] = 0
        end
        for k = 1:NC
            _ob_shift_plane!(codeb, (k - 1) * wpbv, peb, offs[k], nwb)
            _ob_finish!(codeb, (k - 1) * wpbv, nwb, nwv, r)
        end
        generate_carrier_signs!(
            sinw,
            cosw,
            len,
            cps_car;
            phase = carphase + cps_car * blk_off,
        )
        _ob_zeropad!(sinw, 0, nwb, nwv)
        _ob_zeropad!(cosw, 0, nwb, nwv)
        base = signal_start_sample + blk_off
        if band_shared
            for j = 1:M
                _ob_realign_meas!(
                    mrb,
                    mib,
                    (j - 1) * wpbv,
                    band,
                    (j - 1) * bandstride,
                    base,
                    nwb,
                )
                _ob_finish!(mrb, (j - 1) * wpbv, nwb, nwv, r)
                _ob_finish!(mib, (j - 1) * wpbv, nwb, nwv, r)
            end
        else
            for j = 1:M
                _ob_pack_meas!(
                    mrb,
                    mib,
                    (j - 1) * wpbv,
                    signal,
                    j,
                    p_sig,
                    (j - 1) * num_rows,
                    base,
                    len,
                    nwb,
                    r,
                )
                _ob_finish!(mrb, (j - 1) * wpbv, nwb, nwv, r)
                _ob_finish!(mib, (j - 1) * wpbv, nwb, nwv, r)
            end
        end
        ci = 1
        while ci <= nwv
            CCv = vload(SIMD.Vec{_OB_VW,UInt64}, cosw, ci)
            CSv = vload(SIMD.Vec{_OB_VW,UInt64}, sinw, ci)
            for j = 1:M
                MR = vload(SIMD.Vec{_OB_VW,UInt64}, mrb, (j - 1) * wpbv + ci)
                MI = vload(SIMD.Vec{_OB_VW,UInt64}, mib, (j - 1) * wpbv + ci)
                prc = MR ⊻ CCv
                pis = MI ⊻ CSv
                pic = MI ⊻ CCv
                prs = MR ⊻ CSv
                for k = 1:NC
                    CW = vload(SIMD.Vec{_OB_VW,UInt64}, codeb, (k - 1) * wpbv + ci)
                    A[j, k] += Int64(sum(count_ones(CW ⊻ prc)))
                    B[j, k] += Int64(sum(count_ones(CW ⊻ pis)))
                    C[j, k] += Int64(sum(count_ones(CW ⊻ pic)))
                    E[j, k] += Int64(sum(count_ones(CW ⊻ prs)))
                end
            end
            ci += _OB_VW
        end
        blk_off += len
    end

    if M == 1
        return [
            complex(
                Float64(2 * N - 2 * (A[1, k] + B[1, k])),
                Float64(2 * (E[1, k] - C[1, k])),
            ) for k = 1:NC
        ]
    else
        return [
            SVector{M,ComplexF64}(
                ntuple(
                    j -> complex(
                        Float64(2 * N - 2 * (A[j, k] + B[j, k])),
                        Float64(2 * (E[j, k] - C[j, k])),
                    ),
                    M,
                ),
            ) for k = 1:NC
        ]
    end
end

# ── Multi-signal-per-sat tile-share ───────────────────────────────────────────
# Signals on one carrier (e.g. GPS L1 C/A + L1C-D + L1C-P) share the carrier and
# measurement planes, so per block those are built once and each signal correlates its own
# code against them; bit-identical to correlating each signal alone. `@generated` over
# (M, the signals tuple); `SVector` shifts only, like the other tile-shares.
@generated function _onebit_hybrid_blocked_multi!(
    dc::_OneBitDC,
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
    N = length(all_shifts.parameters)
    NCs = [length(all_shifts.parameters[i]) for i = 1:N]
    maxNC = maximum(NCs)

    # Per-signal shift/span/offset locals and per-(signal, antenna, tap) A/B (→ I) and C/E
    # (→ Q) popcount accumulators.
    setup = Expr(:block)
    maxspan = Expr(:call, :max)
    for i = 1:N
        push!(setup.args, :($(Symbol("sh_$i")) = all_shifts[$i]))
        push!(setup.args, :($(Symbol("mins_$i")) = Int(minimum($(Symbol("sh_$i"))))))
        push!(
            setup.args,
            :(
                $(Symbol("span_$i")) =
                    Int(maximum($(Symbol("sh_$i")))) - $(Symbol("mins_$i"))
            ),
        )
        push!(
            setup.args,
            :(
                $(Symbol("cpsc_$i")) =
                    Float64(upreferred(code_freqs[$i] / Hz)) / sampling_freq
            ),
        )
        push!(setup.args, :($(Symbol("cph0_$i")) = Float64(code_phases[$i])))
        for k = 1:NCs[i]
            push!(
                setup.args,
                :(
                    $(Symbol("off_$(i)_$k")) =
                        Int($(Symbol("sh_$i"))[$k]) - $(Symbol("mins_$i"))
                ),
            )
        end
        push!(maxspan.args, Symbol("span_$i"))
        for j = 1:M, k = 1:NCs[i]
            for sfx in ("A", "B", "C", "E")
                push!(
                    setup.args,
                    :($(Symbol("$(sfx)_$(i)_$(j)_$k")) = zero(SIMD.Vec{_OB_VW,UInt64})),
                )
            end
        end
    end

    # CBOC gate per signal (see the single-signal kernel).
    gates = Expr(:block)
    for i = 1:N
        push!(
            gates.args,
            :(
                get_modulation(signal_types[$i]) isa GNSSSignals.CBOC && throw(
                    ArgumentError(
                        string(
                            "OneBitDownconvertAndCorrelator supports binary (±1) codes ",
                            "only (BPSK, BOC, TMBOC); got ",
                            get_signal_name(signal_types[$i]),
                            " with ",
                            typeof(get_modulation(signal_types[$i])),
                            " modulation. Bit-wise correlation keeps only the code sign, ",
                            "so it cannot represent CBOC — a multi-level, amplitude-",
                            "carrying code.",
                        ),
                    ),
                )
            ),
        )
    end

    # Measurement planes per antenna (band realign or direct pack), once per block.
    realign = Expr(:block)
    directpack = Expr(:block)
    for j = 1:M
        push!(
            realign.args,
            quote
                _ob_realign_meas!(
                    mrb,
                    mib,
                    $(j - 1) * wpbv,
                    band,
                    $(j - 1) * bandstride,
                    base,
                    nwb,
                )
                _ob_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                _ob_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
            end,
        )
        push!(
            directpack.args,
            quote
                _ob_pack_meas!(
                    mrb,
                    mib,
                    $(j - 1) * wpbv,
                    signal,
                    $j,
                    p_sig,
                    $(j - 1) * num_rows,
                    base,
                    len,
                    nwb,
                    r,
                )
                _ob_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                _ob_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
            end,
        )
    end

    # Per signal: code fill + pack, then XOR + popcount against the shared planes.
    function signal_block(i)
        b = Expr(:block)
        push!(b.args, :(nwe = cld(len + $(Symbol("span_$i")), 64)))
        push!(
            b.args,
            :(gen_code!(
                view(extb, 1:(len+$(Symbol("span_$i")))),
                signal_types[$i],
                prn,
                sampling_frequency,
                code_freqs[$i],
                $(Symbol("cph0_$i")) + $(Symbol("cpsc_$i")) * blk_off,
                $(Symbol("mins_$i")),
            )),
        )
        push!(
            b.args,
            :(_ob_pack_code!(peb, 0, extb, 0, nwe, (len + $(Symbol("span_$i"))) & 63)),
        )
        push!(b.args, :(
            for _pw = 1:(($(Symbol("span_$i"))>>6)+1)   # zero the funnel-shift read pad
                peb[nwe+_pw] = 0
            end
        ))
        for k = 1:NCs[i]
            push!(
                b.args,
                :(_ob_shift_plane!(
                    codeb,
                    $(k - 1) * wpbv,
                    peb,
                    $(Symbol("off_$(i)_$k")),
                    nwb,
                )),
            )
            push!(b.args, :(_ob_finish!(codeb, $(k - 1) * wpbv, nwb, nwv, r)))
        end
        corr = Expr(:block)
        for j = 1:M
            push!(
                corr.args,
                :(MR = vload(SIMD.Vec{_OB_VW,UInt64}, mrb, $(j - 1) * wpbv + ci)),
            )
            push!(
                corr.args,
                :(MI = vload(SIMD.Vec{_OB_VW,UInt64}, mib, $(j - 1) * wpbv + ci)),
            )
            push!(corr.args, :(prc = MR ⊻ CCv))
            push!(corr.args, :(pis = MI ⊻ CSv))
            push!(corr.args, :(pic = MI ⊻ CCv))
            push!(corr.args, :(prs = MR ⊻ CSv))
            for k = 1:NCs[i]
                push!(
                    corr.args,
                    :(CW = vload(SIMD.Vec{_OB_VW,UInt64}, codeb, $(k - 1) * wpbv + ci)),
                )
                push!(corr.args, :($(Symbol("A_$(i)_$(j)_$k")) += count_ones(CW ⊻ prc)))
                push!(corr.args, :($(Symbol("B_$(i)_$(j)_$k")) += count_ones(CW ⊻ pis)))
                push!(corr.args, :($(Symbol("C_$(i)_$(j)_$k")) += count_ones(CW ⊻ pic)))
                push!(corr.args, :($(Symbol("E_$(i)_$(j)_$k")) += count_ones(CW ⊻ prs)))
            end
        end
        push!(b.args, quote
            ci = 1
            while ci <= nwv
                CCv = vload(SIMD.Vec{_OB_VW,UInt64}, cosw, ci)
                CSv = vload(SIMD.Vec{_OB_VW,UInt64}, sinw, ci)
                $corr
                ci += _OB_VW
            end
        end)
        b
    end
    sigs = Expr(:block)
    for i = 1:N
        push!(sigs.args, signal_block(i))
    end

    # Finalize: per signal, a Vector of NCᵢ tap sums Iₓ = 2N − 2(A+B), Qₓ = 2(E − C). Unlike
    # at the noise reference, the multi-signal `_correlate_signals` lets the compiler elide
    # the `Vector` (0 B for a two-signal satellite).
    function tapval(i, k)
        A(j) = s(Symbol("A_$(i)_$(j)_$k"))
        B(j) = s(Symbol("B_$(i)_$(j)_$k"))
        C(j) = s(Symbol("C_$(i)_$(j)_$k"))
        E(j) = s(Symbol("E_$(i)_$(j)_$k"))
        s(sym) = :(Int64(sum($sym)))
        val(j) = :(complex(
            Float64(2 * N_samp - 2 * ($(A(j)) + $(B(j)))),
            Float64(2 * ($(E(j)) - $(C(j)))),
        ))
        M == 1 ? val(1) : :(SVector{$M,ComplexF64}(tuple($([val(j) for j = 1:M]...))))
    end
    finals = Expr(:tuple)
    for i = 1:N
        taps = [tapval(i, k) for k = 1:NCs[i]]
        push!(
            finals.args,
            M == 1 ? :(ComplexF64[$(taps...)]) : :(SVector{M,ComplexF64}[$(taps...)]),
        )
    end

    quote
        $gates
        N_samp = Int(num_samples)
        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        cps_car = carrier_freq / sampling_freq
        carphase = Float64(carrier_phase)
        num_rows = size(signal, 1)
        p_sig = Ptr{Int16}(pointer(signal))
        band = dc.band                              # band-shared measurement sign planes (>1 sat)
        band_shared = !isempty(band.mrband)
        bandstride = cld(num_rows, 64) + 2
        $setup
        maxspan_v = $maxspan

        bufs = _ob_scratch(dc)
        blk = dc.blk
        wpb = cld(blk, 64)
        wpbv = cld(wpb, _OB_VW) * _OB_VW
        pewords = cld(blk + maxspan_v, 64) + (maxspan_v >> 6) + 1 + _OB_VW   # + funnel-shift pad
        length(bufs.extb) < blk + maxspan_v + 64 && resize!(bufs.extb, blk + maxspan_v + 64)
        length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
        length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
        length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
        length(bufs.codeb) < $maxNC * wpbv && resize!(bufs.codeb, $maxNC * wpbv)
        length(bufs.mrb) < $M * wpbv && resize!(bufs.mrb, $M * wpbv)
        length(bufs.mib) < $M * wpbv && resize!(bufs.mib, $M * wpbv)
        extb = bufs.extb;
        peb = bufs.peb;
        sinw = bufs.sinw;
        cosw = bufs.cosw
        codeb = bufs.codeb;
        mrb = bufs.mrb;
        mib = bufs.mib

        blk_off = 0
        @inbounds while blk_off < N_samp
            len = min(blk, N_samp - blk_off)
            nwb = cld(len, 64)
            nwv = cld(nwb, _OB_VW) * _OB_VW
            r = len & 63
            # carrier sign planes — shared across signals, antennas and taps
            generate_carrier_signs!(
                sinw,
                cosw,
                len,
                cps_car;
                phase = carphase + cps_car * blk_off,
            )
            _ob_zeropad!(sinw, 0, nwb, nwv)
            _ob_zeropad!(cosw, 0, nwb, nwv)
            # measurement sign planes — one pack per sat/block, shared across signals
            base = signal_start_sample + blk_off
            if band_shared
                $realign
            else
                $directpack
            end
            # each signal: its own code, correlated against the shared planes
            $sigs
            blk_off += len
        end
        $finals
    end
end

# ── Correlate / plumbing ──────────────────────────────────────────────────────
# The despread primitive on this backend's kernel; see `_despread_one_signal!` in
# downconvert_and_correlate_cpu.jl (the noise reference must use it: popcount counts and
# float sums are on different scales). `code_replica_size` is ignored, the kernel packs
# its own code plane. `use_band_cache` selects the band-packed measurement planes; the
# noise reference passes `false`, as the planes belong to the group whose loop it rides,
# possibly another band.
@inline _despread_one_signal!(
    dc::_OneBitDC,
    correlator,
    samples,
    signal_type,
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
    use_band_cache::Bool = true,
) = update_accumulator(
    correlator,
    get_accumulators(correlator) .+ _onebit_hybrid_blocked!(
        dc,
        samples,
        _ob_num_ants_val(correlator),
        signal_type,
        prn,
        sample_shifts,
        code_phase,
        carrier_phase,
        code_frequency,
        carrier_frequency,
        sampling_frequency,
        start_sample,
        num_samples,
        use_band_cache,
    ),
)

# Several signals per sat: one carrier + measurement downconvert shared through the
# tile-share kernel. Returns the per-signal `(new_correlator, is_integration_completed)`
# tuples.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal,TrackedSignal,Vararg{TrackedSignal}},
    per_signal_completed::Tuple,
    dc::_OneBitDC,
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
        _signal_replica_params(
            head,
            code_doppler,
            code_phase,
            sampling_frequency,
            num_samples_signal,
        )
    end
    correlators = map(s -> s.correlator, signals)
    signal_types = map(s -> s.signal, signals)
    all_shifts = map(p -> p.sample_shifts, params)
    code_phases = map(p -> p.signal_code_phase, params)
    code_freqs = map(p -> p.code_frequency, params)
    new_accs = _onebit_hybrid_blocked_multi!(
        dc,
        signal,
        _ob_num_ants_val(correlators[1]),
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

# The per-sat chunk loop is the generic `_update_tracked_sat_correlator`; the backend
# boundary is `_correlate_signals` above.

@inline function _dc_one_group!(
    g::SignalGroup,
    dc::_OneBitDC,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
    stop_before_partial::Bool,
    samples_unchanged::Bool,
    noise_items::_NoiseWork,
    n_noise::Int,
)
    vals = g.satellites.values
    isempty(vals) && n_noise == 0 && return nothing
    m = measurements[get_band_id(g.band)]
    _check_sample_type(dc, m)
    # Pack the band's measurement planes once for >1 sat (see `OneBitBandBuffers`); for one
    # sat the direct pack is cheaper, so empty the band buffer to select it. Noise items
    # never read the planes, so only the satellite count decides.
    samples = m.samples
    if samples_unchanged
        # Same buffer as the previous call with this dc: the band state is still valid.
    elseif length(vals) > 1
        nsamp = get_num_samples(m)
        M = samples isa AbstractMatrix ? size(samples, 2) : 1
        num_rows = samples isa AbstractMatrix ? size(samples, 1) : length(samples)
        GC.@preserve samples _ob_pack_band!(
            dc.band,
            samples,
            nsamp,
            M,
            Ptr{Int16}(pointer(samples)),
            num_rows,
        )
    else
        resize!(dc.band.mrband, 0)
        resize!(dc.band.miband, 0)
    end
    num_samples = get_num_samples(m)
    chunk_last_sample =
        _chunk_last_sample(chunk_duration, chunk_index, m.sampling_frequency, num_samples)
    _dc_group_loop!(
        dc,
        vals,
        noise_items,
        n_noise,
        m.samples,
        num_samples,
        chunk_last_sample,
        m.sampling_frequency,
        m.intermediate_frequency,
        stop_before_partial,
    )
end

@inline _threading(::OneBitThreadedDownconvertAndCorrelator) = _BatchLoop()

# The public `downconvert_and_correlate(!)` entry points are the backend-agnostic ones of
# downconvert_and_correlate_cpu.jl.

# Reject non-`Complex{Int16}` sample buffers up front (12-bit ADC contract). A hook rather
# than inline in `_dc_one_group!`, so the noise work, resolved before the group loop,
# throws the same error.
@inline _check_sample_type(::_OneBitDC, m) =
    eltype(m.samples) === Complex{Int16} || throw(
        ArgumentError(
            string(
                "OneBitDownconvertAndCorrelator requires `Complex{Int16}` measurement ",
                "samples (12-bit ADC); got element type ",
                eltype(m.samples),
                ". Use a CPU(Threaded)DownconvertAndCorrelator for floating-point samples.",
            ),
        ),
    )
