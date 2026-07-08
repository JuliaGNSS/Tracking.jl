# Two-bit (sign + magnitude) bit-wise downconvert + correlate backend.
#
# The middle point between the one-bit and Int16 backends: like one-bit it packs
# everything into bit planes and correlates with XOR + popcount, but it keeps a second
# MAGNITUDE bit for the measurement AND the carrier, so both are 4-level `{±1, ±3}`
# quantities. The code stays 1-bit (±1), exact for binary modulations (BPSK/BOC/TMBOC).
#
#   value = s·(1 + 2·b)     s = sign bit;  b = magnitude bit    ∈ {±1, ±3}
#
#   measurement:  b ⇔ |component| ≥ `threshold` (ADC counts; ≈1σ → the classic
#                 near-optimal 4-level quantiser, ≈0.55 dB loss vs float against
#                 ≈1.96 dB for 1-bit)
#   carrier:      SinCosLUT v3.3's `generate_carrier_signs_mags!` — sign + magnitude
#                 bit planes straight off the UInt32 NCO, magnitude threshold on the
#                 45° octant boundary (|sin| ≥ √2/2). The weight pair is the consumer's
#                 choice; at that threshold ±1/±3 and ±1/±2 have the SAME correlation
#                 loss (≈0.25 dB vs a float carrier, against ≈0.91 dB for 1-bit), and
#                 ±1/±3 makes the product weight below strictly cheaper. Both bits are
#                 exact functions of the phase accumulator — no table, no rounding.
#
# A product's SIGN is still an XOR of sign bits; only the WEIGHT changes. With both
# factors `(1+2G)(1+2H)` ∈ {1, 3, 9} the weight has the two-plane expansion
#
#   (1+2G)(1+2H) = 1 + 2·(G ⊻ H) + 8·(G & H)  =  1 + 2·X + 8·A
#
# so each of the four product sign planes (S = codeₓ ⊻ sm ⊻ sc; the one-bit A/B/C/E
# planes mr·cos, mi·sin, mi·cos, mr·sin) costs THREE popcounts — one unmasked (the
# one-bit term) and two masked:
#
#   Σ codeₓ·m·c = (N − 2·pc S) + 2·(pc X − 2·pc(S & X)) + 8·(pc A − 2·pc(S & A))
#
# `pc X + 4·pc A` is a per-(antenna, mask-combo) block constant (4 combos: the pairings
# of measurement mag {Gr, Gi} with carrier mag {Hc, Hs}), accumulated once per chunk and
# shared across taps. Per (antenna, tap) the kernel keeps FOUR `Vec{VW,UInt64}` lane
# accumulators — the same count as one-bit — by merging each tap's I = (mr·cos)+(mi·sin)
# and Q = (mi·cos)−(mr·sin) pairs:
#
#   IS += pc S_A + pc S_B                 IF += F_A + F_B      F_L = pc(S_L & X_L)
#   QS += pc S_C − pc S_E                 QF += F_C − F_E            + 4·pc(S_L & A_L)
#
# The Q accumulators SUBTRACT with wrap-around UInt64 lane arithmetic (exact mod 2⁶⁴;
# reinterpreted signed at finalize — |lane totals| ≤ 5·N/8 keep far away from 2⁶³), so no
# plane needs complementing and zero pad words stay contribution-free. Finalize:
#
#   Iₓ = 2N − 2·IS + 2·(GXA_ac + GXA_bs) − 4·IF        GXA = Σ pc X + 4·pc A of the
#   Qₓ =    − 2·QS + 2·(GXA_cc − GXA_es) − 4·QF        plane's (meas, carrier) combo
#
# With all magnitude bits zero this degenerates exactly to the one-bit backend's result.
# Everything else — strip-mined `blk`-sample blocks, pack-once + funnel-shifted code
# planes, band-shared measurement packing (>1 sat), the multi-signal-per-sat
# carrier+measurement tile-share, the dynamic (runtime tap-count) fallback — mirrors the
# one-bit backend (see downconvert_and_correlate_onebit.jl for the machinery rationale).
#
# The 2-bit measurement + 2-bit carrier recover ≈2 dB of correlation SNR over the one-bit
# backend (≈0.8 dB total quantisation loss vs Float32, against ≈2.9 dB), for ≈3 popcounts
# per plane instead of 1. Accumulators finalize to `ComplexF64` (M=1) /
# `SVector{M,ComplexF64}` (M>1); downstream consumers are ratio/normalised, so the coarse
# amplitude scale is immaterial.

import SinCosLUT
import VectorizationBase
using SinCosLUT: generate_carrier_signs_mags!

# Default strip-mine block length (samples); a multiple of 64 so every block packs a
# whole number of UInt64 words.
const _TWOBIT_BLK = 8192

# UInt64 lanes per SIMD popcount step (see the one-bit backend's `_OB_VW`).
const _TB_VW = 8

# ── sign / magnitude mask helpers ─────────────────────────────────────────────
# Pack the sign bit of 64 lanes into a UInt64 (bit j ⇔ lane j < 0). Same as the one-bit
# `_ob_mm`; duplicated so the two backends stay independent files.
@inline _tb_mm(v::SIMD.Vec{64,Int8}) = Base.llvmcall(
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
@inline _tb_mm(v::SIMD.Vec{64,Int16}) = Base.llvmcall(
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
# Magnitude mask: bit j ⇔ |lane j| ≥ thr, via one unsigned-range compare: with
# `off = thr − 1` and `lim = 2·thr − 1` (as bit patterns), `|v| ≤ thr−1 ⇔
# wrap(v + off) ≤ 2·thr−2` unsigned, so the large plane is `wrap(v + off) ≥ᵤ lim`.
# Exact for EVERY Int16 lane (including −32768, where an `abs`-based test wraps)
# and every `thr ∈ [1, 32767]`.
@inline _tb_mag_mm(v::SIMD.Vec{64,Int16}, off::Int16, lim::Int16) = Base.llvmcall(
    (
        """
define i64 @entry(<64 x i16> %v, i16 %off, i16 %lim) #0 {
  %o0 = insertelement <64 x i16> poison, i16 %off, i32 0
  %ov = shufflevector <64 x i16> %o0, <64 x i16> poison, <64 x i32> zeroinitializer
  %l0 = insertelement <64 x i16> poison, i16 %lim, i32 0
  %lv = shufflevector <64 x i16> %l0, <64 x i16> poison, <64 x i32> zeroinitializer
  %a = add <64 x i16> %v, %ov
  %c = icmp uge <64 x i16> %a, %lv
  %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""",
        "entry",
    ),
    UInt64,
    Tuple{NTuple{64,Base.VecElement{Int16}},Int16,Int16},
    v.data,
    off,
    lim,
)
@inline _tb_mag_off(thr::Int16) = thr - Int16(1)
@inline _tb_mag_lim(thr::Int16) = unsafe_trunc(Int16, 2 * Int32(thr) - Int32(1))
@inline _tb_masklast(w::UInt64, r::Int) = r == 0 ? w : w & ((UInt64(1) << r) - UInt64(1))

# ── x86 fast measurement pack (BMI2) ──────────────────────────────────────────
# The generic pack deinterleaves 64 complex lanes with shufflevector (a pile of
# port-limited permutes once legalised). On x86 with BMI2 it is cheaper to keep the data
# INTERLEAVED — extract each 32-lane register's sign / magnitude masks directly (the
# bits come out re/im-interleaved) — and deinterleave the BITS with scalar `pext`
# (1/cycle, off the vector ports): two pexts per plane word.
@inline _tb_mm32(v::SIMD.Vec{32,Int16}) = Base.llvmcall(
    (
        """
define i32 @entry(<32 x i16> %v) #0 { %c = icmp slt <32 x i16> %v, zeroinitializer
  %m = bitcast <32 x i1> %c to i32 ret i32 %m } attributes #0={alwaysinline}""",
        "entry",
    ),
    UInt32,
    Tuple{NTuple{32,Base.VecElement{Int16}}},
    v.data,
)
@inline _tb_mag_mm32(v::SIMD.Vec{32,Int16}, off::Int16, lim::Int16) = Base.llvmcall(
    (
        """
define i32 @entry(<32 x i16> %v, i16 %off, i16 %lim) #0 {
  %o0 = insertelement <32 x i16> poison, i16 %off, i32 0
  %ov = shufflevector <32 x i16> %o0, <32 x i16> poison, <32 x i32> zeroinitializer
  %l0 = insertelement <32 x i16> poison, i16 %lim, i32 0
  %lv = shufflevector <32 x i16> %l0, <32 x i16> poison, <32 x i32> zeroinitializer
  %a = add <32 x i16> %v, %ov
  %c = icmp uge <32 x i16> %a, %lv
  %m = bitcast <32 x i1> %c to i32 ret i32 %m } attributes #0={alwaysinline}""",
        "entry",
    ),
    UInt32,
    Tuple{NTuple{32,Base.VecElement{Int16}},Int16,Int16},
    v.data,
    off,
    lim,
)
@inline _tb_pext(v::UInt64, m::UInt64) = Base.llvmcall(
    (
        """
declare i64 @llvm.x86.bmi.pext.64(i64, i64)
define i64 @entry(i64 %v, i64 %m) #0 {
  %r = call i64 @llvm.x86.bmi.pext.64(i64 %v, i64 %m)
  ret i64 %r } attributes #0={alwaysinline}""",
        "entry",
    ),
    UInt64,
    Tuple{UInt64,UInt64},
    v,
    m,
)
const _TB_RE_BITS = 0x5555555555555555          # even (re) bit positions
const _TB_IM_BITS = 0xaaaaaaaaaaaaaaaa          # odd (im) bit positions
# `pext` needs BMI2; VectorizationBase resolves the feature at precompile time, so the
# `if _TB_FAST_PACK` below folds away and non-x86 / pre-BMI2 hosts never compile it.
const _TB_FAST_PACK = Bool(VectorizationBase.has_feature(Val(:x86_64_bmi2)))

# One 64-complex-sample word: (sign(re), sign(im), |re| ≥ thr, |im| ≥ thr) plane words,
# from four interleaved 32-lane loads.
@inline function _tb_pack_word_fast(
    p_sig::Ptr{Int16},
    byte_off::Int,
    off::Int16,
    lim::Int16,
)
    v0 = SIMD.vload(SIMD.Vec{32,Int16}, Ptr{Int16}(p_sig + byte_off))
    v1 = SIMD.vload(SIMD.Vec{32,Int16}, Ptr{Int16}(p_sig + byte_off + 64))
    v2 = SIMD.vload(SIMD.Vec{32,Int16}, Ptr{Int16}(p_sig + byte_off + 128))
    v3 = SIMD.vload(SIMD.Vec{32,Int16}, Ptr{Int16}(p_sig + byte_off + 192))
    s01 = UInt64(_tb_mm32(v0)) | (UInt64(_tb_mm32(v1)) << 32)
    s23 = UInt64(_tb_mm32(v2)) | (UInt64(_tb_mm32(v3)) << 32)
    g01 = UInt64(_tb_mag_mm32(v0, off, lim)) | (UInt64(_tb_mag_mm32(v1, off, lim)) << 32)
    g23 = UInt64(_tb_mag_mm32(v2, off, lim)) | (UInt64(_tb_mag_mm32(v3, off, lim)) << 32)
    (
        _tb_pext(s01, _TB_RE_BITS) | (_tb_pext(s23, _TB_RE_BITS) << 32),
        _tb_pext(s01, _TB_IM_BITS) | (_tb_pext(s23, _TB_IM_BITS) << 32),
        _tb_pext(g01, _TB_RE_BITS) | (_tb_pext(g23, _TB_RE_BITS) << 32),
        _tb_pext(g01, _TB_IM_BITS) | (_tb_pext(g23, _TB_IM_BITS) << 32),
    )
end

# ── code plane helpers (code is 1-bit; identical to the one-bit backend) ──────
# Pack tap code sign plane: bit for output sample n = sign(extb[byteoff+n]). `extb` is
# padded ≥ byteoff+nwb·64 so the last (masked) word can be read whole.
@inline function _tb_pack_code!(
    codeb::Vector{UInt64},
    koff::Int,
    extb::Vector{Int8},
    byteoff::Int,
    nwb::Int,
    r::Int,
)
    @inbounds for w = 0:(nwb-1)
        wd = _tb_mm(SIMD.vload(SIMD.Vec{64,Int8}, extb, byteoff + (w << 6) + 1))
        codeb[koff+w+1] = w == nwb - 1 ? _tb_masklast(wd, r) : wd
    end
    nothing
end

# Derive each tap plane from the prompt-extended plane by a funnel bit-shift, for ANY
# `off ≥ 0` (whole-word shift `wsh = off ÷ 64` + sub-word shift). The caller must size /
# zero `pe` so `pe[nwb+wsh+1]` is a valid (zero) pad word (see the per-kernel `pepad`).
@inline function _tb_shift_plane!(
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

# Mask the last valid word to `r` bits (if partial) and zero the VW-pad words [nwb+1, nwv],
# so a whole-VW-chunk correlate sees zeros (which contribute nothing) past the real samples.
@inline function _tb_finish!(buf::Vector{UInt64}, off::Int, nwb::Int, nwv::Int, r::Int)
    @inbounds (r != 0) && (buf[off+nwb] &= (UInt64(1) << r) - UInt64(1))
    @inbounds for w = (nwb+1):nwv
        buf[off+w] = 0
    end
    nothing
end
# The packer already masked the last word; just zero the VW pad.
@inline function _tb_zeropad!(buf::Vector{UInt64}, off::Int, nwb::Int, nwv::Int)
    @inbounds for w = (nwb+1):nwv
        buf[off+w] = 0
    end
    nothing
end

# ── measurement packing (sign + magnitude planes) ─────────────────────────────
# Per-sat direct pack (single-sat group): sign(real)→mrb, |real| ≥ thr→gmrb, likewise
# imag → mib/gmib; full 64-sample words via a deinterleaving vector load, scalar tail.
@inline function _tb_pack_meas!(
    mrb::Vector{UInt64},
    mib::Vector{UInt64},
    gmrb::Vector{UInt64},
    gmib::Vector{UInt64},
    joff::Int,
    signal,
    j::Int,
    p_sig::Ptr{Int16},
    colbase::Int,
    base::Int,
    len::Int,
    nwb::Int,
    r::Int,
    thr::Int16,
)
    off = _tb_mag_off(thr)
    lim = _tb_mag_lim(thr)
    full = fld(len, 64)
    @inbounds for w = 0:(full-1)
        byte_off = (colbase + base + (w << 6) - 1) * 2 * sizeof(Int16)
        if _TB_FAST_PACK
            wr, wi, wgr, wgi = _tb_pack_word_fast(p_sig, byte_off, off, lim)
            mrb[joff+w+1] = wr
            mib[joff+w+1] = wi
            gmrb[joff+w+1] = wgr
            gmib[joff+w+1] = wgi
        else
            re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
            mrb[joff+w+1] = _tb_mm(re)
            mib[joff+w+1] = _tb_mm(im)
            gmrb[joff+w+1] = _tb_mag_mm(re, off, lim)
            gmib[joff+w+1] = _tb_mag_mm(im, off, lim)
        end
    end
    if r != 0
        wr = zero(UInt64);
        wi = zero(UInt64)
        wgr = zero(UInt64);
        wgi = zero(UInt64)
        @inbounds for i = 0:(r-1)
            sig = signal[base+(full<<6)+i, j]
            re = real(sig);
            im = imag(sig)
            (re < 0) && (wr |= UInt64(1) << i)
            (im < 0) && (wi |= UInt64(1) << i)
            (abs(Int(re)) >= thr) && (wgr |= UInt64(1) << i)
            (abs(Int(im)) >= thr) && (wgi |= UInt64(1) << i)
        end
        mrb[joff+full+1] = wr;
        mib[joff+full+1] = wi
        gmrb[joff+full+1] = wgr;
        gmib[joff+full+1] = wgi
    end
    nothing
end

# Band-shared measurement planes (see the one-bit `OneBitBandBuffers`): the sign AND
# magnitude planes are identical for every satellite on a band, so pack them once per
# group and funnel-realign per sat.
struct TwoBitBandBuffers
    mrband::Vector{UInt64}
    miband::Vector{UInt64}
    gmrband::Vector{UInt64}
    gmiband::Vector{UInt64}
end
TwoBitBandBuffers() = TwoBitBandBuffers(UInt64[], UInt64[], UInt64[], UInt64[])

# The full-word body of the band pack: words `wlo:whi` (0-based) of antenna `j`'s planes.
@inline function _tb_pack_band_words!(
    band::TwoBitBandBuffers,
    p_sig::Ptr{Int16},
    off::Int,
    colbase::Int,
    wlo::Int,
    whi::Int,
    moff::Int16,
    mlim::Int16,
)
    @inbounds for w = wlo:whi
        byte_off = (colbase + w * 64) * 2 * sizeof(Int16)
        if _TB_FAST_PACK
            wr, wi, wgr, wgi = _tb_pack_word_fast(p_sig, byte_off, moff, mlim)
            band.mrband[off+w+1] = wr
            band.miband[off+w+1] = wi
            band.gmrband[off+w+1] = wgr
            band.gmiband[off+w+1] = wgi
        else
            re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
            band.mrband[off+w+1] = _tb_mm(re)
            band.miband[off+w+1] = _tb_mm(im)
            band.gmrband[off+w+1] = _tb_mag_mm(re, moff, mlim)
            band.gmiband[off+w+1] = _tb_mag_mm(im, moff, mlim)
        end
    end
    nothing
end

# Below this many full words the whole pack runs serially — a `@batch` spawn costs more
# than it saves on a small capture.
const _TB_BAND_PAR_MIN = 512

# Pack the band's `M` measurement sign + magnitude planes over all `num_samples` samples.
# Plane `j` occupies `stride` words at offset `(j-1)*stride`; the last two words of each
# plane are the zero funnel pad. This pack runs ONCE per group on the critical path
# before the per-sat `@batch`, so for the threaded backend (`threaded = true`) and a
# large enough capture the full-word body is itself split across threads (each word is
# independent: it reads its own 64 samples and writes its own plane words).
function _tb_pack_band!(
    band::TwoBitBandBuffers,
    signal,
    num_samples::Int,
    M::Int,
    p_sig::Ptr{Int16},
    num_rows::Int,
    thr::Int16,
    threaded::Bool,
)
    nwb = cld(num_samples, 64)
    stride = nwb + 2                                # 2 funnel-pad words (realign reads past nwb)
    length(band.mrband) < M * stride && resize!(band.mrband, M * stride)
    length(band.miband) < M * stride && resize!(band.miband, M * stride)
    length(band.gmrband) < M * stride && resize!(band.gmrband, M * stride)
    length(band.gmiband) < M * stride && resize!(band.gmiband, M * stride)
    moff = _tb_mag_off(thr)
    mlim = _tb_mag_lim(thr)
    r = num_samples & 63
    full = fld(num_samples, 64)
    if threaded && M * full >= _TB_BAND_PAR_MIN
        nch = min(Threads.nthreads(), cld(M * full, _TB_BAND_PAR_MIN))
        # Bundle the loop parameters into ONE isbits capture: Polyester passes captured
        # variables through a fixed-size task buffer, and past a handful of captures it
        # falls back to a true-closure `@cfunction` — unsupported on aarch64 (NEON).
        prm = (;
            per = cld(full, nch),
            full,
            M,
            stride,
            num_rows,
            moff,
            mlim,
        )
        @batch for c = 1:nch
            wlo = (c - 1) * prm.per
            whi = min(c * prm.per, prm.full) - 1
            for j = 1:prm.M
                _tb_pack_band_words!(
                    band,
                    p_sig,
                    (j - 1) * prm.stride,
                    (j - 1) * prm.num_rows,
                    wlo,
                    whi,
                    prm.moff,
                    prm.mlim,
                )
            end
        end
    else
        for j = 1:M
            _tb_pack_band_words!(
                band,
                p_sig,
                (j - 1) * stride,
                (j - 1) * num_rows,
                0,
                full - 1,
                moff,
                mlim,
            )
        end
    end
    @inbounds for j = 1:M
        off = (j - 1) * stride
        if r != 0
            wr = zero(UInt64);
            wi = zero(UInt64)
            wgr = zero(UInt64);
            wgi = zero(UInt64)
            for i = 0:(r-1)
                sig = signal[full*64+i+1, j]
                re = real(sig);
                im = imag(sig)
                (re < 0) && (wr |= UInt64(1) << i)
                (im < 0) && (wi |= UInt64(1) << i)
                (abs(Int(re)) >= thr) && (wgr |= UInt64(1) << i)
                (abs(Int(im)) >= thr) && (wgi |= UInt64(1) << i)
            end
            band.mrband[off+full+1] = wr;
            band.miband[off+full+1] = wi
            band.gmrband[off+full+1] = wgr;
            band.gmiband[off+full+1] = wgi
        end
        band.mrband[off+nwb+1] = 0;
        band.miband[off+nwb+1] = 0     # funnel pad
        band.gmrband[off+nwb+1] = 0;
        band.gmiband[off+nwb+1] = 0
        band.mrband[off+nwb+2] = 0;
        band.miband[off+nwb+2] = 0
        band.gmrband[off+nwb+2] = 0;
        band.gmiband[off+nwb+2] = 0
    end
    nothing
end

# Funnel-realign antenna `j`'s block (`len` samples starting at 1-based absolute sample
# `base`) from the shared band planes into the thread-local block planes at word-offset
# `moff`.
@inline function _tb_realign_meas!(
    mrb::Vector{UInt64},
    mib::Vector{UInt64},
    gmrb::Vector{UInt64},
    gmib::Vector{UInt64},
    moff::Int,
    band::TwoBitBandBuffers,
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
            gmrb[moff+w] = band.gmrband[jbase+sw+w]
            gmib[moff+w] = band.gmiband[jbase+sw+w]
        end
    else
        hi = 64 - bit
        for w = 1:nwb
            mrb[moff+w] =
                (band.mrband[jbase+sw+w] >> bit) | (band.mrband[jbase+sw+w+1] << hi)
            mib[moff+w] =
                (band.miband[jbase+sw+w] >> bit) | (band.miband[jbase+sw+w+1] << hi)
            gmrb[moff+w] =
                (band.gmrband[jbase+sw+w] >> bit) | (band.gmrband[jbase+sw+w+1] << hi)
            gmib[moff+w] =
                (band.gmiband[jbase+sw+w] >> bit) | (band.gmiband[jbase+sw+w+1] << hi)
        end
    end
    nothing
end

# Scratch: the one-bit planes plus the measurement magnitude planes (`gmrb`/`gmib`) and
# the carrier magnitude planes (`hsinw`/`hcosw`). Grown lazily and reused, so a hoisted
# backend is allocation-free in steady state. Plane word-strides are padded to a multiple
# of `_TB_VW` (`wpbv`).
struct TwoBitScratchBuffers
    extb::Vector{Int8}
    peb::Vector{UInt64}
    sinw::Vector{UInt64}
    cosw::Vector{UInt64}
    hsinw::Vector{UInt64}
    hcosw::Vector{UInt64}
    codeb::Vector{UInt64}
    mrb::Vector{UInt64}
    mib::Vector{UInt64}
    gmrb::Vector{UInt64}
    gmib::Vector{UInt64}
end
TwoBitScratchBuffers() = TwoBitScratchBuffers(
    Int8[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
)

"""
$(SIGNATURES)

Two-bit (sign + magnitude) bit-wise CPU downconvert + correlate backend
(single-threaded). Opt-in alternative to [`CPUDownconvertAndCorrelator`] for
`Complex{Int16}` sample buffers, sitting between the
[`OneBitDownconvertAndCorrelator`](@ref) and
[`Int16DownconvertAndCorrelator`](@ref): measurement and carrier are 4-level
`{±1, ±3}` (2-bit sign + magnitude, the carrier bit planes straight off
SinCosLUT's NCO), the code stays 1-bit, correlated with XOR + masked popcount.
Recovers ≈2 dB of correlation SNR over the one-bit backend (≈0.8 dB total
quantisation loss vs Float32) while keeping the bit-wise speed advantage over
the integer/float backends. Construct **once outside** the `track!` loop and
pass it via the `downconvert_and_correlator` keyword for an allocation-free
steady state.

`threshold` is the measurement magnitude split point in ADC counts
(`|component| ≥ threshold` ⇒ the `±3` level); set it near 1σ of your front
end's input for the near-optimal 4-level quantiser. Downstream consumers are
ratio-normalised, so the absolute scale is immaterial.
"""
struct TwoBitDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::TwoBitScratchBuffers
    band::TwoBitBandBuffers
    blk::Int
    threshold::Int16
end

"""
$(SIGNATURES)

Multi-threaded two-bit bit-wise backend. One `TwoBitScratchBuffers` per thread
(indexed by `Threads.threadid()` inside `@batch`). See
[`TwoBitDownconvertAndCorrelator`](@ref).
"""
struct TwoBitThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{TwoBitScratchBuffers}
    band::TwoBitBandBuffers
    blk::Int
    threshold::Int16
end

function _tb_check_threshold(threshold::Integer)
    (1 <= threshold <= typemax(Int16)) || throw(
        ArgumentError(
            "threshold must be in [1, $(typemax(Int16))] ADC counts; got $threshold",
        ),
    )
    Int16(threshold)
end

TwoBitDownconvertAndCorrelator(; threshold::Integer = 512, blk::Integer = _TWOBIT_BLK) =
    TwoBitDownconvertAndCorrelator(
        TwoBitScratchBuffers(),
        TwoBitBandBuffers(),
        Int(blk),
        _tb_check_threshold(threshold),
    )

TwoBitThreadedDownconvertAndCorrelator(;
    threshold::Integer = 512,
    blk::Integer = _TWOBIT_BLK,
) = TwoBitThreadedDownconvertAndCorrelator(
    [TwoBitScratchBuffers() for _ = 1:Threads.maxthreadid()],
    TwoBitBandBuffers(),
    Int(blk),
    _tb_check_threshold(threshold),
)

const _TwoBitDC =
    Union{TwoBitDownconvertAndCorrelator,TwoBitThreadedDownconvertAndCorrelator}

@inline _tb_scratch(dc::TwoBitDownconvertAndCorrelator) = dc.buffers
@inline _tb_scratch(dc::TwoBitThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

@inline _tb_num_ants_val(::AbstractCorrelator{M}) where {M} = NumAnts{M}()

# The bit-wise CBOC gate (see the one-bit backend for the rationale: the code plane keeps
# only the sign, so a multi-level amplitude-carrying code cannot be represented).
@inline function _tb_check_modulation(signal_type)
    get_modulation(signal_type) isa GNSSSignals.CBOC && throw(
        ArgumentError(
            string(
                "TwoBitDownconvertAndCorrelator supports binary (±1) codes only ",
                "(BPSK, BOC, TMBOC); got ",
                typeof(signal_type),
                " with ",
                typeof(get_modulation(signal_type)),
                " modulation. Bit-wise correlation keeps only the code sign, so it ",
                "cannot represent CBOC — a multi-level, amplitude-carrying code.",
            ),
        ),
    )
    nothing
end

# ── The two-bit hybrid-blocked kernel ─────────────────────────────────────────
# Returns this integration's correlation contribution: `Vector` of `NC` complex sums (one
# per tap) — `ComplexF64` (M=1) or `SVector{M,ComplexF64}` (M>1) — to be added to the
# correlator's running accumulators. `@generated` over (NC, M): the per-(antenna, tap)
# IS/IF/QS/QF and per-antenna GXA accumulators live in named locals and unroll.
@generated function _twobit_hybrid_blocked!(
    dc::_TwoBitDC,
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
    zc = :(zero(SIMD.Vec{_TB_VW,UInt64}))

    # Accumulator declarations: IS/IF (→ I) and QS/QF (→ Q, wrap-around signed) per
    # (antenna j, tap k), plus the per-antenna mask-total accumulators GXA_(combo)_j.
    init = Expr(:block)
    for k = 1:NC
        push!(init.args, :($(Symbol("off_$k")) = Int(sample_shifts[$k]) - min_shift))
    end
    for j = 1:M
        for c in ("ac", "bs", "cc", "es")
            push!(init.args, :($(Symbol("GXA_$(c)_$j")) = $zc))
        end
        for k = 1:NC, s in ("IS", "IF", "QS", "QF")
            push!(init.args, :($(Symbol("$(s)_$(j)_$k")) = $zc))
        end
    end

    # Per-VW-chunk unrolled accumulate: per antenna the shared sample⊻carrier sign
    # products and the four (measurement mag, carrier mag) mask combos, then per tap the
    # three popcounts per plane merged into the four I/Q accumulators.
    body = Expr(:block)
    for j = 1:M
        push!(body.args, :(MR = vload(SIMD.Vec{_TB_VW,UInt64}, mrb, $(j - 1) * wpbv + i)))
        push!(body.args, :(MI = vload(SIMD.Vec{_TB_VW,UInt64}, mib, $(j - 1) * wpbv + i)))
        push!(body.args, :(GR = vload(SIMD.Vec{_TB_VW,UInt64}, gmrb, $(j - 1) * wpbv + i)))
        push!(body.args, :(GI = vload(SIMD.Vec{_TB_VW,UInt64}, gmib, $(j - 1) * wpbv + i)))
        push!(body.args, :(prc = MR ⊻ CCv))
        push!(body.args, :(pis = MI ⊻ CSv))
        push!(body.args, :(pic = MI ⊻ CCv))
        push!(body.args, :(prs = MR ⊻ CSv))
        push!(body.args, :(Xac = GR ⊻ HCv))
        push!(body.args, :(Aac = GR & HCv))
        push!(body.args, :(Xbs = GI ⊻ HSv))
        push!(body.args, :(Abs = GI & HSv))
        push!(body.args, :(Xcc = GI ⊻ HCv))
        push!(body.args, :(Acc = GI & HCv))
        push!(body.args, :(Xes = GR ⊻ HSv))
        push!(body.args, :(Aes = GR & HSv))
        for (c, X, A) in
            (("ac", :Xac, :Aac), ("bs", :Xbs, :Abs), ("cc", :Xcc, :Acc), ("es", :Xes, :Aes))
            push!(
                body.args,
                :($(Symbol("GXA_$(c)_$j")) += count_ones($X) + (count_ones($A) << 2)),
            )
        end
        for k = 1:NC
            push!(
                body.args,
                :(CW = vload(SIMD.Vec{_TB_VW,UInt64}, codeb, $(k - 1) * wpbv + i)),
            )
            push!(body.args, :(SA = CW ⊻ prc))
            push!(body.args, :(SB = CW ⊻ pis))
            push!(body.args, :(SC = CW ⊻ pic))
            push!(body.args, :(SE = CW ⊻ prs))
            push!(body.args, :($(Symbol("IS_$(j)_$k")) += count_ones(SA) + count_ones(SB)))
            push!(
                body.args,
                :(
                    $(Symbol("IF_$(j)_$k")) +=
                        (count_ones(SA & Xac) + (count_ones(SA & Aac) << 2)) +
                        (count_ones(SB & Xbs) + (count_ones(SB & Abs) << 2))
                ),
            )
            push!(body.args, :($(Symbol("QS_$(j)_$k")) += count_ones(SC) - count_ones(SE)))
            push!(
                body.args,
                :(
                    $(Symbol("QF_$(j)_$k")) +=
                        (count_ones(SC & Xcc) + (count_ones(SC & Acc) << 2)) -
                        (count_ones(SE & Xes) + (count_ones(SE & Aes) << 2))
                ),
            )
        end
    end

    # Finalize: reduce lanes (wrap-around exact, reinterpret signed), then
    #   Iₓ = 2N − 2·IS + 2·(GXA_ac + GXA_bs) − 4·IF
    #   Qₓ =    − 2·QS + 2·(GXA_cc − GXA_es) − 4·QF
    s(sym) = :(sum($sym) % Int64)
    function tapval(k)
        function cplx(j)
            :(complex(
                Float64(
                    2 * N - 2 * $(s(Symbol("IS_$(j)_$k"))) +
                    2 * ($(s(Symbol("GXA_ac_$j"))) + $(s(Symbol("GXA_bs_$j")))) -
                    4 * $(s(Symbol("IF_$(j)_$k"))),
                ),
                Float64(
                    -2 * $(s(Symbol("QS_$(j)_$k"))) +
                    2 * ($(s(Symbol("GXA_cc_$j"))) - $(s(Symbol("GXA_es_$j")))) -
                    4 * $(s(Symbol("QF_$(j)_$k"))),
                ),
            ))
        end
        if M == 1
            cplx(1)
        else
            :(SVector{$M,ComplexF64}(tuple($([cplx(j) for j = 1:M]...))))
        end
    end
    ret =
        M == 1 ? :(ComplexF64[$([tapval(k) for k = 1:NC]...)]) :
        :(SVector{M,ComplexF64}[$([tapval(k) for k = 1:NC]...)])

    quote
        _tb_check_modulation(signal_type)
        N = Int(num_samples)
        min_shift = minimum(sample_shifts)
        span = maximum(sample_shifts) - min_shift
        thr = dc.threshold
        $init

        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        cps_car = carrier_freq / sampling_freq
        carphase = Float64(carrier_phase)
        cps_code = Float64(upreferred(code_frequency / Hz)) / sampling_freq
        code_phase0 = Float64(code_phase)
        num_rows = size(signal, 1)
        p_sig = Ptr{Int16}(pointer(signal))
        band = dc.band                              # band-shared measurement planes (>1 sat)
        band_shared = !isempty(band.mrband)         # pre-packed once per group when it pays off
        bandstride = cld(num_rows, 64) + 2          # per-antenna plane stride in `band`

        bufs = _tb_scratch(dc)
        blk = dc.blk
        wpb = cld(blk, 64)                          # words per block plane
        wpbv = cld(wpb, _TB_VW) * _TB_VW            # …padded to a whole number of VW chunks
        pepad = (span >> 6) + 1                     # zero pad words the funnel shift reads past nwe
        pewords = cld(blk + span, 64) + pepad + _TB_VW   # extended prompt plane + funnel pad
        length(bufs.extb) < blk + span + 64 && resize!(bufs.extb, blk + span + 64)
        length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
        length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
        length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
        length(bufs.hsinw) < wpbv && resize!(bufs.hsinw, wpbv)
        length(bufs.hcosw) < wpbv && resize!(bufs.hcosw, wpbv)
        length(bufs.codeb) < $NC * wpbv && resize!(bufs.codeb, $NC * wpbv)
        length(bufs.mrb) < $M * wpbv && resize!(bufs.mrb, $M * wpbv)
        length(bufs.mib) < $M * wpbv && resize!(bufs.mib, $M * wpbv)
        length(bufs.gmrb) < $M * wpbv && resize!(bufs.gmrb, $M * wpbv)
        length(bufs.gmib) < $M * wpbv && resize!(bufs.gmib, $M * wpbv)
        extb = bufs.extb;
        peb = bufs.peb;
        sinw = bufs.sinw;
        cosw = bufs.cosw
        hsinw = bufs.hsinw;
        hcosw = bufs.hcosw
        codeb = bufs.codeb;
        mrb = bufs.mrb;
        mib = bufs.mib
        gmrb = bufs.gmrb;
        gmib = bufs.gmib

        blk_off = 0
        @inbounds while blk_off < N
            len = min(blk, N - blk_off)
            nwb = cld(len, 64)
            nwv = cld(nwb, _TB_VW) * _TB_VW
            r = len & 63
            # code (Int8 ±1) for samples [min_shift, len+span); pack the extended plane
            # ONCE, then derive each tap by a funnel bit-shift.
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
            _tb_pack_code!(peb, 0, extb, 0, nwe, (len + span) & 63)
            for _pw = 1:pepad                           # zero the words the funnel shift reads past nwe
                peb[nwe+_pw] = 0
            end
            $(Expr(
                :block,
                (
                    quote
                        _tb_shift_plane!(
                            codeb,
                            $(k - 1) * wpbv,
                            peb,
                            $(Symbol("off_$k")),
                            nwb,
                        )
                        _tb_finish!(codeb, $(k - 1) * wpbv, nwb, nwv, r)
                    end for k = 1:NC
                )...,
            ))
            # carrier sign + magnitude planes (shared across antennas & taps) — the four
            # bit planes straight off SinCosLUT's 2-bit NCO.
            generate_carrier_signs_mags!(
                sinw,
                cosw,
                hsinw,
                hcosw,
                len,
                cps_car;
                phase = carphase + cps_car * blk_off,
            )
            _tb_zeropad!(sinw, 0, nwb, nwv)
            _tb_zeropad!(cosw, 0, nwb, nwv)
            _tb_zeropad!(hsinw, 0, nwb, nwv)
            _tb_zeropad!(hcosw, 0, nwb, nwv)
            # measurement sign + magnitude planes, per antenna. >1 sat: funnel-realign
            # from the band-shared planes. 1 sat: pack directly (no realign copy).
            base = signal_start_sample + blk_off
            if band_shared
                $(Expr(
                    :block,
                    (
                        quote
                            _tb_realign_meas!(
                                mrb,
                                mib,
                                gmrb,
                                gmib,
                                $(j - 1) * wpbv,
                                band,
                                $(j - 1) * bandstride,
                                base,
                                nwb,
                            )
                            _tb_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                            _tb_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
                            _tb_finish!(gmrb, $(j - 1) * wpbv, nwb, nwv, r)
                            _tb_finish!(gmib, $(j - 1) * wpbv, nwb, nwv, r)
                        end for j = 1:M
                    )...,
                ))
            else
                $(Expr(
                    :block,
                    (
                        quote
                            _tb_pack_meas!(
                                mrb,
                                mib,
                                gmrb,
                                gmib,
                                $(j - 1) * wpbv,
                                signal,
                                $j,
                                p_sig,
                                $(j - 1) * num_rows,
                                base,
                                len,
                                nwb,
                                r,
                                thr,
                            )
                            _tb_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                            _tb_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
                            _tb_finish!(gmrb, $(j - 1) * wpbv, nwb, nwv, r)
                            _tb_finish!(gmib, $(j - 1) * wpbv, nwb, nwv, r)
                        end for j = 1:M
                    )...,
                ))
            end
            # XOR + masked popcount accumulate over whole VW chunks
            i = 1
            while i <= nwv
                CCv = vload(SIMD.Vec{_TB_VW,UInt64}, cosw, i)
                CSv = vload(SIMD.Vec{_TB_VW,UInt64}, sinw, i)
                HCv = vload(SIMD.Vec{_TB_VW,UInt64}, hcosw, i)
                HSv = vload(SIMD.Vec{_TB_VW,UInt64}, hsinw, i)
                $body
                i += _TB_VW
            end
            blk_off += len
        end
        $ret
    end
end

# Dynamic (runtime tap count) fallback: correlators whose sample shifts are a
# runtime-sized `AbstractVector`. Mirrors the one-bit dynamic method: same block
# pipeline, but loops taps/antennas at runtime and horizontally sums each popcount chunk
# into per-(antenna, tap) Int64 totals. Not the hot path, but bit-identical to the
# `@generated` kernel. `SVector` shifts dispatch to the `@generated` method above, so
# EPL/VEPL keep the fast unrolled path.
function _twobit_hybrid_blocked!(
    dc::_TwoBitDC,
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
    _tb_check_modulation(signal_type)
    NC = length(sample_shifts)
    N = Int(num_samples)
    min_shift = Int(minimum(sample_shifts))
    span = Int(maximum(sample_shifts)) - min_shift
    offs = [Int(sample_shifts[k]) - min_shift for k = 1:NC]
    thr = dc.threshold

    sampling_freq = Float64(upreferred(sampling_frequency / Hz))
    carrier_freq = Float64(upreferred(carrier_frequency / Hz))
    cps_car = carrier_freq / sampling_freq
    carphase = Float64(carrier_phase)
    cps_code = Float64(upreferred(code_frequency / Hz)) / sampling_freq
    code_phase0 = Float64(code_phase)
    num_rows = size(signal, 1)
    p_sig = Ptr{Int16}(pointer(signal))
    band = dc.band
    band_shared = !isempty(band.mrband)
    bandstride = cld(num_rows, 64) + 2

    bufs = _tb_scratch(dc)
    blk = dc.blk
    wpb = cld(blk, 64)
    wpbv = cld(wpb, _TB_VW) * _TB_VW
    pepad = (span >> 6) + 1
    pewords = cld(blk + span, 64) + pepad + _TB_VW
    length(bufs.extb) < blk + span + 64 && resize!(bufs.extb, blk + span + 64)
    length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
    length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
    length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
    length(bufs.hsinw) < wpbv && resize!(bufs.hsinw, wpbv)
    length(bufs.hcosw) < wpbv && resize!(bufs.hcosw, wpbv)
    length(bufs.codeb) < NC * wpbv && resize!(bufs.codeb, NC * wpbv)
    length(bufs.mrb) < M * wpbv && resize!(bufs.mrb, M * wpbv)
    length(bufs.mib) < M * wpbv && resize!(bufs.mib, M * wpbv)
    length(bufs.gmrb) < M * wpbv && resize!(bufs.gmrb, M * wpbv)
    length(bufs.gmib) < M * wpbv && resize!(bufs.gmib, M * wpbv)
    extb = bufs.extb
    peb = bufs.peb
    sinw = bufs.sinw
    cosw = bufs.cosw
    hsinw = bufs.hsinw
    hcosw = bufs.hcosw
    codeb = bufs.codeb
    mrb = bufs.mrb
    mib = bufs.mib
    gmrb = bufs.gmrb
    gmib = bufs.gmib

    IS = zeros(Int64, M, NC)
    IF = zeros(Int64, M, NC)
    QS = zeros(Int64, M, NC)
    QF = zeros(Int64, M, NC)
    GXA = zeros(Int64, M, 4)                       # combos: 1 = ac, 2 = bs, 3 = cc, 4 = es

    blk_off = 0
    @inbounds while blk_off < N
        len = min(blk, N - blk_off)
        nwb = cld(len, 64)
        nwv = cld(nwb, _TB_VW) * _TB_VW
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
        _tb_pack_code!(peb, 0, extb, 0, nwe, (len + span) & 63)
        for _pw = 1:pepad
            peb[nwe+_pw] = 0
        end
        for k = 1:NC
            _tb_shift_plane!(codeb, (k - 1) * wpbv, peb, offs[k], nwb)
            _tb_finish!(codeb, (k - 1) * wpbv, nwb, nwv, r)
        end
        generate_carrier_signs_mags!(
            sinw,
            cosw,
            hsinw,
            hcosw,
            len,
            cps_car;
            phase = carphase + cps_car * blk_off,
        )
        _tb_zeropad!(sinw, 0, nwb, nwv)
        _tb_zeropad!(cosw, 0, nwb, nwv)
        _tb_zeropad!(hsinw, 0, nwb, nwv)
        _tb_zeropad!(hcosw, 0, nwb, nwv)
        base = signal_start_sample + blk_off
        if band_shared
            for j = 1:M
                _tb_realign_meas!(
                    mrb,
                    mib,
                    gmrb,
                    gmib,
                    (j - 1) * wpbv,
                    band,
                    (j - 1) * bandstride,
                    base,
                    nwb,
                )
                _tb_finish!(mrb, (j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(mib, (j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmrb, (j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmib, (j - 1) * wpbv, nwb, nwv, r)
            end
        else
            for j = 1:M
                _tb_pack_meas!(
                    mrb,
                    mib,
                    gmrb,
                    gmib,
                    (j - 1) * wpbv,
                    signal,
                    j,
                    p_sig,
                    (j - 1) * num_rows,
                    base,
                    len,
                    nwb,
                    r,
                    thr,
                )
                _tb_finish!(mrb, (j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(mib, (j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmrb, (j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmib, (j - 1) * wpbv, nwb, nwv, r)
            end
        end
        ci = 1
        while ci <= nwv
            CCv = vload(SIMD.Vec{_TB_VW,UInt64}, cosw, ci)
            CSv = vload(SIMD.Vec{_TB_VW,UInt64}, sinw, ci)
            HCv = vload(SIMD.Vec{_TB_VW,UInt64}, hcosw, ci)
            HSv = vload(SIMD.Vec{_TB_VW,UInt64}, hsinw, ci)
            for j = 1:M
                MR = vload(SIMD.Vec{_TB_VW,UInt64}, mrb, (j - 1) * wpbv + ci)
                MI = vload(SIMD.Vec{_TB_VW,UInt64}, mib, (j - 1) * wpbv + ci)
                GR = vload(SIMD.Vec{_TB_VW,UInt64}, gmrb, (j - 1) * wpbv + ci)
                GI = vload(SIMD.Vec{_TB_VW,UInt64}, gmib, (j - 1) * wpbv + ci)
                prc = MR ⊻ CCv
                pis = MI ⊻ CSv
                pic = MI ⊻ CCv
                prs = MR ⊻ CSv
                Xac = GR ⊻ HCv
                Aac = GR & HCv
                Xbs = GI ⊻ HSv
                Abs = GI & HSv
                Xcc = GI ⊻ HCv
                Acc = GI & HCv
                Xes = GR ⊻ HSv
                Aes = GR & HSv
                GXA[j, 1] += sum(count_ones(Xac) + (count_ones(Aac) << 2)) % Int64
                GXA[j, 2] += sum(count_ones(Xbs) + (count_ones(Abs) << 2)) % Int64
                GXA[j, 3] += sum(count_ones(Xcc) + (count_ones(Acc) << 2)) % Int64
                GXA[j, 4] += sum(count_ones(Xes) + (count_ones(Aes) << 2)) % Int64
                for k = 1:NC
                    CW = vload(SIMD.Vec{_TB_VW,UInt64}, codeb, (k - 1) * wpbv + ci)
                    SA = CW ⊻ prc
                    SB = CW ⊻ pis
                    SC = CW ⊻ pic
                    SE = CW ⊻ prs
                    IS[j, k] += sum(count_ones(SA) + count_ones(SB)) % Int64
                    IF[j, k] +=
                        sum(
                            (count_ones(SA & Xac) + (count_ones(SA & Aac) << 2)) +
                            (count_ones(SB & Xbs) + (count_ones(SB & Abs) << 2)),
                        ) % Int64
                    QS[j, k] += sum(count_ones(SC)) % Int64 - sum(count_ones(SE)) % Int64
                    QF[j, k] +=
                        sum(count_ones(SC & Xcc) + (count_ones(SC & Acc) << 2)) % Int64 -
                        sum(count_ones(SE & Xes) + (count_ones(SE & Aes) << 2)) % Int64
                end
            end
            ci += _TB_VW
        end
        blk_off += len
    end

    tap(j, k) = complex(
        Float64(2 * N - 2 * IS[j, k] + 2 * (GXA[j, 1] + GXA[j, 2]) - 4 * IF[j, k]),
        Float64(-2 * QS[j, k] + 2 * (GXA[j, 3] - GXA[j, 4]) - 4 * QF[j, k]),
    )
    if M == 1
        return [tap(1, k) for k = 1:NC]
    else
        return [SVector{M,ComplexF64}(ntuple(j -> tap(j, k), M)) for k = 1:NC]
    end
end

# ── Multi-signal-per-sat tile-share ───────────────────────────────────────────
# A satellite carrying several signals on one carrier shares the carrier and the
# measurement — and therefore the whole carrier wipe-off AND the magnitude mask totals
# (they involve only measurement and carrier planes, never the code). So per strip-mine
# block the carrier planes are generated ONCE, the measurement planes packed/realigned
# ONCE (per antenna), and the per-antenna GXA mask totals accumulated ONCE (during the
# first signal's pass); each signal then correlates its own code planes against them.
# Bit-identical to correlating each signal alone. Mirrors the one-bit
# `_onebit_hybrid_blocked_multi!`; static tap counts (SVector shifts) only, like the
# other tile-shares.
@generated function _twobit_hybrid_blocked_multi!(
    dc::_TwoBitDC,
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
    zc = :(zero(SIMD.Vec{_TB_VW,UInt64}))

    # Per-signal shift/span/offset locals + per-(signal, antenna, tap) IS/IF/QS/QF and
    # per-antenna GXA accumulators (measurement/carrier only — sat-shared, NOT per signal).
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
        for j = 1:M, k = 1:NCs[i], s in ("IS", "IF", "QS", "QF")
            push!(setup.args, :($(Symbol("$(s)_$(i)_$(j)_$k")) = $zc))
        end
    end
    for j = 1:M, c in ("ac", "bs", "cc", "es")
        push!(setup.args, :($(Symbol("GXA_$(c)_$j")) = $zc))
    end

    gates = Expr(:block)
    for i = 1:N
        push!(gates.args, :(_tb_check_modulation(signal_types[$i])))
    end

    # Per antenna, measurement packing: band-shared funnel-realign (>1 sat) vs direct
    # pack (1 sat). Filled ONCE per block, shared across the sat's signals.
    realign = Expr(:block)
    directpack = Expr(:block)
    for j = 1:M
        push!(
            realign.args,
            quote
                _tb_realign_meas!(
                    mrb,
                    mib,
                    gmrb,
                    gmib,
                    $(j - 1) * wpbv,
                    band,
                    $(j - 1) * bandstride,
                    base,
                    nwb,
                )
                _tb_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmrb, $(j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmib, $(j - 1) * wpbv, nwb, nwv, r)
            end,
        )
        push!(
            directpack.args,
            quote
                _tb_pack_meas!(
                    mrb,
                    mib,
                    gmrb,
                    gmib,
                    $(j - 1) * wpbv,
                    signal,
                    $j,
                    p_sig,
                    $(j - 1) * num_rows,
                    base,
                    len,
                    nwb,
                    r,
                    thr,
                )
                _tb_finish!(mrb, $(j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(mib, $(j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmrb, $(j - 1) * wpbv, nwb, nwv, r)
                _tb_finish!(gmib, $(j - 1) * wpbv, nwb, nwv, r)
            end,
        )
    end

    # Per-signal, within a block: one-shot code fill + pack into `codeb`, then the masked
    # popcount correlate reading the shared carrier/measurement planes. The sat-shared
    # GXA totals accumulate during signal 1's pass only.
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
            :(_tb_pack_code!(peb, 0, extb, 0, nwe, (len + $(Symbol("span_$i"))) & 63)),
        )
        push!(b.args, :(
            for _pw = 1:(($(Symbol("span_$i"))>>6)+1)   # zero the funnel-shift read pad
                peb[nwe+_pw] = 0
            end
        ))
        for k = 1:NCs[i]
            push!(
                b.args,
                :(_tb_shift_plane!(
                    codeb,
                    $(k - 1) * wpbv,
                    peb,
                    $(Symbol("off_$(i)_$k")),
                    nwb,
                )),
            )
            push!(b.args, :(_tb_finish!(codeb, $(k - 1) * wpbv, nwb, nwv, r)))
        end
        corr = Expr(:block)
        for j = 1:M
            push!(
                corr.args,
                :(MR = vload(SIMD.Vec{_TB_VW,UInt64}, mrb, $(j - 1) * wpbv + ci)),
            )
            push!(
                corr.args,
                :(MI = vload(SIMD.Vec{_TB_VW,UInt64}, mib, $(j - 1) * wpbv + ci)),
            )
            push!(
                corr.args,
                :(GR = vload(SIMD.Vec{_TB_VW,UInt64}, gmrb, $(j - 1) * wpbv + ci)),
            )
            push!(
                corr.args,
                :(GI = vload(SIMD.Vec{_TB_VW,UInt64}, gmib, $(j - 1) * wpbv + ci)),
            )
            push!(corr.args, :(prc = MR ⊻ CCv))
            push!(corr.args, :(pis = MI ⊻ CSv))
            push!(corr.args, :(pic = MI ⊻ CCv))
            push!(corr.args, :(prs = MR ⊻ CSv))
            push!(corr.args, :(Xac = GR ⊻ HCv))
            push!(corr.args, :(Aac = GR & HCv))
            push!(corr.args, :(Xbs = GI ⊻ HSv))
            push!(corr.args, :(Abs = GI & HSv))
            push!(corr.args, :(Xcc = GI ⊻ HCv))
            push!(corr.args, :(Acc = GI & HCv))
            push!(corr.args, :(Xes = GR ⊻ HSv))
            push!(corr.args, :(Aes = GR & HSv))
            if i == 1
                for (c, X, A) in (
                    ("ac", :Xac, :Aac),
                    ("bs", :Xbs, :Abs),
                    ("cc", :Xcc, :Acc),
                    ("es", :Xes, :Aes),
                )
                    push!(
                        corr.args,
                        :(
                            $(Symbol("GXA_$(c)_$j")) +=
                                count_ones($X) + (count_ones($A) << 2)
                        ),
                    )
                end
            end
            for k = 1:NCs[i]
                push!(
                    corr.args,
                    :(CW = vload(SIMD.Vec{_TB_VW,UInt64}, codeb, $(k - 1) * wpbv + ci)),
                )
                push!(corr.args, :(SA = CW ⊻ prc))
                push!(corr.args, :(SB = CW ⊻ pis))
                push!(corr.args, :(SC = CW ⊻ pic))
                push!(corr.args, :(SE = CW ⊻ prs))
                push!(
                    corr.args,
                    :($(Symbol("IS_$(i)_$(j)_$k")) += count_ones(SA) + count_ones(SB)),
                )
                push!(
                    corr.args,
                    :(
                        $(Symbol("IF_$(i)_$(j)_$k")) +=
                            (count_ones(SA & Xac) + (count_ones(SA & Aac) << 2)) +
                            (count_ones(SB & Xbs) + (count_ones(SB & Abs) << 2))
                    ),
                )
                push!(
                    corr.args,
                    :($(Symbol("QS_$(i)_$(j)_$k")) += count_ones(SC) - count_ones(SE)),
                )
                push!(
                    corr.args,
                    :(
                        $(Symbol("QF_$(i)_$(j)_$k")) +=
                            (count_ones(SC & Xcc) + (count_ones(SC & Acc) << 2)) -
                            (count_ones(SE & Xes) + (count_ones(SE & Aes) << 2))
                    ),
                )
            end
        end
        push!(b.args, quote
            ci = 1
            while ci <= nwv
                CCv = vload(SIMD.Vec{_TB_VW,UInt64}, cosw, ci)
                CSv = vload(SIMD.Vec{_TB_VW,UInt64}, sinw, ci)
                HCv = vload(SIMD.Vec{_TB_VW,UInt64}, hcosw, ci)
                HSv = vload(SIMD.Vec{_TB_VW,UInt64}, hsinw, ci)
                $corr
                ci += _TB_VW
            end
        end)
        b
    end
    sigs = Expr(:block)
    for i = 1:N
        push!(sigs.args, signal_block(i))
    end

    # Finalize: per signal, a Vector of NCᵢ tap sums (matching the single-signal kernel's
    # return so `_correlate_signals` is shared).
    s(sym) = :(sum($sym) % Int64)
    function tapval(i, k)
        function cplx(j)
            :(complex(
                Float64(
                    2 * N_samp - 2 * $(s(Symbol("IS_$(i)_$(j)_$k"))) +
                    2 * ($(s(Symbol("GXA_ac_$j"))) + $(s(Symbol("GXA_bs_$j")))) -
                    4 * $(s(Symbol("IF_$(i)_$(j)_$k"))),
                ),
                Float64(
                    -2 * $(s(Symbol("QS_$(i)_$(j)_$k"))) +
                    2 * ($(s(Symbol("GXA_cc_$j"))) - $(s(Symbol("GXA_es_$j")))) -
                    4 * $(s(Symbol("QF_$(i)_$(j)_$k"))),
                ),
            ))
        end
        M == 1 ? cplx(1) : :(SVector{$M,ComplexF64}(tuple($([cplx(j) for j = 1:M]...))))
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
        thr = dc.threshold
        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        cps_car = carrier_freq / sampling_freq
        carphase = Float64(carrier_phase)
        num_rows = size(signal, 1)
        p_sig = Ptr{Int16}(pointer(signal))
        band = dc.band                              # band-shared measurement planes (>1 sat)
        band_shared = !isempty(band.mrband)
        bandstride = cld(num_rows, 64) + 2
        $setup
        maxspan_v = $maxspan

        bufs = _tb_scratch(dc)
        blk = dc.blk
        wpb = cld(blk, 64)
        wpbv = cld(wpb, _TB_VW) * _TB_VW
        pewords = cld(blk + maxspan_v, 64) + (maxspan_v >> 6) + 1 + _TB_VW   # + funnel-shift pad
        length(bufs.extb) < blk + maxspan_v + 64 && resize!(bufs.extb, blk + maxspan_v + 64)
        length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
        length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
        length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
        length(bufs.hsinw) < wpbv && resize!(bufs.hsinw, wpbv)
        length(bufs.hcosw) < wpbv && resize!(bufs.hcosw, wpbv)
        length(bufs.codeb) < $maxNC * wpbv && resize!(bufs.codeb, $maxNC * wpbv)
        length(bufs.mrb) < $M * wpbv && resize!(bufs.mrb, $M * wpbv)
        length(bufs.mib) < $M * wpbv && resize!(bufs.mib, $M * wpbv)
        length(bufs.gmrb) < $M * wpbv && resize!(bufs.gmrb, $M * wpbv)
        length(bufs.gmib) < $M * wpbv && resize!(bufs.gmib, $M * wpbv)
        extb = bufs.extb;
        peb = bufs.peb;
        sinw = bufs.sinw;
        cosw = bufs.cosw
        hsinw = bufs.hsinw;
        hcosw = bufs.hcosw
        codeb = bufs.codeb;
        mrb = bufs.mrb;
        mib = bufs.mib
        gmrb = bufs.gmrb;
        gmib = bufs.gmib

        blk_off = 0
        @inbounds while blk_off < N_samp
            len = min(blk, N_samp - blk_off)
            nwb = cld(len, 64)
            nwv = cld(nwb, _TB_VW) * _TB_VW
            r = len & 63
            # carrier sign + magnitude planes — shared across signals, antennas and taps
            generate_carrier_signs_mags!(
                sinw,
                cosw,
                hsinw,
                hcosw,
                len,
                cps_car;
                phase = carphase + cps_car * blk_off,
            )
            _tb_zeropad!(sinw, 0, nwb, nwv)
            _tb_zeropad!(cosw, 0, nwb, nwv)
            _tb_zeropad!(hsinw, 0, nwb, nwv)
            _tb_zeropad!(hcosw, 0, nwb, nwv)
            # measurement planes — one pack per sat/block, shared across signals
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

# ── Correlate / plumbing (mirrors the one-bit backend) ────────────────────────
# Single signal per sat: one kernel call.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal},
    per_signal_completed::Tuple{Bool},
    dc::_TwoBitDC,
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
    p = _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    )
    new_acc = _twobit_hybrid_blocked!(
        dc,
        signal,
        _tb_num_ants_val(head.correlator),
        head.signal,
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
    new_corr =
        update_accumulator(head.correlator, get_accumulators(head.correlator) .+ new_acc)
    ((new_corr, per_signal_completed[1]),)
end

# Multiple signals per sat: share one carrier + measurement downconvert across the
# sat's signals via the tile-share kernel. Returns the per-signal
# `(new_correlator, is_integration_completed)` tuples.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal,TrackedSignal,Vararg{TrackedSignal}},
    per_signal_completed::Tuple,
    dc::_TwoBitDC,
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
    new_accs = _twobit_hybrid_blocked_multi!(
        dc,
        signal,
        _tb_num_ants_val(correlators[1]),
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

function _update_tracked_sat_correlator(
    sat::TrackedSat,
    dc::_TwoBitDC,
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
    samples_to_integrate == 0 && return sat
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

@inline function _dc_one_group!(
    g::SignalGroup,
    dc::_TwoBitDC,
    measurements::BandMeasurements,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    m = measurements[get_band_id(g.band)]
    eltype(m.samples) === Complex{Int16} || throw(
        ArgumentError(
            string(
                "TwoBitDownconvertAndCorrelator requires `Complex{Int16}` measurement ",
                "samples (12-bit ADC); got element type ",
                eltype(m.samples),
                ". Use a CPU(Threaded)DownconvertAndCorrelator for floating-point samples.",
            ),
        ),
    )
    # Pack the band's measurement planes ONCE, shared across sats — but only for >1 sat:
    # for a single sat the shared pack + per-sat realign copy is slower than packing that
    # one sat directly, so empty the band buffer to select the kernel's per-sat path.
    samples = m.samples
    if length(vals) > 1
        nsamp = get_num_samples(m)
        M = samples isa AbstractMatrix ? size(samples, 2) : 1
        num_rows = samples isa AbstractMatrix ? size(samples, 1) : length(samples)
        GC.@preserve samples _tb_pack_band!(
            dc.band,
            samples,
            nsamp,
            M,
            Ptr{Int16}(pointer(samples)),
            num_rows,
            dc.threshold,
            dc isa TwoBitThreadedDownconvertAndCorrelator,
        )
    else
        resize!(dc.band.mrband, 0)
        resize!(dc.band.miband, 0)
        resize!(dc.band.gmrband, 0)
        resize!(dc.band.gmiband, 0)
    end
    _dc_group_loop!(
        dc,
        vals,
        m.samples,
        get_num_samples(m),
        m.sampling_frequency,
        m.intermediate_frequency,
    )
end

@inline function _dc_group_loop!(
    dc::TwoBitDownconvertAndCorrelator,
    vals,
    args::Vararg{Any,4},
)
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

@inline function _dc_group_loop!(
    dc::TwoBitThreadedDownconvertAndCorrelator,
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

Downconvert and correlate all satellites with the two-bit bit-wise backend.
"""
function downconvert_and_correlate(
    dc::_TwoBitDC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    downconvert_and_correlate!(dc, measurements, new_track_state)
end

"""
$(SIGNATURES)

In-place two-bit downconvert and correlate. Returns the same `track_state`.
"""
function downconvert_and_correlate!(
    dc::_TwoBitDC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    _foreach_group!(_dc_one_group!, track_state.groups, dc, measurements)
    return track_state
end
