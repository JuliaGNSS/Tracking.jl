# Two-bit (sign+magnitude) bit-wise downconvert + correlate backend.
#
# The middle point between the one-bit and Int16 backends. Like one-bit it packs everything into
# bit planes and correlates with XOR + popcount, but it keeps a second bit — a MAGNITUDE bit — for
# the measurement (always) and, when `CB == 2`, for the carrier, giving the classic near-optimal
# 4-level `{±1, ±3}` quantiser. The code stays 1-bit (±1), which is exact for binary modulations.
#
#   value = s·(1 + 2·b)          s = sign bit;  b = magnitude bit (1 ⇔ |sample| ≥ threshold)
#         ∈ {±1, ±3}
#
# Because code and (for CB=1) carrier stay ±1, a product's SIGN is still an XOR of sign bits; only
# the WEIGHT changes. For the `mᵣ·cos` term of tap x (sign plane S = codeₓ ⊻ smr ⊻ scos), with
# measurement magnitude plane Gmr and (CB=2) carrier magnitude plane Hcos:
#
#   Σ codeₓ·mᵣ·cos = (N − 2·pc S)                       ← one-bit result
#                  + 2·(pc Gmr  − 2·pc(Gmr & S))        ← measurement magnitude correction
#                  + 2·(pc Hcos − 2·pc(Hcos & S))       ← carrier magnitude correction  (CB=2)
#                  + 4·(pc(Gmr&Hcos) − 2·pc(Gmr&Hcos&S))                                (CB=2)
#
# The four complex sign planes are the one-bit A/B/C/E planes (mr·cos, mi·sin, mi·cos, mr·sin);
# Iₓ = (mr·cos)+(mi·sin), Qₓ = (mi·cos)−(mr·sin). `pc Gmr`, `pc Hcos`, `pc(Gmr&Hcos)` are per-block
# constants (accumulated once per chunk); the per-tap cost is +1 (CB=1) or +3 (CB=2) masked
# popcounts per sign plane over one-bit. Strip-mined into `blk`-sample blocks with L1-resident
# scratch reused across blocks, mirroring the one-bit / integer hybrid-blocked backends.
#
# 2-bit `{±1,±3}` quantisation costs ≈0.55 dB vs float (against ≈1.96 dB for 1-bit): measurement
# 2-bit + 1-bit carrier recovers ≈1.4 dB over the one-bit backend, CB=2 a further ≈0.6 dB.
# Accumulators are `Int64`, converted to `ComplexF64` (M=1) / `SVector{M,ComplexF64}` (M>1) at
# finalize; every downstream consumer is ratio/normalised, so the coarse amplitude is immaterial.
# Scope mirrors one-bit: static tap counts, any antenna M, band-shared measurement for >1 sat.

import SinCosLUT
using SinCosLUT: generate_carrier_signs!, SinCosTable, carrier_engine

# Default strip-mine block length (samples); a multiple of 64 so every block packs whole UInt64s.
const _TWOBIT_BLK = 8192
# UInt64 lanes per SIMD popcount step (see the one-bit backend's `_OB_VW`).
const _TB_VW = 8

# The 2-bit carrier (CB=2) reuses SinCosLUT's Int8 carrier LUT (the same generator the Int16
# backend drives via `_int16_fill_carrier!`) at **amplitude 3**, so a table lane holds
# `round(3·sin) ∈ {0,±1,±2,±3}`. Its sign is the carrier sign plane and `|lane| ≥ 2` is the
# magnitude plane — since `|round(3·sin)| ≥ 2 ⇔ |sin| ≥ 2/3`, exactly the `{±1,±3}` crossover.
# Both planes thus come from one proven source and are self-consistent, no bespoke phase math.
const _TB_CARR_AMP = 3          # Int8 carrier table amplitude → values {0,±1,±2,±3}
const _TB_CARR_THR = Int8(2)    # |value| ≥ 2 ⇔ |sin|,|cos| ≥ 2/3 of peak (the ±3 level)

# ── sign / magnitude mask helpers ─────────────────────────────────────────────
# Pack the sign bit of 64 lanes into a UInt64 (bit j ⇔ lane j < 0). Same as the one-bit `_ob_mm`;
# duplicated here so the two backends stay independent files. Int16 (measurement) + Int8 (carrier).
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
# Magnitude mask: bit j ⇔ |lane j| ≥ thr. `|v| − thr` has its sign bit CLEAR when large, so the
# large plane is the complement of the sign-mask of `|v| − thr` (measurement components ≤ 2047
# counts and carrier lanes ≤ 3 keep `|v|−thr` in range). The caller masks the last (partial) word.
@inline _tb_mag_mm(v::SIMD.Vec{64,Int16}, thr::SIMD.Vec{64,Int16}) = ~_tb_mm(abs(v) - thr)
@inline _tb_mag_mm(v::SIMD.Vec{64,Int8}, thr::SIMD.Vec{64,Int8}) = ~_tb_mm(abs(v) - thr)
@inline _tb_masklast(w::UInt64, r::Int) = r == 0 ? w : w & ((UInt64(1) << r) - UInt64(1))

# ── code plane helpers (code is 1-bit; identical to the one-bit backend) ──────
@inline function _tb_pack_code!(
    codeb::Vector{UInt64},
    koff::Int,
    extb::Vector{Int8},
    byteoff::Int,
    nwb::Int,
    r::Int,
)
    @inbounds for w = 0:(nwb-1)
        wd = Base.llvmcall(
            (
                """
define i64 @entry(<64 x i8> %v) #0 { %c = icmp slt <64 x i8> %v, zeroinitializer
  %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""",
                "entry",
            ),
            UInt64,
            Tuple{NTuple{64,Base.VecElement{Int8}}},
            SIMD.vload(SIMD.Vec{64,Int8}, extb, byteoff + (w << 6) + 1).data,
        )
        codeb[koff+w+1] = w == nwb - 1 ? _tb_masklast(wd, r) : wd
    end
    nothing
end

# Derive Early/Late tap planes from the prompt-extended plane by a funnel bit-shift (see one-bit).
@inline function _tb_shift_plane!(
    dst::Vector{UInt64},
    doff::Int,
    pe::Vector{UInt64},
    off::Int,
    nwb::Int,
)
    if off == 0
        @inbounds for w = 1:nwb
            dst[doff+w] = pe[w]
        end
    else
        lo = off
        hi = 64 - off
        @inbounds for w = 1:nwb
            dst[doff+w] = (pe[w] >> lo) | (pe[w+1] << hi)
        end
    end
    nothing
end

# Mask the last valid word to `r` bits (if partial) and zero the VW-pad words [nwb+1, nwv].
@inline function _tb_finish!(buf::Vector{UInt64}, off::Int, nwb::Int, nwv::Int, r::Int)
    @inbounds (r != 0) && (buf[off+nwb] &= (UInt64(1) << r) - UInt64(1))
    @inbounds for w = (nwb+1):nwv
        buf[off+w] = 0
    end
    nothing
end
@inline function _tb_zeropad!(buf::Vector{UInt64}, off::Int, nwb::Int, nwv::Int)
    @inbounds for w = (nwb+1):nwv
        buf[off+w] = 0
    end
    nothing
end

# ── carrier planes from the Int8 LUT (CB=2) ───────────────────────────────────
# Pack the block's Int8 carrier (`csb`/`ccb`, amplitude 3, filled by `_int16_fill_carrier!`) into
# sign planes (`sinw`/`cosw`: bit j ⇔ lane j < 0) and magnitude planes (`hsinw`/`hcosw`: bit j ⇔
# |lane j| ≥ 2), 64 lanes per UInt64 word. The last (partial) word is masked to `r` bits so the
# per-block magnitude-count popcount excludes the invalid high lanes; the VW pad is zeroed by the
# caller. Sign and magnitude come from the same table lane, so they are self-consistent.
@inline function _tb_pack_carrier!(
    sinw::Vector{UInt64},
    cosw::Vector{UInt64},
    hsinw::Vector{UInt64},
    hcosw::Vector{UInt64},
    csb::Vector{Int8},
    ccb::Vector{Int8},
    len::Int,
)
    r = len & 63
    full = fld(len, 64)
    thrv = SIMD.Vec{64,Int8}(_TB_CARR_THR)
    @inbounds for w = 0:(full-1)
        sv = SIMD.vload(SIMD.Vec{64,Int8}, csb, (w << 6) + 1)
        cv = SIMD.vload(SIMD.Vec{64,Int8}, ccb, (w << 6) + 1)
        sinw[w+1] = _tb_mm(sv)
        cosw[w+1] = _tb_mm(cv)
        hsinw[w+1] = _tb_mag_mm(sv, thrv)
        hcosw[w+1] = _tb_mag_mm(cv, thrv)
    end
    if r != 0
        sv = SIMD.vload(SIMD.Vec{64,Int8}, csb, (full << 6) + 1)
        cv = SIMD.vload(SIMD.Vec{64,Int8}, ccb, (full << 6) + 1)
        sinw[full+1] = _tb_masklast(_tb_mm(sv), r)
        cosw[full+1] = _tb_masklast(_tb_mm(cv), r)
        hsinw[full+1] = _tb_masklast(_tb_mag_mm(sv, thrv), r)
        hcosw[full+1] = _tb_masklast(_tb_mag_mm(cv, thrv), r)
    end
    nothing
end

# ── measurement packing (sign + magnitude) ────────────────────────────────────
# Per-sat direct pack (single-sat group): sign(real)→mrb, |real|≥thr→gmrb, likewise imag.
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
    thrv = SIMD.Vec{64,Int16}(thr)
    full = fld(len, 64)
    @inbounds for w = 0:(full-1)
        byte_off = (colbase + base + (w << 6) - 1) * 2 * sizeof(Int16)
        re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
        mrb[joff+w+1] = _tb_mm(re)
        mib[joff+w+1] = _tb_mm(im)
        gmrb[joff+w+1] = _tb_mag_mm(re, thrv)
        gmib[joff+w+1] = _tb_mag_mm(im, thrv)
    end
    if r != 0
        wr = zero(UInt64)
        wi = zero(UInt64)
        wgr = zero(UInt64)
        wgi = zero(UInt64)
        @inbounds for i = 0:(r-1)
            sig = signal[base+(full<<6)+i, j]
            re = real(sig)
            im = imag(sig)
            (re < 0) && (wr |= UInt64(1) << i)
            (im < 0) && (wi |= UInt64(1) << i)
            (abs(re) >= thr) && (wgr |= UInt64(1) << i)
            (abs(im) >= thr) && (wgi |= UInt64(1) << i)
        end
        mrb[joff+full+1] = wr
        mib[joff+full+1] = wi
        gmrb[joff+full+1] = wgr
        gmib[joff+full+1] = wgi
    end
    nothing
end

# Band-shared measurement (see one-bit `OneBitBandBuffers`): sign + magnitude planes packed once
# per group, funnel-shared across sats.
struct TwoBitBandBuffers
    mrband::Vector{UInt64}
    miband::Vector{UInt64}
    gmrband::Vector{UInt64}
    gmiband::Vector{UInt64}
end
TwoBitBandBuffers() = TwoBitBandBuffers(UInt64[], UInt64[], UInt64[], UInt64[])

function _tb_pack_band!(
    band::TwoBitBandBuffers,
    signal,
    num_samples::Int,
    M::Int,
    p_sig::Ptr{Int16},
    num_rows::Int,
    thr::Int16,
)
    nwb = cld(num_samples, 64)
    stride = nwb + 2                                # 2 funnel-pad words (realign reads past nwb)
    length(band.mrband) < M * stride && resize!(band.mrband, M * stride)
    length(band.miband) < M * stride && resize!(band.miband, M * stride)
    length(band.gmrband) < M * stride && resize!(band.gmrband, M * stride)
    length(band.gmiband) < M * stride && resize!(band.gmiband, M * stride)
    thrv = SIMD.Vec{64,Int16}(thr)
    r = num_samples & 63
    full = fld(num_samples, 64)
    @inbounds for j = 1:M
        off = (j - 1) * stride
        colbase = (j - 1) * num_rows
        for w = 0:(full-1)
            byte_off = (colbase + w * 64) * 2 * sizeof(Int16)
            re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
            band.mrband[off+w+1] = _tb_mm(re)
            band.miband[off+w+1] = _tb_mm(im)
            band.gmrband[off+w+1] = _tb_mag_mm(re, thrv)
            band.gmiband[off+w+1] = _tb_mag_mm(im, thrv)
        end
        if r != 0
            wr = zero(UInt64)
            wi = zero(UInt64)
            wgr = zero(UInt64)
            wgi = zero(UInt64)
            for i = 0:(r-1)
                sig = signal[full*64+i+1, j]
                re = real(sig)
                im = imag(sig)
                (re < 0) && (wr |= UInt64(1) << i)
                (im < 0) && (wi |= UInt64(1) << i)
                (abs(re) >= thr) && (wgr |= UInt64(1) << i)
                (abs(im) >= thr) && (wgi |= UInt64(1) << i)
            end
            band.mrband[off+full+1] = wr
            band.miband[off+full+1] = wi
            band.gmrband[off+full+1] = wgr
            band.gmiband[off+full+1] = wgi
        end
        for pad in (band.mrband, band.miband, band.gmrband, band.gmiband)
            pad[off+nwb+1] = 0
            pad[off+nwb+2] = 0
        end
    end
    nothing
end

# Funnel-realign antenna j's block from the shared band planes into thread-local block planes.
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
    b0 = base - 1
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

# ── scratch ────────────────────────────────────────────────────────────────
# One-bit's planes plus the measurement magnitude planes (gmrb/gmib) and, for CB=2, the carrier
# magnitude planes (hsinw/hcosw) and the Int8 carrier fill scratch (csb/ccb) — all empty and unused
# for CB=1. Grown lazily, reused across blocks.
struct TwoBitScratchBuffers
    extb::Vector{Int8}
    peb::Vector{UInt64}
    sinw::Vector{UInt64}
    cosw::Vector{UInt64}
    hsinw::Vector{UInt64}
    hcosw::Vector{UInt64}
    csb::Vector{Int8}
    ccb::Vector{Int8}
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
    Int8[],
    Int8[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
    UInt64[],
)

"""
$(SIGNATURES)

Two-bit (sign+magnitude) bit-wise CPU downconvert + correlate backend
(single-threaded). Opt-in alternative to [`OneBitDownconvertAndCorrelator`] for
`Complex{Int16}` sample buffers: the measurement is 2-bit `{±1,±3}`, the carrier
is `carrier_bits`-bit (1 or 2), the code stays 1-bit, correlated with
XOR + masked popcount. Better SNR than the one-bit backend (≈1.4 dB from the
2-bit measurement, a further ≈0.6 dB with `carrier_bits = 2`), in the same speed
class as the Int16 backend (≈2–3× faster than Float32) with `carrier_bits = 1`.
Construct **once outside** the `track!` loop and pass it via the
`downconvert_and_correlator` keyword for an allocation-free steady state.

`threshold` is the measurement magnitude split point in ADC counts (`|sample| ≥ threshold` ⇒ the `±3` level); set it near 1σ of your front end's input for the
near-optimal 4-level quantiser. Downstream consumers are ratio-normalised, so
the absolute scale is immaterial. With `carrier_bits = 2` the carrier is drawn
from SinCosLUT's Int8 LUT (amplitude 3) and thresholded to `{±1,±3}`; `steps`
sets that table's phase resolution.
"""
struct TwoBitDownconvertAndCorrelator{CB,TBL} <: AbstractDownconvertAndCorrelator
    buffers::TwoBitScratchBuffers
    band::TwoBitBandBuffers
    table::TBL
    blk::Int
    threshold::Int16
end

"""
$(SIGNATURES)

Multi-threaded two-bit bit-wise backend. One `TwoBitScratchBuffers` per thread
(indexed by `Threads.threadid()` inside `@batch`); the carrier `table` is
immutable and shared. See [`TwoBitDownconvertAndCorrelator`](@ref).
"""
struct TwoBitThreadedDownconvertAndCorrelator{CB,TBL} <: AbstractDownconvertAndCorrelator
    buffers::Vector{TwoBitScratchBuffers}
    band::TwoBitBandBuffers
    table::TBL
    blk::Int
    threshold::Int16
end

function _tb_check_carrier_bits(carrier_bits::Integer)
    (carrier_bits == 1 || carrier_bits == 2) || throw(
        ArgumentError(
            "carrier_bits must be 1 (sign-only carrier) or 2 (sign+magnitude); got $carrier_bits",
        ),
    )
    Int(carrier_bits)
end

# Int8 carrier table for the CB=2 sign+magnitude carrier (amplitude 3 → {0,±1,±2,±3} lanes).
_tb_carrier_table(steps::Integer) =
    SinCosTable(Int8; steps = Int(steps), amplitude = _TB_CARR_AMP)

function TwoBitDownconvertAndCorrelator(;
    carrier_bits::Integer = 1,
    threshold::Integer = 512,
    steps::Integer = 64,
    blk::Integer = _TWOBIT_BLK,
)
    CB = _tb_check_carrier_bits(carrier_bits)
    table = _tb_carrier_table(steps)
    TwoBitDownconvertAndCorrelator{CB,typeof(table)}(
        TwoBitScratchBuffers(),
        TwoBitBandBuffers(),
        table,
        Int(blk),
        Int16(threshold),
    )
end

function TwoBitThreadedDownconvertAndCorrelator(;
    carrier_bits::Integer = 1,
    threshold::Integer = 512,
    steps::Integer = 64,
    blk::Integer = _TWOBIT_BLK,
)
    CB = _tb_check_carrier_bits(carrier_bits)
    table = _tb_carrier_table(steps)
    TwoBitThreadedDownconvertAndCorrelator{CB,typeof(table)}(
        [TwoBitScratchBuffers() for _ = 1:Threads.maxthreadid()],
        TwoBitBandBuffers(),
        table,
        Int(blk),
        Int16(threshold),
    )
end

const _TwoBitDC =
    Union{TwoBitDownconvertAndCorrelator,TwoBitThreadedDownconvertAndCorrelator}

@inline _tb_scratch(dc::TwoBitDownconvertAndCorrelator) = dc.buffers
@inline _tb_scratch(dc::TwoBitThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

@inline _tb_num_ants_val(::AbstractCorrelator{M}) where {M} = NumAnts{M}()

# ── the two-bit hybrid-blocked kernel ─────────────────────────────────────────
# `@generated` over (NC, M) from the argument types and CB from `dc`'s concrete type. Emits, per
# (antenna j, tap k), the masked-popcount accumulators for the four sign planes A/B/C/E, plus (CB=2)
# their carrier-magnitude and both-magnitude counterparts, and per-block magnitude-count constants.
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
    CB = dc.parameters[1]::Int
    zc = :(zero(SIMD.Vec{_TB_VW,UInt64}))

    # Accumulator declarations.
    init = Expr(:block)
    for k = 1:NC
        push!(init.args, :($(Symbol("off_$k")) = Int(sample_shifts[$k]) - min_shift))
    end
    for j = 1:M
        push!(init.args, :($(Symbol("GmrN_$j")) = $zc))
        push!(init.args, :($(Symbol("GmiN_$j")) = $zc))
        if CB == 2
            push!(init.args, :($(Symbol("GHac_$j")) = $zc))  # gmr & hcos (A: mr·cos)
            push!(init.args, :($(Symbol("GHbs_$j")) = $zc))  # gmi & hsin (B: mi·sin)
            push!(init.args, :($(Symbol("GHcc_$j")) = $zc))  # gmi & hcos (C: mi·cos)
            push!(init.args, :($(Symbol("GHes_$j")) = $zc))  # gmr & hsin (E: mr·sin)
        end
        for k = 1:NC, L in ("A", "B", "C", "E")
            push!(init.args, :($(Symbol("$(L)_$(j)_$k")) = $zc))
            push!(init.args, :($(Symbol("$(L)g_$(j)_$k")) = $zc))
            if CB == 2
                push!(init.args, :($(Symbol("$(L)h_$(j)_$k")) = $zc))
                push!(init.args, :($(Symbol("$(L)gh_$(j)_$k")) = $zc))
            end
        end
    end
    if CB == 2
        push!(init.args, :(HcosN = $zc))
        push!(init.args, :(HsinN = $zc))
    end

    # Per-VW-chunk accumulate.
    body = Expr(:block)
    if CB == 2
        push!(body.args, :(Hc = vload(SIMD.Vec{_TB_VW,UInt64}, hcosw, i)))
        push!(body.args, :(Hs = vload(SIMD.Vec{_TB_VW,UInt64}, hsinw, i)))
        push!(body.args, :(HcosN += count_ones(Hc)))
        push!(body.args, :(HsinN += count_ones(Hs)))
    end
    for j = 1:M
        push!(body.args, :(MR = vload(SIMD.Vec{_TB_VW,UInt64}, mrb, $(j - 1) * wpbv + i)))
        push!(body.args, :(MI = vload(SIMD.Vec{_TB_VW,UInt64}, mib, $(j - 1) * wpbv + i)))
        push!(body.args, :(GR = vload(SIMD.Vec{_TB_VW,UInt64}, gmrb, $(j - 1) * wpbv + i)))
        push!(body.args, :(GI = vload(SIMD.Vec{_TB_VW,UInt64}, gmib, $(j - 1) * wpbv + i)))
        push!(body.args, :($(Symbol("GmrN_$j")) += count_ones(GR)))
        push!(body.args, :($(Symbol("GmiN_$j")) += count_ones(GI)))
        push!(body.args, :(prc = MR ⊻ CCv))
        push!(body.args, :(pis = MI ⊻ CSv))
        push!(body.args, :(pic = MI ⊻ CCv))
        push!(body.args, :(prs = MR ⊻ CSv))
        if CB == 2
            push!(body.args, :($(Symbol("GHac_$j")) += count_ones(GR & Hc)))
            push!(body.args, :($(Symbol("GHbs_$j")) += count_ones(GI & Hs)))
            push!(body.args, :($(Symbol("GHcc_$j")) += count_ones(GI & Hc)))
            push!(body.args, :($(Symbol("GHes_$j")) += count_ones(GR & Hs)))
        end
        # (sign-plane var, measurement-mag Vec, carrier-mag Vec)
        planes = (
            ("A", :prc, :GR, :Hc),
            ("B", :pis, :GI, :Hs),
            ("C", :pic, :GI, :Hc),
            ("E", :prs, :GR, :Hs),
        )
        for k = 1:NC
            push!(
                body.args,
                :(CW = vload(SIMD.Vec{_TB_VW,UInt64}, codeb, $(k - 1) * wpbv + i)),
            )
            for (L, sp, gv, hv) in planes
                s = Symbol("S$L")
                push!(body.args, :($s = CW ⊻ $sp))
                push!(body.args, :($(Symbol("$(L)_$(j)_$k")) += count_ones($s)))
                push!(body.args, :($(Symbol("$(L)g_$(j)_$k")) += count_ones($gv & $s)))
                if CB == 2
                    push!(body.args, :($(Symbol("$(L)h_$(j)_$k")) += count_ones($hv & $s)))
                    push!(
                        body.args,
                        :($(Symbol("$(L)gh_$(j)_$k")) += count_ones($gv & $hv & $s)),
                    )
                end
            end
        end
    end

    # Finalize.
    s(sym) = :(Int64(sum($sym)))
    # Σ code·(measurement·carrier) for sign-plane L at (j,k), including magnitude corrections.
    function sigma(L, j, k, gN, HN, GHN)
        base = :(N - 2 * $(s(Symbol("$(L)_$(j)_$k"))))
        gterm = :(2 * ($(s(gN)) - 2 * $(s(Symbol("$(L)g_$(j)_$k")))))
        if CB == 2
            hterm = :(2 * ($(s(HN)) - 2 * $(s(Symbol("$(L)h_$(j)_$k")))))
            ghterm = :(4 * ($(s(GHN)) - 2 * $(s(Symbol("$(L)gh_$(j)_$k")))))
            :($base + $gterm + $hterm + $ghterm)
        else
            :($base + $gterm)
        end
    end
    function tapval(k)
        cplx(j) = quote
            mrcos = $(sigma("A", j, k, Symbol("GmrN_$j"), :HcosN, Symbol("GHac_$j")))
            misin = $(sigma("B", j, k, Symbol("GmiN_$j"), :HsinN, Symbol("GHbs_$j")))
            micos = $(sigma("C", j, k, Symbol("GmiN_$j"), :HcosN, Symbol("GHcc_$j")))
            mrsin = $(sigma("E", j, k, Symbol("GmrN_$j"), :HsinN, Symbol("GHes_$j")))
            complex(Float64(mrcos + misin), Float64(micos - mrsin))
        end
        if M == 1
            :(
                let
                    $(cplx(1))
                end
            )
        else
            ant = [:(
                let
                    $(cplx(j))
                end
            ) for j = 1:M]
            :(SVector{$M,ComplexF64}(tuple($(ant...))))
        end
    end
    ret =
        M == 1 ? :(ComplexF64[$([tapval(k) for k = 1:NC]...)]) :
        :(SVector{M,ComplexF64}[$([tapval(k) for k = 1:NC]...)])

    quote
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
        band = dc.band
        band_shared = !isempty(band.mrband)
        bandstride = cld(num_rows, 64) + 2

        bufs = _tb_scratch(dc)
        blk = dc.blk
        wpb = cld(blk, 64)
        wpbv = cld(wpb, _TB_VW) * _TB_VW
        pewords = cld(blk + span, 64) + _TB_VW
        length(bufs.extb) < blk + span + 64 && resize!(bufs.extb, blk + span + 64)
        length(bufs.peb) < pewords && resize!(bufs.peb, pewords)
        length(bufs.sinw) < wpbv && resize!(bufs.sinw, wpbv)
        length(bufs.cosw) < wpbv && resize!(bufs.cosw, wpbv)
        length(bufs.codeb) < $NC * wpbv && resize!(bufs.codeb, $NC * wpbv)
        length(bufs.mrb) < $M * wpbv && resize!(bufs.mrb, $M * wpbv)
        length(bufs.mib) < $M * wpbv && resize!(bufs.mib, $M * wpbv)
        length(bufs.gmrb) < $M * wpbv && resize!(bufs.gmrb, $M * wpbv)
        length(bufs.gmib) < $M * wpbv && resize!(bufs.gmib, $M * wpbv)
        $(CB == 2 ? quote
            ncar = (cld(blk, _INT16_W) + 4) * _INT16_W   # Int8 carrier fill (4-way tail room)
            length(bufs.hsinw) < wpbv && resize!(bufs.hsinw, wpbv)
            length(bufs.hcosw) < wpbv && resize!(bufs.hcosw, wpbv)
            length(bufs.csb) < ncar && resize!(bufs.csb, ncar)
            length(bufs.ccb) < ncar && resize!(bufs.ccb, ncar)
            csb = bufs.csb
            ccb = bufs.ccb
            reng = carrier_engine(dc.table, cps_car)
        end : :())
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

        blk_off = 0
        @inbounds while blk_off < N
            len = min(blk, N - blk_off)
            nwb = cld(len, 64)
            nwv = cld(nwb, _TB_VW) * _TB_VW
            r = len & 63
            # code (1-bit): pack extended prompt plane once, derive taps by funnel shift.
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
            peb[nwe+1] = 0
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
            # carrier planes. CB=1: SinCosLUT's 1-bit NCO (sign only, fast square-wave path).
            # CB=2: fill the Int8 amplitude-3 carrier and pack sign + magnitude from the same lanes.
            $(
                CB == 1 ?
                quote
                    generate_carrier_signs!(
                        sinw,
                        cosw,
                        len,
                        cps_car;
                        phase = carphase + cps_car * blk_off,
                    )
                    _tb_zeropad!(sinw, 0, nwb, nwv)
                    _tb_zeropad!(cosw, 0, nwb, nwv)
                end :
                quote
                    _int16_fill_carrier!(
                        csb,
                        ccb,
                        reng,
                        blk_off,
                        len,
                        carphase,
                        Val(_INT16_W),
                    )
                    _tb_pack_carrier!(sinw, cosw, hsinw, hcosw, csb, ccb, len)
                    _tb_zeropad!(sinw, 0, nwb, nwv)
                    _tb_zeropad!(cosw, 0, nwb, nwv)
                    _tb_zeropad!(hsinw, 0, nwb, nwv)
                    _tb_zeropad!(hcosw, 0, nwb, nwv)
                end
            )
            # measurement sign + magnitude planes.
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
            # XOR + masked popcount accumulate over whole VW chunks.
            i = 1
            while i <= nwv
                CCv = vload(SIMD.Vec{_TB_VW,UInt64}, cosw, i)
                CSv = vload(SIMD.Vec{_TB_VW,UInt64}, sinw, i)
                $body
                i += _TB_VW
            end
            blk_off += len
        end
        $ret
    end
end

# ── correlate / plumbing (mirrors the one-bit backend) ────────────────────────
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
    new_corrs = map(signals) do head
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
        update_accumulator(head.correlator, get_accumulators(head.correlator) .+ new_acc)
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
    m = measurements[band_key(g.band)]
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
