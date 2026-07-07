# One-bit (hard-limited) bit-wise downconvert + correlate backend.
#
# The bit-wise counterpart of the Float32 fused / integer paths. It takes the same
# `Complex{Int16}` (12-bit ADC) sample buffers, but keeps only the SIGN of every
# operand (measurement, carrier, code), so the carrier wipe-off collapses to XOR
# and the tap accumulate to popcount:
#
#   value +1 → bit 0, −1 → bit 1 (signbit); a ±1 product is an XOR of sign bits,
#   and Σ over N samples of a ±1 sequence is  N − 2·popcount.
#
#   DI = mᵣ·cos + mᵢ·sin,  DQ = mᵢ·cos − mᵣ·sin      (all ±1 in 1-bit form) ⇒
#   Iₓ = Σ codeₓ·DI = (N − 2·pc(codeₓ⊻mᵣ⊻cos)) + (N − 2·pc(codeₓ⊻mᵢ⊻sin))
#   Qₓ = Σ codeₓ·DQ = (N − 2·pc(codeₓ⊻mᵢ⊻cos)) − (N − 2·pc(codeₓ⊻mᵣ⊻sin))
#
# Per block: pack the measurement sign planes (per antenna), generate the carrier
# sign planes straight off the NCO's top bit (sign(sin)=MSB(acc), sign(cos)=MSB(acc+¼)),
# pack the code sign plane per tap (movemask of the GNSSSignals Int8 ±1 replica), then
# XOR + `count_ones` accumulate. Strip-mined into `blk`-sample blocks with L1-resident
# scratch reused across blocks, mirroring the integer hybrid-blocked backend.
#
# 1-bit quantisation costs ≈2 dB of SNR (plus ≈1 dB for the 1-bit carrier), the classic
# trade for bit-wise speed and 1-bit memory bandwidth. Accumulators are converted to
# `ComplexF64` (M=1) / `SVector{M,ComplexF64}` (M>1) at finalize; every downstream
# consumer (discriminators, C/N0, bit buffer) is ratio/normalised, so the absolute 1-bit
# scale is immaterial. Scope mirrors the integer backend: static tap counts (with a
# runtime/dynamic tap-count fallback), any antenna M. A multi-signal sat shares one
# carrier + measurement downconvert across its signals (tile-share).

import SinCosLUT
using SinCosLUT: generate_carrier_signs!

# Default strip-mine block length (samples); a multiple of 64 so every block packs a
# whole number of UInt64 words.
const _ONEBIT_BLK = 8192

# UInt64 lanes per SIMD popcount step (`count_ones(::Vec{VW,UInt64})` → VPOPCNTQ on AVX-512,
# a movemask/popcount sequence elsewhere). Block planes are padded to a multiple of VW words
# and the pad is zeroed, so the correlate loop runs whole VW chunks with no scalar tail.
const _OB_VW = 8

# ── sign-mask helpers: pack the sign bit of 64 lanes into a UInt64 (bit j ⇔ lane j < 0) ──
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
@inline _ob_masklast(w::UInt64, r::Int) = r == 0 ? w : w & ((UInt64(1) << r) - UInt64(1))

# Pack tap code sign plane: bit for output sample n = sign(extb[byteoff+n]). `extb` is
# padded ≥ byteoff+nwb·64 so the last (masked) word can be read whole.
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

# Derive each tap plane from the prompt-extended plane by a funnel bit-shift, instead of a
# separate movemask per tap: `dst[w] = pe[w+off]` in bit terms, for ANY `off ≥ 0`. One movemask
# of the whole code block then feeds all NC taps. `off` splits into a whole-word shift
# `wsh = off ÷ 64` and a sub-word shift `bit = off mod 64`; a tap output word combines source
# words `pe[w+wsh]` and `pe[w+wsh+1]`. The caller must size/zero `pe` so `pe[nwb+wsh+1]` is a
# valid (zero) pad word (see the per-kernel `pepad`).
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

# Mask the last valid word to `r` bits (if partial) and zero the VW-pad words [nwb+1, nwv], so a
# whole-VW-chunk correlate sees zeros (which contribute nothing) past the real samples.
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

# Per-sat measurement sign packing (used only for a single-sat group, where band-sharing +
# realign would add a copy the direct pack avoids). sign(real)→mrb, sign(imag)→mib; full
# 64-sample words via a deinterleaving vector load, scalar tail for the last < 64.
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

# Scratch: block code buffer `extb` (Int8), the packed prompt-extended code plane `peb`, carrier
# sign planes, per-tap code sign planes (`codeb`, NC·wpbv), per-antenna measurement sign planes
# (`mrb`/`mib`, M·wpbv). Grown lazily and reused, so a hoisted backend is allocation-free in steady
# state. Plane word-strides are padded to a multiple of `_OB_VW` (`wpbv`).
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

# `sign(measurement)` is identical for every satellite on a band, so pack the whole band's
# measurement sign planes ONCE per group (`M` planes each, aligned to sample 1, one pad word for
# the funnel realign) and share them across sats. Each sat then funnel-shifts its slice into its
# thread-local block buffer instead of re-deinterleaving/movemasking the Int16 samples. Filled
# serially in `_dc_one_group!` before the per-sat (`@batch`) loop, then read-only across threads.
struct OneBitBandBuffers
    mrband::Vector{UInt64}
    miband::Vector{UInt64}
end
OneBitBandBuffers() = OneBitBandBuffers(UInt64[], UInt64[])

# Pack the band's `M` measurement sign planes over all `num_samples` samples. Plane `j` occupies
# `stride` words at offset `(j-1)*stride`; the last word of each plane is the zero funnel pad.
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

# Funnel-realign antenna `j`'s block (`len` samples starting at 1-based absolute sample `base`)
# from the shared band plane into the thread-local `mrb`/`mib` at word-offset `moff`.
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
`Complex{Int16}` sample buffers: it 1-bit-quantises the measurement, carrier and
code and correlates with XOR + popcount. Construct **once outside** the `track!`
loop and pass it via the `downconvert_and_correlator` keyword for an
allocation-free steady state. 1-bit quantisation trades ≈2–3 dB of SNR for
bit-wise speed; downstream consumers are ratio-normalised, so the coarse
amplitude is immaterial.
"""
struct OneBitDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::OneBitScratchBuffers
    band::OneBitBandBuffers
    blk::Int
end

"""
$(SIGNATURES)

Multi-threaded one-bit bit-wise backend. One `OneBitScratchBuffers` per thread
(indexed by `Threads.threadid()` inside `@batch`). See
[`OneBitDownconvertAndCorrelator`](@ref).
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
# Returns this integration's correlation contribution: `Vector` of `NC` complex
# sums (one per tap) — `ComplexF64` (M=1) or `SVector{M,ComplexF64}` (M>1) — to be
# added to the correlator's running accumulators. `@generated` over (NC, M): the
# per-(antenna, tap) popcount accumulators live in named locals and unroll.
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
    ret =
        M == 1 ? :(ComplexF64[$([tapval(k) for k = 1:NC]...)]) :
        :(SVector{M,ComplexF64}[$([tapval(k) for k = 1:NC]...)])

    quote
        # Bit-wise correlation keeps only the code's SIGN, so it needs a binary (±1,
        # two-level) code — BPSK, BOC and TMBOC all qualify. CBOC is the one GNSS
        # modulation that carries an amplitude (a multi-level weighted sum of two BOCs),
        # so its sign discards information. Gate on the modulation, not the code element
        # type: newer GNSSSignals code generation quantises CBOC to an *integer* replica,
        # so a `get_code_type <: Integer` test no longer excludes it.
        get_modulation(signal_type) isa GNSSSignals.CBOC && throw(
            ArgumentError(
                string(
                    "OneBitDownconvertAndCorrelator supports binary (±1) codes only ",
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
        band_shared = !isempty(band.mrband)         # pre-packed once per group when it pays off
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
            # code (Int8 ±1) for samples [min_shift, len+span); pack the extended plane ONCE, then
            # derive each tap by a funnel bit-shift.
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
            # carrier sign planes (shared across antennas & taps) — SinCosLUT's 1-bit NCO
            generate_carrier_signs!(
                sinw,
                cosw,
                len,
                cps_car;
                phase = carphase + cps_car * blk_off,
            )
            _ob_zeropad!(sinw, 0, nwb, nwv)
            _ob_zeropad!(cosw, 0, nwb, nwv)
            # measurement sign planes, per antenna. >1 sat: funnel-realign from the band-shared
            # planes (packed once per group). 1 sat: pack directly (no realign copy).
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

# Dynamic (runtime tap count) fallback: correlators whose sample shifts are a
# runtime-sized `AbstractVector` (e.g. a `Vector`-accumulator correlator — issue
# #126 (b)). Mirrors the integer backend's `_int16_hybrid_blocked!` AbstractVector
# method: same block pipeline (pack the extended code plane once, funnel-shift each
# tap, carrier + measurement sign planes, XOR + popcount), but loops taps/antennas
# at runtime and horizontally sums each popcount chunk into per-(antenna, tap) Int64
# totals. Not the hot path (it allocates the offset/total scratch and sums per chunk),
# but bit-identical to the `@generated` kernel: `sum` over chunks of `sum(count_ones)`
# equals `sum(count_ones)` over all chunks. Returns `Vector{ComplexF64}` (M=1) or
# `Vector{SVector{M,ComplexF64}}` (M>1). `SVector` shifts are more specific and
# dispatch to the `@generated` method above, so EPL/VEPL keep the fast unrolled path.
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
) where {M}
    get_modulation(signal_type) isa GNSSSignals.CBOC && throw(
        ArgumentError(
            string(
                "OneBitDownconvertAndCorrelator supports binary (±1) codes only ",
                "(BPSK, BOC, TMBOC); got ",
                typeof(signal_type),
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
    band_shared = !isempty(band.mrband)
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
# A satellite carrying several signals on one carrier (e.g. GPS L1 C/A + L1C-D +
# L1C-P) shares the carrier and the measurement — and therefore the whole carrier
# wipe-off. So per strip-mine block we generate the carrier sign planes ONCE and
# pack/realign the measurement sign planes ONCE (per antenna), then correlate each
# signal's own code against them: one carrier+measurement pass per sat, not per
# signal. Bit-identical to correlating each signal alone — the carrier and
# measurement planes are sat-shared, so a signal sees exactly the planes the
# single-signal kernel would build. Mirrors the integer `_int16_hybrid_blocked_multi!`
# and the Float32 `downconvert_and_correlate_fused_tuple!`, adapted to the 1-bit
# XOR+popcount pipeline. `@generated` over (M, the signals tuple): each signal's tap
# count NCᵢ and the M antenna passes unroll; the per-(signal, antenna, tap) popcount
# accumulators live in named locals. Static tap counts (SVector shifts) and any
# antenna M — the scope of the other tile-shares; dynamic (AbstractVector) tap
# counts are not supported here (nor by the integer/Float32 tile-shares).
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

    # Per-signal shift/span/offset locals + per-(signal, antenna, tap) A/B (→I) and
    # C/E (→Q) popcount accumulators.
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

    # CBOC gate, per signal (see the single-signal kernel): bit-wise correlation
    # keeps only the code sign, so a multi-level amplitude-carrying code is rejected.
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
                            typeof(signal_types[$i]),
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

    # Per antenna, measurement packing: band-shared funnel-realign (>1 sat) vs direct
    # pack (1 sat). Filled ONCE per block, shared across the sat's signals.
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

    # Per-signal, within a block: one-shot code fill + pack into `codeb`, then the
    # XOR+popcount correlate reading the shared carrier/measurement planes and
    # accumulating into this signal's A/B/C/E counters.
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

    # Finalize: per signal, a Vector of NCᵢ tap sums — Iₓ = 2N − 2(A+B), Qₓ = 2(E − C)
    # (matching the single-signal kernel's return so `_correlate_signals` is shared).
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

# ── Correlate / plumbing (mirrors the CPU/Int16 backends) ─────────────────────
# Single signal per sat: one kernel call.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal},
    per_signal_completed::Tuple{Bool},
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
    head = signals[1]
    p = _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    )
    new_acc = _onebit_hybrid_blocked!(
        dc,
        signal,
        _ob_num_ants_val(head.correlator),
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

function _update_tracked_sat_correlator(
    sat::TrackedSat,
    dc::_OneBitDC,
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
    dc::_OneBitDC,
    measurements::BandMeasurements,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    m = measurements[get_band_id(g.band)]
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
    # Pack the band's measurement sign planes ONCE, shared across sats — but only for >1 sat:
    # for a single sat the shared pack + per-sat realign copy is slower than packing that one sat
    # directly, so empty the band buffer to select the kernel's per-sat path.
    samples = m.samples
    if length(vals) > 1
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
    dc::OneBitDownconvertAndCorrelator,
    vals,
    args::Vararg{Any,4},
)
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

@inline function _dc_group_loop!(
    dc::OneBitThreadedDownconvertAndCorrelator,
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

Downconvert and correlate all satellites with the one-bit bit-wise backend.
"""
function downconvert_and_correlate(
    dc::_OneBitDC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    downconvert_and_correlate!(dc, measurements, new_track_state)
end

"""
$(SIGNATURES)

In-place one-bit downconvert and correlate. Returns the same `track_state`.
"""
function downconvert_and_correlate!(
    dc::_OneBitDC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    _foreach_group!(_dc_one_group!, track_state.groups, dc, measurements)
    return track_state
end
