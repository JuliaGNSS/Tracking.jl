"""
$(SIGNATURES)

Downconvert a signal that has only a single dimension, e.g. like
a single antenna channel.
This is for the case that the signal is of type StructArray.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},1},
    signal::StructArray{Complex{ST},1},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    s_re = signal.re
    s_im = signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1
        ds_re[i] = s_re[i] * c_re[i] + s_im[i] * c_im[i]
        ds_im[i] = s_im[i] * c_re[i] - s_re[i] * c_im[i]
    end
    downconverted_signal
end

"""
$(SIGNATURES)

Downconvert a signal that has multiple dimensions, e.g. like
multiple antenna channels.
This is for the case that the signal is of type StructArray.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},2},
    signal::StructArray{Complex{ST},2},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    s_re = signal.re
    s_im = signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1, j = 1:size(s_re, 2)
        ds_re[i, j] = s_re[i, j] * c_re[i] + s_im[i, j] * c_im[i]
        ds_im[i, j] = s_im[i, j] * c_re[i] - s_re[i, j] * c_im[i]
    end
    downconverted_signal
end

"""
$(SIGNATURES)

Downconvert a signal that has only a single dimension, e.g. like
a single antenna channel.
This is for the case that the signal is a regular complex array.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},1},
    signal::AbstractArray{Complex{ST},1},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    signal_real = reinterpret(reshape, ST, signal)
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1
        ds_re[i] = signal_real[1, i] * c_re[i] + signal_real[2, i] * c_im[i]
        ds_im[i] = signal_real[2, i] * c_re[i] - signal_real[1, i] * c_im[i]
    end
    downconverted_signal
end

"""
$(SIGNATURES)

Downconvert a signal that has multiple dimensions, e.g. like
multiple antenna channels.
This is for the case that the signal is a regular complex array.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},2},
    signal::AbstractArray{Complex{ST},2},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    signal_real = reinterpret(reshape, ST, signal)
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1, j = 1:size(ds_re, 2)
        ds_re[i, j] = signal_real[1, i, j] * c_re[i] + signal_real[2, i, j] * c_im[i]
        ds_im[i, j] = signal_real[2, i, j] * c_re[i] - signal_real[1, i, j] * c_im[i]
    end
    downconverted_signal
end

# Fused carrier generation and downconversion: generates the carrier
# on-the-fly during downconversion, eliminating the carrier replica buffer.
# Uses FastSinCos SIMD for fast sincos computation.

using Base.Cartesian: @nexprs
using VectorizationBase: pick_vector_width

# Optimal SIMD width for type T, determined at precompile time from CPU features.
_simd_width(::Type{T}) where {T} = Int(pick_vector_width(T))

# Convert Vec{N,S} to Vec{N,T}. No-op when S == T.
@inline _to_vec(::Type{SIMD.Vec{N,T}}, v::SIMD.Vec{N,T}) where {N,T} = v
@inline _to_vec(::Type{SIMD.Vec{N,T}}, v::SIMD.Vec{N,S}) where {N,T,S} =
    SIMD.Vec{N,T}(ntuple(k -> T(v[k]), N))

# SIMD deinterleave load: load N interleaved complex pairs [re1,im1,re2,im2,...]
# and separate into (re_vec, im_vec) using shufflevector.
@inline @generated function _deinterleave_load(::Type{SIMD.Vec{N,T}}, p::Ptr{ST}, byte_offset::Int) where {N,T,ST}
    re_idx = ntuple(k -> 2(k - 1), N)
    im_idx = ntuple(k -> 2(k - 1) + 1, N)
    quote
        base = Ptr{ST}(p + byte_offset)
        lo = vload(SIMD.Vec{$N,ST}, base)
        hi = vload(SIMD.Vec{$N,ST}, base + $N * sizeof(ST))
        re_raw = shufflevector(lo, hi, Val{$re_idx}())
        im_raw = shufflevector(lo, hi, Val{$im_idx}())
        _to_vec(SIMD.Vec{$N,T}, re_raw), _to_vec(SIMD.Vec{$N,T}, im_raw)
    end
end

# Fused carrier generation + downconversion for interleaved complex array signals.
# Handles both single-antenna (NANT=1) and multi-antenna cases; the compiler
# unrolls the `for j in 1:NANT` loop when NANT is a compile-time constant.
function downconvert!(
    downconverted_signal::StructArray{Complex{T}},
    signal::AbstractArray{Complex{ST}},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples,
    ::NumAnts{NANT},
) where {T,ST,NANT}
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    carrier_freq = T(upreferred(carrier_frequency / Hz))
    sampling_freq = T(upreferred(sampling_frequency / Hz))
    phase0 = T(start_phase)
    two_pi = T(2π)
    freq_ratio = carrier_freq / sampling_freq
    W = _simd_width(T)
    @nexprs 4 u -> off_u = SIMD.Vec{W,T}(ntuple(k -> T(k - 1 + (u - 1) * W), W))
    num_samples_signal = size(signal, 1)
    p_sig = Ptr{ST}(pointer(signal))
    sizeof_ST = sizeof(ST)
    col_stride_d = size(ds_re, 1)
    p_dre = pointer(ds_re)
    p_dim = pointer(ds_im)
    sizeof_T = sizeof(T)
    sig_col_bytes = num_samples_signal * 2 * sizeof_ST
    i = start_sample
    last = start_sample + num_samples - 1
    @inbounds while i + 4W - 1 <= last
        base_idx = T(i - start_sample)
        @nexprs 4 u -> p_u = two_pi * ((base_idx + off_u) * freq_ratio + phase0)
        @nexprs 4 u -> (ci_u, cr_u) = fast_sincos_u100k(p_u)
        row_byte_off = (i - 1) * 2 * sizeof_ST
        for j in 1:NANT
            col_byte_off = (j - 1) * sig_col_bytes + row_byte_off
            @nexprs 4 u -> (sr_u, si_u) = _deinterleave_load(SIMD.Vec{W,T}, p_sig, col_byte_off + (u - 1) * W * 2 * sizeof_ST)
            d_off = ((j - 1) * col_stride_d + i - 1) * sizeof_T
            @nexprs 4 u -> vstore(sr_u * cr_u + si_u * ci_u, p_dre + d_off + (u - 1) * W * sizeof_T)
            @nexprs 4 u -> vstore(si_u * cr_u - sr_u * ci_u, p_dim + d_off + (u - 1) * W * sizeof_T)
        end
        i += 4W
    end
    @inbounds while i + W - 1 <= last
        base_idx = T(i - start_sample)
        phase = two_pi * ((base_idx + off_1) * freq_ratio + phase0)
        c_im, c_re = fast_sincos_u100k(phase)
        row_byte_off = (i - 1) * 2 * sizeof_ST
        for j in 1:NANT
            sr, si = _deinterleave_load(SIMD.Vec{W,T}, p_sig, (j - 1) * sig_col_bytes + row_byte_off)
            d_off = ((j - 1) * col_stride_d + i - 1) * sizeof_T
            vstore(sr * c_re + si * c_im, p_dre + d_off)
            vstore(si * c_re - sr * c_im, p_dim + d_off)
        end
        i += W
    end
    @inbounds while i <= last
        ph = two_pi * (T(i - start_sample) * freq_ratio + phase0)
        c_im_s, c_re_s = sincos(ph)
        for j in 1:NANT
            sig = signal[i, j]
            ds_re[i, j] = T(real(sig)) * c_re_s + T(imag(sig)) * c_im_s
            ds_im[i, j] = T(imag(sig)) * c_re_s - T(real(sig)) * c_im_s
        end
        i += 1
    end
    downconverted_signal
end
