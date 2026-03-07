# SIMD helpers shared by the fused downconvert+correlate kernel.

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
