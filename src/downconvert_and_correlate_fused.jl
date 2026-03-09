# Fused downconvert + correlate SIMD kernel.
# Generates the carrier on-the-fly, downconverts, and accumulates against the
# code replica in a single pass — keeping downconverted samples in registers
# instead of writing them to an intermediate buffer.

using Base.Cartesian: @nexprs

# Horizontal sum: reduce SIMD vector lanes to a scalar.
@inline _hsum(v::SIMD.Vec) = sum(Tuple(v))

# Build SIMD offset vector for unroll index u (1-based).
@inline function _make_offset(::Type{T}, W::Int, u::Int) where {T}
    SIMD.Vec{W,T}(ntuple(k -> T(k - 1 + (u - 1) * W), W))
end

# Stack-friendly mutable copy / immutable conversion for accumulators.
# SVector → MVector (stack-allocated), Vector → copy (heap).
@inline _mutable_copy(v::SVector) = MVector(v)
@inline _mutable_copy(v::Vector) = copy(v)
@inline _to_immutable(v::MVector) = SVector(v)
@inline _to_immutable(v::Vector) = v

# Add a single antenna's correlation result to accumulator element.
# M==1: scalar complex; M>1: SVector with one antenna updated.
@inline _add_antenna(acc::Complex, prev::Complex, j::Int, val::ComplexF64) = prev + val
@inline function _add_antenna(acc::SVector{M}, prev::SVector{M}, j::Int, val::ComplexF64) where {M}
    Base.setindex(acc, acc[j] + val, j)
end

"""
    downconvert_and_correlate_fused!(
        correlator, signal, code_replica, sample_shifts,
        carrier_frequency, sampling_frequency, carrier_phase,
        start_sample, num_samples,
    )

Fused downconvert-and-correlate: generates the carrier replica on-the-fly,
downconverts the signal, and immediately accumulates the result against
shifted code replicas.  Returns an updated correlator of the same type.

Both antennas (M) and taps (NC) are fully unrolled at compile time via
`@generated`, so all accumulators live in named local variables (registers).
"""
@generated function downconvert_and_correlate_fused!(
    correlator::AbstractCorrelator{M},
    signal::AbstractArray{Complex{ST}},
    code_replica,
    sample_shifts::SVector{NC},
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample::Integer,
    num_samples::Integer,
) where {M,ST,NC}
    W = _simd_width(Float32)  # Compute at generation time, embed as literal
    # ── Accumulator init: acc_re_j_k, acc_im_j_k (SIMD), s_re_j_k, s_im_j_k (scalar) ──
    acc_init = Expr(:block)
    for j in 1:M, k in 1:NC
        push!(acc_init.args, :($(Symbol("acc_re_$(j)_$(k)")) = z))
        push!(acc_init.args, :($(Symbol("acc_im_$(j)_$(k)")) = z))
        push!(acc_init.args, :($(Symbol("s_re_$(j)_$(k)")) = 0.0))
        push!(acc_init.args, :($(Symbol("s_im_$(j)_$(k)")) = 0.0))
    end

    # ── Precompute per-tap shifted code pointers ──
    shift_init = Expr(:block)
    for k in 1:NC
        push!(shift_init.args, :($(Symbol("p_code_$(k)")) = p_code + (sample_shifts[$k] - min_shift) * sizeof_CT))
    end

    # ── Main 4x-unrolled accumulate block ──
    main_accum = Expr(:block)
    for j in 1:M
        push!(main_accum.args, :(col_byte_off = $(j - 1) * sig_col_bytes + row_byte_off))
        push!(main_accum.args, :(@nexprs 4 u -> (sr_u, si_u) = _deinterleave_load(
            SIMD.Vec{W,T}, p_sig,
            col_byte_off + (u - 1) * W * 2 * sizeof_ST,
        )))
        push!(main_accum.args, :(@nexprs 4 u -> dre_u = sr_u * cr_u + si_u * ci_u))
        push!(main_accum.args, :(@nexprs 4 u -> dim_u = si_u * cr_u - sr_u * ci_u))
        for k in 1:NC
            are = Symbol("acc_re_$(j)_$(k)")
            aim = Symbol("acc_im_$(j)_$(k)")
            pck = Symbol("p_code_$(k)")
            push!(main_accum.args, quote
                @nexprs 4 u -> begin
                    code_u = vload(SIMD.Vec{W,CT}, $pck + (i - 1 + (u - 1) * W) * sizeof_CT)
                    code_f_u = _to_vec(SIMD.Vec{W,T}, code_u)
                    $are = muladd(dre_u, code_f_u, $are)
                    $aim = muladd(dim_u, code_f_u, $aim)
                end
            end)
        end
    end

    # ── 1x cleanup accumulate block ──
    cleanup_accum = Expr(:block)
    for j in 1:M
        sr_j = Symbol("sr_c_$(j)")
        si_j = Symbol("si_c_$(j)")
        dre_j = Symbol("dre_c_$(j)")
        dim_j = Symbol("dim_c_$(j)")
        push!(cleanup_accum.args, :(($sr_j, $si_j) = _deinterleave_load(
            SIMD.Vec{W,T}, p_sig,
            $(j - 1) * sig_col_bytes + row_byte_off,
        )))
        push!(cleanup_accum.args, :($dre_j = $sr_j * c_re + $si_j * c_im))
        push!(cleanup_accum.args, :($dim_j = $si_j * c_re - $sr_j * c_im))
        for k in 1:NC
            are = Symbol("acc_re_$(j)_$(k)")
            aim = Symbol("acc_im_$(j)_$(k)")
            cf_jk = Symbol("code_f_c_$(j)_$(k)")
            pck = Symbol("p_code_$(k)")
            push!(cleanup_accum.args, :($cf_jk = _to_vec(SIMD.Vec{W,T}, vload(SIMD.Vec{W,CT}, $pck + (i - 1) * sizeof_CT))))
            push!(cleanup_accum.args, :($are = muladd($dre_j, $cf_jk, $are)))
            push!(cleanup_accum.args, :($aim = muladd($dim_j, $cf_jk, $aim)))
        end
    end

    # ── Scalar remainder accumulate block ──
    scalar_accum = Expr(:block)
    for j in 1:M
        sig_j = Symbol("sig_s_$(j)")
        sr_j = Symbol("sr_s_$(j)")
        si_j = Symbol("si_s_$(j)")
        dre_j = Symbol("dre_s_$(j)")
        dim_j = Symbol("dim_s_$(j)")
        push!(scalar_accum.args, :($sig_j = signal[i, $j]))
        push!(scalar_accum.args, :($sr_j = T(real($sig_j))))
        push!(scalar_accum.args, :($si_j = T(imag($sig_j))))
        push!(scalar_accum.args, :($dre_j = $sr_j * c_re_s + $si_j * c_im_s))
        push!(scalar_accum.args, :($dim_j = $si_j * c_re_s - $sr_j * c_im_s))
        for k in 1:NC
            sre = Symbol("s_re_$(j)_$(k)")
            sim = Symbol("s_im_$(j)_$(k)")
            cv_jk = Symbol("cv_s_$(j)_$(k)")
            pck = Symbol("p_code_$(k)")
            push!(scalar_accum.args, :($cv_jk = T(unsafe_load(Ptr{CT}($pck + (i - 1) * sizeof_CT)))))
            push!(scalar_accum.args, :($sre += Float64($dre_j * $cv_jk)))
            push!(scalar_accum.args, :($sim += Float64($dim_j * $cv_jk)))
        end
    end

    # ── Final reduction: build result tuple directly ──
    if M == 1
        tap_exprs = [quote
            complex(
                Float64(_hsum($(Symbol("acc_re_1_$(k)")))) + $(Symbol("s_re_1_$(k)")),
                Float64(_hsum($(Symbol("acc_im_1_$(k)")))) + $(Symbol("s_im_1_$(k)")),
            )
        end for k in 1:NC]
    else
        tap_exprs = [begin
            ant_exprs = [quote
                complex(
                    Float64(_hsum($(Symbol("acc_re_$(j)_$(k)")))) + $(Symbol("s_re_$(j)_$(k)")),
                    Float64(_hsum($(Symbol("acc_im_$(j)_$(k)")))) + $(Symbol("s_im_$(j)_$(k)")),
                )
            end for j in 1:M]
            :(SVector(tuple($(ant_exprs...))))
        end for k in 1:NC]
    end

    result_expr = quote
        accumulators_result = SVector(tuple($(tap_exprs...)))
        prev = get_accumulators(correlator)
        new_acc = map(+, prev, accumulators_result)
        update_accumulator(correlator, new_acc)
    end

    quote
        T = Float32
        CT = eltype(code_replica)
        sizeof_CT = sizeof(CT)
        sizeof_ST = sizeof(ST)
        W = $W
        min_shift = minimum(sample_shifts)

        carrier_freq = T(upreferred(carrier_frequency / Hz))
        sampling_freq = T(upreferred(sampling_frequency / Hz))
        phase0 = T(carrier_phase)
        two_pi = T(2π)
        freq_ratio = carrier_freq / sampling_freq

        num_samples_signal = size(signal, 1)
        p_sig = Ptr{ST}(pointer(signal))
        sig_col_bytes = num_samples_signal * 2 * sizeof_ST
        p_code = Ptr{CT}(pointer(code_replica))

        z = zero(SIMD.Vec{W,T})
        $acc_init
        $shift_init

        off_1 = _make_offset(T, W, 1)
        off_2 = _make_offset(T, W, 2)
        off_3 = _make_offset(T, W, 3)
        off_4 = _make_offset(T, W, 4)

        # Phase via muladd: p_u = base_phase + init_u where
        #   base_phase = 2π·freq_ratio·(i - start_sample)  (scalar, broadcast)
        #   init_u     = 2π·(off_u·freq_ratio + phase0)    (precomputed per-vector)
        # One muladd per vector, no drift, no per-iteration multiply chain.
        two_pi_fr = SIMD.Vec{W,T}(two_pi * freq_ratio)
        init_1 = two_pi * (off_1 * freq_ratio + phase0)
        init_2 = two_pi * (off_2 * freq_ratio + phase0)
        init_3 = two_pi * (off_3 * freq_ratio + phase0)
        init_4 = two_pi * (off_4 * freq_ratio + phase0)

        i = start_sample
        last = start_sample + num_samples - 1

        @inbounds while i + 4W - 1 <= last
            base_phase = SIMD.Vec{W,T}(T(i - start_sample))
            p_1 = muladd(base_phase, two_pi_fr, init_1)
            p_2 = muladd(base_phase, two_pi_fr, init_2)
            p_3 = muladd(base_phase, two_pi_fr, init_3)
            p_4 = muladd(base_phase, two_pi_fr, init_4)
            ci_1, cr_1 = fast_sincos_u100k(p_1)
            ci_2, cr_2 = fast_sincos_u100k(p_2)
            ci_3, cr_3 = fast_sincos_u100k(p_3)
            ci_4, cr_4 = fast_sincos_u100k(p_4)
            row_byte_off = (i - 1) * 2 * sizeof_ST
            $main_accum
            i += 4W
        end

        @inbounds while i + W - 1 <= last
            base_phase = SIMD.Vec{W,T}(T(i - start_sample))
            phase = muladd(base_phase, two_pi_fr, init_1)
            c_im, c_re = fast_sincos_u100k(phase)
            row_byte_off = (i - 1) * 2 * sizeof_ST
            $cleanup_accum
            i += W
        end

        @inbounds while i <= last
            ph = two_pi * (T(i - start_sample) * freq_ratio + phase0)
            c_im_s, c_re_s = sincos(ph)
            $scalar_accum
            i += 1
        end

        $result_expr
    end
end


# ── Fallback for dynamic-length sample_shifts (AbstractVector) ────────
# Fused downconvert + correlate using stack-allocated SoA tile buffers
# (via @no_escape) and @simd for the tap accumulation loop.
# No dependency on @avx / LoopVectorization.
function downconvert_and_correlate_fused!(
    correlator::AbstractCorrelator{M},
    signal::AbstractArray{Complex{ST}},
    code_replica,
    sample_shifts::AbstractVector,
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample::Integer,
    num_samples::Integer,
) where {M,ST}
    T = Float32
    sizeof_ST = sizeof(ST)
    W = _simd_width(T)
    num_taps = length(sample_shifts)
    min_shift = minimum(sample_shifts)

    carrier_freq = T(upreferred(carrier_frequency / Hz))
    phase0 = T(carrier_phase)
    two_pi = T(2π)
    freq_ratio = carrier_freq / T(upreferred(sampling_frequency / Hz))

    num_samples_signal = size(signal, 1)
    p_sig = Ptr{ST}(pointer(signal))
    sig_col_bytes = num_samples_signal * 2 * sizeof_ST

    off_1 = _make_offset(T, W, 1)
    two_pi_fr = SIMD.Vec{W,T}(two_pi * freq_ratio)
    init_1 = two_pi * (off_1 * freq_ratio + phase0)

    last = start_sample + num_samples - 1

    @no_escape begin
        # Flat buffer: M antennas × num_samples, SoA layout
        tile_re = @alloc(T, num_samples * M)
        tile_im = @alloc(T, num_samples * M)

        # Downconvert each antenna into its tile slice
        i = start_sample
        idx = 1
        @inbounds while i + W - 1 <= last
            base_phase = SIMD.Vec{W,T}(T(i - start_sample))
            phase = muladd(base_phase, two_pi_fr, init_1)
            ci, cr = fast_sincos_u100k(phase)
            row_byte_off = (i - 1) * 2 * sizeof_ST
            for j in 1:M
                sr, si = _deinterleave_load(
                    SIMD.Vec{W,T}, p_sig,
                    (j - 1) * sig_col_bytes + row_byte_off,
                )
                dre = sr * cr + si * ci
                dim = si * cr - sr * ci
                off = ((j - 1) * num_samples + idx - 1) * sizeof(T)
                vstore(dre, pointer(tile_re) + off)
                vstore(dim, pointer(tile_im) + off)
            end
            i += W
            idx += W
        end
        @inbounds while i <= last
            ph = two_pi * (T(i - start_sample) * freq_ratio + phase0)
            c_im_s, c_re_s = sincos(ph)
            for j in 1:M
                sig = signal[i, j]
                tile_re[(j-1)*num_samples + idx] = T(real(sig)) * c_re_s + T(imag(sig)) * c_im_s
                tile_im[(j-1)*num_samples + idx] = T(imag(sig)) * c_re_s - T(real(sig)) * c_im_s
            end
            i += 1
            idx += 1
        end

        # Correlate: tap-outer, antenna-inner with @simd
        prev = get_accumulators(correlator)
        new_acc = _mutable_copy(prev)
        @inbounds for k in 1:num_taps
            shift_offset = sample_shifts[k] - min_shift
            for j in 1:M
                acc_r = zero(T)
                acc_i = zero(T)
                ant_off = (j - 1) * num_samples
                @simd for n in 1:num_samples
                    c = code_replica[n + shift_offset]
                    acc_r += tile_re[ant_off + n] * c
                    acc_i += tile_im[ant_off + n] * c
                end
                corr_val = complex(Float64(acc_r), Float64(acc_i))
                new_acc[k] = _add_antenna(new_acc[k], prev[k], j, corr_val)
            end
        end

        update_accumulator(correlator, _to_immutable(new_acc))
    end
end

