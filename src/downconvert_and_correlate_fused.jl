# Fused downconvert + correlate SIMD kernel.
# Generates the carrier on-the-fly, downconverts, and accumulates against the
# code replica in a single pass — keeping downconverted samples in registers
# instead of writing them to an intermediate buffer.

using Base.Cartesian: @nexprs

# Horizontal sum: reduce SIMD vector lanes to a scalar.
@inline _hsum(v::SIMD.Vec) = sum(Tuple(v))

"""
    downconvert_and_correlate_fused!(
        correlator, signal, code_replica, sample_shifts,
        carrier_frequency, sampling_frequency, carrier_phase,
        start_sample, num_samples,
    )

Fused downconvert-and-correlate: generates the carrier replica on-the-fly,
downconverts the signal, and immediately accumulates the result against
shifted code replicas.  Returns an updated correlator of the same type.

Works for both single-antenna (`M=1`) and multi-antenna (`M>1`) systems.
"""
function downconvert_and_correlate_fused!(
    correlator::AbstractCorrelator{M},
    signal::AbstractArray{Complex{ST}},
    code_replica,
    sample_shifts,
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample::Integer,
    num_samples::Integer,
) where {M,ST}
    T = Float32  # computation type (matches existing downconvert buffer)
    CT = eltype(code_replica)
    sizeof_CT = sizeof(CT)
    sizeof_ST = sizeof(ST)
    W = _simd_width(T)
    num_taps = length(sample_shifts)
    min_shift = minimum(sample_shifts)

    carrier_freq = T(upreferred(carrier_frequency / Hz))
    sampling_freq = T(upreferred(sampling_frequency / Hz))
    phase0 = T(carrier_phase)
    two_pi = T(2π)
    freq_ratio = carrier_freq / sampling_freq

    # SIMD lane offsets for 4x-unrolled loop
    @nexprs 4 u -> off_u = SIMD.Vec{W,T}(ntuple(k -> T(k - 1 + (u - 1) * W), W))

    # Signal pointer setup
    num_samples_signal = size(signal, 1)
    p_sig = Ptr{ST}(pointer(signal))
    sig_col_bytes = num_samples_signal * 2 * sizeof_ST

    # Code replica pointer
    p_code = Ptr{CT}(pointer(code_replica))

    # SIMD accumulators: M antennas × num_taps correlator taps
    acc_re = [zero(SIMD.Vec{W,T}) for _ in 1:M, _ in 1:num_taps]
    acc_im = [zero(SIMD.Vec{W,T}) for _ in 1:M, _ in 1:num_taps]

    # Scalar accumulators for the remainder loop
    s_acc_re = zeros(Float64, M, num_taps)
    s_acc_im = zeros(Float64, M, num_taps)

    i = start_sample
    last = start_sample + num_samples - 1

    # ── Main 4×W unrolled loop ──────────────────────────────────────────
    @inbounds while i + 4W - 1 <= last
        base_idx = T(i - start_sample)
        @nexprs 4 u -> p_u = two_pi * ((base_idx + off_u) * freq_ratio + phase0)
        @nexprs 4 u -> (ci_u, cr_u) = fast_sincos_u100k(p_u)
        row_byte_off = (i - 1) * 2 * sizeof_ST
        for j in 1:M
            col_byte_off = (j - 1) * sig_col_bytes + row_byte_off
            @nexprs 4 u -> (sr_u, si_u) = _deinterleave_load(
                SIMD.Vec{W,T}, p_sig,
                col_byte_off + (u - 1) * W * 2 * sizeof_ST,
            )
            @nexprs 4 u -> dre_u = sr_u * cr_u + si_u * ci_u
            @nexprs 4 u -> dim_u = si_u * cr_u - sr_u * ci_u
            for k in 1:num_taps
                shift = sample_shifts[k] - min_shift
                @nexprs 4 u -> begin
                    code_u = vload(
                        SIMD.Vec{W,CT},
                        p_code + (i - 1 + shift + (u - 1) * W) * sizeof_CT,
                    )
                    code_f_u = _to_vec(SIMD.Vec{W,T}, code_u)
                    acc_re[j, k] = muladd(dre_u, code_f_u, acc_re[j, k])
                    acc_im[j, k] = muladd(dim_u, code_f_u, acc_im[j, k])
                end
            end
        end
        i += 4W
    end

    # ── 1×W cleanup loop ───────────────────────────────────────────────
    @inbounds while i + W - 1 <= last
        base_idx = T(i - start_sample)
        phase = two_pi * ((base_idx + off_1) * freq_ratio + phase0)
        c_im, c_re = fast_sincos_u100k(phase)
        row_byte_off = (i - 1) * 2 * sizeof_ST
        for j in 1:M
            sr, si = _deinterleave_load(
                SIMD.Vec{W,T}, p_sig,
                (j - 1) * sig_col_bytes + row_byte_off,
            )
            dre = sr * c_re + si * c_im
            dim = si * c_re - sr * c_im
            for k in 1:num_taps
                shift = sample_shifts[k] - min_shift
                code_v = vload(
                    SIMD.Vec{W,CT},
                    p_code + (i - 1 + shift) * sizeof_CT,
                )
                code_f = _to_vec(SIMD.Vec{W,T}, code_v)
                acc_re[j, k] = muladd(dre, code_f, acc_re[j, k])
                acc_im[j, k] = muladd(dim, code_f, acc_im[j, k])
            end
        end
        i += W
    end

    # ── Scalar remainder loop ──────────────────────────────────────────
    @inbounds while i <= last
        ph = two_pi * (T(i - start_sample) * freq_ratio + phase0)
        c_im_s, c_re_s = sincos(ph)
        for j in 1:M
            sig = signal[i, j]
            sr_s = T(real(sig))
            si_s = T(imag(sig))
            dre_s = sr_s * c_re_s + si_s * c_im_s
            dim_s = si_s * c_re_s - sr_s * c_im_s
            for k in 1:num_taps
                shift = sample_shifts[k] - min_shift
                c_val = T(code_replica[i + shift])
                s_acc_re[j, k] += Float64(dre_s * c_val)
                s_acc_im[j, k] += Float64(dim_s * c_val)
            end
        end
        i += 1
    end

    # ── Final reduction ────────────────────────────────────────────────
    _build_result(correlator, acc_re, acc_im, s_acc_re, s_acc_im, Val(M), num_taps)
end

# Single-antenna specialisation: signal is a Vector
function downconvert_and_correlate_fused!(
    correlator::AbstractCorrelator{1},
    signal::AbstractVector{Complex{ST}},
    code_replica,
    sample_shifts,
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample::Integer,
    num_samples::Integer,
) where {ST}
    # Reshape to 2-D (N×1) so the generic kernel handles it uniformly
    sig2d = reshape(signal, :, 1)
    downconvert_and_correlate_fused!(
        correlator, sig2d, code_replica, sample_shifts,
        carrier_frequency, sampling_frequency, carrier_phase,
        start_sample, num_samples,
    )
end

# ── Result builders ────────────────────────────────────────────────────

# Single-antenna (M == 1): accumulators are SVector{NC, ComplexF64}
function _build_result(
    correlator::AbstractCorrelator{1},
    acc_re, acc_im, s_acc_re, s_acc_im, ::Val{1}, num_taps,
)
    prev = get_accumulators(correlator)
    accumulators_result = ntuple(num_taps) do k
        re = Float64(_hsum(acc_re[1, k])) + s_acc_re[1, k]
        im = Float64(_hsum(acc_im[1, k])) + s_acc_im[1, k]
        complex(re, im)
    end
    new_acc = map(+, prev, SVector(accumulators_result))
    update_accumulator(correlator, new_acc)
end

# Multi-antenna (M > 1): accumulators are SVector{NC, SVector{M, ComplexF64}}
function _build_result(
    correlator::AbstractCorrelator{M},
    acc_re, acc_im, s_acc_re, s_acc_im, ::Val{M}, num_taps,
) where {M}
    prev = get_accumulators(correlator)
    accumulators_result = ntuple(num_taps) do k
        SVector{M,ComplexF64}(ntuple(M) do j
            re = Float64(_hsum(acc_re[j, k])) + s_acc_re[j, k]
            im = Float64(_hsum(acc_im[j, k])) + s_acc_im[j, k]
            complex(re, im)
        end)
    end
    new_acc = map(+, prev, SVector(accumulators_result))
    update_accumulator(correlator, new_acc)
end
