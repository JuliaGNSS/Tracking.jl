# Fused downconvert + correlate SIMD kernel.
# Generates the carrier on-the-fly, downconverts, and accumulates against the
# code replica in a single pass — keeping downconverted samples in registers
# instead of writing them to an intermediate buffer.

using Base.Cartesian: @nexprs

# Horizontal sum: reduce SIMD vector lanes to a scalar.
@inline _hsum(v::SIMD.Vec) = sum(Tuple(v))

# Build SIMD offset vector for unroll index u (1-based).
# Uses @generated to produce a literal tuple, avoiding ntuple heap allocation on AVX-512.
@inline @generated function _make_offset(::Type{T}, ::Val{W}, ::Val{U}) where {T,W,U}
    elems = [:(T($(k - 1 + (U - 1) * W))) for k in 1:W]
    :(SIMD.Vec{$W,T}(($(elems...),)))
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
    # Block-accumulation period: flush the Float32 lane accumulators into the
    # Float64 running totals every FLUSH_EVERY main-loop iterations. Each main
    # iteration adds 4 terms per lane (the 4× unroll), so a lane accumulates
    # 4·FLUSH_EVERY = 128 Float32 terms before being flushed — small enough that
    # the Float32 partial sum stays accurate (width-difference stays ~1e-8,
    # far below the ~1e-7 floor that matters for tracking), large enough that
    # the Float64 flush (amortised over 4·FLUSH_EVERY·W samples) is negligible.
    FLUSH_EVERY = 32
    # ── Accumulators ──
    # The correlator sum is the sum of ~10^5-10^6 per-sample products. Summing
    # those in Float32 loses precision that *scales with the SIMD width*: each
    # of the W lanes accumulates a different strided subset of the terms, and
    # the per-lane Float32 rounding grows with the partial-sum magnitude, so the
    # final hsum — and the closed tracking loop downstream — depended on the
    # host's vector width (AVX2 width-8 vs AVX-512 width-16, a >100 Hz Doppler
    # swing; issue #152).
    #
    # Fix: two-level (blocked) accumulation. The hot loop keeps the fast
    # full-width Float32 FMA into per-lane block accumulators `f_re/f_im`; every
    # FLUSH_EVERY iterations those are widened and added into Float64 running
    # totals `acc_re/acc_im` and reset to zero. Bounding the Float32 partial sum
    # this way makes the result width-independent to ~1e-9 (vs ~4e-6 unblocked)
    # while keeping Float32 throughput — the product is exact for the ±1 code
    # replica, and only the cheap flush runs in Float64. `s_re/s_im` hold the
    # scalar-remainder tail in Float64.
    acc_init = Expr(:block)
    for j in 1:M, k in 1:NC
        push!(acc_init.args, :($(Symbol("acc_re_$(j)_$(k)")) = z64))
        push!(acc_init.args, :($(Symbol("acc_im_$(j)_$(k)")) = z64))
        push!(acc_init.args, :($(Symbol("f_re_$(j)_$(k)")) = z32))
        push!(acc_init.args, :($(Symbol("f_im_$(j)_$(k)")) = z32))
        push!(acc_init.args, :($(Symbol("s_re_$(j)_$(k)")) = 0.0))
        push!(acc_init.args, :($(Symbol("s_im_$(j)_$(k)")) = 0.0))
    end

    # ── Flush block: acc_{re,im} += widen(f_{re,im}); reset f_{re,im} = 0 ──
    flush_block = Expr(:block)
    for j in 1:M, k in 1:NC
        are = Symbol("acc_re_$(j)_$(k)"); aim = Symbol("acc_im_$(j)_$(k)")
        fre = Symbol("f_re_$(j)_$(k)"); fim = Symbol("f_im_$(j)_$(k)")
        push!(flush_block.args, :($are += _to_vec(SIMD.Vec{W,Float64}, $fre)))
        push!(flush_block.args, :($aim += _to_vec(SIMD.Vec{W,Float64}, $fim)))
        push!(flush_block.args, :($fre = z32))
        push!(flush_block.args, :($fim = z32))
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
            fre = Symbol("f_re_$(j)_$(k)")
            fim = Symbol("f_im_$(j)_$(k)")
            pck = Symbol("p_code_$(k)")
            push!(main_accum.args, quote
                @nexprs 4 u -> begin
                    code_u = vload(SIMD.Vec{W,CT}, $pck + (i - 1 + (u - 1) * W) * sizeof_CT)
                    code_f_u = _to_vec(SIMD.Vec{W,T}, code_u)
                    $fre = muladd(dre_u, code_f_u, $fre)
                    $fim = muladd(dim_u, code_f_u, $fim)
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
            push!(cleanup_accum.args, :($are += _to_vec(SIMD.Vec{W,Float64}, $dre_j * $cf_jk)))
            push!(cleanup_accum.args, :($aim += _to_vec(SIMD.Vec{W,Float64}, $dim_j * $cf_jk)))
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

        z32 = zero(SIMD.Vec{W,T})
        z64 = zero(SIMD.Vec{W,Float64})
        $acc_init
        $shift_init

        off_1 = _make_offset(T, Val{W}(), Val{1}())
        off_2 = _make_offset(T, Val{W}(), Val{2}())
        off_3 = _make_offset(T, Val{W}(), Val{3}())
        off_4 = _make_offset(T, Val{W}(), Val{4}())

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

        blk = 0
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
            blk += 1
            if blk == $FLUSH_EVERY
                $flush_block
                blk = 0
            end
        end
        # Final flush of the partial Float32 block into the Float64 totals.
        $flush_block

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


# ── Shared downconvert-into-tile phase ────────────────────────────────
# Generates the carrier on-the-fly and downconverts `signal` into the SoA
# tile (one `num_samples`-long slice per antenna):
# `tile_re/tile_im[(j-1)*num_samples + idx]`. Shared by the dynamic-shifts
# fused kernel and the tuple tile-share kernel below so the two phase-1
# loops cannot drift (issue #133). `@generated` so the M antenna slices
# unroll via `@nexprs`, exactly as both kernels did before.
@generated function _downconvert_into_tile!(
    tile_re,
    tile_im,
    signal::AbstractArray{Complex{ST}},
    ::NumAnts{M},
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample::Integer,
    num_samples::Integer,
) where {ST,M}
    quote
        T = Float32
        sizeof_ST = sizeof(ST)
        W = $(_simd_width(Float32))

        carrier_freq = T(upreferred(carrier_frequency / Hz))
        phase0 = T(carrier_phase)
        two_pi = T(2π)
        freq_ratio = carrier_freq / T(upreferred(sampling_frequency / Hz))

        num_samples_signal = size(signal, 1)
        p_sig = Ptr{ST}(pointer(signal))
        sig_col_bytes = num_samples_signal * 2 * sizeof_ST

        off_1 = _make_offset(T, Val{W}(), Val{1}())
        two_pi_fr = SIMD.Vec{W,T}(two_pi * freq_ratio)
        init_1 = two_pi * (off_1 * freq_ratio + phase0)

        last = start_sample + num_samples - 1

        i = start_sample
        idx = 1
        @inbounds while i + W - 1 <= last
            base_phase = SIMD.Vec{W,T}(T(i - start_sample))
            phase = muladd(base_phase, two_pi_fr, init_1)
            ci, cr = fast_sincos_u100k(phase)
            row_byte_off = (i - 1) * 2 * sizeof_ST
            Base.Cartesian.@nexprs $M j -> begin
                sr_j, si_j = _deinterleave_load(
                    SIMD.Vec{W,T}, p_sig,
                    (j - 1) * sig_col_bytes + row_byte_off,
                )
                dre_j = sr_j * cr + si_j * ci
                dim_j = si_j * cr - sr_j * ci
                off_j = ((j - 1) * num_samples + idx - 1) * sizeof(T)
                vstore(dre_j, pointer(tile_re) + off_j)
                vstore(dim_j, pointer(tile_im) + off_j)
            end
            i += W
            idx += W
        end
        @inbounds while i <= last
            ph = two_pi * (T(i - start_sample) * freq_ratio + phase0)
            c_im_s, c_re_s = sincos(ph)
            Base.Cartesian.@nexprs $M j -> begin
                sig_j = signal[i, j]
                tile_re[(j - 1) * num_samples + idx] =
                    T(real(sig_j)) * c_re_s + T(imag(sig_j)) * c_im_s
                tile_im[(j - 1) * num_samples + idx] =
                    T(imag(sig_j)) * c_re_s - T(real(sig_j)) * c_im_s
            end
            i += 1
            idx += 1
        end
        return nothing
    end
end

# ── Overload for dynamic-length sample_shifts (AbstractVector) ────────
# Fired when the caller passes a runtime-sized `AbstractVector` of
# shifts; the common EPL/VEPL correlators use `SVector{NC}` and
# dispatch to the `@generated` overload above (which has no tile
# buffers). This overload requires the caller to supply SoA tile
# buffers (`tile_re`, `tile_im`) of length `num_samples * M` and
# element type `Float32` — typically pulled from a long-lived per-
# thread scratch pool so the kernel itself stays allocation-free. Any
# `AbstractVector{Float32}` compatible with `pointer()` and `vstore`
# works (`Vector{Float32}` for direct callers, `ScratchView{Float32}`
# for the threaded `track!` hot path).
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
    tile_re,
    tile_im,
) where {M,ST}
    T = Float32
    num_taps = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    CHUNK = 512  # block-accumulation length; see issue #152

    # Downconvert each antenna into its tile slice
    _downconvert_into_tile!(
        tile_re, tile_im, signal, NumAnts{M}(),
        carrier_frequency, sampling_frequency, carrier_phase,
        start_sample, num_samples,
    )

    # Correlate: tap-outer, antenna-inner with @simd. The code replica is
    # written at absolute index `start_sample` by `gen_code_replica!`, so
    # reads must be offset by `start_sample - 1` (like the in-register and
    # tuple tile-share kernels).
    prev = get_accumulators(correlator)
    new_acc = _mutable_copy(prev)
    @inbounds for k in 1:num_taps
        shift_offset = start_sample - 1 + sample_shifts[k] - min_shift
        for j in 1:M
            # Block accumulation: fast Float32 @simd inner sum over CHUNK-sample
            # blocks, flushed into Float64 totals so the result is accurate and
            # width-independent without paying full-Float64 throughput (#152).
            acc_r = zero(Float64)
            acc_i = zero(Float64)
            ant_off = (j - 1) * num_samples
            n0 = 1
            while n0 <= num_samples
                n1 = ifelse(n0 + CHUNK - 1 < num_samples, n0 + CHUNK - 1, num_samples)
                fr = zero(Float32)
                fi = zero(Float32)
                @simd for n in n0:n1
                    c = code_replica[n + shift_offset]
                    fr += tile_re[ant_off + n] * c
                    fi += tile_im[ant_off + n] * c
                end
                acc_r += Float64(fr)
                acc_i += Float64(fi)
                n0 = n1 + 1
            end
            corr_val = complex(acc_r, acc_i)
            new_acc[k] = _add_antenna(new_acc[k], prev[k], j, corr_val)
        end
    end

    update_accumulator(correlator, _to_immutable(new_acc))
end

# ── Tuple-of-correlators tile-share fused kernel ──────────────────────
# For multi-signal-per-satellite tracking: one downconvert into the
# `tile_re` / `tile_im` SoA tile (one slice per antenna), followed by M
# antenna-outer correlate passes that accumulate into N correlators
# (each with its own code replica and `sample_shifts`). All N correlators
# must agree on the antenna count M.
#
# At N=2 this beats fused-N-times by ~37%; at N=3 by ~51% — see the
# multi-signal-tracking design doc in docs/plans. The single-signal path
# (N=1) is intentionally NOT routed here — it still uses the in-register
# static-shifts kernel above, which is ~24% faster at N=1 because the
# downconverted samples never leave registers.
#
# The downconvert phase mirrors the dynamic-shifts kernel above. The
# correlate phase is `@generated`-unrolled into M sample-outer passes,
# one per antenna; each pass keeps only N·NC accumulators live so the
# inner loop fits in 16 ymm registers (AVX2) without spilling. At M=1
# this collapses to a single pass (same shape as before); at M=2..4 it
# wins 15-31% over a single fused pass over M·N·NC accumulators.
@generated function downconvert_and_correlate_fused_tuple!(
    correlators::Tuple{AbstractCorrelator{M}, Vararg{AbstractCorrelator{M}, NM1}},
    signal::AbstractArray{Complex{ST}},
    code_replicas::Tuple{Any, Vararg{Any, NM1}},
    all_sample_shifts::Tuple{SVector, Vararg{SVector, NM1}},
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample::Integer,
    num_samples::Integer,
    tile_re,
    tile_im,
) where {M, ST, NM1}
    N = NM1 + 1
    # Per-signal NC from each SVector's `length` (a compile-time constant
    # available on the type itself, via the StaticArray Size interface).
    NC_per_signal = Int[length(all_sample_shifts.parameters[i]) for i in 1:N]

    # Block-accumulation chunk length (issue #152). The inner `@simd` reduction
    # runs in Float32 (fast, full-width), but a plain Float32 sum over the whole
    # ~10^5-10^6-sample integration loses precision that depends on the SIMD
    # width. Summing each CHUNK-sample block in Float32 and flushing into a
    # Float64 total keeps the Float32 partial sum short enough to stay accurate
    # and width-independent (~1e-8) while preserving Float32 throughput.
    CHUNK = 512

    # Correlate phase: emit M independent passes over `n`, one per antenna.
    # Each pass keeps only N·NC accumulators live (vs M·N·NC for a single
    # fused pass), which avoids register spilling on 16-ymm AVX2 targets
    # and gives 15-31% speedups at M=2..4 (measured against the
    # single-fused-pass variant during the multi-signal kernel work).
    # Each antenna's tile slice is streamed exactly once (same as the
    # single-pass shape); the code replicas are re-read M times across
    # passes, which costs nothing since each replica fits in L1.
    correlate_passes = Expr(:block)
    for j in 1:M
        # Per-pass accumulator init: ar_j_i_k / ai_j_i_k for THIS antenna.
        # We keep the global naming so `finalize` resolves to these locals.
        pass_init = Expr(:block)
        block_init = Expr(:block)
        flush_block = Expr(:block)
        for i in 1:N, k in 1:NC_per_signal[i]
            ar = Symbol("ar_$(j)_$(i)_$(k)"); ai = Symbol("ai_$(j)_$(i)_$(k)")
            far = Symbol("far_$(j)_$(i)_$(k)"); fai = Symbol("fai_$(j)_$(i)_$(k)")
            # Float64 running totals + Float32 per-block accumulators (issue #152).
            push!(pass_init.args, :($ar = zero(Float64)))
            push!(pass_init.args, :($ai = zero(Float64)))
            push!(block_init.args, :($far = zero(Float32)))
            push!(block_init.args, :($fai = zero(Float32)))
            push!(flush_block.args, :($ar += Float64($far)))
            push!(flush_block.args, :($ai += Float64($fai)))
        end
        # Per-pass locals: code_replica pointers + per-tap shift offsets.
        pass_locals = Expr(:block)
        for i in 1:N
            push!(pass_locals.args, :($(Symbol("cr_$(i)")) = code_replicas[$i]))
            push!(pass_locals.args,
                :($(Symbol("ms_$(i)")) = minimum(all_sample_shifts[$i])))
            for k in 1:NC_per_signal[i]
                push!(pass_locals.args, :(
                    $(Symbol("sh_$(i)_$(k)")) =
                        all_sample_shifts[$i][$k] - $(Symbol("ms_$(i)"))
                ))
            end
        end
        # Inner per-sample body for antenna `j`: one tile load + N·NC FMAs.
        pass_body = Expr(:block)
        push!(pass_body.args, :(tr = tile_re[$(j - 1) * num_samples + n]))
        push!(pass_body.args, :(ti = tile_im[$(j - 1) * num_samples + n]))
        for i in 1:N
            cr = Symbol("cr_$(i)")
            for k in 1:NC_per_signal[i]
                sh = Symbol("sh_$(i)_$(k)")
                # `gen_code_replica!` writes its first sample at index
                # `start_sample` (absolute), so the prompt code for the
                # tile's sample `n` lives at `start_sample - 1 + n` (+tap
                # shift). Indexing from `n` alone silently reads the wrong
                # code for every integration whose `start_sample != 1`
                # (i.e. all but the first per chunk) — see the in-register
                # kernel above, which reads at the same absolute offset.
                push!(pass_body.args, :(c = Float32($cr[start_sample - 1 + n + $sh])))
                push!(pass_body.args, :($(Symbol("far_$(j)_$(i)_$(k)")) =
                    muladd(tr, c, $(Symbol("far_$(j)_$(i)_$(k)")))))
                push!(pass_body.args, :($(Symbol("fai_$(j)_$(i)_$(k)")) =
                    muladd(ti, c, $(Symbol("fai_$(j)_$(i)_$(k)")))))
            end
        end
        push!(correlate_passes.args, quote
            $pass_init
            $pass_locals
            n0 = 1
            @inbounds while n0 <= num_samples
                n1 = ifelse(n0 + $CHUNK - 1 < num_samples, n0 + $CHUNK - 1, num_samples)
                $block_init
                @simd for n in n0:n1
                    $pass_body
                end
                $flush_block
                n0 = n1 + 1
            end
        end)
    end

    # Build the per-tap accumulator expression: scalar Complex for M==1,
    # SVector{M, ComplexF64} for M>1.
    function tap_accum_expr(i, k)
        if M == 1
            return quote
                complex(
                    Float64($(Symbol("ar_1_$(i)_$(k)"))),
                    Float64($(Symbol("ai_1_$(i)_$(k)"))),
                )
            end
        else
            ant_exprs = [
                quote
                    complex(
                        Float64($(Symbol("ar_$(j)_$(i)_$(k)"))),
                        Float64($(Symbol("ai_$(j)_$(i)_$(k)"))),
                    )
                end
                for j in 1:M
            ]
            return :(SVector{$M}(tuple($(ant_exprs...))))
        end
    end

    # Finalize: write each accumulator back into its correlator.
    finalize = Expr(:block)
    push!(finalize.args, :(updated = ()))
    for i in 1:N
        nc = NC_per_signal[i]
        tap_exprs = [tap_accum_expr(i, k) for k in 1:nc]
        push!(finalize.args, :(corr_i = correlators[$i]))
        push!(finalize.args, :(prev_i = get_accumulators(corr_i)))
        push!(finalize.args, :(new_accs_i = SVector{$nc}(($(tap_exprs...),))))
        push!(finalize.args, :(new_corr_i = update_accumulator(corr_i, prev_i .+ new_accs_i)))
        push!(finalize.args, :(updated = (updated..., new_corr_i)))
    end

    quote
        # ── Phase 1: downconvert into the shared tile (M antennas). ──
        _downconvert_into_tile!(
            tile_re, tile_im, signal, NumAnts{$M}(),
            carrier_frequency, sampling_frequency, carrier_phase,
            start_sample, num_samples,
        )

        # ── Phase 2: antenna-outer fused correlate, N·NC accumulators per pass. ──
        $correlate_passes

        # ── Finalize: rebuild N updated correlators. ──
        $finalize
        return updated
    end
end

