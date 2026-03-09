using KernelAbstractions: @kernel, @index, get_backend, @Const, @uniform, @localmem, @synchronize, @groupsize, @private, synchronize
using GNSSSignals: LOC, BOC, CBOC, BOCsin, BOCcos, get_subcarrier_code, get_modulation

"""
    KADownconvertAndCorrelator

GPU-accelerated downconvert-and-correlate using KernelAbstractions.

Expanded code tables (with subcarrier baked in) live permanently on GPU.
Only a tiny param buffer is copied per call.

Uses Int64 fixed-point arithmetic with 32 fractional bits for essentially
exact subcarrier indexing (zero quantization errors across 100K samples).
"""
struct KADownconvertAndCorrelator{CT,AR<:AbstractArray{ComplexF64,3},AP<:AbstractMatrix{Float32}} <:
       AbstractDownconvertAndCorrelator
    expanded_codes_gpu::CT        # Tuple of GPU Float32 matrices: (code_length*sub_per_chip) × num_prns
    sub_per_chip::Vector{Int32}   # sub-phases per chip for each system (1 for LOC, 12 for CBOC, etc.)
    system_indices::Dict{UInt64,Int32}  # objectid(system) → system_idx in tuple
    results_gpu::AR
    results_cpu::Array{ComplexF64,3}
    sat_params_gpu::AP            # header_rows × max_sats Float32 (tiny param buffer)
    sat_params_cpu::Matrix{Float32}
    num_threads::Int
    max_sats::Int
    max_taps::Int
    header_rows::Int
end

"""
    _get_sub_per_chip(modulation)

Number of subcarrier sub-phases per code chip.
LOC = 1, BOCsin(m,n) = 2m, BOCcos(m,n) = 2m, CBOC = lcm(2m1, 2m2).
"""
_get_sub_per_chip(::LOC) = 1
_get_sub_per_chip(mod::BOCsin) = 2 * mod.m
_get_sub_per_chip(mod::BOCcos) = 2 * mod.m
function _get_sub_per_chip(mod::CBOC)
    lcm(2 * mod.boc1.m, 2 * mod.boc2.m)
end

"""
    _build_expanded_code_table(codes, modulation, code_length)

Build expanded code table with subcarrier baked in.
Entry [(i-1)*sub_per_chip + j, prn] = codes[i, prn] * subcarrier((j-1)/sub_per_chip).
For LOC (sub_per_chip=1), this is just Float32.(codes).
"""
function _build_expanded_code_table(codes::AbstractMatrix, modulation::LOC, code_length::Int)
    Float32.(codes)
end

function _build_expanded_code_table(codes::AbstractMatrix, modulation, code_length::Int)
    spc = _get_sub_per_chip(modulation)
    num_prns = size(codes, 2)
    expanded = Matrix{Float32}(undef, code_length * spc, num_prns)
    for prn in 1:num_prns
        for i in 1:code_length
            chip_val = Float32(codes[i, prn])
            for j in 1:spc
                sc_phase = (j - 1) / spc
                sc_val = Float32(get_subcarrier_code(modulation, sc_phase))
                expanded[(i - 1) * spc + j, prn] = chip_val * sc_val
            end
        end
    end
    expanded
end

"""
    KADownconvertAndCorrelator(systems, ArrayType; kwargs...)

Construct a GPU downconvert-and-correlator using Int64 fixed-point
arithmetic with 32 fractional bits for exact subcarrier indexing.
"""
function KADownconvertAndCorrelator(
    systems,
    ArrayType::Type{<:AbstractArray};
    num_threads::Int = 256,
    max_ants::Int = 1,
    max_taps::Int = 5,
    max_sats::Int = 32,
)
    expanded_codes_gpu = Tuple(
        let mod = get_modulation(typeof(sys)),
            codes = get_codes(sys),
            code_length = get_code_length(sys)
            ArrayType{Float32}(_build_expanded_code_table(codes, mod, code_length))
        end
        for sys in systems
    )
    sub_per_chip_vec = Int32[_get_sub_per_chip(get_modulation(typeof(sys))) for sys in systems]
    system_indices = Dict{UInt64,Int32}(
        objectid(sys) => Int32(i) for (i, sys) in enumerate(systems)
    )

    # Header layout per column (all Float32; integer values bit-packed via reinterpret):
    #
    #   Row 1: carrier_freq_ratio       Float32
    #   Row 2: carrier_phase            Float32
    #   Row 3: num_samples              Int32 reinterpreted
    #   Row 4: start_sample             Int32 reinterpreted
    #   Row 5: prn                      Int32 reinterpreted
    #   Rows 6-7: code_phase_fixed      Int64 as 2×Float32
    #   Rows 8-9: delta_code_phase_fixed Int64 as 2×Float32
    #   Rows 10..9+2*max_taps: tap_shift_fixed  Int64 as 2×Float32 each
    header_rows = 9 + 2 * max_taps

    results_gpu = ArrayType{ComplexF64}(undef, max_sats, max_ants, max_taps)
    results_cpu = Array{ComplexF64}(undef, max_sats, max_ants, max_taps)
    sat_params_gpu = ArrayType{Float32}(undef, header_rows, max_sats)
    sat_params_cpu = Matrix{Float32}(undef, header_rows, max_sats)
    KADownconvertAndCorrelator{typeof(expanded_codes_gpu),typeof(results_gpu),typeof(sat_params_gpu)}(
        expanded_codes_gpu,
        sub_per_chip_vec,
        system_indices,
        results_gpu,
        results_cpu,
        sat_params_gpu,
        sat_params_cpu,
        num_threads,
        max_sats,
        max_taps,
        header_rows,
    )
end

# Helper: store Int32 as Float32 bits (bit-exact round-trip)
@inline _f32_from_i32(x::Int32) = reinterpret(Float32, x)
# Helper: recover Int32 from Float32 bits
@inline _i32_from_f32(x::Float32) = reinterpret(Int32, x)

# Helper: store Int64 as two Float32 values (lo, hi halves via UInt32)
@inline function _f32_pair_from_i64(x::Int64)
    bits = reinterpret(UInt64, x)
    lo = reinterpret(Float32, UInt32(bits & 0xFFFFFFFF))
    hi = reinterpret(Float32, UInt32((bits >> 32) & 0xFFFFFFFF))
    lo, hi
end
# Helper: recover Int64 from two Float32 values
@inline function _i64_from_f32_pair(lo::Float32, hi::Float32)
    lo_bits = UInt64(reinterpret(UInt32, lo))
    hi_bits = UInt64(reinterpret(UInt32, hi))
    reinterpret(Int64, lo_bits | (hi_bits << 32))
end

# Fixed-point constants: 32 fractional bits
const FIXED_POINT_BITS = 32
const FIXED_POINT_SCALE = Int64(1) << FIXED_POINT_BITS  # 4294967296

# ============================================================
# Unified kernel: handles single and multi-antenna via num_ants
# ============================================================
@kernel function ka_dc_kernel!(
    results,
    @Const(signal),
    @Const(expanded_codes),
    ::Val{num_taps},
    ::Val{combined},
    expanded_code_length::Int32,
    num_ants::Int32,
    @Const(sat_params),
    result_offset::Int32,
) where {num_taps, combined}
    N = @uniform @groupsize()[1]
    local_tid = @index(Local, Linear)
    group_idx = @index(Group, Linear)

    shmem_re = @localmem Float64 (combined ? num_taps : 1, 256)
    shmem_im = @localmem Float64 (combined ? num_taps : 1, 256)
    sat_fp = @localmem Float32 (2,)
    sat_ip = @localmem Int32 (3,)
    sat_phase = @localmem Int64 (2,)
    sat_sh = @localmem Int64 (num_taps,)

    if local_tid == 1
        col = result_offset + Int32(group_idx)
        @inbounds sat_fp[1] = sat_params[1, col]
        @inbounds sat_fp[2] = sat_params[2, col]
        @inbounds sat_ip[1] = _i32_from_f32(sat_params[3, col])
        @inbounds sat_ip[2] = _i32_from_f32(sat_params[4, col])
        @inbounds sat_ip[3] = _i32_from_f32(sat_params[5, col])
        @inbounds sat_phase[1] = _i64_from_f32_pair(sat_params[6, col], sat_params[7, col])
        @inbounds sat_phase[2] = _i64_from_f32_pair(sat_params[8, col], sat_params[9, col])
        @inbounds for t in 1:num_taps
            sat_sh[t] = _i64_from_f32_pair(sat_params[9 + 2*(t-1) + 1, col], sat_params[9 + 2*(t-1) + 2, col])
        end
    end
    @synchronize

    acc_re = @private Float32 (num_taps,)
    acc_im = @private Float32 (num_taps,)

    for ant in Int32(1):num_ants
        for tap in 1:num_taps
            @inbounds acc_re[tap] = zero(Float32)
            @inbounds acc_im[tap] = zero(Float32)
        end

        @inbounds begin
            initial_code_fixed = sat_phase[1]
            delta_code_fixed = sat_phase[2]
            prn = sat_ip[3]

            ecl_shifted = Int64(expanded_code_length) << Int64(FIXED_POINT_BITS)
            delta_stride = delta_code_fixed * Int64(N)

            code_phase = mod((Int64(local_tid) - Int64(1)) * delta_code_fixed + initial_code_fixed, ecl_shifted)

            @fastmath begin
                step_phase = Float32(6.283185307179586) * sat_fp[1] * Float32(N)
                s_step, c_step = sincos(step_phase)
                init_phase = Float32(6.283185307179586) * (Float32(Int32(local_tid) - Int32(1)) * sat_fp[1] + sat_fp[2])
                s, c = sincos(init_phase)
            end

            sample = Int32(local_tid)
            while sample <= sat_ip[1]
                sig = signal[sample + sat_ip[2] - Int32(1), ant]
                sig_re = real(sig)
                sig_im = imag(sig)
                dc_re = sig_re * c + sig_im * s
                dc_im = sig_im * c - sig_re * s

                for tap in 1:num_taps
                    tap_phase = code_phase + sat_sh[tap]
                    tap_phase -= (tap_phase >= ecl_shifted) * ecl_shifted
                    tap_phase += (tap_phase < Int64(0)) * ecl_shifted
                    idx = Int32(tap_phase >> Int64(FIXED_POINT_BITS)) + Int32(1)
                    code_val = expanded_codes[idx, prn]
                    acc_re[tap] += dc_re * code_val
                    acc_im[tap] += dc_im * code_val
                end

                code_phase += delta_stride
                code_phase -= (code_phase >= ecl_shifted) * ecl_shifted

                @fastmath begin
                    c_new = c * c_step - s * s_step
                    s_new = s * c_step + c * s_step
                    c = c_new
                    s = s_new
                end
                sample += Int32(N)
            end
        end

        # Reduction
        if combined
            @inbounds for tap in 1:num_taps
                shmem_re[tap, local_tid] = Float64(acc_re[tap])
                shmem_im[tap, local_tid] = Float64(acc_im[tap])
            end
            for half in (128, 64, 32, 16, 8, 4, 2)
                @synchronize
                if local_tid <= half
                    @inbounds for tap in 1:num_taps
                        shmem_re[tap, local_tid] += shmem_re[tap, local_tid + half]
                        shmem_im[tap, local_tid] += shmem_im[tap, local_tid + half]
                    end
                end
            end
            @synchronize
            if local_tid == 1
                col = result_offset + Int32(group_idx)
                @inbounds for tap in 1:num_taps
                    results[col, ant, tap] = complex(
                        shmem_re[tap, 1] + shmem_re[tap, 2],
                        shmem_im[tap, 1] + shmem_im[tap, 2],
                    )
                end
            end
        else
            for tap in 1:num_taps
                @inbounds shmem_re[1, local_tid] = Float64(acc_re[tap])
                @inbounds shmem_im[1, local_tid] = Float64(acc_im[tap])
                for half in (128, 64, 32, 16, 8, 4, 2)
                    @synchronize
                    if local_tid <= half
                        @inbounds shmem_re[1, local_tid] += shmem_re[1, local_tid + half]
                        @inbounds shmem_im[1, local_tid] += shmem_im[1, local_tid + half]
                    end
                end
                @synchronize
                if local_tid == 1
                    col = result_offset + Int32(group_idx)
                    @inbounds results[col, ant, tap] = complex(
                        shmem_re[1, 1] + shmem_re[1, 2],
                        shmem_im[1, 1] + shmem_im[1, 2],
                    )
                end
                @synchronize
            end
        end
    end
end

# ============================================================
# CPU-side result update helpers
# ============================================================
function ka_update_from_results(
    correlator::AbstractCorrelator{1},
    results,
    batch_idx,
)
    accumulators = get_accumulators(correlator)
    new_accums = map(enumerate(accumulators)) do (tap, prev)
        @inbounds prev + results[batch_idx, 1, tap]
    end
    update_accumulator(correlator, new_accums)
end

function ka_update_from_results(
    correlator::AbstractCorrelator{M},
    results,
    batch_idx,
) where {M}
    accumulators = get_accumulators(correlator)
    new_accums = map(enumerate(accumulators)) do (tap, prev)
        ant_vals = ntuple(M) do ant
            @inbounds results[batch_idx, ant, tap]
        end
        prev + SVector{M}(ant_vals)
    end
    update_accumulator(correlator, new_accums)
end

# ============================================================
# Param packing (Int64 expanded units)
# ============================================================
function _pack_sat_params!(
    sat_params_cpu, batch_idx,
    carrier_freq_ratio, carrier_phase,
    signal_samples_to_integrate, signal_start_sample, prn,
    code_freq_ratio, norm_phase, code_length, spc,
    sample_shifts, num_taps,
)
    expanded_freq_ratio = code_freq_ratio * spc
    expanded_norm_phase = norm_phase * spc

    delta_code_fixed = round(Int64, expanded_freq_ratio * FIXED_POINT_SCALE)
    code_phase_fixed = round(Int64, expanded_norm_phase * FIXED_POINT_SCALE)

    @inbounds sat_params_cpu[1, batch_idx] = carrier_freq_ratio
    @inbounds sat_params_cpu[2, batch_idx] = carrier_phase
    @inbounds sat_params_cpu[3, batch_idx] = _f32_from_i32(Int32(signal_samples_to_integrate))
    @inbounds sat_params_cpu[4, batch_idx] = _f32_from_i32(Int32(signal_start_sample))
    @inbounds sat_params_cpu[5, batch_idx] = _f32_from_i32(Int32(prn))

    lo, hi = _f32_pair_from_i64(code_phase_fixed)
    @inbounds sat_params_cpu[6, batch_idx] = lo
    @inbounds sat_params_cpu[7, batch_idx] = hi
    lo, hi = _f32_pair_from_i64(delta_code_fixed)
    @inbounds sat_params_cpu[8, batch_idx] = lo
    @inbounds sat_params_cpu[9, batch_idx] = hi

    for t in 1:num_taps
        tap_shift_expanded = Float64(sample_shifts[t]) * expanded_freq_ratio
        tap_shift_fixed = round(Int64, tap_shift_expanded * FIXED_POINT_SCALE)
        lo, hi = _f32_pair_from_i64(tap_shift_fixed)
        @inbounds sat_params_cpu[9 + 2*(t-1) + 1, batch_idx] = lo
        @inbounds sat_params_cpu[9 + 2*(t-1) + 2, batch_idx] = hi
    end
end

# ============================================================
# Main dispatch: downconvert_and_correlate
# ============================================================
function downconvert_and_correlate(
    dc::KADownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    backend = get_backend(first(dc.expanded_codes_gpu))

    batch_idx = 0
    M = get_num_ants(track_state)
    header_rows = dc.header_rows

    system_data = map(track_state.multiple_system_sats_state) do system_sats_state
        system = system_sats_state.system
        sys_idx = dc.system_indices[objectid(system)]
        states = system_sats_state.states
        code_length = get_code_length(system)
        spc = dc.sub_per_chip[sys_idx]

        batch_offset = batch_idx

        sat_infos = map(states) do sat_state
            signal_samples_to_integrate, is_integration_completed =
                calc_signal_samples_to_integrate(
                    system,
                    sat_state.signal_start_sample,
                    sampling_frequency,
                    sat_state.code_doppler,
                    sat_state.code_phase,
                    preferred_num_code_blocks_to_integrate,
                    has_bit_or_secondary_code_been_found(sat_state),
                    num_samples_signal,
                )
            carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
            code_frequency = sat_state.code_doppler + get_code_frequency(system)
            sample_shifts = get_correlator_sample_shifts(
                sat_state.correlator,
                sampling_frequency,
                code_frequency,
            )
            (;
                signal_samples_to_integrate,
                is_integration_completed,
                carrier_frequency,
                code_frequency,
                sample_shifts,
            )
        end

        key_to_batch = Dict{eltype(keys(states)),Int}()
        system_num_taps = 0

        for key in keys(states)
            info = sat_infos[key]
            info.signal_samples_to_integrate == 0 && continue

            batch_idx += 1
            key_to_batch[key] = batch_idx
            sat_state = states[key]
            num_taps = length(info.sample_shifts)
            system_num_taps = num_taps

            code_freq_ratio = Float64(upreferred(info.code_frequency / Hz) / upreferred(sampling_frequency / Hz))
            norm_phase = mod(Float64(sat_state.code_phase), code_length)

            carrier_freq_ratio = Float32(
                upreferred(info.carrier_frequency / Hz) /
                upreferred(sampling_frequency / Hz),
            )

            _pack_sat_params!(
                dc.sat_params_cpu, batch_idx,
                carrier_freq_ratio, Float32(sat_state.carrier_phase),
                info.signal_samples_to_integrate, sat_state.signal_start_sample, sat_state.prn,
                code_freq_ratio, norm_phase, code_length, spc,
                info.sample_shifts, num_taps,
            )
        end

        num_active = batch_idx - batch_offset

        (; system_sats_state, system, sys_idx, sat_infos, key_to_batch, batch_offset, num_active, system_num_taps,
           code_length, spc)
    end

    num_total_active = batch_idx

    if num_total_active == 0
        return track_state
    end

    # Tiny copy to GPU
    n = num_total_active
    elems = header_rows * n
    copyto!(vec(dc.sat_params_gpu), 1, vec(dc.sat_params_cpu), 1, elems)

    num_threads = dc.num_threads

    # Per-system kernel launches
    for sd in system_data
        sd.num_active == 0 && continue

        codes = dc.expanded_codes_gpu[sd.sys_idx]
        combined = Val(sd.system_num_taps <= 8)
        expanded_code_length = Int32(sd.code_length * sd.spc)

        ka_dc_kernel!(backend, (num_threads,))(
            dc.results_gpu, signal, codes,
            Val(sd.system_num_taps), combined,
            expanded_code_length, Int32(M),
            dc.sat_params_gpu, Int32(sd.batch_offset);
            ndrange = num_threads * sd.num_active,
        )
    end

    # Synchronize + copy results
    synchronize(backend)
    copyto!(dc.results_cpu, dc.results_gpu)

    # Distribute results back
    new_multiple_system_sats_state = map(system_data) do sd
        states = sd.system_sats_state.states
        new_sat_states = map(pairs(states)) do (key, sat_state)
            batch_pos = get(sd.key_to_batch, key, nothing)
            if batch_pos === nothing
                return sat_state
            end
            info = sd.sat_infos[key]
            new_correlator = ka_update_from_results(
                sat_state.correlator,
                dc.results_cpu,
                batch_pos,
            )::typeof(sat_state.correlator)
            return update(
                sd.system,
                sat_state,
                info.signal_samples_to_integrate,
                intermediate_frequency,
                sampling_frequency,
                new_correlator,
                info.is_integration_completed,
            )
        end
        return SystemSatsState(sd.system_sats_state, new_sat_states)
    end

    return TrackState(
        track_state;
        multiple_system_sats_state = new_multiple_system_sats_state,
    )
end
