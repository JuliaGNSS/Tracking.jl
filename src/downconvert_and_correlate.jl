#=
# This is currently slower than splitting the loop.
# See https://github.com/JuliaSIMD/LoopVectorization.jl/issues/284
function downconvert_and_correlate(
    signal::StructArray{Complex{T}},
    correlator::C,
    code,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples
) where {T, C <: AbstractCorrelator}
    s_re = signal.re; s_im = signal.im
    accumulators = zero_accumulators(get_accumulators(correlator), signal)
    a_re = real.(accumulators)
    a_im = imag.(accumulators)
    @avx for i = start_sample:start_sample + num_samples - 1
        c_im, c_re = sincos(T(2π) * ((i - start_sample) * T(upreferred(carrier_frequency / Hz)) / T(upreferred(sampling_frequency / Hz)) + T(start_phase)))
        d_re = s_re[i] * c_re + s_im[i] * c_im
        d_im = s_im[i] * c_re - s_re[i] * c_im
        for j = 1:length(a_re)
            sample_shift = correlator_sample_shifts[j] - correlator_sample_shifts[1]
            a_re[j] += d_re * code[i + sample_shift]
            a_im[j] += d_im * code[i + sample_shift]
        end
    end
    accumulators_result = complex.(a_re, a_im)
    C(map(+, get_accumulators(correlator), accumulators_result))
end
=#

# CUDA Kernel 
function downconvert_and_correlate_kernel(
    res_re,
    res_im,
    signal_re,
    signal_im,
    carrier_re,
    carrier_im,
    codes,
    code_frequency,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    start_code_phase,
    carrier_phase,
    code_length,
    prn,
    num_samples,
    num_ants,
    num_corrs
)   
    cache = @cuDynamicSharedMem(Float32, (2 * blockDim().x, num_ants, num_corrs))   
    sample_idx   = 1 + ((blockIdx().x - 1) * blockDim().x + (threadIdx().x - 1))
    antenna_idx  = 1 + ((blockIdx().y - 1) * blockDim().y + (threadIdx().y - 1))
    corr_idx     = 1 + ((blockIdx().z - 1) * blockDim().z + (threadIdx().z - 1))
    iq_offset = blockDim().x
    cache_index = threadIdx().x - 1 

    code_phase = accum_re = accum_im = dw_re = dw_im = 0.0f0
    mod_floor_code_phase = Int(0)

    if sample_idx <= num_samples && antenna_idx <= num_ants && corr_idx <= num_corrs
        # generate carrier
        carrier_im[sample_idx], carrier_re[sample_idx] = CUDA.sincos(2π * ((sample_idx - 1) * carrier_frequency / sampling_frequency + carrier_phase))
    
        # downconvert with the conjugate of the carrier
        dw_re = signal_re[sample_idx, antenna_idx] * carrier_re[sample_idx] + signal_im[sample_idx, antenna_idx] * carrier_im[sample_idx]
        dw_im = signal_im[sample_idx, antenna_idx] * carrier_re[sample_idx] - signal_re[sample_idx, antenna_idx] * carrier_im[sample_idx]

        # calculate the code phase
        code_phase = code_frequency / sampling_frequency * ((sample_idx - 1) + correlator_sample_shifts[corr_idx]) + start_code_phase

        # wrap the code phase around the code length e.g. phase = 1024 -> modfloorphase = 1
        mod_floor_code_phase = 1 + mod(floor(Int32, code_phase), code_length)

        # multiply elementwise with the code
        accum_re += codes[mod_floor_code_phase, prn] * dw_re
        accum_im += codes[mod_floor_code_phase, prn] * dw_im
    end

    cache[1 + cache_index + 0 * iq_offset, antenna_idx, corr_idx] = accum_re
    cache[1 + cache_index + 1 * iq_offset, antenna_idx, corr_idx] = accum_im

    ## Reduction
    # wait until all the accumulators have done writing the results to the cache
    sync_threads()

    i::Int = blockDim().x ÷ 2
    @inbounds while i != 0
        if cache_index < i
            cache[1 + cache_index + 0 * iq_offset, antenna_idx, corr_idx] += cache[1 + cache_index + 0 * iq_offset + i, antenna_idx, corr_idx]
            cache[1 + cache_index + 1 * iq_offset, antenna_idx, corr_idx] += cache[1 + cache_index + 1 * iq_offset + i, antenna_idx, corr_idx]
        end
        sync_threads()
        i ÷= 2
    end
    
    if (threadIdx().x - 1) == 0
        res_re[blockIdx().x, antenna_idx, corr_idx] += cache[1 + 0 * iq_offset, antenna_idx, corr_idx]
        res_im[blockIdx().x, antenna_idx, corr_idx] += cache[1 + 1 * iq_offset, antenna_idx, corr_idx]
    end
    return nothing
end

function downconvert_and_correlate_kernel_wrapper(
    system,
    signal,
    correlator,
    code_replica,
    code_phase,
    carrier_replica,
    carrier_phase,
    downconverted_signal,
    code_frequency,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn
)
    num_corrs = length(correlator_sample_shifts)
    num_ants = size(signal, 2)
    num_samples = size(signal, 1)
    block_dim_z = num_corrs
    block_dim_y = num_ants
    # keep num_corrs and num_ants in seperate dimensions, truncate num_samples accordingly to fit
    block_dim_x = prevpow(2, 1024 ÷ block_dim_y ÷ block_dim_z)
    threads = (block_dim_x, block_dim_y, block_dim_z)
    blocks = cld(size(signal, 1), block_dim_x)
    res_re = CUDA.zeros(Float32, blocks, block_dim_y, block_dim_z)
    res_im = CUDA.zeros(Float32, blocks, block_dim_y, block_dim_z)
    shmem_size = sizeof(ComplexF32)*block_dim_x*block_dim_y*block_dim_z
    @cuda threads=threads blocks=blocks shmem=shmem_size downconvert_and_correlate_kernel(
        res_re, 
        res_im, 
        signal.re, 
        signal.im, 
        carrier_replica.carrier.re, 
        carrier_replica.carrier.im, 
        system.codes,
        Float32(code_frequency),
        correlator_sample_shifts,
        Float32(carrier_frequency),
        Float32(sampling_frequency),
        Float32(code_phase),
        Float32(carrier_phase),
        size(system.codes, 1),
        prn,
        num_samples, 
        num_ants,
        num_corrs
    )
    return sum(res_re .+ 1im*res_im, dims=1)
end