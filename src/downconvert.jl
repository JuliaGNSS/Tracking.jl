function downconvert!(
    downconverted_signal::StructArray{Complex{T}, 1},
    signal::StructArray{Complex{ST}, 1},
    carrier_replica::StructArray{Complex{T}, 1},
    start_sample::Integer,
    num_samples::Integer
) where {T, ST}
    ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
    s_re = signal.re; s_im = signal.im
    c_re = carrier_replica.re; c_im = carrier_replica.im
    @avx for i = start_sample:num_samples + start_sample - 1
        ds_re[i] = s_re[i] * c_re[i] + s_im[i] * c_im[i]
        ds_im[i] = s_im[i] * c_re[i] - s_re[i] * c_im[i]
    end
    downconverted_signal
end

# This function creates the carrier and downconverts
# directly. It is a little bit faster than first creating the
# replica and then downconverting. This is not in use at the moment
# to be consistent (see next function).
function downconvert!(
    downconverted_signal::StructArray{Complex{T}},
    signal::StructArray{Complex{TS}},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples
) where {T, TS}
    ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
    s_re = signal.re; s_im = signal.im
    carrier_freq = upreferred(carrier_frequency / Hz)
    sampling_freq = upreferred(sampling_frequency / Hz)
    @avx for i = start_sample:start_sample + num_samples - 1
        c_im, c_re = sincos(T(2π) * ((i - start_sample) * T(carrier_freq) / T(sampling_freq) + T(start_phase)))
        ds_re[i] = s_re[i] * c_re + s_im[i] * c_im
        ds_im[i] = s_im[i] * c_re - s_re[i] * c_im
    end
    downconverted_signal
end

# Same as above but for the multiple antenna case. It is faster
# for less than 4 antennas, but slower otherwise. This is not in
# use due to this circumstances.
function downconvert!(
    downconverted_signal::StructArray{Complex{T}, 2},
    signal::StructArray{Complex{TS}, 2},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples
) where {N, T, TS}
    ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
    s_re = signal.re; s_im = signal.im
    carrier_freq = upreferred(carrier_frequency / Hz)
    sampling_freq = upreferred(sampling_frequency / Hz)
    @avx for i = start_sample:start_sample + num_samples - 1
        c_im, c_re = sincos(T(2π) * ((i - start_sample) * T(carrier_freq) / T(sampling_freq) + T(start_phase)))
        for j = 1:size(s_re, 2)
            ds_re[i,j] = s_re[i,j] * c_re + s_im[i,j] * c_im
            ds_im[i,j] = s_im[i,j] * c_re - s_re[i,j] * c_im
        end
    end
    downconverted_signal
end

# StructArray GPU downconvert function Ants = 1
function downconvert!(
    downconverted_signal_re::CuVector{T},
    downconverted_signal_im::CuVector{T},
    carrier_re::CuVector{T},
    carrier_im::CuVector{T},
    signal_re::CuVector{T},
    signal_im::CuVector{T},
    start_sample::Integer,
    num_samples_left::Integer
) where {T <: AbstractFloat}
    sample_range = start_sample:start_sample + num_samples_left - 1
    @. @views downconverted_signal_re[sample_range] = signal_re[sample_range] * carrier_re[sample_range] + signal_im[sample_range] * carrier_im[sample_range]
    @. @views downconverted_signal_im[sample_range] = signal_im[sample_range] * carrier_re[sample_range] - signal_re[sample_range] * carrier_im[sample_range]
end

# StructArray GPU downconvert function Ants > 1
# function downconvert!(
#     downconverted_signal_re::CuVector{T},
#     downconverted_signal_im::CuVector{T},
#     carrier_re::CuVector{T},
#     carrier_im::CuVector{T},
#     signal_re::CuVector{T},
#     signal_im::CuVector{T},
#     start_sample::Integer,
#     num_samples_left::Integer
# ) where T <: AbstractFloat
#     idxs = start_sample:start_sample + num_samples_left - 1
#     @. @views downconverted_signal_re[idxs] = signal_re[idxs] * carrier_re[idxs] + signal_im[idxs] * carrier_im[idxs]
#     @. @views downconverted_signal_im[idxs] = signal_im[idxs] * carrier_re[idxs] - signal_re[idxs] * carrier_im[idxs]
# end

# CuArray GPU downconvert function Ants = 1
function downconvert!(
    downconverted_signal::CuVector{Complex{T}},
    signal::CuVector{Complex{T}},
    carrier::CuVector{Complex{T}},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat
    sample_range = start_sample:num_samples_left + start_sample - 1
    downconverted_signal[sample_range] .= @views signal[sample_range] .* conj.(carrier[sample_range])
    return downconverted_signal
end

# CuArray GPU downconvert function Ants > 1
function downconvert!(
    downconverted_signal::CuMatrix{Complex{T}},
    carrier::CuVector{Complex{T}},
    signal::CuMatrix{Complex{T}},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat
    @. @views downconverted_signal[start_sample:num_samples_left + start_sample - 1] =
        signal[start_sample:num_samples_left + start_sample - 1] * 
        conj(carrier[start_sample:num_samples_left + start_sample - 1])
end

function downconvert!(
    downconverted_signal::StructArray{Complex{T}, 2},
    signal::StructArray{Complex{ST}, 2},
    carrier_replica::StructArray{Complex{T}, 1},
    start_sample::Integer,
    num_samples::Integer
) where {T, ST}
    ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
    s_re = signal.re; s_im = signal.im
    c_re = carrier_replica.re; c_im = carrier_replica.im
    @avx for i = start_sample:num_samples + start_sample - 1, j = 1:size(s_re, 2)
        ds_re[i, j] = s_re[i, j] * c_re[i] + s_im[i, j] * c_im[i]
        ds_im[i, j] = s_im[i, j] * c_re[i] - s_re[i, j] * c_im[i]
    end
    downconverted_signal
end

@static if VERSION >= v"1.6"
    function downconvert!(
        downconverted_signal::StructArray{Complex{T}, 1},
        signal::AbstractArray{Complex{ST}, 1},
        carrier_replica::StructArray{Complex{T}, 1},
        start_sample::Integer,
        num_samples::Integer
    ) where {T, ST}
        signal_real = reinterpret(reshape, ST, signal)
        ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
        c_re = carrier_replica.re; c_im = carrier_replica.im
        @avx for i = start_sample:num_samples + start_sample - 1
            ds_re[i] = signal_real[1, i] * c_re[i] + signal_real[2, i] * c_im[i]
            ds_im[i] = signal_real[2, i] * c_re[i] - signal_real[1, i] * c_im[i]
        end
        downconverted_signal
    end

    function downconvert!(
        downconverted_signal::StructArray{Complex{T}, 2},
        signal::AbstractArray{Complex{ST}, 2},
        carrier_replica::StructArray{Complex{T}, 1},
        start_sample::Integer,
        num_samples::Integer
    ) where {T, ST}
        signal_real = reinterpret(reshape, ST, signal)
        ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
        c_re = carrier_replica.re; c_im = carrier_replica.im
        @avx for i = start_sample:num_samples + start_sample - 1, j = 1:size(signal_real, 3)
            ds_re[i, j] = signal_real[1, i, j] * c_re[i] + signal_real[2, i, j] * c_im[i]
            ds_im[i, j] = signal_real[2, i, j] * c_re[i] - signal_real[1, i, j] * c_im[i]
        end
        downconverted_signal
    end
else
    function downconvert!(
        downconverted_signal,
        signal::AbstractMatrix,
        carrier_replica,
        start_sample::Integer,
        num_samples::Integer
    )    
        sample_range = start_sample:start_sample + num_samples - 1
        downconverted_signal[sample_range,:] .= @view(signal[sample_range,:]) .* conj.(@view(carrier_replica[sample_range]))
        downconverted_signal
    end

    function downconvert!(
        downconverted_signal,
        signal::AbstractVector,
        carrier_replica,
        start_sample::Integer,
        num_samples::Integer
    )    
        sample_range = start_sample:start_sample + num_samples - 1
        downconverted_signal[sample_range,:] .= @view(signal[sample_range]) .* conj.(@view(carrier_replica[sample_range]))
        downconverted_signal
    end
end