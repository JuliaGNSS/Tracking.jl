function downconvert!(
    downconverted_signal_re::AbstractVector,
    downconverted_signal_im::AbstractVector,
    carrier_re,
    carrier_im,
    signal_re::AbstractVector,
    signal_im::AbstractVector,
    start_sample::Integer,
    num_samples_left::Integer
)
    @avx unroll = 3 for i = start_sample:num_samples_left + start_sample - 1
        downconverted_signal_re[i] = signal_re[i] * carrier_re[i] +
            signal_im[i] * carrier_im[i]
        downconverted_signal_im[i] = signal_im[i] * carrier_re[i] -
            signal_re[i] * carrier_im[i]
    end
end

function downconvert!(
    downconverted_signal_re::AbstractMatrix,
    downconverted_signal_im::AbstractMatrix,
    carrier_re,
    carrier_im,
    signal_re::AbstractMatrix,
    signal_im::AbstractMatrix,
    start_sample::Integer,
    num_samples_left::Integer
)
    @avx unroll = 3 for i = start_sample:num_samples_left + start_sample - 1, j = 1:size(signal_re, 2)
        # Calculate signal * carrier'
        downconverted_signal_re[i, j] = signal_re[i, j] * carrier_re[i] +
            signal_im[i, j] * carrier_im[i]
        downconverted_signal_im[i, j] = signal_im[i, j] * carrier_re[i] -
            signal_re[i, j] * carrier_im[i]
    end
end

# StructArray GPU downconvert function Ants = 1
function downconvert!(
    downconverted_signal_re::CuMatrix{Complex{T}},
    downconverted_signal_im::CuMatrix{Complex{T}},
    carrier_re::CuVector{Complex{T}},
    carrier_im::CuVector{Complex{T}},
    signal_re::CuMatrix{Complex{T}},
    signal_im::CuMatrix{Complex{T}},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat
    idxs = start_sample:start_sample + num_samples_left - 1
    @. @views downconverted_signal_re[idxs] = signal_re[idxs] * carrier_re[idxs] + signal_im[idxs] * carrier_im[idxs]
    @. @views downconverted_signal_im[idxs] = signal_im[idxs] * carrier_re[idxs] - signal_re[idxs] * carrier_im[idxs]
end

# StructArray GPU downconvert function Ants > 1
function downconvert!(
    downconverted_signal_re::CuVector{T},
    downconverted_signal_im::CuVector{T},
    carrier_re::CuVector{T},
    carrier_im::CuVector{T},
    signal_re::CuVector{T},
    signal_im::CuVector{T},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat
    idxs = start_sample:start_sample + num_samples_left - 1
    @. @views downconverted_signal_re[idxs] = signal_re[idxs] * carrier_re[idxs] + signal_im[idxs] * carrier_im[idxs]
    @. @views downconverted_signal_im[idxs] = signal_im[idxs] * carrier_re[idxs] - signal_re[idxs] * carrier_im[idxs]
end

# CuArray GPU downconvert function Ants = 1
function downconvert!(
    downconverted_signal::CuVector{Complex{T}},
    carrier::CuVector{Complex{T}},
    signal::CuVector{Complex{T}},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat    
    @. @views downconverted_signal[start_sample:num_samples_left + start_sample - 1] =
        signal[start_sample:num_samples_left + start_sample - 1] * 
        conj(carrier[start_sample:num_samples_left + start_sample - 1])
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
    downconverted_signal,
    signal,
    carrier_replica,
    start_sample::Integer,
    num_samples_left::Integer
)
    downconvert!(
        downconverted_signal.re,
        downconverted_signal.im,
        carrier_replica.re,
        carrier_replica.im,
        signal.re,
        signal.im,
        start_sample,
        num_samples_left
    )
    downconverted_signal
end
