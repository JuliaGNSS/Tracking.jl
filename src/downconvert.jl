function downconvert!(
    downconverted_signal_re::Vector{T},
    downconverted_signal_im::Vector{T},
    carrier_re::Vector{T},
    carrier_im::Vector{T},
    signal_re::Vector{T},
    signal_im::Vector{T},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat
    @avx unroll = 3 for i = start_sample:num_samples_left + start_sample - 1
        downconverted_signal_re[i] = signal_re[i] * carrier_re[i] +
            signal_im[i] * carrier_im[i]
        downconverted_signal_im[i] = signal_im[i] * carrier_re[i] -
            signal_re[i] * carrier_im[i]
    end
end

function downconvert!(
    downconverted_signal_re::Matrix{T},
    downconverted_signal_im::Matrix{T},
    carrier_re::Vector{T},
    carrier_im::Vector{T},
    signal_re::Matrix{T},
    signal_im::Matrix{T},
    start_sample::Integer,
    num_samples_left::Integer
) where T <: AbstractFloat
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
    downconverted_signal::StructArray,
    signal::StructArray,
    carrier_replica::StructArray,
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
