# TODO: This function might not be needed! Without @avx it can be merged with
# the two dimensional case
function downconvert!(
    downconverted_signal::StructArray{Complex{T}, 1},
    signal::StructArray{Complex{T}, 1},
    carrier_replica::StructArray{Complex{T}, 1},
    start_sample::Integer,
    num_samples::Integer
) where T
    ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
    s_re = signal.re; s_im = signal.im
    c_re = carrier_replica.re; c_im = carrier_replica.im
    @avx for i = start_sample:num_samples + start_sample - 1
        ds_re[i] = s_re[i] * c_re[i] + s_im[i] * c_im[i]
        ds_im[i] = s_im[i] * c_re[i] - s_re[i] * c_im[i]
    end
    downconverted_signal
end

function downconvert!(
    downconverted_signal::StructArray{Complex{T}, 2},
    signal::StructArray{Complex{T}, 2},
    carrier_replica::StructArray{Complex{T}, 1},
    start_sample::Integer,
    num_samples::Integer
) where T
    ds_re = downconverted_signal.re; ds_im = downconverted_signal.im
    s_re = signal.re; s_im = signal.im
    c_re = carrier_replica.re; c_im = carrier_replica.im
    @avx for i = start_sample:num_samples + start_sample - 1, j = 1:size(s_re, 2)
        ds_re[i, j] = s_re[i, j] * c_re[i] + s_im[i, j] * c_im[i]
        ds_im[i, j] = s_im[i, j] * c_re[i] - s_re[i, j] * c_im[i]
    end
    downconverted_signal
end

function downconvert!(
    downconverted_signal,
    signal,
    carrier_replica,
    start_sample::Integer,
    num_samples::Integer
)    
    sample_range = start_sample:start_sample + num_samples - 1
    downconverted_signal[sample_range] .= @view(signal[sample_range]) .* conj.(@view(carrier_replica[sample_range]))
    downconverted_signal
end
