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
    @avx for i = start_sample:num_samples_left + start_sample - 1
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
    @avx for i = start_sample:num_samples_left + start_sample - 1, j = 1:size(signal_re, 2)
        # Calculate signal * carrier'
        downconverted_signal_re[i, j] = signal_re[i, j] * carrier_re[i] +
            signal_im[i, j] * carrier_im[i]
        downconverted_signal_im[i, j] = signal_im[i, j] * carrier_re[i] -
            signal_re[i, j] * carrier_im[i]
    end
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
