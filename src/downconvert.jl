"""
$(SIGNATURES)

Downconvert a signal that has only a single dimension, e.g. like
a single antenna channel.
This is for the case that the signal is of type StructArray.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},1},
    signal::StructArray{Complex{ST},1},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    s_re = signal.re
    s_im = signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1
        ds_re[i] = s_re[i] * c_re[i] + s_im[i] * c_im[i]
        ds_im[i] = s_im[i] * c_re[i] - s_re[i] * c_im[i]
    end
    downconverted_signal
end

"""
$(SIGNATURES)

Downconvert a signal that has multiple dimensions, e.g. like
multiple antenna channels.
This is for the case that the signal is of type StructArray.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},2},
    signal::StructArray{Complex{ST},2},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    s_re = signal.re
    s_im = signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1, j = 1:size(s_re, 2)
        ds_re[i, j] = s_re[i, j] * c_re[i] + s_im[i, j] * c_im[i]
        ds_im[i, j] = s_im[i, j] * c_re[i] - s_re[i, j] * c_im[i]
    end
    downconverted_signal
end

"""
$(SIGNATURES)

Downconvert a signal that has only a single dimension, e.g. like
a single antenna channel.
This is for the case that the signal is a regular complex array.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},1},
    signal::AbstractArray{Complex{ST},1},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    signal_real = reinterpret(reshape, ST, signal)
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1
        ds_re[i] = signal_real[1, i] * c_re[i] + signal_real[2, i] * c_im[i]
        ds_im[i] = signal_real[2, i] * c_re[i] - signal_real[1, i] * c_im[i]
    end
    downconverted_signal
end

"""
$(SIGNATURES)

Downconvert a signal that has multiple dimensions, e.g. like
multiple antenna channels.
This is for the case that the signal is a regular complex array.
"""
function downconvert!(
    downconverted_signal::StructArray{Complex{T},2},
    signal::AbstractArray{Complex{ST},2},
    carrier_replica::StructArray{Complex{T},1},
    start_sample::Integer,
    num_samples::Integer,
) where {T,ST}
    signal_real = reinterpret(reshape, ST, signal)
    ds_re = downconverted_signal.re
    ds_im = downconverted_signal.im
    c_re = carrier_replica.re
    c_im = carrier_replica.im
    @avx for i = start_sample:num_samples+start_sample-1, j = 1:size(ds_re, 2)
        ds_re[i, j] = signal_real[1, i, j] * c_re[i] + signal_real[2, i, j] * c_im[i]
        ds_im[i, j] = signal_real[2, i, j] * c_re[i] - signal_real[1, i, j] * c_im[i]
    end
    downconverted_signal
end

# Fused carrier generation and downconversion: generates the carrier
# on-the-fly during downconversion, eliminating the carrier replica buffer.
# Uses @generated to unroll the antenna dimension at compile time, which allows
# @avx to vectorize purely along the sample dimension while computing sincos
# once per sample.

# StructArray signal, any number of antennas
function downconvert!(
    downconverted_signal::StructArray{Complex{T}},
    signal::StructArray{Complex{TS}},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples,
    ::NumAnts{NANT},
) where {T,TS,NANT}
    _fused_downconvert_unrolled!(
        downconverted_signal.re,
        downconverted_signal.im,
        signal.re,
        signal.im,
        T(upreferred(carrier_frequency / Hz)),
        T(upreferred(sampling_frequency / Hz)),
        T(start_phase),
        T(2π),
        start_sample,
        num_samples,
        NumAnts{NANT}(),
    )
    downconverted_signal
end

# Regular complex array signal, any number of antennas
function downconvert!(
    downconverted_signal::StructArray{Complex{T}},
    signal::AbstractArray{Complex{ST}},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples,
    ::NumAnts{NANT},
) where {T,ST,NANT}
    signal_real = reinterpret(reshape, ST, signal)
    _fused_downconvert_reinterp_unrolled!(
        downconverted_signal.re,
        downconverted_signal.im,
        signal_real,
        T(upreferred(carrier_frequency / Hz)),
        T(upreferred(sampling_frequency / Hz)),
        T(start_phase),
        T(2π),
        start_sample,
        num_samples,
        NumAnts{NANT}(),
    )
    downconverted_signal
end

@generated function _fused_downconvert_unrolled!(
    ds_re, ds_im, s_re, s_im, carrier_freq, sampling_freq, start_phase,
    two_pi, start_sample, num_samples, ::NumAnts{NANT},
) where {NANT}
    body_lines = Expr[]
    for j in 1:NANT
        push!(body_lines, :(ds_re[i, $j] = s_re[i, $j] * c_re + s_im[i, $j] * c_im))
        push!(body_lines, :(ds_im[i, $j] = s_im[i, $j] * c_re - s_re[i, $j] * c_im))
    end
    quote
        @avx for i = start_sample:start_sample+num_samples-1
            c_im, c_re = sincos(
                two_pi *
                ((i - start_sample) * carrier_freq / sampling_freq + start_phase),
            )
            $(body_lines...)
        end
    end
end

@generated function _fused_downconvert_reinterp_unrolled!(
    ds_re, ds_im, signal_real, carrier_freq, sampling_freq, start_phase,
    two_pi, start_sample, num_samples, ::NumAnts{NANT},
) where {NANT}
    body_lines = Expr[]
    if NANT == 1
        push!(body_lines, :(ds_re[i] = signal_real[1, i] * c_re + signal_real[2, i] * c_im))
        push!(body_lines, :(ds_im[i] = signal_real[2, i] * c_re - signal_real[1, i] * c_im))
    else
        for j in 1:NANT
            push!(body_lines, :(ds_re[i, $j] = signal_real[1, i, $j] * c_re + signal_real[2, i, $j] * c_im))
            push!(body_lines, :(ds_im[i, $j] = signal_real[2, i, $j] * c_re - signal_real[1, i, $j] * c_im))
        end
    end
    quote
        @avx for i = start_sample:start_sample+num_samples-1
            c_im, c_re = sincos(
                two_pi *
                ((i - start_sample) * carrier_freq / sampling_freq + start_phase),
            )
            $(body_lines...)
        end
    end
end
