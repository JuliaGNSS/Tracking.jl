"""
$(SIGNATURES)

Fixed point CPU StructArray carrier replica generation
"""
function gen_carrier_replica!(
    carrier_replica::StructArray{Complex{T}},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples
) where T
    c_re = carrier_replica.re; c_im = carrier_replica.im
    carrier_freq = upreferred(carrier_frequency / Hz)
    sampling_freq = upreferred(sampling_frequency / Hz)
    @avx for i in 0:num_samples - 1
        c_im[i + start_sample], c_re[i + start_sample] =
            sincos(T(2π) * (i * T(carrier_freq) / T(sampling_freq) + T(start_phase)))
    end
    carrier_replica
end

"""
$(SIGNATURES)

Floating point CPU StructArray carrier generation
"""
function gen_carrier_replica!(
    carrier_replica::StructArray{Complex{T},1,NamedTuple{(:re, :im),Tuple{Array{T,1},Array{T,1}}},Int64},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    carrier_amplitude_power::Val{N},
    start_sample,
    num_samples
) where {
    T <: AbstractFloat,
    N
}
    sample_range = start_sample:num_samples + start_sample - 1
    @views @. carrier_replica.re[sample_range] = 2pi * (sample_range) * carrier_frequency / sampling_frequency + start_phase
    @. carrier_replica.im[sample_range] = sin(carrier_replica.re[sample_range])
    @. carrier_replica.re[sample_range] = cos(carrier_replica.re[sample_range])
    return carrier_replica
end

"""
$(SIGNATURES)

Updates the carrier phase.
"""
function update_carrier_phase(
    num_samples,
    carrier_frequency,
    sampling_frequency,
    start_phase,
)
    phase = num_samples * carrier_frequency / sampling_frequency + start_phase
    mod(phase + 0.5, 1) - 0.5
end

"""
$(SIGNATURES)

Calculates the current carrier frequency.
"""
function get_current_carrier_frequency(intermediate_frequency, carrier_doppler)
    carrier_doppler + intermediate_frequency
end
