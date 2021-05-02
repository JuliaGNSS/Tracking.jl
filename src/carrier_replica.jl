function gen_carrier_replica!(
    carrier_replica::StructArray{Complex{T}},
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples
) where T
    c_re = carrier_replica.re; c_im = carrier_replica.im
    carr_freq = T(upreferred(carrier_frequency / Hz))
    sample_freq = T(upreferred(sampling_frequency / Hz))
    twopi = T(2π)
    phase = T(start_phase)
    @avx for i in 0:num_samples - 1
        c_im[i + start_sample], c_re[i + start_sample] =
            sincos(twopi * (i * carr_freq / sample_freq + phase))
    end
    carrier_replica
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
