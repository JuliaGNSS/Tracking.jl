function gen_carrier_replica!(
    carrier_replica::StructArray{Complex{Int16}},
    carrier_frequency,
    sample_frequency,
    start_phase,
    carrier_amplitude_power::Val{N},
    start_sample,
    num_samples
) where N
    fpcarrier!(
        carrier_replica,
        carrier_frequency,
        sample_frequency,
        start_phase,
        start_sample = start_sample,
        num_samples = num_samples,
        bits = carrier_amplitude_power
    )
    carrier_replica
end

"""
$(SIGNATURES)

GPU carrier generation
"""
function gen_carrier_replica!(
    carrier_replica::CuArray{Complex{Float16}},
    carrier_frequency,
    sample_frequency,
    start_phase,
    carrier_amplitude_power::Val{N},
    start_sample,
    num_samples
) where N 
    # Here fpcarrier!() from GNSSSignals is called
    # TODO multimpe dispatch fpcarrier!() in GNSSSignals.jl
    # TODO call the fpcarrier!() described above
    carrier_replica = exp.(1im * (2pi * (1:num_samples) * carrier_frequency / sample_frequency + start_phase))
    return carrier_replica
end

"""
$(SIGNATURES)

Updates the carrier phase.
"""
function update_carrier_phase(
    num_samples,
    carrier_frequency,
    sample_frequency,
    start_phase,
    carrier_amplitude_power::Val{N}
) where N
    n = N + 2
    fixed_point = 32 - n - 2
    delta = floor(Int32, carrier_frequency * 1 << (fixed_point + n) / sample_frequency)
    fixed_point_start_phase = floor(Int32, start_phase * 1 << (fixed_point + n))
    phase_fixed_point = delta * num_samples + fixed_point_start_phase
    mod(
        phase_fixed_point / 1 << (fixed_point + n) + 0.5,
        1
    ) - 0.5
end

"""
$(SIGNATURES)

Calculates the current carrier frequency.
"""
function get_current_carrier_frequency(intermediate_frequency, carrier_doppler)
    carrier_doppler + intermediate_frequency
end
