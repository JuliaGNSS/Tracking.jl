function gen_carrier_replica!(
    carrier_sin,
    carrier_cos,
    carrier_frequency,
    sample_frequency,
    start_phase::AbstractFloat,
    signal_start_sample::Integer,
    num_samples_left::Integer
)
    # Phase in Int32 and sin and cos in Int16
    n = get_quadrant_size_power(zero(Int16)) + 2
    fixed_point = 32 - n - 2
    delta = floor(Int32, carrier_frequency * 1 << (fixed_point + n) / sample_frequency)
    fixed_point_start_phase = floor(Int32, start_phase * 1 << (fixed_point + n))
    fixed_point_phase = fixed_point_start_phase - delta
    @inbounds for i = signal_start_sample:num_samples_left + signal_start_sample - 1
        fixed_point_phase = fixed_point_phase + delta
        carrier_sin[i], carrier_cos[i] = fpsincos(Int16(fixed_point_phase >> fixed_point))
    end
end

function gen_carrier_replica!(
    carrier_replica::StructArray{Complex{Int16}},
    carrier_frequency,
    sample_frequency,
    start_phase,
    signal_start_sample,
    num_samples_left
)
    gen_carrier_replica!(
        carrier_replica.im,
        carrier_replica.re,
        carrier_frequency,
        sample_frequency,
        start_phase,
        signal_start_sample,
        num_samples_left
    )
    carrier_replica
end

"""
$(SIGNATURES)

Updates the carrier phase.
"""
function update_carrier_phase(
    num_samples,
    carrier_frequency,
    sample_frequency,
    start_phase
)
    n = get_quadrant_size_power(zero(Int16)) + 2
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
