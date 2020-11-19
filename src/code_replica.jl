function gen_code_replica!(
    code_replica,
    ::Type{S},
    code_frequency,
    sampling_frequency,
    start_code_phase::AbstractFloat,
    start_sample::Integer,
    num_samples::Integer,
    correlator_sample_shifts::SVector,
    prn::Integer
) where S <: AbstractGNSSSystem
    early_sample_shift = correlator_sample_shifts[end]
    late_sample_shift  = correlator_sample_shifts[1]
    early_late_samples = early_sample_shift - late_sample_shift
    # code_replica[start_sample:num_samples + early_late_samples + start_sample - 1] .= GNSSSignals.get_code.(S, (0 + late_sample_shift:num_samples - 1 + early_sample_shift) .* code_frequency ./ sampling_frequency .+ start_code_phase, prn)
    # code_replica
    fixed_point = sizeof(Int) * 8 - 1 - min_bits_for_code_length(S)
    delta = floor(Int, code_frequency * 1 << fixed_point / sampling_frequency)
    modded_start_code_phase = mod(
        start_code_phase,
        get_code_length(S) * get_secondary_code_length(S)
    )
    fixed_point_start_code_phase = floor(Int, modded_start_code_phase * 1 << fixed_point)
    @inbounds for i = start_sample:num_samples + early_late_samples + start_sample - 1
        fixed_point_code_phase = (i + late_sample_shift - start_sample) * delta +
            fixed_point_start_code_phase
        code_index = fixed_point_code_phase >> fixed_point
        code_replica[i] = get_code_unsafe(S, code_index, prn)
    end
    code_replica
end

"""
$(SIGNATURES)

Updates the code phase.
"""
function update_code_phase(
    ::Type{S},
    num_samples,
    code_frequency,
    sampling_frequency,
    start_code_phase,
    secondary_code_or_bit_found
) where S <: AbstractGNSSSystem
    if get_data_frequency(S) == 0Hz
        secondary_code_or_bit_length = get_secondary_code_length(S)
    else
        secondary_code_or_bit_length =
            Int(get_code_frequency(S) / (get_data_frequency(S) * get_code_length(S)))
    end
    code_length = get_code_length(S) *
        (secondary_code_or_bit_found ? secondary_code_or_bit_length : 1)
    mod(code_frequency * num_samples / sampling_frequency + start_code_phase, code_length)
#    fixed_point = sizeof(Int) * 8 - 1 - min_bits_for_code_length(S)
#    delta = floor(Int, code_frequency * 1 << fixed_point / sampling_frequency)
#    fixed_point_start_phase = floor(Int, start_code_phase * 1 << fixed_point)
#    phase_fixed_point = delta * num_samples + fixed_point_start_phase
#    mod(phase_fixed_point / 1 << fixed_point, code_length)
end

"""
$(SIGNATURES)

Calculates the current code frequency.
"""
function get_current_code_frequency(::Type{S}, code_doppler) where S <: AbstractGNSSSystem
    code_doppler + get_code_frequency(S)
end
