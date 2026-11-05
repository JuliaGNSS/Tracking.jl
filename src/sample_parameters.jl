"""
$(SIGNATURES)

The number of chips left from `current_code_phase` to the end of an integration of
`num_code_periods` code periods.
"""
function calc_num_chips_left_to_integrate(
    signal::AbstractGNSSSignal,
    num_code_periods::Int,
    current_code_phase,
)
    max_phase = num_code_periods * get_code_length(signal)
    current_phase_mod_max_phase = mod(current_code_phase, max_phase)
    return max_phase - current_phase_mod_max_phase
end

"""
$(SIGNATURES)

The number of samples left from `current_code_phase` to the end of an integration of
`num_code_periods` code periods at the code Doppler `current_code_doppler`, rounded up.
"""
function calc_num_samples_left_to_integrate(
    signal::AbstractGNSSSignal,
    num_code_periods::Int,
    sampling_frequency,
    current_code_doppler,
    current_code_phase,
)
    num_chips_left_to_integrate =
        calc_num_chips_left_to_integrate(signal, num_code_periods, current_code_phase)
    code_frequency = get_chip_rate(signal) + current_code_doppler
    ceil(Int, num_chips_left_to_integrate * sampling_frequency / code_frequency)
end
