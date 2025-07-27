"""
$(SIGNATURES)

Returns the appropiate number of code blocks to integrate. It will be just a single
code block when the secondary code or bit isn't found. The maximum number of code
blocks to integrate is limited by the bit shift.
"""
function calc_num_code_blocks_to_integrate(
    system::AbstractGNSS,
    preferred_num_code_blocks::Int,
    secondary_code_or_bit_found::Bool,
)
    ifelse(
        secondary_code_or_bit_found,
        min(
            Int(
                get_code_frequency(system) /
                (get_code_length(system) * get_data_frequency(system)),
            ),
            preferred_num_code_blocks,
        ),
        1,
    )
end

"""
$(SIGNATURES)

Calculates the number of chips to integrate.
"""
function calc_num_chips_to_integrate(
    system::AbstractGNSS,
    num_code_blocks::Int,
    current_code_phase,
)
    max_phase = num_code_blocks * get_code_length(system)
    current_phase_mod_max_phase = mod(current_code_phase, max_phase)
    return max_phase - current_phase_mod_max_phase
end

"""
$(SIGNATURES)

Calculates the number of samples to integrate.
"""
function calc_num_samples_left_to_integrate(
    system::AbstractGNSS,
    num_code_blocks::Int,
    sampling_frequency,
    current_code_doppler,
    current_code_phase,
)
    chips_to_integrate =
        calc_num_chips_to_integrate(system, num_code_blocks, current_code_phase)
    code_frequency = get_code_frequency(system) + current_code_doppler
    ceil(Int, chips_to_integrate * sampling_frequency / code_frequency)
end

function calc_signal_samples_to_integrate(
    system::AbstractGNSS,
    signal_start_sample::Int,
    sampling_frequency,
    code_doppler,
    code_phase,
    preferred_num_code_blocks_to_integrate::Int,
    secondary_code_or_bit_found::Bool,
    num_samples_signal::Int,
)
    num_code_blocks_to_integrate = calc_num_code_blocks_to_integrate(
        system,
        preferred_num_code_blocks_to_integrate,
        secondary_code_or_bit_found,
    )
    samples_to_integrate_until_code_end = calc_num_samples_left_to_integrate(
        system,
        num_code_blocks_to_integrate,
        sampling_frequency,
        code_doppler,
        code_phase,
    )
    signal_samples_left = num_samples_signal - signal_start_sample + 1
    samples_to_integrate = min(samples_to_integrate_until_code_end, signal_samples_left)
    return samples_to_integrate, samples_to_integrate == samples_to_integrate_until_code_end
end