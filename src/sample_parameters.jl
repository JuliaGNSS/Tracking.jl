"""
$(SIGNATURES)

Returns the appropiate number of code blocks to integrate. It will be just a single
code block when the secondary code or bit isn't found. The maximum number of code
blocks to integrate is limited by the bit shift.

For pilot signals with `data_frequency == 0` (e.g. GPS L1C-P), longer coherent
integration would be bounded by the secondary-code period instead, but secondary-
code sync isn't wired up yet — clamp to one code block until it lands.
"""
function calc_num_code_blocks_to_integrate(
    signal::AbstractGNSSSignal,
    preferred_num_code_blocks::Int,
    secondary_code_or_bit_found::Bool,
)
    secondary_code_or_bit_found || return 1
    data_freq = get_data_frequency(signal)
    iszero(data_freq) && return 1
    min(
        Int(get_code_frequency(signal) / (get_code_length(signal) * data_freq)),
        preferred_num_code_blocks,
    )
end

"""
$(SIGNATURES)

Calculates the number of chips to integrate.
"""
function calc_num_chips_to_integrate(
    signal::AbstractGNSSSignal,
    num_code_blocks::Int,
    current_code_phase,
)
    max_phase = num_code_blocks * get_code_length(signal)
    current_phase_mod_max_phase = mod(current_code_phase, max_phase)
    return max_phase - current_phase_mod_max_phase
end

"""
$(SIGNATURES)

Calculates the number of samples to integrate.
"""
function calc_num_samples_left_to_integrate(
    signal::AbstractGNSSSignal,
    num_code_blocks::Int,
    sampling_frequency,
    current_code_doppler,
    current_code_phase,
)
    chips_to_integrate =
        calc_num_chips_to_integrate(signal, num_code_blocks, current_code_phase)
    code_frequency = get_code_frequency(signal) + current_code_doppler
    ceil(Int, chips_to_integrate * sampling_frequency / code_frequency)
end

"""
$(SIGNATURES)

Calculate the number of signal samples to integrate and whether the integration
will complete a full code period. Returns a tuple of (samples_to_integrate,
is_integration_completed). Takes into account available signal samples and
code block boundaries.
"""
function calc_signal_samples_to_integrate(
    signal::AbstractGNSSSignal,
    signal_start_sample::Int,
    sampling_frequency,
    code_doppler,
    code_phase,
    preferred_num_code_blocks_to_integrate::Int,
    secondary_code_or_bit_found::Bool,
    num_samples_signal::Int,
)
    num_code_blocks_to_integrate = calc_num_code_blocks_to_integrate(
        signal,
        preferred_num_code_blocks_to_integrate,
        secondary_code_or_bit_found,
    )
    samples_to_integrate_until_code_end = calc_num_samples_left_to_integrate(
        signal,
        num_code_blocks_to_integrate,
        sampling_frequency,
        code_doppler,
        code_phase,
    )
    signal_samples_left = num_samples_signal - signal_start_sample + 1
    samples_to_integrate = min(samples_to_integrate_until_code_end, signal_samples_left)
    return samples_to_integrate, samples_to_integrate == samples_to_integrate_until_code_end
end