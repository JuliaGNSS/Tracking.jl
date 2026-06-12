"""
$(SIGNATURES)

Returns the appropriate number of code blocks to integrate. It will be just a
single code block as long as the secondary code or bit hasn't been found. Once
found, the coherent integration is capped by one full symbol period — the
data-bit period for data-bearing signals, or the secondary-code period for
pilot signals (`data_frequency == 0`, e.g. GPS L1C-P; `set_preferred_num_code_blocks_to_integrate!`
used to be a silent no-op for these, issue #134).

The number of code blocks is then clamped to the largest divisor of that symbol
period not exceeding the preferred value: an integration length that straddled
a symbol boundary would keep the post-sync accumulator from ever hitting the
boundary, silently stalling bit/secondary emission (issue #128). Preferred
values are additionally validated when a [`TrackedSignal`](@ref) is built, so
this clamp only matters for direct callers.
"""
function calc_num_code_blocks_to_integrate(
    signal::AbstractGNSSSignal,
    preferred_num_code_blocks::Int,
    secondary_code_or_bit_found::Bool,
)
    secondary_code_or_bit_found || return 1
    data_freq = get_data_frequency(signal)
    # One full symbol = one data bit for data-bearing signals, or one
    # secondary-code period for pilots (`data_frequency == 0`).
    num_code_blocks_that_form_a_symbol = iszero(data_freq) ?
        get_secondary_code_length(signal) :
        Int(get_code_frequency(signal) / (get_code_length(signal) * data_freq))
    num_code_blocks = clamp(preferred_num_code_blocks, 1, num_code_blocks_that_form_a_symbol)
    while num_code_blocks_that_form_a_symbol % num_code_blocks != 0
        num_code_blocks -= 1
    end
    num_code_blocks
end

"""
$(SIGNATURES)

Number of primary code blocks to credit the bit-buffer accumulator with for a
just-completed integration.

Pre-sync this is always 1: the bit/secondary-code search shifts exactly one
prompt sign into its window per call, and a non-block-aligned start can make
the first integration a *fractional* block — reporting anything but 1 would
trip the single-block invariant in `_buffer_find_bit`.

Post-sync it is the number of whole blocks the integration actually covered,
recovered from its sample count. Integrations then end on a secondary-/data-bit
boundary so this divides cleanly (the `round` absorbs the sub-chip Doppler
residual in `code_frequency`). This can differ from the *intended* count
returned by [`calc_num_code_blocks_to_integrate`](@ref): the first post-sync
integration is truncated to `secondary_code_length − secondary_phase` blocks so
it lands on the data-bit boundary, even though the preferred length is the full
bit period. Crediting the accumulator with the intended length would misalign
the decoded bits (issue #125).
"""
@inline function calc_num_code_blocks_for_bit_buffer(
    signal::AbstractGNSSSignal,
    integrated_samples::Integer,
    sampling_frequency,
    secondary_code_or_bit_found::Bool,
)
    secondary_code_or_bit_found || return 1
    round(
        Int,
        integrated_samples * get_code_frequency(signal) /
        (get_code_length(signal) * sampling_frequency),
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