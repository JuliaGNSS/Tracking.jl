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

struct SampleParams
    signal_samples_to_integrate::Int
    signal_start_sample::Int
    samples_to_integrate_until_code_end::Int
    num_code_blocks_to_integrate::Int
end

function SampleParams(
    system::AbstractGNSS,
    preferred_num_code_blocks_to_integrate::Int,
    secondary_code_or_bit_found::Bool,
)
    num_code_blocks_to_integrate = calc_num_code_blocks_to_integrate(
        system,
        preferred_num_code_blocks_to_integrate,
        secondary_code_or_bit_found,
    )
    SampleParams(0, 1, 0, num_code_blocks_to_integrate)
end

function SampleParams(
    system::AbstractGNSS,
    prev_sample_params::SampleParams,
    num_samples_signal::Int,
    sampling_frequency,
    code_doppler,
    code_phase,
    preferred_num_code_blocks_to_integrate::Int,
    secondary_code_or_bit_found::Bool,
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
    signal_start_sample =
        prev_sample_params.signal_start_sample +
        prev_sample_params.signal_samples_to_integrate
    signal_samples_left = num_samples_signal - signal_start_sample + 1
    signal_samples_to_integrate =
        min(samples_to_integrate_until_code_end, signal_samples_left)
    SampleParams(
        signal_samples_to_integrate,
        signal_start_sample,
        samples_to_integrate_until_code_end,
        num_code_blocks_to_integrate,
    )
end

function init_sample_params(
    system_sats_states::TupleLike{<:NTuple{N,SystemSatsState}},
    preferred_num_code_blocks_to_integrate::Int,
) where {N}
    map(system_sats_states) do states
        map(states.states) do state
            SampleParams(
                states.system,
                preferred_num_code_blocks_to_integrate,
                found(state.sc_bit_detector),
            )
        end
    end
end

function calc_sample_params(
    system_sats_states::TupleLike{<:NTuple{N,SystemSatsState}},
    prev_system_sats_sample_params,
    num_samples_signal,
    sampling_frequency,
    preferred_num_code_blocks_to_integrate::Int,
) where {N}
    map(
        system_sats_states,
        prev_system_sats_sample_params,
    ) do states, prev_sats_sample_params
        map(states.states, prev_sats_sample_params) do state, prev_sample_params
            SampleParams(
                states.system,
                prev_sample_params,
                num_samples_signal,
                sampling_frequency,
                state.code_doppler,
                state.code_phase,
                preferred_num_code_blocks_to_integrate,
                found(state.sc_bit_detector),
            )
        end
    end
end
