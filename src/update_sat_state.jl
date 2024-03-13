function update(
    system_sats_states::TupleLike{<:NTuple{N,SystemSatsState}},
    sampling_frequency,
    intermediate_frequency,
    new_correlators::TupleLike{<:NTuple{N,Vector{<:AbstractCorrelator}}},
    params::TupleLike{<:NTuple{N,Vector{SampleParams}}},
    num_samples_signal,
) where {N}
    map(
        new_correlators,
        params,
        system_sats_states,
    ) do system_correlators, system_params, system_sats
        sats_state = map(
            system_correlators,
            system_params,
            system_sats.states,
        ) do new_correlator, sat_params, sat_state
            carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
            code_frequency =
                sat_state.code_doppler + get_code_frequency(system_sats.system)
            carrier_phase = update_carrier_phase(
                sat_params.signal_samples_to_integrate,
                carrier_frequency,
                sampling_frequency,
                sat_state.carrier_phase,
            )
            code_phase = update_code_phase(
                system_sats.system,
                sat_params.signal_samples_to_integrate,
                code_frequency,
                sampling_frequency,
                sat_state.code_phase,
                found(sat_state.sc_bit_detector),
            )
            sample_of_last_fully_integrated_correlator =
                sat_state.sample_of_last_fully_integrated_correlator -
                (sat_params.signal_start_sample == 1 ? num_samples_signal : 0)
            SatState(
                sat_state.prn,
                code_phase,
                sat_state.code_doppler,
                carrier_phase,
                sat_state.carrier_doppler,
                sat_state.integrated_samples + sat_params.signal_samples_to_integrate,
                new_correlator,
                sat_state.last_fully_integrated_correlator,
                sat_state.last_fully_integrated_filtered_prompt,
                sample_of_last_fully_integrated_correlator,
                sat_state.sc_bit_detector,
                sat_state.prompts_buffer,
                sat_state.bit_buffer,
            )
        end
        SystemSatsState(system_sats.system, sats_state)
    end
end

function update(
    system_sats_state::TupleLike{<:NTuple{N,SystemSatsState}},
    dopplers_and_filtered_prompts::TupleLike{
        <:NTuple{N,Vector{<:Union{Nothing,DopplersAndFilteredPrompt}}},
    },
    sat_sample_params::TupleLike{<:NTuple{N,Vector{SampleParams}}},
) where {N}
    map(
        system_sats_state,
        dopplers_and_filtered_prompts,
        sat_sample_params,
    ) do system_sats, dopplers_and_filtered_prompts_per_system, sat_sample_params_per_system
        sats_state = map(
            system_sats.states,
            dopplers_and_filtered_prompts_per_system,
            sat_sample_params_per_system,
        ) do state, dopplers_and_filtered_prompt, sample_params
            return if dopplers_and_filtered_prompt.filtered_prompt != complex(0.0, 0.0)
                filtered_prompt = dopplers_and_filtered_prompt.filtered_prompt
                carrier_doppler = dopplers_and_filtered_prompt.carrier_doppler
                code_doppler = dopplers_and_filtered_prompt.code_doppler
                prompts_buffer = update(state.prompts_buffer, filtered_prompt)
                bit_buffer = buffer(
                    system_sats.system,
                    state.bit_buffer,
                    sample_params.num_code_blocks_to_integrate,
                    found(state.sc_bit_detector),
                    filtered_prompt,
                )
                sc_bit_detector =
                    find(system_sats.system, state.sc_bit_detector, filtered_prompt)
                correlator = zero(state.correlator)
                integrated_samples = 0

                new_state = SatState(
                    state.prn,
                    state.code_phase,
                    code_doppler,
                    state.carrier_phase,
                    carrier_doppler,
                    integrated_samples,
                    correlator,
                    state.correlator,
                    dopplers_and_filtered_prompt.filtered_prompt,
                    sample_params.signal_start_sample +
                    sample_params.signal_samples_to_integrate - 1,
                    sc_bit_detector,
                    prompts_buffer,
                    bit_buffer,
                )
                new_state
            else
                state
            end
        end
        SystemSatsState(system_sats.system, sats_state)
    end
end
