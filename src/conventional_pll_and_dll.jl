@kwdef struct SatConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
end

function SatConventionalPLLAndDLL(sat_state::SatState)
    SatConventionalPLLAndDLL(
        sat_state.carrier_doppler,
        sat_state.code_doppler,
        ThirdOrderBilinearLF(),
        SecondOrderBilinearLF(),
    )
end

function SatConventionalPLLAndDLL(
    sat_conventional_pll_and_dll::SatConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL{CA,CO}(
        sat_conventional_pll_and_dll.init_carrier_doppler,
        sat_conventional_pll_and_dll.init_code_doppler,
        isnothing(carrier_loop_filter) ? sat_conventional_pll_and_dll.carrier_loop_filter :
        carrier_loop_filter,
        isnothing(code_loop_filter) ? sat_conventional_pll_and_dll.code_loop_filter :
        code_loop_filter,
    )
end

@kwdef struct ConventionalPLLAndDLL{N,I,CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
              AbstractDopplerEstimator{N,I}
    states::MultipleSystemSatType{N,I,SatConventionalPLLAndDLL{CA,CO}}
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
end

function ConventionalPLLAndDLL(
    multiple_system_sats_state::MultipleSystemSatsState;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
)
    states = map(multiple_system_sats_state) do system_sats_state
        map(SatConventionalPLLAndDLL, system_sats_state.states)
    end
    ConventionalPLLAndDLL(states, carrier_loop_filter_bandwidth, code_loop_filter_bandwidth)
end

function ConventionalPLLAndDLL(
    conventional_pll_and_dll::ConventionalPLLAndDLL{N,I,CA,CO},
    states::MultipleSystemSatType{N,I,SatConventionalPLLAndDLL{CA,CO}},
) where {N,I,CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{N,I,CA,CO}(
        states,
        conventional_pll_and_dll.carrier_loop_filter_bandwidth,
        conventional_pll_and_dll.code_loop_filter_bandwidth,
    )
end

# Be careful when calling this.
# It might lead to types that are inferred at runtime?!
# Tested with 1.11.6
function ConventionalPLLAndDLL(
    conventional_pll_and_dll::ConventionalPLLAndDLL{N,I,CA,CO};
    states::Maybe{MultipleSystemSatType{N,I,SatConventionalPLLAndDLL{CA,CO}}} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {N,I,CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{N,I,CA,CO}(
        isnothing(states) ? conventional_pll_and_dll.states : states,
        isnothing(carrier_loop_filter_bandwidth) ?
        conventional_pll_and_dll.carrier_loop_filter_bandwidth :
        carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ?
        conventional_pll_and_dll.code_loop_filter_bandwidth : code_loop_filter_bandwidth,
    )
end

function merge_sats(
    pll_and_dll::ConventionalPLLAndDLL{N,I,CA,CO},
    system_idx,
    sats_state::Dictionary{I,<:SatState},
) where {N,I,CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    new_sats_state = map(sats_state) do sat_state
        SatConventionalPLLAndDLL{CA,CO}(
            sat_state.carrier_doppler,
            sat_state.code_doppler,
            constructorof(CA)(),
            constructorof(CO)(),
        )
    end
    @set pll_and_dll.states[system_idx] =
        merge(pll_and_dll.states[system_idx], new_sats_state)
end

function filter_out_sats(
    pll_and_dll::ConventionalPLLAndDLL{N,I,CA,CO},
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {N,I,CA,CO}
    filtered_pll_and_dlls = map(
        last,
        filter(((id,),) -> !in(id, identifiers), pairs(pll_and_dll.states[system_idx])),
    )
    @set pll_and_dll.states[system_idx] = filtered_pll_and_dlls
end

"""
$(SIGNATURES)

Aid dopplers. That is velocity aiding for the carrier doppler and carrier aiding
for the code doppler.
"""
function aid_dopplers(
    system::AbstractGNSS,
    init_carrier_doppler,
    init_code_doppler,
    carrier_freq_update,
    code_freq_update,
)
    carrier_doppler = carrier_freq_update
    code_doppler =
        code_freq_update + carrier_doppler * get_code_center_frequency_ratio(system)
    init_carrier_doppler + carrier_doppler, init_code_doppler + code_doppler
end

"""
$(SIGNATURES)

Estimate Dopplers and filter prompts for all satellites where the correlation has reached
the end of the code or multiples of that. This function uses the
conventional PLL and DLL implementation to estimate Dopplers for
carrier and code. Those Doppler estimations will be used to create the next
replicas to downconvert and decode the incoming signal. In addition to the
Doppler estimation it will also filter the prompt with the configured 
post correlation filter.
In the case that the that the correlation hasn't reached the end, e.g. in the case
the incoming signal did not provide enough samples, it will return struct with
zeroed values.
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:MultipleSystemSatsState,<:ConventionalPLLAndDLL},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    new_multiple_system_states = map(
        track_state.multiple_system_sats_state,
        track_state.doppler_estimator.states,
    ) do system_sats_state, pll_and_dll_states
        new_states = map(
            system_sats_state.states,
            pll_and_dll_states,
        ) do sat_state, pll_and_dll_state
            if !is_zero(sat_state.correlator) ||
               is_zero(sat_state.last_fully_integrated_correlator) ||
               sat_state.integrated_samples == 0
                return (estimator = pll_and_dll_state, state = sat_state)
            end
            integrated_code_blocks = calc_num_code_blocks_to_integrate(
                system_sats_state.system,
                preferred_num_code_blocks_to_integrate,
                has_bit_or_secondary_code_been_found(sat_state.bit_buffer),
            )

            normalized_correlator = normalize(
                sat_state.last_fully_integrated_correlator,
                sat_state.integrated_samples,
            )
            post_corr_filter = update(
                sat_state.post_corr_filter,
                get_prompt(normalized_correlator),
            )
            filtered_correlator = apply(post_corr_filter, normalized_correlator)
            prompt = get_prompt(filtered_correlator)
            cn0_estimator = update(get_cn0_estimator(sat_state), prompt)
            bit_buffer = buffer(
                system_sats_state.system,
                sat_state.bit_buffer,
                integrated_code_blocks,
                prompt,
            )
            integration_time = sat_state.integrated_samples / sampling_frequency
            pll_discriminator = pll_disc(system_sats_state.system, filtered_correlator)
            dll_discriminator = dll_disc(
                system_sats_state.system,
                filtered_correlator,
                sat_state.code_doppler,
                sampling_frequency,
            )
            carrier_freq_update, carrier_loop_filter = filter_loop(
                pll_and_dll_state.carrier_loop_filter,
                pll_discriminator,
                integration_time,
                track_state.doppler_estimator.carrier_loop_filter_bandwidth,
            )
            code_freq_update, code_loop_filter = filter_loop(
                pll_and_dll_state.code_loop_filter,
                dll_discriminator,
                integration_time,
                track_state.doppler_estimator.code_loop_filter_bandwidth,
            )
            carrier_doppler, code_doppler = aid_dopplers(
                system_sats_state.system,
                pll_and_dll_state.init_carrier_doppler,
                pll_and_dll_state.init_code_doppler,
                carrier_freq_update,
                code_freq_update,
            )
            doppler_estimator = SatConventionalPLLAndDLL(
                pll_and_dll_state;
                carrier_loop_filter,
                code_loop_filter,
            )
            return (
                estimator = doppler_estimator,
                state = SatState(
                    sat_state;
                    carrier_doppler,
                    code_doppler,
                    integrated_samples = 0,
                    last_fully_integrated_filtered_prompt = prompt,
                    bit_buffer,
                    cn0_estimator,
                    post_corr_filter,
                ),
            )
        end

        return (
            estimators = map(x -> x.estimator, new_states),
            states = SystemSatsState(system_sats_state, map(x -> x.state, new_states)),
        )
    end
    new_doppler_estimators = map(x -> x.estimators, new_multiple_system_states)
    new_doppler_estimator =
        ConventionalPLLAndDLL(track_state.doppler_estimator, new_doppler_estimators)
    new_multiple_system_sats_state = map(x -> x.states, new_multiple_system_states)
    return TrackState(track_state, new_multiple_system_sats_state, new_doppler_estimator)
end