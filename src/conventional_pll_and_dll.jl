abstract type AbstractDopplerEstimator end

struct ConventionalPLLAndDLL{
    CA<:AbstractLoopFilter,
    CO<:AbstractLoopFilter,
    F<:AbstractPostCorrFilter,
}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA
    code_loop_filter::CO
    carrier_loop_filter_bandwidth::typeof(1.0Hz)
    code_loop_filter_bandwidth::typeof(1.0Hz)
    post_corr_filter::F
end

function ConventionalPLLAndDLL(
    sat_state::SatState,
    carrier_loop_filter = ThirdOrderBilinearLF(),
    code_loop_filter = SecondOrderBilinearLF(),
    carrier_loop_filter_bandwidth = 18Hz,
    code_loop_filter_bandwidth = 1Hz,
    post_corr_filter = DefaultPostCorrFilter(),
)
    ConventionalPLLAndDLL(
        sat_state.carrier_doppler,
        sat_state.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        float(carrier_loop_filter_bandwidth),
        float(code_loop_filter_bandwidth),
        post_corr_filter,
    )
end

struct ConventionalPLLsAndDLLs{
    N,
    PD<:TupleLike{<:NTuple{N,Vector{<:ConventionalPLLAndDLL}}},
} <: AbstractDopplerEstimator
    plls_and_dlls::PD
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
    doppler_estimator::ConventionalPLLsAndDLLs,
    system_sats_state::TupleLike{<:NTuple{N,SystemSatsState}},
    sat_sample_params::TupleLike{<:NTuple{N,Vector{SampleParams}}},
    sampling_frequency,
    min_integration_time,
) where {N}
    results = map(
        doppler_estimator.plls_and_dlls,
        system_sats_state,
        sat_sample_params,
    ) do estimator_per_system, system_sats, sat_sample_params_per_system
        map(
            estimator_per_system,
            system_sats.states,
            sat_sample_params_per_system,
        ) do estimator, state, sample_params
            integration_time = state.integrated_samples / sampling_frequency
            return if sample_params.signal_samples_to_integrate ==
                      sample_params.samples_to_integrate_until_code_end &&
                      integration_time >= min_integration_time
                correlator = normalize(state.correlator, state.integrated_samples)
                post_corr_filter =
                    update(estimator.post_corr_filter, get_prompt(correlator))
                filtered_correlator = apply(post_corr_filter, correlator)
                pll_discriminator = pll_disc(system_sats.system, filtered_correlator)
                dll_discriminator = dll_disc(
                    system_sats.system,
                    filtered_correlator,
                    (state.code_doppler + get_code_frequency(system_sats.system)) /
                    sampling_frequency,
                )
                carrier_freq_update, carrier_loop_filter = filter_loop(
                    estimator.carrier_loop_filter,
                    pll_discriminator,
                    integration_time,
                    estimator.carrier_loop_filter_bandwidth,
                )
                code_freq_update, code_loop_filter = filter_loop(
                    estimator.code_loop_filter,
                    dll_discriminator,
                    integration_time,
                    estimator.code_loop_filter_bandwidth,
                )
                carrier_doppler, code_doppler = aid_dopplers(
                    system_sats.system,
                    estimator.init_carrier_doppler,
                    estimator.init_code_doppler,
                    carrier_freq_update,
                    code_freq_update,
                )
                dopplers_and_filtered_prompt = DopplersAndFilteredPrompt(
                    carrier_doppler,
                    code_doppler,
                    get_prompt(filtered_correlator),
                )
                new_estimator = ConventionalPLLAndDLL(
                    estimator.init_carrier_doppler,
                    estimator.init_code_doppler,
                    carrier_loop_filter,
                    code_loop_filter,
                    estimator.carrier_loop_filter_bandwidth,
                    estimator.code_loop_filter_bandwidth,
                    post_corr_filter,
                )
                (
                    estimator = new_estimator,
                    dopplers_and_filtered_prompt = dopplers_and_filtered_prompt,
                )
            else
                (
                    estimator = estimator,
                    dopplers_and_filtered_prompt = DopplersAndFilteredPrompt(
                        0.0Hz,
                        0.0Hz,
                        complex(0.0, 0.0),
                    ),
                )
            end
        end
    end
    new_doppler_estimator =
        ConventionalPLLsAndDLLs(map(x -> map(y -> y.estimator, x), results))
    new_dopplers_and_filtered_prompts =
        map(x -> map(y -> y.dopplers_and_filtered_prompt, x), results)
    new_doppler_estimator, new_dopplers_and_filtered_prompts
end
