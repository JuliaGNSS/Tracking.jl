struct SystemConventionalPLLAndDLL <: AbstractSystemDopplerEstimator end
struct TrackingConventionalPLLAndDLL <: AbstractTrackingDopplerEstimator end

struct ConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
       AbstractSatDopplerEstimator
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA
    code_loop_filter::CO
    carrier_loop_filter_bandwidth::typeof(1.0Hz)
    code_loop_filter_bandwidth::typeof(1.0Hz)
end

function ConventionalPLLAndDLL(
    carrier_doppler,
    code_doppler;
    carrier_loop_filter = ThirdOrderBilinearLF(),
    code_loop_filter = SecondOrderBilinearLF(),
    carrier_loop_filter_bandwidth = 18Hz,
    code_loop_filter_bandwidth = 1Hz,
)
    ConventionalPLLAndDLL(
        float(carrier_doppler),
        float(code_doppler),
        carrier_loop_filter,
        code_loop_filter,
        float(carrier_loop_filter_bandwidth),
        float(code_loop_filter_bandwidth),
    )
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
    track_state::TrackState{
        <:MultipleSystemSatsState{
            N,
            <:AbstractGNSS,
            <:SatState{
                <:AbstractCorrelator,
                <:AbstractPostCorrFilter,
                <:ConventionalPLLAndDLL,
            },
            <:SystemConventionalPLLAndDLL,
        },
        <:TrackingConventionalPLLAndDLL,
    },
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
) where {N}
    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            new_sat_states = map(system_sats_state.states) do sat_state
                if !is_zero(sat_state.correlator) ||
                   is_zero(sat_state.last_fully_integrated_correlator) ||
                   sat_state.integrated_samples == 0
                    return sat_state
                end
                integrated_code_blocks = calc_num_code_blocks_to_integrate(
                    system_sats_state.system,
                    preferred_num_code_blocks_to_integrate,
                    found(sat_state.sc_bit_detector),
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
                    found(sat_state.sc_bit_detector), # TODO: should the detector be updated before this?
                    prompt,
                )
                sc_bit_detector =
                    find(system_sats_state.system, sat_state.sc_bit_detector, prompt)
                integration_time = sat_state.integrated_samples / sampling_frequency
                pll_discriminator =
                    pll_disc(system_sats_state.system, filtered_correlator)
                dll_discriminator = dll_disc(
                    system_sats_state.system,
                    filtered_correlator,
                    (
                        sat_state.code_doppler +
                        get_code_frequency(system_sats_state.system)
                    ) / sampling_frequency,
                )
                carrier_freq_update, carrier_loop_filter = filter_loop(
                    sat_state.doppler_estimator.carrier_loop_filter,
                    pll_discriminator,
                    integration_time,
                    sat_state.doppler_estimator.carrier_loop_filter_bandwidth,
                )
                code_freq_update, code_loop_filter = filter_loop(
                    sat_state.doppler_estimator.code_loop_filter,
                    dll_discriminator,
                    integration_time,
                    sat_state.doppler_estimator.code_loop_filter_bandwidth,
                )
                carrier_doppler, code_doppler = aid_dopplers(
                    system_sats_state.system,
                    sat_state.doppler_estimator.init_carrier_doppler,
                    sat_state.doppler_estimator.init_code_doppler,
                    carrier_freq_update,
                    code_freq_update,
                )
                new_doppler_estimator = ConventionalPLLAndDLL(
                    sat_state.doppler_estimator.init_carrier_doppler,
                    sat_state.doppler_estimator.init_code_doppler,
                    carrier_loop_filter,
                    code_loop_filter,
                    sat_state.doppler_estimator.carrier_loop_filter_bandwidth,
                    sat_state.doppler_estimator.code_loop_filter_bandwidth,
                )
                return SatState(
                    sat_state;
                    carrier_doppler,
                    code_doppler,
                    integrated_samples = 0,
                    doppler_estimator = new_doppler_estimator,
                    last_fully_integrated_filtered_prompt = prompt,
                    bit_buffer,
                    sc_bit_detector,
                    cn0_estimator,
                    post_corr_filter,
                )
            end
            return SystemSatsState(system_sats_state, new_sat_states)
        end
    return TrackState(track_state, new_multiple_system_sats_state)
end