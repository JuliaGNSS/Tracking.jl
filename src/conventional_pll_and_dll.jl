"""
Per-satellite state for the conventional PLL and DLL Doppler estimator.
Holds initial Doppler values and loop filter states.
"""
@kwdef struct SatConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
end

function SatConventionalPLLAndDLL(
    sat_state::SatState,
    carrier_loop_filter::CA,
    code_loop_filter::CO;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL(
        sat_state.carrier_doppler,
        sat_state.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

function SatConventionalPLLAndDLL(
    sat_conventional_pll_and_dll::SatConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL{CA,CO}(
        sat_conventional_pll_and_dll.init_carrier_doppler,
        sat_conventional_pll_and_dll.init_code_doppler,
        isnothing(carrier_loop_filter) ? sat_conventional_pll_and_dll.carrier_loop_filter :
        carrier_loop_filter,
        isnothing(code_loop_filter) ? sat_conventional_pll_and_dll.code_loop_filter :
        code_loop_filter,
        isnothing(carrier_loop_filter_bandwidth) ?
        sat_conventional_pll_and_dll.carrier_loop_filter_bandwidth :
        carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ?
        sat_conventional_pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Conventional Phase-Locked Loop (PLL) and Delay-Locked Loop (DLL) Doppler
estimator. Configuration-only — per-satellite state lives in each
[`TrackedSat`](@ref) wrapper, produced via [`init_estimator_state`](@ref).

Type parameters `CA` and `CO` select the carrier and code loop filter types;
the bandwidth fields configure the default loop bandwidths used when seeding
new satellites.
"""
struct ConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
       AbstractDopplerEstimator
    carrier_loop_filter_bandwidth::typeof(1.0Hz)
    code_loop_filter_bandwidth::typeof(1.0Hz)
end

function ConventionalPLLAndDLL(
    ::Type{CA} = ThirdOrderBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{CA,CO}(
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Create a ConventionalPLLAndDLL with FLL-assisted carrier tracking. This is the
default Doppler estimator used by TrackState. Uses a ThirdOrderAssistedBilinearLF
for the carrier loop filter which combines PLL and FLL discriminators for
improved tracking under high dynamics.
"""
function ConventionalAssistedPLLAndDLL(
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL(
        ThirdOrderAssistedBilinearLF,
        CO;
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Kwarg-update constructor for tweaking bandwidths in place.
function ConventionalPLLAndDLL(
    pll_and_dll::ConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ?
        pll_and_dll.code_loop_filter_bandwidth : code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Build the per-satellite estimator state stored in a [`TrackedSat`](@ref) for a
satellite tracked under [`ConventionalPLLAndDLL`](@ref).
"""
function init_estimator_state(
    estimator::ConventionalPLLAndDLL{CA,CO},
    sat_state::SatState,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    carrier_loop_filter = constructorof(CA)()
    code_loop_filter = constructorof(CO)()
    SatConventionalPLLAndDLL(
        sat_state.carrier_doppler,
        sat_state.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        estimator.carrier_loop_filter_bandwidth,
        estimator.code_loop_filter_bandwidth,
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

# Per-sat update for the conventional PLL/DLL estimator. Pure: takes a
# TrackedSat and returns the updated TrackedSat. Shared by the immutable
# `estimate_dopplers_and_filter_prompt` and the in-place
# `estimate_dopplers_and_filter_prompt!` so the two cannot drift.
function _update_tracked_sat_doppler(
    tracked_sat::TrackedSat,
    system,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    sat_state = tracked_sat.sat_state
    pll_and_dll_state = tracked_sat.estimator_state
    if !sat_state.is_integration_completed || sat_state.integrated_samples == 0
        return tracked_sat
    end
    integrated_code_blocks = calc_num_code_blocks_to_integrate(
        system,
        preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(sat_state.bit_buffer),
    )

    normalized_correlator = normalize(sat_state.correlator, sat_state.integrated_samples)
    post_corr_filter =
        update(sat_state.post_corr_filter, get_prompt(normalized_correlator))
    filtered_correlator = apply(post_corr_filter, normalized_correlator)
    prompt = get_prompt(filtered_correlator)
    push!(sat_state.filtered_prompts, prompt)
    cn0_estimator = update(get_cn0_estimator(sat_state), prompt)
    bit_buffer =
        buffer(system, sat_state.bit_buffer, integrated_code_blocks, prompt)
    integration_time = sat_state.integrated_samples / sampling_frequency

    carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
        system,
        pll_and_dll_state.carrier_loop_filter,
        filtered_correlator,
        get_last_fully_integrated_filtered_prompt(sat_state),
        integration_time,
        pll_and_dll_state.carrier_loop_filter_bandwidth,
    )
    code_freq_update, code_loop_filter = calculate_code_frequency_update(
        system,
        pll_and_dll_state.code_loop_filter,
        filtered_correlator,
        sat_state.code_doppler,
        sampling_frequency,
        integration_time,
        pll_and_dll_state.code_loop_filter_bandwidth,
    )
    carrier_doppler, code_doppler = aid_dopplers(
        system,
        pll_and_dll_state.init_carrier_doppler,
        pll_and_dll_state.init_code_doppler,
        carrier_freq_update,
        code_freq_update,
    )
    new_estimator_state = SatConventionalPLLAndDLL(
        pll_and_dll_state;
        carrier_loop_filter,
        code_loop_filter,
    )
    new_sat_state = SatState(
        sat_state;
        carrier_doppler,
        code_doppler,
        integrated_samples = 0,
        is_integration_completed = false,
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        correlator = zero(sat_state.correlator),
        last_fully_integrated_correlator = sat_state.correlator,
    )
    TrackedSat(new_sat_state, new_estimator_state)
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
    # Detach slot vectors from the input, then delegate to the in-place
    # form. The per-sat doppler update is identical between the two —
    # only the storage ownership differs.
    new_track_state = TrackState(
        track_state;
        multiple_system_sats_state =
            _copy_slot_vectors(track_state.multiple_system_sats_state),
    )
    estimate_dopplers_and_filter_prompt!(
        new_track_state, preferred_num_code_blocks_to_integrate, sampling_frequency,
    )
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref). Walks each
system's `Vector{TrackedSat}` backing storage and overwrites slots with the
new immutable `TrackedSat` value. Returns the same `track_state` object —
allocation-free in steady state when [`track!`](@ref)'s preconditions are met.
"""
# Per-system body for the doppler estimator. Pulled out so
# `_foreach_system!` can call it without boxing when the system tuple
# is heterogeneous (e.g. GPS L1 + Galileo E1B).
@inline function _est_one_system!(
    sss::SystemSatsState, preferred_num_code_blocks_to_integrate, sampling_frequency,
)
    system = sss.system
    vals = sss.states.values
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_doppler(
            vals[i],
            system,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )
    end
    return nothing
end

function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:MultipleSystemSatsState,<:ConventionalPLLAndDLL},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    _foreach_system!(
        _est_one_system!, track_state.multiple_system_sats_state,
        preferred_num_code_blocks_to_integrate, sampling_frequency,
    )
    return track_state
end

function calculate_carrier_frequency_update(
    system::AbstractGNSS,
    carrier_loop_filter::ThirdOrderAssistedBilinearLF,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(system, correlator)
    fll_discriminator = fll_disc(system, correlator, previous_prompt, integration_time)
    filter_loop(
        carrier_loop_filter,
        (pll_discriminator, fll_discriminator),
        integration_time,
        loop_bandwidth,
    )
end

function calculate_carrier_frequency_update(
    system::AbstractGNSS,
    carrier_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(system, correlator)
    filter_loop(carrier_loop_filter, pll_discriminator, integration_time, loop_bandwidth)
end

function calculate_code_frequency_update(
    system::AbstractGNSS,
    code_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    code_doppler,
    sampling_frequency,
    integration_time,
    loop_bandwidth,
)
    dll_discriminator = dll_disc(system, correlator, code_doppler, sampling_frequency)
    filter_loop(code_loop_filter, dll_discriminator, integration_time, loop_bandwidth)
end