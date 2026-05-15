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
    sat::TrackedSat,
    carrier_loop_filter::CA,
    code_loop_filter::CO;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
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
    sat::TrackedSat,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    carrier_loop_filter = constructorof(CA)()
    code_loop_filter = constructorof(CO)()
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
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
    system::AbstractGNSSSignal,
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
    sat::TrackedSat,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    # Walk all signals. For each one whose integration completed this
    # iteration, normalize/filter its prompt, advance CN0 and bit buffer,
    # and move its correlator to `last_fully_integrated_*`. Additionally,
    # for `signals[1]` (the Doppler source), run PLL/DLL and update the
    # sat-shared carrier/code Doppler.
    pll_and_dll_state = sat.doppler_estimator_state
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    new_head, new_doppler_estimator_state, new_carrier_doppler, new_code_doppler =
        _process_doppler_source_signal(
            head,
            sat,
            pll_and_dll_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )

    new_tail = _process_passenger_signals(
        tail_signals,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    TrackedSat(
        sat;
        carrier_doppler = new_carrier_doppler,
        code_doppler = new_code_doppler,
        signals = (new_head, new_tail...),
        doppler_estimator_state = new_doppler_estimator_state,
    )
end

# Process the Doppler-source signal: if its integration completed, run the
# PLL/DLL plus prompt filter / CN0 / bit-buffer update and return new
# values for carrier_doppler, code_doppler, and doppler_estimator_state.
# Otherwise return unchanged values.
@inline function _process_doppler_source_signal(
    signal::TrackedSignal,
    sat::TrackedSat,
    pll_and_dll_state::SatConventionalPLLAndDLL,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    if !signal.is_integration_completed || signal.integrated_samples == 0
        return signal, pll_and_dll_state, sat.carrier_doppler, sat.code_doppler
    end
    system = signal.signal
    integrated_code_blocks = calc_num_code_blocks_to_integrate(
        system,
        preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(signal.bit_buffer),
    )

    normalized_correlator = normalize(signal.correlator, signal.integrated_samples)
    post_corr_filter = update(signal.post_corr_filter, get_prompt(normalized_correlator))
    filtered_correlator = apply(post_corr_filter, normalized_correlator)
    prompt = get_prompt(filtered_correlator)
    push!(signal.filtered_prompts, prompt)
    cn0_estimator = update(get_cn0_estimator(signal), prompt)
    bit_buffer = buffer(system, signal.bit_buffer, integrated_code_blocks, prompt)
    integration_time = signal.integrated_samples / sampling_frequency

    carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
        system,
        pll_and_dll_state.carrier_loop_filter,
        filtered_correlator,
        get_last_fully_integrated_filtered_prompt(signal),
        integration_time,
        pll_and_dll_state.carrier_loop_filter_bandwidth,
    )
    code_freq_update, code_loop_filter = calculate_code_frequency_update(
        system,
        pll_and_dll_state.code_loop_filter,
        filtered_correlator,
        sat.code_doppler,
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
    new_doppler_estimator_state = SatConventionalPLLAndDLL(
        pll_and_dll_state;
        carrier_loop_filter,
        code_loop_filter,
    )
    new_signal = TrackedSignal(
        signal;
        integrated_samples = 0,
        is_integration_completed = false,
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        correlator = zero(signal.correlator),
        last_fully_integrated_correlator = signal.correlator,
    )
    return new_signal, new_doppler_estimator_state, carrier_doppler, code_doppler
end

# Process the non-Doppler-source signals: per-signal prompt filter, CN0
# update, bit-buffer advance, correlator hand-off. No loop-filter work.
# Walks the tuple recursively to keep type-stability and avoid boxing.
@inline _process_passenger_signals(
    ::Tuple{}, _, _,
) = ()
@inline function _process_passenger_signals(
    signals::Tuple,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    head = first(signals)
    new_head = _process_one_passenger_signal(
        head,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )
    (new_head, _process_passenger_signals(
        Base.tail(signals),
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )...)
end

@inline function _process_one_passenger_signal(
    signal::TrackedSignal,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    if !signal.is_integration_completed || signal.integrated_samples == 0
        return signal
    end
    system = signal.signal
    integrated_code_blocks = calc_num_code_blocks_to_integrate(
        system,
        preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(signal.bit_buffer),
    )
    normalized_correlator = normalize(signal.correlator, signal.integrated_samples)
    post_corr_filter = update(signal.post_corr_filter, get_prompt(normalized_correlator))
    filtered_correlator = apply(post_corr_filter, normalized_correlator)
    prompt = get_prompt(filtered_correlator)
    push!(signal.filtered_prompts, prompt)
    cn0_estimator = update(get_cn0_estimator(signal), prompt)
    bit_buffer = buffer(system, signal.bit_buffer, integrated_code_blocks, prompt)
    TrackedSignal(
        signal;
        integrated_samples = 0,
        is_integration_completed = false,
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        correlator = zero(signal.correlator),
        last_fully_integrated_correlator = signal.correlator,
    )
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
    track_state::TrackState{<:TrackedSystems,<:ConventionalPLLAndDLL},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    # Detach slot vectors from the input, then delegate to the in-place
    # form. The per-sat doppler update is identical between the two —
    # only the storage ownership differs.
    new_track_state = TrackState(
        track_state;
        tracked_systems =
            _copy_slot_vectors(track_state.tracked_systems),
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
# is heterogeneous (e.g. GPS L1 + Galileo E1B). The per-signal system
# type is recovered from each signal inside `_update_tracked_sat_doppler`.
@inline function _est_one_system!(
    sats::Dictionary{<:Any,<:TrackedSat},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    vals = sats.values
    isempty(vals) && return nothing
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_doppler(
            vals[i],
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )
    end
    return nothing
end

function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:TrackedSystems,<:ConventionalPLLAndDLL},
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    _foreach_system!(
        _est_one_system!, track_state.tracked_systems,
        preferred_num_code_blocks_to_integrate, sampling_frequency,
    )
    return track_state
end

function calculate_carrier_frequency_update(
    system::AbstractGNSSSignal,
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
    system::AbstractGNSSSignal,
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
    system::AbstractGNSSSignal,
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