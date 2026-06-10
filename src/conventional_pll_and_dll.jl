"""
$(SIGNATURES)

Recommended carrier-loop-filter bandwidth for `signal`'s primary integration
period. Sized so that the PLL time-bandwidth product `BL * T` lands at
about 0.018 (≈10× margin from the 0.18 stability edge of the bilinear
third-order filter). Used by [`TrackState(; signal=…)`](@ref) when the
user doesn't pass an explicit `doppler_estimator`.

Override by defining a method for your signal type, or by constructing
[`ConventionalAssistedPLLAndDLL`](@ref) yourself with explicit
`carrier_loop_filter_bandwidth =` / `code_loop_filter_bandwidth =` kwargs.

```julia
T = get_code_length(signal) / get_code_frequency(signal)   # primary period
BL = 0.018 / T                                              # this default
```

For L1 C/A (T = 1 ms) this returns 18 Hz — matching the historical
hand-picked default. For L1C-D / L1C-P / L5I (T = 10 ms) it returns 1.8 Hz,
which is the well-inside-stability value the multi-signal flagship use case
needs.
"""
function default_carrier_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    # T = primary_code_period_seconds. The estimator's bandwidth fields are
    # typed `typeof(1.0Hz)`, so explicitly land on Hz (otherwise `1/s`
    # propagates and trips the typed field assignment).
    primary_period = get_code_length(signal) / get_code_frequency(signal)
    uconvert(Hz, 0.018 / primary_period)
end

"""
$(SIGNATURES)

Recommended code-loop-filter (DLL) bandwidth for `signal`'s primary
integration period. Picks 1/18 of [`default_carrier_loop_filter_bandwidth`](@ref)
— the historical 18:1 carrier:code-bandwidth ratio that gives the DLL good
noise rejection without lagging the PLL.

For L1 C/A this returns 1 Hz; for L1C-D / L1C-P / L5I it returns 0.1 Hz.
Override by defining a method for your signal type.
"""
function default_code_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    default_carrier_loop_filter_bandwidth(signal) / 18
end

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
new satellites. The literal defaults on this constructor are 18 Hz / 1 Hz,
sized for GPS L1 C/A's 1 ms integration; when [`TrackState`](@ref) builds
this estimator implicitly it picks per-signal bandwidths via
[`default_carrier_loop_filter_bandwidth`](@ref) and
[`default_code_loop_filter_bandwidth`](@ref) so longer-period signals
(L1C-D, L1C-P, Galileo E1B) get appropriately tighter loops.
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
    signal::AbstractGNSSSignal,
    init_carrier_doppler,
    init_code_doppler,
    carrier_freq_update,
    code_freq_update,
)
    carrier_doppler = carrier_freq_update
    code_doppler =
        code_freq_update + carrier_doppler * get_code_center_frequency_ratio(signal)
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
    # for `signals[1]` (the estimator-driver signal), run PLL/DLL and
    # update the sat-shared carrier/code Doppler.
    pll_and_dll_state = sat.doppler_estimator_state
    head = first(sat.signals)
    tail_signals = Base.tail(sat.signals)

    new_head, new_doppler_estimator_state, new_carrier_doppler, new_code_doppler =
        _process_estimator_driver_signal(
            head,
            sat,
            pll_and_dll_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
        )

    new_tail = _process_passenger_signals(
        tail_signals,
        sat.prn,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )

    # Phase-snap fallback chain. Picks the synced signal with the
    # longest `(primary × secondary)` code length, and uses its
    # secondary-code phase to anchor `sat.code_phase` to the right
    # secondary-chip window.
    #
    # This is a *one-time* anchoring applied only on the iteration a
    # signal transitions `found == false → true`. The snap aligns
    # `code_phase` to a primary-code boundary within the secondary-code
    # window, which is only valid when `code_phase` actually sits on
    # such a boundary — i.e. just after a completed integration, which
    # is exactly when sync is detected. Re-running it on later
    # iterations would clobber the within-primary-block phase that
    # `update` advances during a partial (chunk-bounded) integration,
    # pinning `code_phase` to the boundary and wedging the satellite
    # into a state where no chunk can ever complete the block again
    # (see issue #117). After sync, `update`'s `mod(…, current_code_wrap)`
    # maintains the secondary alignment on its own.
    new_signals = (new_head, new_tail...)
    snapped_code_phase = if _any_signal_just_synced(sat.signals, new_signals)
        _snap_code_phase_from_synced_signal(new_signals, sat.code_phase)
    else
        sat.code_phase
    end

    TrackedSat(
        sat;
        code_phase = snapped_code_phase,
        carrier_doppler = new_carrier_doppler,
        code_doppler = new_code_doppler,
        signals = new_signals,
        doppler_estimator_state = new_doppler_estimator_state,
    )
end

# Process the estimator-driver signal (signals[1]): if its integration
# completed, run the PLL/DLL plus prompt filter / CN0 / bit-buffer update
# and return new values for carrier_doppler, code_doppler, and
# doppler_estimator_state. Otherwise return unchanged values. This is
# where ConventionalPLLAndDLL hard-codes the "signals[1] drives the loop
# filter" rule — a custom AbstractDopplerEstimator may use any/all
# signals' state.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    pll_and_dll_state::SatConventionalPLLAndDLL,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    if !tracked_signal.is_integration_completed || tracked_signal.integrated_samples == 0
        return tracked_signal, pll_and_dll_state, sat.carrier_doppler, sat.code_doppler
    end
    signal = tracked_signal.signal
    integrated_code_blocks = calc_num_code_blocks_to_integrate(
        signal,
        preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    )

    normalized_correlator = normalize(tracked_signal.correlator, tracked_signal.integrated_samples)
    post_corr_filter = update(tracked_signal.post_corr_filter, get_prompt(normalized_correlator))
    filtered_correlator = apply(post_corr_filter, normalized_correlator)
    prompt = get_prompt(filtered_correlator)
    push!(tracked_signal.filtered_prompts, prompt)
    cn0_estimator = update(get_cn0_estimator(tracked_signal), prompt)
    bit_buffer = buffer(signal, sat.prn, tracked_signal.bit_buffer, integrated_code_blocks, prompt)
    integration_time = tracked_signal.integrated_samples / sampling_frequency

    carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
        signal,
        pll_and_dll_state.carrier_loop_filter,
        filtered_correlator,
        get_last_fully_integrated_filtered_prompt(tracked_signal),
        integration_time,
        pll_and_dll_state.carrier_loop_filter_bandwidth,
    )
    code_freq_update, code_loop_filter = calculate_code_frequency_update(
        signal,
        pll_and_dll_state.code_loop_filter,
        filtered_correlator,
        sat.code_doppler,
        sampling_frequency,
        integration_time,
        pll_and_dll_state.code_loop_filter_bandwidth,
    )
    carrier_doppler, code_doppler = aid_dopplers(
        signal,
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
        tracked_signal;
        integrated_samples = 0,
        is_integration_completed = false,
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        correlator = zero(tracked_signal.correlator),
        last_fully_integrated_correlator = tracked_signal.correlator,
    )
    return new_signal, new_doppler_estimator_state, carrier_doppler, code_doppler
end

# Process the non-driver signals (signals[2:end]): per-signal prompt
# filter, CN0 update, bit-buffer advance, correlator hand-off. No loop-
# filter work. Walks the tuple recursively to keep type-stability and
# avoid boxing.
@inline _process_passenger_signals(
    ::Tuple{}, _, _,
) = ()
@inline _process_passenger_signals(
    ::Tuple{}, prn::Integer, _, _,
) = ()
@inline function _process_passenger_signals(
    signals::Tuple,
    prn::Integer,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    head = first(signals)
    new_head = _process_one_passenger_signal(
        head,
        prn,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )
    (new_head, _process_passenger_signals(
        Base.tail(signals),
        prn,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
    )...)
end

@inline function _process_one_passenger_signal(
    tracked_signal::TrackedSignal,
    prn::Integer,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
)
    if !tracked_signal.is_integration_completed || tracked_signal.integrated_samples == 0
        return tracked_signal
    end
    signal = tracked_signal.signal
    integrated_code_blocks = calc_num_code_blocks_to_integrate(
        signal,
        preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    )
    normalized_correlator = normalize(tracked_signal.correlator, tracked_signal.integrated_samples)
    post_corr_filter = update(tracked_signal.post_corr_filter, get_prompt(normalized_correlator))
    filtered_correlator = apply(post_corr_filter, normalized_correlator)
    prompt = get_prompt(filtered_correlator)
    push!(tracked_signal.filtered_prompts, prompt)
    cn0_estimator = update(get_cn0_estimator(tracked_signal), prompt)
    bit_buffer = buffer(signal, prn, tracked_signal.bit_buffer, integrated_code_blocks, prompt)
    TrackedSignal(
        tracked_signal;
        integrated_samples = 0,
        is_integration_completed = false,
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        correlator = zero(tracked_signal.correlator),
        last_fully_integrated_correlator = tracked_signal.correlator,
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
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    measurements::Measurements,
    preferred_num_code_blocks_to_integrate,
)
    # Detach slot vectors from the input, then delegate to the in-place
    # form. The per-sat doppler update is identical between the two —
    # only the storage ownership differs.
    new_track_state = TrackState(
        track_state;
        groups = _copy_groups_slot_vectors(track_state.groups),
    )
    estimate_dopplers_and_filter_prompt!(
        new_track_state, measurements, preferred_num_code_blocks_to_integrate,
    )
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref). Walks each
group's `Vector{TrackedSat}` backing storage and overwrites slots with the
new immutable `TrackedSat` value. Returns the same `track_state` object —
allocation-free in steady state when [`track!`](@ref)'s preconditions are met.
"""
# Per-group body for the doppler estimator. Pulled out so
# `_foreach_group!` can call it without boxing when the groups tuple
# is heterogeneous (e.g. GPS L1 + Galileo E1B). The per-signal type
# is recovered from each signal inside `_update_tracked_sat_doppler`.
# Routes to this group's band's `Measurement` for sampling frequency.
@inline function _est_one_group!(
    g::SignalGroup,
    measurements::Measurements,
    preferred_num_code_blocks_to_integrate,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    sampling_frequency = measurements[band_key(g.band)].sampling_frequency
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
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    measurements::Measurements,
    preferred_num_code_blocks_to_integrate,
)
    _foreach_group!(
        _est_one_group!, track_state.groups,
        measurements, preferred_num_code_blocks_to_integrate,
    )
    return track_state
end

function calculate_carrier_frequency_update(
    signal::AbstractGNSSSignal,
    carrier_loop_filter::ThirdOrderAssistedBilinearLF,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(signal, correlator)
    fll_discriminator = fll_disc(signal, correlator, previous_prompt, integration_time)
    filter_loop(
        carrier_loop_filter,
        (pll_discriminator, fll_discriminator),
        integration_time,
        loop_bandwidth,
    )
end

function calculate_carrier_frequency_update(
    signal::AbstractGNSSSignal,
    carrier_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(signal, correlator)
    filter_loop(carrier_loop_filter, pll_discriminator, integration_time, loop_bandwidth)
end

function calculate_code_frequency_update(
    signal::AbstractGNSSSignal,
    code_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    code_doppler,
    sampling_frequency,
    integration_time,
    loop_bandwidth,
)
    dll_discriminator = dll_disc(signal, correlator, code_doppler, sampling_frequency)
    filter_loop(code_loop_filter, dll_discriminator, integration_time, loop_bandwidth)
end