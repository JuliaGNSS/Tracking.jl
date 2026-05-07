"""
$(SIGNATURES)

Holds the state of a single satellite being tracked. Contains code and carrier
phase/doppler information, correlator state, CN0 estimator, and bit buffer.
"""
struct SatState{C<:AbstractCorrelator,PCF<:AbstractPostCorrFilter}
    prn::Int
    code_phase::Float64
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    carrier_doppler::typeof(1.0Hz)
    integrated_samples::Int
    signal_start_sample::Int
    is_integration_completed::Bool
    correlator::C
    last_fully_integrated_correlator::C
    last_fully_integrated_filtered_prompt::ComplexF64
    cn0_estimator::MomentsCN0Estimator
    bit_buffer::BitBuffer
    post_corr_filter::PCF
    filtered_prompts::Vector{ComplexF64}
end

"""
$(SIGNATURES)

Get the PRN (Pseudo-Random Noise) number of the satellite.
"""
get_prn(s::SatState) = s.prn
get_num_ants(s::SatState{<:AbstractCorrelator{M}}) where {M} = M

"""
$(SIGNATURES)

Get the current code phase in chips.
"""
get_code_phase(s::SatState) = s.code_phase

"""
$(SIGNATURES)

Get the current code Doppler frequency.
"""
get_code_doppler(s::SatState) = s.code_doppler

"""
$(SIGNATURES)

Get the current carrier phase in radians.
"""
get_carrier_phase(s::SatState) = s.carrier_phase * 2π

"""
$(SIGNATURES)

Get the current carrier Doppler frequency.
"""
get_carrier_doppler(s::SatState) = s.carrier_doppler

"""
$(SIGNATURES)

Get the number of samples that have been integrated so far.
"""
get_integrated_samples(s::SatState) = s.integrated_samples

"""
$(SIGNATURES)

Get the starting sample index in the signal for the next integration.
"""
get_signal_start_sample(s::SatState) = s.signal_start_sample

"""
$(SIGNATURES)

Get the current correlator state.
"""
get_correlator(s::SatState) = s.correlator

"""
$(SIGNATURES)

Get the correlator from the last fully completed integration period.
"""
get_last_fully_integrated_correlator(s::SatState) = s.last_fully_integrated_correlator

"""
$(SIGNATURES)

Get the filtered prompt value from the last fully completed integration period.
"""
get_last_fully_integrated_filtered_prompt(s::SatState) =
    s.last_fully_integrated_filtered_prompt

"""
$(SIGNATURES)

Get all filtered prompt correlation results that were produced during the most
recent `track` call. The vector is reset at the start of each `track` call and
appended to for every completed integration, so its length equals the number of
correlations that completed within that call.
"""
get_filtered_prompts(s::SatState) = s.filtered_prompts
get_post_corr_filter(s::SatState) = s.post_corr_filter
get_cn0_estimator(s::SatState) = s.cn0_estimator

"""
$(SIGNATURES)

Get the bit buffer containing decoded navigation bits.
"""
get_bit_buffer(s::SatState) = s.bit_buffer

"""
$(SIGNATURES)

Get the decoded navigation bits as an integer.
"""
get_bits(s::SatState) = get_bits(get_bit_buffer(s))

"""
$(SIGNATURES)

Get the number of decoded navigation bits.
"""
get_num_bits(s::SatState) = length(get_bit_buffer(s))

"""
$(SIGNATURES)

Check if the bit or secondary code synchronization has been achieved.
"""
has_bit_or_secondary_code_been_found(s::SatState) =
    has_bit_or_secondary_code_been_found(get_bit_buffer(s))

function SatState(
    system::AbstractGNSS,
    prn::Int,
    code_phase,
    carrier_doppler;
    num_ants::NumAnts = NumAnts(1),
    correlator = get_default_correlator(system, num_ants),
    carrier_phase = 0.0,
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(system),
    num_prompts_for_cn0_estimation::Int = 100,
    post_corr_filter::AbstractPostCorrFilter = DefaultPostCorrFilter(),
)
    SatState(
        prn,
        float(code_phase),
        float(code_doppler),
        float(carrier_phase) / 2π,
        float(carrier_doppler),
        0,
        1,
        false,
        correlator,
        correlator,
        complex(0.0, 0.0),
        MomentsCN0Estimator(num_prompts_for_cn0_estimation),
        BitBuffer(),
        post_corr_filter,
        ComplexF64[],
    )
end

function SatState(
    sat_state::SatState{C,PCF};
    prn = nothing,
    code_phase = nothing,
    code_doppler = nothing,
    carrier_phase = nothing,
    carrier_doppler = nothing,
    integrated_samples = nothing,
    signal_start_sample = nothing,
    is_integration_completed = nothing,
    correlator::Maybe{C} = nothing,
    last_fully_integrated_correlator = nothing,
    last_fully_integrated_filtered_prompt = nothing,
    cn0_estimator = nothing,
    bit_buffer = nothing,
    post_corr_filter::Maybe{PCF} = nothing,
    filtered_prompts::Maybe{Vector{ComplexF64}} = nothing,
) where {C<:AbstractCorrelator,PCF<:AbstractPostCorrFilter}
    SatState{C,PCF}(
        isnothing(prn) ? sat_state.prn : prn,
        isnothing(code_phase) ? sat_state.code_phase : code_phase,
        isnothing(code_doppler) ? sat_state.code_doppler : code_doppler,
        isnothing(carrier_phase) ? sat_state.carrier_phase : carrier_phase,
        isnothing(carrier_doppler) ? sat_state.carrier_doppler : carrier_doppler,
        isnothing(integrated_samples) ? sat_state.integrated_samples : integrated_samples,
        isnothing(signal_start_sample) ? sat_state.signal_start_sample :
        signal_start_sample,
        isnothing(is_integration_completed) ? sat_state.is_integration_completed :
        is_integration_completed,
        isnothing(correlator) ? sat_state.correlator : correlator,
        isnothing(last_fully_integrated_correlator) ?
        sat_state.last_fully_integrated_correlator : last_fully_integrated_correlator,
        isnothing(last_fully_integrated_filtered_prompt) ?
        sat_state.last_fully_integrated_filtered_prompt :
        last_fully_integrated_filtered_prompt,
        isnothing(cn0_estimator) ? sat_state.cn0_estimator : cn0_estimator,
        isnothing(bit_buffer) ? sat_state.bit_buffer : bit_buffer,
        isnothing(post_corr_filter) ? sat_state.post_corr_filter : post_corr_filter,
        isnothing(filtered_prompts) ? sat_state.filtered_prompts : filtered_prompts,
    )
end

function reset_start_sample_and_bit_buffer(sat_state::SatState)
    empty!(sat_state.filtered_prompts)
    SatState(sat_state; signal_start_sample = 1, bit_buffer = reset(sat_state.bit_buffer))
end

"""
$(SIGNATURES)

Wrapper holding a [`SatState`](@ref) together with the per-satellite state used
by the active doppler estimator. Per-satellite estimator state lives here
rather than in a parallel dictionary on the estimator, so the two are always
updated together and traversed in a single `map`.

`E` is the estimator-specific per-satellite state type. Each
[`AbstractDopplerEstimator`](@ref) provides an
[`init_estimator_state`](@ref) method that produces an `E` value for a given
[`SatState`](@ref); the wrapper is built once per satellite at the
acquisition→tracking handoff and updated in lockstep thereafter.
"""
struct TrackedSat{S<:SatState,E}
    sat_state::S
    estimator_state::E
end

get_sat_state(t::TrackedSat) = t.sat_state
get_estimator_state(t::TrackedSat) = t.estimator_state

function reset_start_sample_and_bit_buffer(t::TrackedSat)
    TrackedSat(reset_start_sample_and_bit_buffer(t.sat_state), t.estimator_state)
end

"""
$(SIGNATURES)

Build the per-satellite state used by `estimator` for `sat_state`. Called once
per satellite when a new sat enters the track set (see [`merge_sats`](@ref)).

A custom doppler estimator must define this method for its
[`AbstractDopplerEstimator`](@ref) subtype.
"""
function init_estimator_state end

"""
$(SIGNATURES)

Optionally update an estimator's cross-satellite or cross-system shared state
when new satellites enter the track set. Called once per [`merge_sats`](@ref)
call with the dictionary of incoming `SatState`s, *after* per-sat seeding via
[`init_estimator_state`](@ref).

The default returns `estimator` unchanged, so estimators with no shared state
need not implement it.

The returned estimator must have the same concrete type as the input —
[`TrackState`](@ref) is parameterized on the estimator type, and changing it
would break inference. For growing shared state, hold the storage in a
resizable container (e.g. `Vector`, `Matrix`) on an otherwise immutable
estimator and `push!`/`resize!` it in place; rebuild the estimator with
[`Setfield.@set`](https://jw3126.github.io/Setfield.jl/stable/) or a copying
constructor when fields need replacing.
"""
update_estimator_on_handoff(estimator::AbstractDopplerEstimator, _new_sats) = estimator

"""
$(SIGNATURES)

Holds the state of multiple satellites for a single GNSS system. Contains the
system definition and a dictionary of [`TrackedSat`](@ref) entries indexed by
identifier.
"""
struct SystemSatsState{S<:AbstractGNSS,T<:TrackedSat,I}
    system::S
    states::Dictionary{I,T}
end

"""
Type alias for a tuple or named tuple of `SystemSatsState` objects, representing
tracking state across multiple GNSS systems.
"""
const MultipleSystemSatsState{N,I,S,T} =
    TupleLike{<:NTuple{N,SystemSatsState{<:S,<:T,<:I}}}

"""
$(SIGNATURES)

Merge already-wrapped tracked satellites into the existing tracking state.
Used internally by [`TrackState`](@ref)'s `merge_sats` after wrapping incoming
`SatState`s with their estimator state — external callers should prefer the
`TrackState`-level method.
"""
function merge_sats(
    multiple_system_sats_state::MultipleSystemSatsState,
    system_idx,
    new_tracked_sats::Dictionary{<:Any,<:TrackedSat},
)
    system_sats_state = get_system_sats_state(multiple_system_sats_state, system_idx)
    @set multiple_system_sats_state[system_idx].states =
        merge(system_sats_state.states, new_tracked_sats)
end

"""
$(SIGNATURES)

Remove satellites with the specified identifiers from the tracking state.
"""
function filter_out_sats(
    multiple_system_sats_state::MultipleSystemSatsState,
    system_idx::Union{Symbol,Integer},
    identifiers,
)
    filtered_sat_states = map(
        last,
        filter(
            ((id,),) -> !in(id, identifiers),
            pairs(multiple_system_sats_state[system_idx].states),
        ),
    )
    @set multiple_system_sats_state[system_idx].states = filtered_sat_states
end

function reset_start_sample_and_bit_buffer(
    multiple_system_sats_state::MultipleSystemSatsState,
)
    map(multiple_system_sats_state) do system_sats_state
        new_sat_states = map(reset_start_sample_and_bit_buffer, system_sats_state.states)
        SystemSatsState(system_sats_state, new_sat_states)
    end
end

# In-place variants: walk Vector{TrackedSat} slots and overwrite each entry
# with the reset value. The vector itself, the Dictionary, and the
# SystemSatsState/MultipleSystemSatsState wrappers are all reused.
function reset_start_sample_and_bit_buffer!(t::TrackedSat)
    TrackedSat(reset_start_sample_and_bit_buffer(t.sat_state), t.estimator_state)
end

function reset_start_sample_and_bit_buffer!(
    multiple_system_sats_state::MultipleSystemSatsState,
)
    for system_sats_state in multiple_system_sats_state
        vals = system_sats_state.states.values
        @inbounds for i in eachindex(vals)
            vals[i] = reset_start_sample_and_bit_buffer!(vals[i])
        end
    end
    return multiple_system_sats_state
end

function to_dictionary(sat_states::Dictionary{I,<:SatState}) where {I}
    sat_states
end

function to_dictionary(sat_states::Vector{<:SatState})
    Dictionary(map(get_prn, sat_states), sat_states)
end

function to_dictionary(sat_state::SatState)
    dictionary((get_prn(sat_state) => sat_state,))
end

function to_dictionary(tracked_sats::Dictionary{I,<:TrackedSat}) where {I}
    tracked_sats
end

function to_dictionary(tracked_sats::Vector{<:TrackedSat})
    Dictionary(map(t -> get_prn(t.sat_state), tracked_sats), tracked_sats)
end

function to_dictionary(t::TrackedSat)
    dictionary((get_prn(t.sat_state) => t,))
end

"""
$(SIGNATURES)

Wrap each `SatState` in a [`TrackedSat`](@ref) by pairing it with the
per-satellite estimator state produced by [`init_estimator_state`](@ref).
"""
function wrap_sats(estimator, sat_states::Dictionary{I,<:SatState}) where {I}
    map(s -> TrackedSat(s, init_estimator_state(estimator, s)), sat_states)
end

function wrap_sats(estimator, sat_states::Vector{<:SatState})
    wrap_sats(estimator, to_dictionary(sat_states))
end

function wrap_sats(estimator, sat_state::SatState)
    wrap_sats(estimator, to_dictionary(sat_state))
end

"""
$(SIGNATURES)

Build a `SystemSatsState` from already-wrapped tracked sats. Accepts a
`Vector{TrackedSat}`, a `Dictionary{I,TrackedSat}`, or a single `TrackedSat`.

To build from raw [`SatState`](@ref) input, supply the doppler estimator —
the convenience method [`SystemSatsState`](@ref)`(estimator, system, sats)`
wraps each sat via [`init_estimator_state`](@ref).
"""
function SystemSatsState(
    system::AbstractGNSS,
    tracked_sats::Union{TrackedSat,Vector{<:TrackedSat},Dictionary{<:Any,<:TrackedSat}},
)
    SystemSatsState(system, to_dictionary(tracked_sats))
end

function SystemSatsState(
    system_sats_state::SystemSatsState,
    states::Dictionary{I,<:TrackedSat},
) where {I}
    SystemSatsState(system_sats_state.system, states)
end

"""
$(SIGNATURES)

Convenience constructor that wraps each `SatState` with its initial estimator
state via [`init_estimator_state`](@ref) before building the
`SystemSatsState`. Useful when assembling multi-system tracking state by
hand.
"""
function SystemSatsState(
    estimator::AbstractDopplerEstimator,
    system::AbstractGNSS,
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{<:Any,<:SatState}},
)
    SystemSatsState(system, wrap_sats(estimator, sat_states))
end

"""
$(SIGNATURES)

Get the GNSS system definition from a SystemSatsState.
"""
get_system(sss::SystemSatsState) = sss.system
get_states(sss::SystemSatsState) = sss.states

"""
$(SIGNATURES)

Get the satellite state for a specific satellite identifier. Unwraps the
internal [`TrackedSat`](@ref) — callers see the [`SatState`](@ref) directly.
"""
get_sat_state(sss::SystemSatsState, identifier) = sss.states[identifier].sat_state
get_sat_state(sss::SystemSatsState) = only(sss.states).sat_state

"""
$(SIGNATURES)

Get the per-satellite estimator state for a specific satellite identifier.
"""
get_estimator_state(sss::SystemSatsState, identifier) =
    sss.states[identifier].estimator_state
get_estimator_state(sss::SystemSatsState) = only(sss.states).estimator_state

function estimate_cn0(sss::SystemSatsState, id...)
    estimate_cn0(sss.system, get_sat_state(sss, id...))
end

function estimate_cn0(system::AbstractGNSS, sat_state::SatState)
    estimate_cn0(
        get_cn0_estimator(sat_state),
        get_code_length(system) / get_code_frequency(system),
    )
end
