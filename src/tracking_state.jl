function TrackState(
    system::AbstractGNSS,
    sat_states;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    wrapped = wrap_sats(doppler_estimator, sat_states)
    tracked_systems = (TrackedSystem(system, wrapped),)
    TrackState(tracked_systems, doppler_estimator)
end

function TrackState(
    tracked_system::TrackedSystem;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    tracked_systems = (tracked_system,)
    TrackState(tracked_systems, doppler_estimator)
end

function TrackState(
    tracked_systems::TrackedSystems;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(),
)
    TrackState(tracked_systems, doppler_estimator)
end

function TrackState(
    track_state::TrackState{S,DE},
    tracked_systems::S,
    doppler_estimator::DE,
) where {S<:TrackedSystems,DE<:AbstractDopplerEstimator}
    TrackState{S,DE}(tracked_systems, doppler_estimator)
end

# Be careful when calling this.
# It might lead to types that are inferred at runtime?!
# Tested with 1.11.6
function TrackState(
    track_state::TrackState{S,DE};
    tracked_systems::Maybe{S} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
) where {S<:TrackedSystems,DE<:AbstractDopplerEstimator}
    TrackState{S,DE}(
        isnothing(tracked_systems) ? track_state.tracked_systems :
        tracked_systems,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
    )
end

function reset_start_sample_and_bit_buffer(track_state::TrackState)
    TrackState(
        track_state;
        tracked_systems = reset_start_sample_and_bit_buffer(
            track_state.tracked_systems,
        ),
    )
end

function reset_start_sample_and_bit_buffer!(track_state::TrackState)
    reset_start_sample_and_bit_buffer!(track_state.tracked_systems)
    return track_state
end

function has_integration_reached_signal_end_for_all_satellites(
    track_state::TrackState,
    num_samples::Int,
)
    target = num_samples + 1
    # NamedTuple unwraps to its underlying Tuple via Tuple(...) — concrete and
    # cheap.
    _all_sats_at(Tuple(track_state.tracked_systems), target)
end

# Recursive tuple-walk: each step has fully concrete types, no closure boxes,
# no allocation. Base case is an empty tuple.
@inline _all_sats_at(::Tuple{}, target::Int) = true
@inline function _all_sats_at(t::Tuple, target::Int)
    sss = first(t)
    @inbounds for tracked_sat in sss.states.values
        tracked_sat.sat_state.signal_start_sample == target || return false
    end
    _all_sats_at(Base.tail(t), target)
end

"""
$(SIGNATURES)

Get the TrackedSystem for a specific GNSS system from a TrackState or
TrackedSystems by index or symbol.
"""
function get_tracked_system(
    track_state::TrackState{<:TrackedSystems{N}},
    system_idx,
) where {N}
    get_tracked_system(track_state.tracked_systems, system_idx)
end

function get_tracked_system(
    tracked_systems::TrackedSystems{N},
    system_idx::Union{Symbol,Integer},
) where {N}
    tracked_systems[system_idx]
end

function get_tracked_system(
    tracked_systems::TrackedSystems{N},
    system_idx::Val{M},
) where {N,M}
    tracked_systems[M]
end

"""
$(SIGNATURES)

Get the dictionary of satellite states for a specific GNSS system.
"""
function get_sat_states(
    tracked_systems::TrackedSystems{N},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_tracked_system(tracked_systems, system_idx).states
end

function get_sat_states(tracked_systems::TrackedSystems{1})
    get_sat_states(tracked_systems, 1)
end

function get_sat_states(
    track_state::TrackState{<:TrackedSystems{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sat_states(track_state.tracked_systems, system_idx)
end

function get_sat_states(track_state::TrackState{<:TrackedSystems{1}})
    get_sat_states(track_state.tracked_systems)
end

function get_system(
    track_state::TrackState{<:TrackedSystems{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_tracked_system(track_state, system_idx).system
end

function get_system(track_state::TrackState{<:TrackedSystems{1}})
    get_system(track_state, 1)
end

function get_sat_state(
    track_state::TrackState{<:TrackedSystems{N}},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
) where {N}
    get_sat_state(get_tracked_system(track_state, system_idx), sat_identifier)
end

function get_sat_state(
    track_state::TrackState{<:TrackedSystems{1}},
    sat_identifier,
)
    get_sat_state(track_state, 1, sat_identifier)
end

function get_sat_state(track_state::TrackState{<:TrackedSystems{1}})
    only(get_sat_states(track_state, 1)).sat_state
end

function estimate_cn0(
    track_state::TrackState{<:TrackedSystems},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
)
    estimate_cn0(get_tracked_system(track_state, system_idx), sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:TrackedSystems{1}}, sat_identifier)
    estimate_cn0(track_state, 1, sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:TrackedSystems{1}})
    estimate_cn0(get_tracked_system(track_state, 1))
end

function merge_sats(
    track_state::TrackState{S,DE},
    system_idx::Union{Symbol,Integer},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{<:Any,<:SatState}},
) where {S<:TrackedSystems,DE<:AbstractDopplerEstimator}
    new_sats_dict = to_dictionary(sat_states)
    wrapped = wrap_sats(track_state.doppler_estimator, new_sats_dict)
    new_estimator =
        update_estimator_on_handoff(track_state.doppler_estimator, new_sats_dict)
    TrackState{S,DE}(
        merge_sats(track_state.tracked_systems, system_idx, wrapped),
        new_estimator,
    )
end

function merge_sats(
    track_state::TrackState{<:TrackedSystems{1}},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{Any,<:SatState}},
)
    merge_sats(track_state, 1, sat_states)
end

function filter_out_sats(
    track_state::TrackState{S,DE},
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {S<:TrackedSystems,DE<:AbstractDopplerEstimator}
    TrackState{S,DE}(
        filter_out_sats(track_state.tracked_systems, system_idx, identifiers),
        track_state.doppler_estimator,
    )
end

function filter_out_sats(track_state::TrackState{<:TrackedSystems{1}}, identifiers)
    filter_out_sats(track_state, 1, identifiers)
end

# Convenient methods
get_prn(s::TrackState, id...) = get_prn(get_sat_state(s, id...))
get_num_ants(s::TrackState, id...) = get_num_ants(get_sat_state(s, id...))
get_code_phase(s::TrackState, id...) = get_code_phase(get_sat_state(s, id...))
get_code_doppler(s::TrackState, id...) = get_code_doppler(get_sat_state(s, id...))
get_carrier_phase(s::TrackState, id...) = get_carrier_phase(get_sat_state(s, id...))
get_carrier_doppler(s::TrackState, id...) = get_carrier_doppler(get_sat_state(s, id...))
get_integrated_samples(s::TrackState, id...) =
    get_integrated_samples(get_sat_state(s, id...))
get_signal_start_sample(s::TrackState, id...) =
    get_signal_start_sample(get_sat_state(s, id...))
get_correlator(s::TrackState, id...) = get_correlator(get_sat_state(s, id...))
get_last_fully_integrated_correlator(s::TrackState, id...) =
    get_last_fully_integrated_correlator(get_sat_state(s, id...))
get_last_fully_integrated_filtered_prompt(s::TrackState, id...) =
    get_last_fully_integrated_filtered_prompt(get_sat_state(s, id...))
get_post_corr_filter(s::TrackState, id...) = get_post_corr_filter(get_sat_state(s, id...))
get_cn0_estimator(s::TrackState, id...) = get_cn0_estimator(get_sat_state(s, id...))
get_bit_buffer(s::TrackState, id...) = get_bit_buffer(get_sat_state(s, id...))
get_bits(s::TrackState, id...) = get_bits(get_sat_state(s, id...))
get_num_bits(s::TrackState, id...) = get_num_bits(get_sat_state(s, id...))
has_bit_or_secondary_code_been_found(s::TrackState, id...) =
    has_bit_or_secondary_code_been_found(get_sat_state(s, id...))