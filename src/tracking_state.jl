function TrackState(
    system::AbstractGNSS,
    sat_states;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL((
        SystemSatsState(system, sat_states),
    )),
)
    multiple_system_sats_state = (SystemSatsState(system, sat_states),)
    TrackState(multiple_system_sats_state, doppler_estimator)
end

function TrackState(
    system_sat_states::SystemSatsState;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL((
        system_sat_states,
    )),
)
    multiple_system_sats_state = (system_sat_states,)
    TrackState(multiple_system_sats_state, doppler_estimator)
end

function TrackState(
    multiple_system_sats_state::MultipleSystemSatsState;
    doppler_estimator::AbstractDopplerEstimator = ConventionalAssistedPLLAndDLL(
        multiple_system_sats_state,
    ),
)
    TrackState(multiple_system_sats_state, doppler_estimator)
end

function TrackState(
    track_state::TrackState{S,DE},
    multiple_system_sats_state::S,
    doppler_estimator::DE,
) where {N,I,S<:MultipleSystemSatsState{N,I},DE<:AbstractDopplerEstimator{N,I}}
    TrackState{S,DE}(multiple_system_sats_state, doppler_estimator)
end

# Be careful when calling this.
# It might lead to types that are inferred at runtime?!
# Tested with 1.11.6
function TrackState(
    track_state::TrackState{S,DE};
    multiple_system_sats_state::Maybe{S} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
) where {N,I,S<:MultipleSystemSatsState{N,I},DE<:AbstractDopplerEstimator{N,I}}
    TrackState{S,DE}(
        isnothing(multiple_system_sats_state) ? track_state.multiple_system_sats_state :
        multiple_system_sats_state,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
    )
end

function reset_start_sample_and_bit_buffer(track_state::TrackState)
    TrackState(
        track_state;
        multiple_system_sats_state = reset_start_sample_and_bit_buffer(
            track_state.multiple_system_sats_state,
        ),
    )
end

function has_integration_reached_signal_end_for_all_satellites(
    track_state::TrackState,
    num_samples::Int,
)
    all(track_state.multiple_system_sats_state) do system_sats_state
        all(system_sats_state.states) do sat_state
            sat_state.signal_start_sample == num_samples + 1
        end
    end
end

"""
$(SIGNATURES)

Get the SystemSatsState for a specific GNSS system from a TrackState or
MultipleSystemSatsState by index or symbol.
"""
function get_system_sats_state(
    track_state::TrackState{<:MultipleSystemSatsState{N}},
    system_idx,
) where {N}
    get_system_sats_state(track_state.multiple_system_sats_state, system_idx)
end

function get_system_sats_state(
    multiple_system_sats_state::MultipleSystemSatsState{N},
    system_idx::Union{Symbol,Integer},
) where {N}
    multiple_system_sats_state[system_idx]
end

function get_system_sats_state(
    multiple_system_sats_state::MultipleSystemSatsState{N},
    system_idx::Val{M},
) where {N,M}
    multiple_system_sats_state[M]
end

"""
$(SIGNATURES)

Get the dictionary of satellite states for a specific GNSS system.
"""
function get_sat_states(
    multiple_system_sats_state::MultipleSystemSatsState{N},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_system_sats_state(multiple_system_sats_state, system_idx).states
end

function get_sat_states(multiple_system_sats_state::MultipleSystemSatsState{1})
    get_sat_states(multiple_system_sats_state, 1)
end

function get_sat_states(
    track_state::TrackState{<:MultipleSystemSatsState{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sat_states(track_state.multiple_system_sats_state, system_idx)
end

function get_sat_states(track_state::TrackState{<:MultipleSystemSatsState{1}})
    get_sat_states(track_state.multiple_system_sats_state)
end

function get_system(
    track_state::TrackState{<:MultipleSystemSatsState{N}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_system_sats_state(track_state, system_idx).system
end

function get_system(track_state::TrackState{<:MultipleSystemSatsState{1}})
    get_system(track_state, 1)
end

function get_sat_state(
    track_state::TrackState{<:MultipleSystemSatsState{N}},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
) where {N}
    get_sat_state(get_system_sats_state(track_state, system_idx), sat_identifier)
end

function get_sat_state(
    track_state::TrackState{<:MultipleSystemSatsState{1}},
    sat_identifier,
)
    get_sat_state(track_state, 1, sat_identifier)
end

function get_sat_state(track_state::TrackState{<:MultipleSystemSatsState{1}})
    only(get_sat_states(track_state, 1))
end

function estimate_cn0(
    track_state::TrackState{<:MultipleSystemSatsState},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
)
    estimate_cn0(get_system_sats_state(track_state, system_idx), sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:MultipleSystemSatsState{1}}, sat_identifier)
    estimate_cn0(track_state, 1, sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:MultipleSystemSatsState{1}})
    estimate_cn0(get_system_sats_state(track_state, 1))
end

function merge_sats(
    track_state::TrackState{S,DE},
    system_idx::Union{Symbol,Integer},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{<:Any,<:SatState}},
) where {S<:MultipleSystemSatsState,DE<:AbstractDopplerEstimator}
    TrackState{S,DE}(
        merge_sats(
            track_state.multiple_system_sats_state,
            system_idx,
            to_dictionary(sat_states),
        ),
        merge_sats(track_state.doppler_estimator, system_idx, to_dictionary(sat_states)),
    )
end

function merge_sats(
    track_state::TrackState{<:MultipleSystemSatsState{1}},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{Any,<:SatState}},
)
    merge_sats(track_state, 1, sat_states)
end

function filter_out_sats(
    track_state::TrackState{S,DE},
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {S<:MultipleSystemSatsState,DE<:AbstractDopplerEstimator}
    TrackState{S,DE}(
        filter_out_sats(track_state.multiple_system_sats_state, system_idx, identifiers),
        filter_out_sats(track_state.doppler_estimator, system_idx, identifiers),
    )
end

function filter_out_sats(track_state::TrackState{<:MultipleSystemSatsState{1}}, identifiers)
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