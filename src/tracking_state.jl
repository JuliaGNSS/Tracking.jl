function TrackState(
    system::AbstractGNSS,
    sat_states;
    num_samples,
    doppler_estimator::AbstractTrackingDopplerEstimator = TrackingConventionalPLLAndDLL(),
    post_process::AbstractTrackingPostProcess = NoTrackingPostProcess(),
)
    multiple_system_sats_state = (SystemSatsState(system, sat_states),)
    TrackState(multiple_system_sats_state; doppler_estimator, post_process, num_samples)
end

function TrackState(
    system_sat_states::SystemSatsState;
    num_samples,
    doppler_estimator::AbstractTrackingDopplerEstimator = TrackingConventionalPLLAndDLL(),
    post_process::AbstractTrackingPostProcess = NoTrackingPostProcess(),
)
    multiple_system_sats_state = (system_sat_states,)
    TrackState(multiple_system_sats_state; doppler_estimator, post_process, num_samples)
end

function TrackState(
    multiple_system_sats_state::MultipleSystemSatsState{N};
    num_samples,
    doppler_estimator::AbstractTrackingDopplerEstimator = TrackingConventionalPLLAndDLL(),
    post_process::AbstractTrackingPostProcess = NoTrackingPostProcess(),
) where {N}
    TrackState(
        initiate_downconvert_and_correlator(multiple_system_sats_state, num_samples),
        doppler_estimator,
        post_process,
        num_samples,
    )
end

function initiate_downconvert_and_correlator(multiple_system_sats_state, num_samples)
    multiple_system_sats_state
end

function initiate_downconvert_and_correlator(
    multiple_system_sats_state::MultipleSystemSatsState{
        N,
        <:AbstractGNSS,
        <:SatState{<:AbstractCorrelator,<:AbstractSatDopplerEstimator,<:Nothing},
        <:AbstractSystemDopplerEstimator,
        <:CPUSystemDownconvertAndCorrelator,
    },
    num_samples::Int,
) where {N}
    return map(multiple_system_sats_state) do system_sats_state
        sat_states = map(system_sats_state.states) do sat_state
            initiate_downconvert_and_correlator(
                system_sats_state.system,
                sat_state,
                num_samples,
            )
        end
        return SystemSatsState(system_sats_state, sat_states)
    end
end

function initiate_downconvert_and_correlator(
    system::AbstractGNSS,
    sat_state::SatState,
    num_samples::Int,
)
    sat_state
end

function initiate_downconvert_and_correlator(
    system::AbstractGNSS,
    sat_state::SatState{<:AbstractCorrelator,<:AbstractSatDopplerEstimator,<:Nothing},
    num_samples::Int,
)
    downconvert_and_correlator =
        CPUSatDownconvertAndCorrelator(system, sat_state.correlator, num_samples)
    return @set sat_state.downconvert_and_correlator = downconvert_and_correlator
end

function TrackState(
    track_state::TrackState,
    multiple_system_sats_state::MultipleSystemSatsState,
)
    TrackState(
        multiple_system_sats_state,
        track_state.doppler_estimator,
        track_state.post_process,
        track_state.num_samples,
    )
end

function TrackState(
    track_state::TrackState,
    multiple_system_sats_state::MultipleSystemSatsState,
    doppler_estimator::AbstractTrackingDopplerEstimator,
)
    TrackState(
        multiple_system_sats_state,
        doppler_estimator,
        track_state.post_process,
        track_state.num_samples,
    )
end

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

function estimate_cn0(
    track_state::TrackState{<:MultipleSystemSatsState{N}},
    system_idx::Union{Symbol,Integer,Val},
    sat_identifier,
) where {N}
    estimate_cn0(get_system_sats_state(track_state, system_idx), sat_identifier)
end

function estimate_cn0(track_state::TrackState{<:MultipleSystemSatsState{1}}, sat_identifier)
    estimate_cn0(track_state, 1, sat_identifier)
end

function merge_sats(
    track_state::TS,
    system_idx::Union{Symbol,Integer},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{Any,<:SatState}},
) where {TS<:TrackState{<:MultipleSystemSatsState{N}} where {N}}
    TS(
        merge_sats(
            track_state.multiple_system_sats_state,
            system_idx,
            to_dictionary(sat_states),
            track_state.num_samples,
        ),
        track_state.doppler_estimator,
        track_state.post_process,
        track_state.num_samples,
    )
end

function merge_sats(
    track_state::TrackState{<:MultipleSystemSatsState{1}},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{Any,<:SatState}},
)
    merge_sats(track_state, 1, sat_states)
end

function filter_out_sats(
    track_state::TS,
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {TS<:TrackState{<:MultipleSystemSatsState{N}} where {N}}
    new_states = map(
        ((id, sat_state),) -> sat_state,
        filter(
            ((id, sat_state),) -> !in(id, identifiers),
            pairs(track_state.multiple_system_sats_state[system_idx].states),
        ),
    )
    new_track_state =
        @set track_state.multiple_system_sats_state[system_idx].states = new_states
    return new_track_state::TS
end

function filter_out_sats(track_state::TrackState{<:MultipleSystemSatsState{1}}, identifiers)
    filter_out_sats(track_state, 1, identifiers)
end
