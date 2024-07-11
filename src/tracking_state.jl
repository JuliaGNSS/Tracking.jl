function TrackState(
    system::AbstractGNSS,
    sat_states;
    num_samples,
    doppler_estimator::AbstractTrackingDopplerEstimator = TrackingConventionalPLLAndDLL(),
    post_process::AbstractTrackingPostProcess = NoTrackingPostProcess(),
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(),
)
    multiple_system_sats_state = (SystemSatsState(system, sat_states),)
    TrackState(
        multiple_system_sats_state;
        doppler_estimator,
        post_process,
        downconvert_and_correlator,
        num_samples,
    )
end

function TrackState(
    system_sat_states::SystemSatsState;
    num_samples,
    doppler_estimator::AbstractTrackingDopplerEstimator = TrackingConventionalPLLAndDLL(),
    post_process::AbstractTrackingPostProcess = NoTrackingPostProcess(),
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(),
)
    multiple_system_sats_state = (system_sat_states,)
    TrackState(
        multiple_system_sats_state;
        doppler_estimator,
        post_process,
        downconvert_and_correlator,
        num_samples,
    )
end

function TrackState(
    multiple_system_sats_state::MultipleSystemSatsState{N};
    num_samples,
    doppler_estimator::AbstractTrackingDopplerEstimator = TrackingConventionalPLLAndDLL(),
    post_process::AbstractTrackingPostProcess = NoTrackingPostProcess(),
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(),
) where {N}
    TrackState(
        initiate_downconvert_and_correlator(
            downconvert_and_correlator,
            multiple_system_sats_state,
            num_samples,
        ),
        doppler_estimator,
        downconvert_and_correlator,
        post_process,
        num_samples,
    )
end

function initiate_downconvert_and_correlator(
    downconvert_and_correlator,
    multiple_system_sats_state,
    num_samples,
)
    multiple_system_sats_state
end

function initiate_downconvert_and_correlator(
    downconvert_and_correlator::CPUDownconvertAndCorrelator,
    multiple_system_sats_state::MultipleSystemSatsState{
        N,
        <:AbstractGNSS,
        <:SatState{<:AbstractCorrelator,<:AbstractSatDopplerEstimator,<:Nothing},
        <:AbstractSystemDopplerEstimator,
        <:Nothing,
    },
    num_samples::Int,
) where {N}
    return map(multiple_system_sats_state) do system_sats_state
        sat_states = map(system_sats_state.states) do sat_state
            sat_downconvert_and_correlator = CPUSatDownconvertAndCorrelator(
                system_sats_state.system,
                sat_state.correlator,
                num_samples,
            )
            return @set sat_state.downconvert_and_correlator =
                sat_downconvert_and_correlator
        end
        return SystemSatsState(
            system_sats_state,
            sat_states;
            downconvert_and_correlator = CPUSystemDownconvertAndCorrelator(),
        )
    end
end

function initiate_downconvert_and_correlator(
    downconvert_and_correlator::GPUDownconvertAndCorrelator,
    multiple_system_sats_state::MultipleSystemSatsState{
        N,
        <:AbstractGNSS,
        <:SatState{<:AbstractCorrelator,<:AbstractSatDopplerEstimator,<:Nothing},
        <:AbstractSystemDopplerEstimator,
        <:Nothing,
    },
    num_samples::Int,
) where {N}
    return map(multiple_system_sats_state) do system_sats_state
        sat_states = map(system_sats_state.states) do sat_state
            sat_downconvert_and_correlator = GPUSatDownconvertAndCorrelator(
                system_sats_state.system,
                sat_state.correlator,
                num_samples,
            )
            return @set sat_state.downconvert_and_correlator =
                sat_downconvert_and_correlator
        end
        return SystemSatsState(
            system_sats_state,
            sat_states;
            downconvert_and_correlator = GPUSystemDownconvertAndCorrelator(
                system_sats_state.system,
            ),
        )
    end
end

function TrackState(
    track_state::TrackState,
    multiple_system_sats_state::MultipleSystemSatsState,
)
    TrackState(
        multiple_system_sats_state,
        track_state.doppler_estimator,
        track_state.downconvert_and_correlator,
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
        track_state.downconvert_and_correlator,
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
            instantiate_downconvert_and_correlator(
                track_state,
                to_dictionary(sat_states),
                track_state.multiple_system_sats_state[system_idx].system,
            ),
            track_state.num_samples,
        ),
        track_state.doppler_estimator,
        track_state.downconvert_and_correlator,
        track_state.post_process,
        track_state.num_samples,
    )
end

function instantiate_downconvert_and_correlator(track_state, sat_states, system)
    sat_states
end

function instantiate_downconvert_and_correlator(
    track_state::TrackState{
        <:MultipleSystemSatsState,
        <:AbstractTrackingDopplerEstimator,
        <:CPUDownconvertAndCorrelator,
    },
    sat_states::Dictionary{
        I,
        <:SatState{<:AbstractCorrelator,<:Maybe{<:AbstractSatDopplerEstimator},Nothing},
    },
    system::AbstractGNSS,
) where {I}
    map(sat_states) do sat_state
        sat_downconvert_and_correlator = CPUSatDownconvertAndCorrelator(
            system,
            sat_state.correlator,
            track_state.num_samples,
        )
        return @set sat_state.downconvert_and_correlator = sat_downconvert_and_correlator
    end
end

function instantiate_downconvert_and_correlator(
    track_state::TrackState{
        <:MultipleSystemSatsState,
        <:AbstractTrackingDopplerEstimator,
        <:GPUDownconvertAndCorrelator,
    },
    sat_states::Dictionary{
        I,
        <:SatState{<:AbstractCorrelator,<:Maybe{<:AbstractSatDopplerEstimator},Nothing},
    },
    system::AbstractGNSS,
) where {I}
    map(sat_states) do sat_state
        sat_downconvert_and_correlator = GPUSatDownconvertAndCorrelator(
            system,
            sat_state.correlator,
            track_state.num_samples,
        )
        return @set sat_state.downconvert_and_correlator = sat_downconvert_and_correlator
    end
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
