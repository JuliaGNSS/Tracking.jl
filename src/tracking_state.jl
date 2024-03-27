struct TrackState{
    S<:TupleLike{<:NTuple{N,SystemSatsState}} where {N},
    DE<:AbstractDopplerEstimator,
    DC<:AbstractDownconvertAndCorrelator,
}
    system_sats_states::S
    doppler_estimator::DE
    downconvert_and_correlator::DC
    num_samples::Int
end

function TrackState(
    system::AbstractGNSS,
    sat_states::Vector{<:SatState};
    num_samples,
    doppler_estimator = map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states),
    downconvert_and_correlator = CPUDownconvertAndCorrelator(
        (SystemSatsState(system, sat_states),),
        num_samples,
    ),
)
    system_sats_stats = (SystemSatsState(system, sat_states),)
    TrackState(
        system_sats_stats,
        ConventionalPLLsAndDLLs((doppler_estimator,)),
        downconvert_and_correlator,
        num_samples,
    )
end

function TrackState(
    system_sats_states::TupleLike{<:NTuple{N,SystemSatsState}};
    num_samples,
    doppler_estimators = ConventionalPLLsAndDLLs(
        map(
            sat_states ->
                map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states.states),
            system_sats_states,
        ),
    ),
    downconvert_and_correlator = CPUDownconvertAndCorrelator(
        system_sats_states,
        num_samples,
    ),
) where {N}
    TrackState(
        system_sats_states,
        doppler_estimators,
        downconvert_and_correlator,
        num_samples,
    )
end

function get_sats_states(
    track_state::TrackState{<:TupleLike{<:NTuple{N,SystemSatsState}}},
    system_idx::Union{Symbol,Integer},
) where {N}
    track_state.system_sats_states[system_idx]
end

function get_sats_states(
    track_state::TrackState{<:TupleLike{<:NTuple{N,SystemSatsState}}},
    system_idx::Val{M},
) where {N,M}
    track_state.system_sats_states[M]
end

function get_sat_states(
    track_state::TrackState{<:TupleLike{<:NTuple{N,SystemSatsState}}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sats_states(track_state, system_idx).states
end

function get_sat_states(track_state::TrackState{<:TupleLike{<:NTuple{1,SystemSatsState}}})
    get_sat_states(track_state, 1)
end

function get_system(
    track_state::TrackState{<:TupleLike{<:NTuple{N,SystemSatsState}}},
    system_idx::Union{Symbol,Integer,Val},
) where {N}
    get_sats_states(track_state, system_idx).system
end

function get_system(track_state::TrackState{<:TupleLike{<:NTuple{1,SystemSatsState}}})
    get_system(track_state, 1)
end

function get_sat_state(
    track_state::TrackState{<:TupleLike{<:NTuple{N,SystemSatsState}}},
    system_idx::Union{Symbol,Integer,Val},
    sat_idx::Integer,
) where {N}
    get_sat_state(get_sats_states(track_state, system_idx), sat_idx)
end

function get_sat_state(
    track_state::TrackState{<:TupleLike{<:NTuple{1,SystemSatsState}}},
    sat_idx::Integer,
)
    get_sat_state(track_state, 1, sat_idx)
end

function estimate_cn0(
    track_state::TrackState{<:TupleLike{<:NTuple{N,SystemSatsState}}},
    system_idx::Union{Symbol,Integer,Val},
    sat_idx::Integer,
) where {N}
    estimate_cn0(get_sats_states(track_state, system_idx), sat_idx)
end

function estimate_cn0(
    track_state::TrackState{<:TupleLike{<:NTuple{1,SystemSatsState}}},
    sat_idx::Integer,
)
    estimate_cn0(track_state, 1, sat_idx)
end

create_buffer(buffers::Vector{<:CPUBuffers}, system, sat_states, num_samples) =
    CPUBuffers(system, sat_states, num_samples)

create_buffer(buffers::Vector{<:GPUBuffers}, system, sat_states, num_samples) =
    GPUBuffers(system, sat_states, num_samples)

function add_sats!(
    track_state::TrackState{
        <:TupleLike{<:NTuple{N,SystemSatsState}},
        <:ConventionalPLLsAndDLLs,
    },
    system_idx::Union{Symbol,Integer},
    system::AbstractGNSS,
    sat_states::Union{SatState,Vector{<:SatState}},
) where {N}
    system_sats_state = track_state.system_sats_states[system_idx]
    pll_and_dlls = track_state.doppler_estimator.plls_and_dlls[system_idx]
    buffers = track_state.downconvert_and_correlator.buffers[system_idx]
    if sat_states isa Vector
        append!(system_sats_state.states, sat_states)
        append!(
            pll_and_dlls,
            map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states),
        )
        append!(
            buffers,
            map(
                sat_state ->
                    create_buffer(buffers, system, sat_state, track_state.num_samples),
                sat_states,
            ),
        )
    else
        push!(system_sats_state.states, sat_states)
        push!(pll_and_dlls, ConventionalPLLAndDLL(sat_states))
        push!(buffers, create_buffer(buffers, system, sat_states, track_state.num_samples))
    end
    track_state
end

function add_sats!(
    track_state::TrackState{
        <:TupleLike{<:NTuple{1,SystemSatsState}},
        <:ConventionalPLLsAndDLLs,
    },
    system::AbstractGNSS,
    sat_states::Union{SatState,Vector{<:SatState}},
)
    add_sats!(track_state, 1, system, sat_states)
end

function remove_sats!(
    track_state::TrackState{
        <:TupleLike{<:NTuple{N,SystemSatsState}},
        <:ConventionalPLLsAndDLLs,
    },
    system_idx::Union{Symbol,Integer},
    prns::Union{Int,Vector{Int}},
) where {N}
    system_sats_state = track_state.system_sats_states[system_idx]
    pll_and_dlls = track_state.doppler_estimator.plls_and_dlls[system_idx]
    buffers = track_state.downconvert_and_correlator.buffers[system_idx]
    idxs_to_delete = findall(state -> state.prn in prns, system_sats_state.states)
    deleteat!(system_sats_state.states, idxs_to_delete)
    deleteat!(pll_and_dlls, idxs_to_delete)
    deleteat!(buffers, idxs_to_delete)
    track_state
end

function remove_sats!(
    track_state::TrackState{
        <:TupleLike{<:NTuple{1,SystemSatsState}},
        <:ConventionalPLLsAndDLLs,
    },
    prns::Union{Int,Vector{Int}},
)
    remove_sats!(track_state, 1, prns)
end
