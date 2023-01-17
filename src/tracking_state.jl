struct TrackState{
    S <: NTuple{N, SystemSatsState} where N,
    DE <: AbstractDopplerEstimator,
    DC <: AbstractDownconvertAndCorrelator,
}
    system_sats_states::S
    doppler_estimator::DE
    downconvert_and_correlator::DC
    num_samples::Int
end

function TrackState(
    system::AbstractGNSS,
    sat_states::Vector{<: SatState};
    num_samples,
    doppler_estimator = map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states),
    downconvert_and_correlator = CPUDownconvertAndCorrelator((SystemSatsState(system, sat_states),), num_samples)
)
    system_sats_stats = (SystemSatsState(system, sat_states),)
    TrackState(
        system_sats_stats,
        ConventionalPLLsAndDLLs((doppler_estimator,)),
        downconvert_and_correlator,
        num_samples
    )
end

function TrackState(
    system_sats_states;
    num_samples,
    doppler_estimators = ConventionalPLLsAndDLLs(map(sat_states -> map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states.states), system_sats_states)),
    downconvert_and_correlator = CPUDownconvertAndCorrelator(system_sats_states, num_samples), 
)
    TrackState(
        system_sats_states,
        doppler_estimators,
        downconvert_and_correlator,
        num_samples
    )
end

function add_sats!(
    track_state::TrackState{S, <: ConventionalPLLsAndDLLs},
    system::AbstractGNSS,
    sat_states::Union{SatState, Vector{<:SatState}},
) where S
    foreach(track_state.system_sats_states, track_state.doppler_estimator.plls_and_dlls, track_state.downconvert_and_correlator.buffers) do system_sats_state, pll_and_dlls, buffers
        if typeof(system) == typeof(system_sats_state.system)
            if sat_states isa Vector
                append!(system_sats_state.states, sat_states)
                append!(pll_and_dlls, map(sat_state -> ConventionalPLLAndDLL(sat_state), sat_states))
                append!(buffers, map(sat_state -> CPUBuffers(system, sat_state, track_state.num_samples), sat_states))
            else
                push!(system_sats_state.states, sat_states)
                push!(pll_and_dlls, ConventionalPLLAndDLL(sat_states))
                push!(buffers, CPUBuffers(system, sat_states, track_state.num_samples))
            end
        end
    end
    track_state
end

function remove_sats!(
    track_state::TrackState{S, <: ConventionalPLLsAndDLLs},
    system::AbstractGNSS,
    prns::Union{Int, Vector{Int}},
) where S
    foreach(track_state.system_sats_states, track_state.doppler_estimator.plls_and_dlls, track_state.downconvert_and_correlator.buffers) do system_sats_state, pll_and_dlls, buffers
        if typeof(system) == typeof(system_sats_state.system)
            idxs_to_delete = findall(state -> state.prn in prns, system_sats_state.states)
            deleteat!(system_sats_state.states, idxs_to_delete)
            deleteat!(pll_and_dlls, idxs_to_delete)
            deleteat!(buffers, idxs_to_delete)
        end
    end
    track_state
end