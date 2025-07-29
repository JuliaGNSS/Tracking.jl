function TrackState(
    system::AbstractGNSS,
    sat_states;
    maximum_expected_sampling_frequency::Maybe{Val} = nothing,
    num_samples,
    doppler_estimator::AbstractDopplerEstimator = ConventionalPLLAndDLL((
        SystemSatsState(system, sat_states),
    )),
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(
        maximum_expected_sampling_frequency,
        (SystemSatsState(system, sat_states),),
        num_samples,
    ),
)
    multiple_system_sats_state = (SystemSatsState(system, sat_states),)
    TrackState(
        multiple_system_sats_state,
        doppler_estimator,
        downconvert_and_correlator,
        num_samples,
    )
end

function TrackState(
    system_sat_states::SystemSatsState;
    maximum_expected_sampling_frequency::Maybe{Val} = nothing,
    num_samples,
    doppler_estimator::AbstractDopplerEstimator = ConventionalPLLAndDLL((
        system_sat_states,
    )),
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(
        maximum_expected_sampling_frequency,
        (system_sat_states,),
        num_samples,
    ),
)
    multiple_system_sats_state = (system_sat_states,)
    TrackState(
        multiple_system_sats_state,
        doppler_estimator,
        downconvert_and_correlator,
        num_samples,
    )
end

function TrackState(
    multiple_system_sats_state::MultipleSystemSatsState;
    maximum_expected_sampling_frequency::Maybe{Val} = nothing,
    num_samples,
    doppler_estimator::AbstractDopplerEstimator = ConventionalPLLAndDLL(
        multiple_system_sats_state,
    ),
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUDownconvertAndCorrelator(
        maximum_expected_sampling_frequency,
        multiple_system_sats_state,
        num_samples,
    ),
)
    TrackState(
        multiple_system_sats_state,
        doppler_estimator,
        downconvert_and_correlator,
        num_samples,
    )
end

function TrackState(
    track_state::TrackState{S,DE,DC},
    multiple_system_sats_state::S,
    doppler_estimator::DE,
) where {
    N,
    I,
    S<:MultipleSystemSatsState{N,I},
    DE<:AbstractDopplerEstimator{N,I},
    DC<:AbstractDownconvertAndCorrelator{N,I},
}
    TrackState{S,DE,DC}(
        multiple_system_sats_state,
        doppler_estimator,
        track_state.downconvert_and_correlator,
        track_state.num_samples,
    )
end

# Be careful when calling this.
# It might lead to types that are inferred at runtime?!
# Tested with 1.11.6
function TrackState(
    track_state::TrackState{S,DE,DC};
    multiple_system_sats_state::Maybe{S} = nothing,
    doppler_estimator::Maybe{DE} = nothing,
    downconvert_and_correlator::Maybe{DC} = nothing,
) where {
    N,
    I,
    S<:MultipleSystemSatsState{N,I},
    DE<:AbstractDopplerEstimator{N,I},
    DC<:AbstractDownconvertAndCorrelator{N,I},
}
    TrackState{S,DE,DC}(
        isnothing(multiple_system_sats_state) ? track_state.multiple_system_sats_state :
        multiple_system_sats_state,
        isnothing(doppler_estimator) ? track_state.doppler_estimator : doppler_estimator,
        isnothing(downconvert_and_correlator) ? track_state.downconvert_and_correlator :
        downconvert_and_correlator,
        track_state.num_samples,
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
    track_state::TrackState{S,DE,DC},
    system_idx::Union{Symbol,Integer},
    sat_states::Union{SatState,Vector{<:SatState},Dictionary{Any,<:SatState}},
) where {
    S<:MultipleSystemSatsState,
    DE<:AbstractDopplerEstimator,
    DC<:AbstractDownconvertAndCorrelator,
}
    TrackState{S,DE,DC}(
        merge_sats(
            track_state.multiple_system_sats_state,
            system_idx,
            to_dictionary(sat_states),
        ),
        merge_sats(track_state.doppler_estimator, system_idx, to_dictionary(sat_states)),
        merge_sats(
            get_system(track_state, system_idx),
            track_state.downconvert_and_correlator,
            system_idx,
            to_dictionary(sat_states),
            track_state.num_samples,
        ),
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
    track_state::TrackState{S,DE,DC},
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {
    S<:MultipleSystemSatsState,
    DE<:AbstractDopplerEstimator,
    DC<:AbstractDownconvertAndCorrelator,
}
    TrackState{S,DE,DC}(
        filter_out_sats(track_state.multiple_system_sats_state, system_idx, identifiers),
        filter_out_sats(track_state.doppler_estimator, system_idx, identifiers),
        filter_out_sats(track_state.downconvert_and_correlator, system_idx, identifiers),
        track_state.num_samples,
    )
end

function filter_out_sats(track_state::TrackState{<:MultipleSystemSatsState{1}}, identifiers)
    filter_out_sats(track_state, 1, identifiers)
end
