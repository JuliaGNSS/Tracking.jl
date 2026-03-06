"""
$(SIGNATURES)

CPU-based implementation of downconversion and correlation. The MESF type
parameter specifies the maximum expected sampling frequency.
"""
struct CPUDownconvertAndCorrelator{MESF,B} <: AbstractDownconvertAndCorrelator
    buffer::B
end

function CPUDownconvertAndCorrelator(
    maximum_expected_sampling_frequency::Val{MESF},
    buffer::B = default_buffer(),
) where {MESF,B}
    CPUDownconvertAndCorrelator{MESF,B}(buffer)
end

"""
$(SIGNATURES)

Multi-threaded CPU downconvert and correlate. Takes systems at construction to
pre-allocate correctly-typed code replica buffers per satellite slot. Spawns
one `Threads.@spawn` task per satellite. Does not use Bumper.

The `max_sats` parameter controls how many satellite slots (and buffer sets)
are pre-allocated. Each slot gets its own code replica buffer.
"""
struct CPUThreadedDownconvertAndCorrelator{MESF,CT<:Tuple} <: AbstractDownconvertAndCorrelator
    max_sats::Int
    max_samples::Int
    system_code_types::Dict{UInt64,Int}  # objectid(system) => index into code_replica_buffers tuple
    code_replica_buffers::CT             # Tuple of Vector{Vector{CodeType}} per system
end

function CPUThreadedDownconvertAndCorrelator(
    systems,
    ::Val{MESF};
    max_sats::Int = 32,
    max_sample_shift::Int = 20,
    max_num_samples::Int = ceil(Int, upreferred(MESF / Hz) * 1e-3),
) where {MESF}
    max_samples = max_num_samples
    code_len = max_samples + 2 * max_sample_shift

    system_code_types = Dict{UInt64,Int}(
        objectid(sys) => i for (i, sys) in enumerate(systems)
    )

    code_replica_buffers = Tuple(
        [Vector{get_code_type(sys)}(undef, code_len) for _ in 1:max_sats]
        for sys in systems
    )

    CPUThreadedDownconvertAndCorrelator{MESF,typeof(code_replica_buffers)}(
        max_sats,
        max_samples,
        system_code_types,
        code_replica_buffers,
    )
end

"""
$(SIGNATURES)

Downconvert and correlate all available satellites on the CPU.
"""
function downconvert_and_correlate(
    downconvert_and_correlator::CPUDownconvertAndCorrelator{MESF},
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
) where {MESF}
    num_samples_signal = get_num_samples(signal)
    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            new_sat_states = map(system_sats_state.states) do sat_state
                signal_samples_to_integrate, is_integration_completed =
                    calc_signal_samples_to_integrate(
                        system_sats_state.system,
                        sat_state.signal_start_sample,
                        sampling_frequency,
                        sat_state.code_doppler,
                        sat_state.code_phase,
                        preferred_num_code_blocks_to_integrate,
                        has_bit_or_secondary_code_been_found(sat_state),
                        num_samples_signal,
                    )
                if signal_samples_to_integrate == 0
                    return sat_state
                end
                carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
                code_frequency =
                    sat_state.code_doppler + get_code_frequency(system_sats_state.system)

                sample_shifts = get_correlator_sample_shifts(
                    sat_state.correlator,
                    sampling_frequency,
                    code_frequency,
                )
                @no_escape downconvert_and_correlator.buffer begin
                    code_replica_buffer = @alloc(
                        get_code_type(system_sats_state.system),
                        num_samples_signal + maximum(sample_shifts) -
                        minimum(sample_shifts)
                    )
                    new_correlator = downconvert_and_correlate!(
                        system_sats_state.system,
                        signal,
                        sat_state.correlator,
                        code_replica_buffer,
                        sat_state.code_phase,
                        sat_state.carrier_phase,
                        code_frequency,
                        carrier_frequency,
                        sampling_frequency,
                        sat_state.signal_start_sample,
                        signal_samples_to_integrate,
                        sat_state.prn,
                        Val{MESF}(),
                    )::typeof(sat_state.correlator)
                end
                return update(
                    system_sats_state.system,
                    sat_state,
                    signal_samples_to_integrate,
                    intermediate_frequency,
                    sampling_frequency,
                    new_correlator,
                    is_integration_completed,
                )
            end
            return SystemSatsState(system_sats_state, new_sat_states)
        end
    return TrackState(
        track_state;
        multiple_system_sats_state = new_multiple_system_sats_state,
    )
end

"""
$(SIGNATURES)

Multi-threaded downconvert and correlate. Spawns one task per satellite across
all systems, using pre-allocated per-satellite code replica buffers and the
fused downconvert+correlate kernel.
"""
function downconvert_and_correlate(
    dc::CPUThreadedDownconvertAndCorrelator{MESF},
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
) where {MESF}
    num_samples_signal = get_num_samples(signal)

    buf_offset = 0
    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            system = system_sats_state.system
            sys_idx = dc.system_code_types[objectid(system)]
            states = system_sats_state.states
            n = length(states)
            new_vals = Vector{valtype(states)}(undef, n)

            Threads.@threads for i in 1:n
                buf_idx = buf_offset + i
                sat_state = states.values[i]

                signal_samples_to_integrate, is_integration_completed =
                    calc_signal_samples_to_integrate(
                        system,
                        sat_state.signal_start_sample,
                        sampling_frequency,
                        sat_state.code_doppler,
                        sat_state.code_phase,
                        preferred_num_code_blocks_to_integrate,
                        has_bit_or_secondary_code_been_found(sat_state),
                        num_samples_signal,
                    )

                if signal_samples_to_integrate == 0
                    new_vals[i] = sat_state
                    continue
                end

                carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
                code_frequency = sat_state.code_doppler + get_code_frequency(system)

                sample_shifts = get_correlator_sample_shifts(
                    sat_state.correlator,
                    sampling_frequency,
                    code_frequency,
                )

                code_replica = dc.code_replica_buffers[sys_idx][buf_idx]
                gen_code_replica!(
                    code_replica,
                    system,
                    code_frequency,
                    sampling_frequency,
                    sat_state.code_phase,
                    sat_state.signal_start_sample,
                    signal_samples_to_integrate,
                    sample_shifts,
                    sat_state.prn,
                    Val{MESF}(),
                )

                new_correlator = downconvert_and_correlate_fused!(
                    sat_state.correlator,
                    signal,
                    code_replica,
                    sample_shifts,
                    carrier_frequency,
                    sampling_frequency,
                    sat_state.carrier_phase,
                    sat_state.signal_start_sample,
                    signal_samples_to_integrate,
                )::typeof(sat_state.correlator)

                new_vals[i] = update(
                    system,
                    sat_state,
                    signal_samples_to_integrate,
                    intermediate_frequency,
                    sampling_frequency,
                    new_correlator,
                    is_integration_completed,
                )
            end

            buf_offset += n
            new_sat_states = Dictionary(keys(states), new_vals)
            SystemSatsState(system_sats_state, new_sat_states)
        end
    TrackState(
        track_state;
        multiple_system_sats_state = new_multiple_system_sats_state,
    )
end

"""
$(SIGNATURES)

Downconvert and correlate a single satellite on the CPU.
"""
function downconvert_and_correlate!(
    system,
    signal,
    correlator::AbstractCorrelator{M},
    code_replica,
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn,
    maximum_expected_sampling_frequency,
) where {M}
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    gen_code_replica!(
        code_replica,
        system,
        code_frequency,
        sampling_frequency,
        code_phase,
        signal_start_sample,
        num_samples_left,
        sample_shifts,
        prn,
        maximum_expected_sampling_frequency,
    )
    downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        signal_start_sample,
        num_samples_left,
    )
end
