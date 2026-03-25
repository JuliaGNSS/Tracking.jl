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

Multi-threaded CPU downconvert and correlate using Bumper.jl for temporary
code-replica allocation. One `SlabBuffer` is pre-allocated per thread so that
concurrent `@batch` iterations never share a buffer.
"""
struct CPUThreadedDownconvertAndCorrelator{MESF} <: AbstractDownconvertAndCorrelator
    buffers::Vector{SlabBuffer}
end

function CPUThreadedDownconvertAndCorrelator(
    ::Val{MESF};
    nthreads::Int = Threads.maxthreadid(),
) where {MESF}
    buffers = [SlabBuffer() for _ = 1:nthreads]
    CPUThreadedDownconvertAndCorrelator{MESF}(buffers)
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

    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            system = system_sats_state.system
            states = system_sats_state.states
            n = length(states)
            new_vals = Vector{valtype(states)}(undef, n)

            @batch for i = 1:n
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

                buf = dc.buffers[Threads.threadid()]
                @no_escape buf begin
                    code_replica = @alloc(
                        get_code_type(system),
                        num_samples_signal + maximum(sample_shifts) -
                        minimum(sample_shifts)
                    )
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
            end

            new_sat_states = Dictionary(keys(states), new_vals)
            SystemSatsState(system_sats_state, new_sat_states)
        end
    TrackState(track_state; multiple_system_sats_state = new_multiple_system_sats_state)
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
