"""
$(SIGNATURES)

CPU-based implementation of downconversion and correlation.
"""
struct CPUDownconvertAndCorrelator{B} <: AbstractDownconvertAndCorrelator
    buffer::B
end

CPUDownconvertAndCorrelator() = CPUDownconvertAndCorrelator(default_buffer())

"""
$(SIGNATURES)

Multi-threaded CPU downconvert and correlate using Bumper.jl for temporary
code-replica allocation. One `SlabBuffer` is pre-allocated per thread so that
concurrent `@batch` iterations never share a buffer.
"""
struct CPUThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{SlabBuffer}
end

function CPUThreadedDownconvertAndCorrelator(;
    nthreads::Int = Threads.maxthreadid(),
)
    buffers = [SlabBuffer() for _ = 1:nthreads]
    CPUThreadedDownconvertAndCorrelator(buffers)
end

# Per-sat downconvert+correlate for the single-threaded CPU backend. Pure:
# returns the updated TrackedSat. Shared by `downconvert_and_correlate` and
# `downconvert_and_correlate!` for `CPUDownconvertAndCorrelator`.
function _update_tracked_sat_correlator(
    tracked_sat::TrackedSat,
    dc::CPUDownconvertAndCorrelator,
    signal,
    system,
    num_samples_signal,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
    intermediate_frequency,
)
    sat_state = tracked_sat.sat_state
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
        return tracked_sat
    end
    carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
    code_frequency = sat_state.code_doppler + get_code_frequency(system)
    sample_shifts = get_correlator_sample_shifts(
        sat_state.correlator,
        sampling_frequency,
        code_frequency,
    )
    new_correlator = @no_escape dc.buffer begin
        code_replica_buffer = @alloc(
            get_code_type(system),
            num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
        )
        downconvert_and_correlate!(
            system,
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
        )::typeof(sat_state.correlator)
    end
    new_sat_state = update(
        system,
        sat_state,
        signal_samples_to_integrate,
        intermediate_frequency,
        sampling_frequency,
        new_correlator,
        is_integration_completed,
    )
    TrackedSat(new_sat_state, tracked_sat.estimator_state)
end

"""
$(SIGNATURES)

Downconvert and correlate all available satellites on the CPU.
"""
function downconvert_and_correlate(
    dc::CPUDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            new_tracked_sats = map(system_sats_state.states) do tracked_sat
                _update_tracked_sat_correlator(
                    tracked_sat,
                    dc,
                    signal,
                    system_sats_state.system,
                    num_samples_signal,
                    preferred_num_code_blocks_to_integrate,
                    sampling_frequency,
                    intermediate_frequency,
                )
            end
            return SystemSatsState(system_sats_state, new_tracked_sats)
        end
    return TrackState(
        track_state;
        multiple_system_sats_state = new_multiple_system_sats_state,
    )
end

"""
$(SIGNATURES)

In-place version: walks each system's `Vector{TrackedSat}` and overwrites
slots in place. Returns the same `track_state`. Allocation-free in steady
state — see [`track!`](@ref).
"""
function downconvert_and_correlate!(
    dc::CPUDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    for system_sats_state in track_state.multiple_system_sats_state
        system = system_sats_state.system
        vals = system_sats_state.states.values
        @inbounds for i in eachindex(vals)
            vals[i] = _update_tracked_sat_correlator(
                vals[i],
                dc,
                signal,
                system,
                num_samples_signal,
                preferred_num_code_blocks_to_integrate,
                sampling_frequency,
                intermediate_frequency,
            )
        end
    end
    return track_state
end

# Per-sat update for the threaded CPU backend (fused kernel + per-thread
# Bumper buffer). Called from a `@batch` loop body; each invocation reads
# its thread's buffer slot independently. Pure: returns the new TrackedSat.
function _update_tracked_sat_correlator_threaded(
    tracked_sat::TrackedSat,
    dc::CPUThreadedDownconvertAndCorrelator,
    signal,
    system,
    num_samples_signal,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
    intermediate_frequency,
)
    sat_state = tracked_sat.sat_state
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
        return tracked_sat
    end
    carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
    code_frequency = sat_state.code_doppler + get_code_frequency(system)
    sample_shifts = get_correlator_sample_shifts(
        sat_state.correlator,
        sampling_frequency,
        code_frequency,
    )
    buf = dc.buffers[Threads.threadid()]
    new_tracked_sat = @no_escape buf begin
        code_replica = @alloc(
            get_code_type(system),
            num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
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
        new_sat_state = update(
            system,
            sat_state,
            signal_samples_to_integrate,
            intermediate_frequency,
            sampling_frequency,
            new_correlator,
            is_integration_completed,
        )
        TrackedSat(new_sat_state, tracked_sat.estimator_state)
    end
    new_tracked_sat
end

"""
$(SIGNATURES)

Multi-threaded downconvert and correlate. Spawns one task per satellite across
all systems, using pre-allocated per-satellite code replica buffers and the
fused downconvert+correlate kernel.
"""
function downconvert_and_correlate(
    dc::CPUThreadedDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            system = system_sats_state.system
            states = system_sats_state.states
            n = length(states)
            new_vals = Vector{valtype(states)}(undef, n)
            @batch for i = 1:n
                @inbounds new_vals[i] = _update_tracked_sat_correlator_threaded(
                    states.values[i],
                    dc,
                    signal,
                    system,
                    num_samples_signal,
                    preferred_num_code_blocks_to_integrate,
                    sampling_frequency,
                    intermediate_frequency,
                )
            end
            new_sat_states = Dictionary(keys(states), new_vals)
            SystemSatsState(system_sats_state, new_sat_states)
        end
    TrackState(track_state; multiple_system_sats_state = new_multiple_system_sats_state)
end

"""
$(SIGNATURES)

In-place version: writes new `TrackedSat` values directly into each system's
existing `Vector{TrackedSat}` backing storage. Different `@batch` iterations
write to disjoint slots, so no synchronization is needed. Returns the same
`track_state`. Allocation-free in steady state — see [`track!`](@ref).
"""
function downconvert_and_correlate!(
    dc::CPUThreadedDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    for system_sats_state in track_state.multiple_system_sats_state
        system = system_sats_state.system
        vals = system_sats_state.states.values
        n = length(vals)
        @batch for i = 1:n
            @inbounds vals[i] = _update_tracked_sat_correlator_threaded(
                vals[i],
                dc,
                signal,
                system,
                num_samples_signal,
                preferred_num_code_blocks_to_integrate,
                sampling_frequency,
                intermediate_frequency,
            )
        end
    end
    return track_state
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
