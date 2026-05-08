# Three independent scratch byte buffers — one per role. The code-replica
# path uses `code_replica`; the fused kernel's `AbstractVector`-shifts
# fallback uses `tile_re` / `tile_im` concurrently. Holding separate
# buffers per role (instead of bump-arena offset math on one) keeps the
# call sites simple. Both CPU backends share this struct: the
# single-threaded one holds a single `ScratchBuffers`, the threaded one
# holds a `Vector{ScratchBuffers}` (one per thread).
struct ScratchBuffers
    code_replica::Vector{UInt8}
    tile_re::Vector{UInt8}
    tile_im::Vector{UInt8}
end

ScratchBuffers() = ScratchBuffers(UInt8[], UInt8[], UInt8[])

# Bitstype pointer-and-length view over a `Vector{UInt8}` slot, used as
# the typed handle the kernels accept. Implements just enough of the
# `AbstractArray` interface for the call sites: linear indexing,
# `length`/`size`, and `pointer` (so SIMD `vstore` works). The struct is
# fully isbits, so building one in `_with_scratch_buffer` is free.
struct ScratchView{T} <: DenseVector{T}
    ptr::Ptr{T}
    len::Int
end
@inline Base.size(v::ScratchView) = (v.len,)
@inline Base.IndexStyle(::Type{<:ScratchView}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::ScratchView, i::Int)
    @boundscheck checkbounds(v, i)
    unsafe_load(v.ptr, i)
end
Base.@propagate_inbounds function Base.setindex!(v::ScratchView, x, i::Int)
    @boundscheck checkbounds(v, i)
    unsafe_store!(v.ptr, x, i)
    v
end
Base.unsafe_convert(::Type{Ptr{T}}, v::ScratchView{T}) where {T} = v.ptr
Base.pointer(v::ScratchView) = v.ptr
Base.elsize(::Type{ScratchView{T}}) where {T} = sizeof(T)

"""
$(SIGNATURES)

CPU-based implementation of downconversion and correlation. Holds one
`ScratchBuffers` — three long-lived `Vector{UInt8}` byte
buffers, one per scratch role (code replica + the fused kernel's two
tile halves). Buffers grow lazily on first use and are reused
thereafter, so a hoisted instance has zero allocations per `track!`
call in steady state.

For real-time loops, construct the correlator **once outside** the
`track!` loop and pass it via the `downconvert_and_correlator` keyword
argument — the default value rebuilds the buffers on every call.
"""
struct CPUDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::ScratchBuffers
end

CPUDownconvertAndCorrelator() = CPUDownconvertAndCorrelator(ScratchBuffers())

"""
$(SIGNATURES)

Multi-threaded CPU downconvert and correlate. Holds one
`ScratchBuffers` per thread, indexed by `Threads.threadid()`
inside `@batch` (which pins each iteration to a fixed thread). Buffers
grow lazily on first use and are reused thereafter, so a hoisted
instance has near-zero allocations per `track!` call in steady state
(Polyester's `@batch` keeps a small irreducible per-call closure
allocation).

For real-time loops, construct the correlator **once outside** the
`track!` loop and pass it via the `downconvert_and_correlator` keyword
argument — the default value rebuilds the per-thread buffers on every
call.
"""
struct CPUThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{ScratchBuffers}
end

CPUThreadedDownconvertAndCorrelator() = CPUThreadedDownconvertAndCorrelator(
    [ScratchBuffers() for _ in 1:Threads.maxthreadid()],
)

# Look up the active `ScratchBuffers` for this thread. Single-threaded
# backend has just one; threaded backend indexes by `Threads.threadid()`
# (stable under Polyester `@batch`).
@inline _scratch_buffers(dc::CPUDownconvertAndCorrelator) = dc.buffers
@inline _scratch_buffers(dc::CPUThreadedDownconvertAndCorrelator) =
    @inbounds dc.buffers[Threads.threadid()]

# Grow the role's byte buffer to fit `n` elements of `T` and hand a
# typed `ScratchView` of exactly that size to `f`. The view is bitstype,
# so building one is free; the underlying `Vector{UInt8}` is reused
# across calls — once the buffer has reached its working size every
# subsequent call is allocation-free. `GC.@preserve` keeps the byte
# vector rooted while `f` runs (the view holds a raw `Ptr`, untracked
# by GC).
@inline function _with_scratch_view(f, buf::Vector{UInt8}, ::Type{T}, n::Int) where {T}
    nbytes = n * sizeof(T)
    length(buf) < nbytes && resize!(buf, nbytes)
    GC.@preserve buf begin
        f(ScratchView{T}(Ptr{T}(pointer(buf)), n))
    end
end

# Convenience wrappers per role.
@inline function _with_code_replica_buffer(
    f,
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    ::Type{T},
    n::Int,
) where {T}
    _with_scratch_view(f, _scratch_buffers(dc).code_replica, T, n)
end

# Per-backend correlation kernel. Single-threaded backend uses the split
# downconvert→correlate kernel; threaded backend pre-generates the code
# replica and uses the fused kernel. Both write into the supplied
# `code_replica` buffer (allocated from the slab) and return a value of
# the same correlator type as `sat_state.correlator`.
@inline function _correlate_with_buffer!(
    ::CPUDownconvertAndCorrelator,
    code_replica,
    system,
    sat_state,
    signal,
    sample_shifts,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_samples_to_integrate,
)
    downconvert_and_correlate!(
        system,
        signal,
        sat_state.correlator,
        code_replica,
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

@inline function _correlate_with_buffer!(
    dc::CPUThreadedDownconvertAndCorrelator,
    code_replica,
    system,
    sat_state,
    signal,
    sample_shifts,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_samples_to_integrate,
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
    _fused_with_tile_scratch!(
        dc,
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
end

# Dispatch helper for the fused kernel that hands SoA tile buffers from
# the calling correlator's `ScratchBuffers` to the kernel — when the
# kernel actually needs them. The `@generated` overload of
# `downconvert_and_correlate_fused!` (for `SVector{NC}` shifts) ignores
# tiles entirely; only the `AbstractVector`-shifts fallback uses them.
# Picking the right path here keeps the fused kernel's existing
# dispatch shape and lets the SoA-tile path be allocation-free.
@inline function _fused_with_tile_scratch!(
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    correlator::AbstractCorrelator{M},
    signal,
    code_replica,
    sample_shifts::SVector,
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample,
    num_samples,
) where {M}
    # Static-shifts overload: no tile buffers needed.
    downconvert_and_correlate_fused!(
        correlator, signal, code_replica, sample_shifts,
        carrier_frequency, sampling_frequency, carrier_phase,
        start_sample, num_samples,
    )::typeof(correlator)
end

@inline function _fused_with_tile_scratch!(
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    correlator::AbstractCorrelator{M},
    signal,
    code_replica,
    sample_shifts::AbstractVector,
    carrier_frequency,
    sampling_frequency,
    carrier_phase,
    start_sample,
    num_samples,
) where {M}
    # Dynamic-shifts fallback: pull SoA tile buffers from the calling
    # correlator's per-thread `ScratchBuffers` (one re-half + one
    # im-half) so the fused kernel doesn't allocate them per call.
    bufs = _scratch_buffers(dc)
    n = num_samples * M
    _with_scratch_view(bufs.tile_re, Float32, n) do tile_re
        _with_scratch_view(bufs.tile_im, Float32, n) do tile_im
            _downconvert_and_correlate_fused_with_tiles!(
                correlator, signal, code_replica, sample_shifts,
                carrier_frequency, sampling_frequency, carrier_phase,
                start_sample, num_samples, tile_re, tile_im,
            )::typeof(correlator)
        end
    end
end

# Per-sat downconvert+correlate. Pure: returns the updated TrackedSat. Shared
# by `downconvert_and_correlate` and `downconvert_and_correlate!` across both
# CPU backends; the per-backend differences (slab buffer source, kernel
# choice) are dispatched via `_slab_buffer` and `_correlate_with_buffer!`.
function _update_tracked_sat_correlator(
    tracked_sat::TrackedSat,
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
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
    code_replica_size =
        num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
    new_correlator = _with_code_replica_buffer(
        dc, get_code_type(system), code_replica_size,
    ) do code_replica
        _correlate_with_buffer!(
            dc,
            code_replica,
            system,
            sat_state,
            signal,
            sample_shifts,
            code_frequency,
            carrier_frequency,
            sampling_frequency,
            signal_samples_to_integrate,
        )
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
                @inbounds new_vals[i] = _update_tracked_sat_correlator(
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
            @inbounds vals[i] = _update_tracked_sat_correlator(
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
