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

# Per-backend per-signal correlation kernel. Single-threaded backend uses
# `downconvert_and_correlate!` (which internally gen-code-replicas then
# runs the fused kernel); threaded backend pre-generates the code replica
# and feeds the fused kernel directly. Both write into the supplied
# `code_replica` buffer (from the calling correlator's `ScratchBuffers`)
# and return a value of the same correlator type as `correlator`.
#
# Args are per-signal scalars: the caller is responsible for picking the
# right `signal_type`, `correlator`, `code_phase` (modded into the signal's
# primary period), `code_frequency`, and `sample_shifts` for this signal.
@inline function _correlate_one_signal!(
    ::CPUDownconvertAndCorrelator,
    code_replica,
    signal_type,
    correlator,
    signal,
    sample_shifts,
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    signal_samples_to_integrate,
    prn,
)
    downconvert_and_correlate!(
        signal_type,
        signal,
        correlator,
        code_replica,
        code_phase,
        carrier_phase,
        code_frequency,
        carrier_frequency,
        sampling_frequency,
        signal_start_sample,
        signal_samples_to_integrate,
        prn,
    )::typeof(correlator)
end

@inline function _correlate_one_signal!(
    dc::CPUThreadedDownconvertAndCorrelator,
    code_replica,
    signal_type,
    correlator,
    signal,
    sample_shifts,
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    signal_samples_to_integrate,
    prn,
)
    gen_code_replica!(
        code_replica,
        signal_type,
        code_frequency,
        sampling_frequency,
        code_phase,
        signal_start_sample,
        signal_samples_to_integrate,
        sample_shifts,
        prn,
    )
    _fused_with_tile_scratch!(
        dc,
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        signal_start_sample,
        signal_samples_to_integrate,
    )::typeof(correlator)
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
            downconvert_and_correlate_fused!(
                correlator, signal, code_replica, sample_shifts,
                carrier_frequency, sampling_frequency, carrier_phase,
                start_sample, num_samples, tile_re, tile_im,
            )::typeof(correlator)
        end
    end
end

# Per-sat downconvert+correlate. Pure: returns the updated TrackedSat. Shared
# by `downconvert_and_correlate` and `downconvert_and_correlate!` across both
# CPU backends; the per-backend differences (kernel choice) are dispatched
# via `_correlate_one_signal!`.
#
# For multi-signal sats, the iteration window is the MIN samples-to-next-
# boundary across all signals (or buffer end, whichever is sooner). Each
# signal in `sat.signals` then runs its own (gen_code_replica +
# fused-kernel) call over that shared window — the carrier/code Doppler
# and start sample are sat-shared; the signal type, correlator, code
# replica buffer, and per-signal code phase differ.
function _update_tracked_sat_correlator(
    sat::TrackedSat,
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    signal,
    num_samples_signal,
    preferred_num_code_blocks_to_integrate,
    sampling_frequency,
    intermediate_frequency,
)
    # MIN samples-to-next-boundary across all signals on this sat, clamped
    # to remaining buffer. The signal-level boundary calc reads each
    # signal's primary-code-length-relative phase via mod(sat.code_phase,
    # get_code_length(signal)).
    samples_to_integrate, per_signal_completed = _calc_min_samples_and_completed(
        sat.signals,
        sat.signal_start_sample,
        sampling_frequency,
        sat.code_doppler,
        sat.code_phase,
        preferred_num_code_blocks_to_integrate,
        num_samples_signal,
    )
    if samples_to_integrate == 0
        return sat
    end
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    new_signals_data = _correlate_all_signals(
        sat.signals,
        per_signal_completed,
        dc,
        signal,
        sat.code_doppler,
        sat.code_phase,
        carrier_frequency,
        sat.carrier_phase,
        sampling_frequency,
        sat.signal_start_sample,
        samples_to_integrate,
        sat.prn,
        num_samples_signal,
    )
    update(
        sat,
        samples_to_integrate,
        intermediate_frequency,
        sampling_frequency,
        new_signals_data,
    )
end

# Compute (samples_to_integrate, per_signal_completed_tuple) via tuple
# recursion. Returns the MIN across all signals' samples-to-next-boundary,
# clamped to remaining buffer samples. Per-signal `completed` flags are
# derived after the MIN is known.
@inline function _calc_min_samples_and_completed(
    signals::Tuple,
    signal_start_sample,
    sampling_frequency,
    code_doppler,
    code_phase,
    preferred_num_code_blocks_to_integrate,
    num_samples_signal,
)
    per_signal_to_boundary = _per_signal_samples_to_boundary(
        signals, signal_start_sample, sampling_frequency, code_doppler, code_phase,
        preferred_num_code_blocks_to_integrate, num_samples_signal,
    )
    samples_to_integrate = _min_of_tuple(per_signal_to_boundary)
    signal_samples_left = num_samples_signal - signal_start_sample + 1
    samples_to_integrate = min(samples_to_integrate, signal_samples_left)
    per_signal_completed = _flag_completed(per_signal_to_boundary, samples_to_integrate)
    return samples_to_integrate, per_signal_completed
end

# For each signal, return its samples-to-next-primary-code-boundary using
# the signal-specific primary-code-relative phase. Tuple recursion keeps
# the heterogeneous walk inline / inference-friendly.
@inline _per_signal_samples_to_boundary(
    ::Tuple{}, _, _, _, _, _, _,
) = ()
@inline function _per_signal_samples_to_boundary(
    signals::Tuple,
    signal_start_sample,
    sampling_frequency,
    code_doppler,
    code_phase,
    preferred_num_code_blocks_to_integrate,
    num_samples_signal,
)
    s = first(signals).signal
    n_blocks = calc_num_code_blocks_to_integrate(
        s,
        preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(first(signals)),
    )
    # The signal's chips-to-next-boundary uses its own primary code length;
    # `mod(code_phase, get_code_length(s))` gives the signal's replica-
    # relative phase.
    per_signal_phase = mod(code_phase, get_code_length(s))
    n = calc_num_samples_left_to_integrate(
        s, n_blocks, sampling_frequency, code_doppler, per_signal_phase,
    )
    (n, _per_signal_samples_to_boundary(
        Base.tail(signals), signal_start_sample, sampling_frequency, code_doppler,
        code_phase, preferred_num_code_blocks_to_integrate, num_samples_signal,
    )...)
end

@inline _min_of_tuple(t::Tuple{Any}) = first(t)
@inline _min_of_tuple(t::Tuple) = min(first(t), _min_of_tuple(Base.tail(t)))

# For each signal, `is_completed = (chosen_samples == samples_to_boundary)`.
@inline _flag_completed(::Tuple{}, _) = ()
@inline _flag_completed(t::Tuple, chosen) =
    (first(t) == chosen, _flag_completed(Base.tail(t), chosen)...)

# For each signal, gen its code replica and run the fused kernel. Returns
# a tuple of `(new_correlator, is_integration_completed)` pairs.
@inline function _correlate_all_signals(
    signals::Tuple,
    per_signal_completed::Tuple,
    dc,
    signal,
    code_doppler,
    code_phase,
    carrier_frequency,
    carrier_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_samples_signal,
)
    head = first(signals)
    s = head.signal
    correlator = head.correlator
    code_frequency = code_doppler + get_code_frequency(s)
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    # Buffer sized for the full signal because the kernel indexes from
    # `signal_start_sample` rather than from sample 1.
    code_replica_size =
        num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
    per_signal_phase = mod(code_phase, get_code_length(s))
    new_corr = _with_code_replica_buffer(
        dc, get_code_type(s), code_replica_size,
    ) do code_replica
        _correlate_one_signal!(
            dc,
            code_replica,
            s,
            correlator,
            signal,
            sample_shifts,
            per_signal_phase,
            carrier_phase,
            code_frequency,
            carrier_frequency,
            sampling_frequency,
            signal_start_sample,
            samples_to_integrate,
            prn,
        )
    end
    head_pair = (new_corr, first(per_signal_completed))
    (head_pair, _correlate_all_signals(
        Base.tail(signals),
        Base.tail(per_signal_completed),
        dc, signal, code_doppler, code_phase, carrier_frequency, carrier_phase,
        sampling_frequency, signal_start_sample, samples_to_integrate, prn,
        num_samples_signal,
    )...)
end

@inline _correlate_all_signals(
    ::Tuple{}, ::Tuple{},
    _, _, _, _, _, _, _, _, _, _, _,
) = ()

"""
$(SIGNATURES)

Downconvert and correlate all available satellites on the CPU.
Returns a new `TrackState` whose slot vectors are detached from the
input — the input `track_state` is left untouched, satisfying
`track`'s immutability contract. Per-call code-replica scratch comes
from the correlator's `ScratchBuffers` so the kernel itself stays
allocation-free; the only per-call allocation is the slot-vector
copy needed for immutability.
"""
function downconvert_and_correlate(
    dc::CPUDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    new_track_state = TrackState(
        track_state;
        tracked_systems =
            _copy_slot_vectors(track_state.tracked_systems),
    )
    downconvert_and_correlate!(
        dc, signal, new_track_state,
        preferred_num_code_blocks_to_integrate, sampling_frequency,
        intermediate_frequency,
    )
end

"""
$(SIGNATURES)

In-place version: walks each system's `Vector{TrackedSat}` and overwrites
slots in place. Returns the same `track_state`. Allocation-free in steady
state — see [`track!`](@ref).
"""
# Per-system body for the single-threaded backend. Pulled out so
# `_foreach_system!` can call it on each per-system satellite dictionary in
# the (possibly heterogeneous) `tracked_systems` tuple without dynamic
# dispatch / boxing.
@inline function _dc_one_system!(
    sats::Dictionary{<:Any,<:TrackedSat}, dc::CPUDownconvertAndCorrelator,
    signal, num_samples_signal, preferred_num_code_blocks_to_integrate,
    sampling_frequency, intermediate_frequency,
)
    vals = sats.values
    isempty(vals) && return nothing
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(
            vals[i],
            dc,
            signal,
            num_samples_signal,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )
    end
    return nothing
end

function downconvert_and_correlate!(
    dc::CPUDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    _foreach_system!(
        _dc_one_system!, track_state.tracked_systems,
        dc, signal, num_samples_signal, preferred_num_code_blocks_to_integrate,
        sampling_frequency, intermediate_frequency,
    )
    return track_state
end

"""
$(SIGNATURES)

Multi-threaded downconvert and correlate. Returns a new `TrackState`
whose slot vectors are detached from the input — the input
`track_state` is left untouched. Per-call scratch comes from the
correlator's per-thread `ScratchBuffers`; the only per-call
allocation is the slot-vector copy needed for immutability.
"""
function downconvert_and_correlate(
    dc::CPUThreadedDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    new_track_state = TrackState(
        track_state;
        tracked_systems =
            _copy_slot_vectors(track_state.tracked_systems),
    )
    downconvert_and_correlate!(
        dc, signal, new_track_state,
        preferred_num_code_blocks_to_integrate, sampling_frequency,
        intermediate_frequency,
    )
end

"""
$(SIGNATURES)

In-place version: writes new `TrackedSat` values directly into each system's
existing `Vector{TrackedSat}` backing storage. Different `@batch` iterations
write to disjoint slots, so no synchronization is needed. Returns the same
`track_state`. Allocation-free in steady state — see [`track!`](@ref).
"""
# Per-system body for the threaded backend. Each `@batch` writes to
# disjoint slots in this system's `Vector{TrackedSat}`. Pulled out so
# `_foreach_system!` can call it without boxing on heterogeneous
# system tuples.
@inline function _dc_one_system_threaded!(
    sats::Dictionary{<:Any,<:TrackedSat}, dc::CPUThreadedDownconvertAndCorrelator,
    signal, num_samples_signal, preferred_num_code_blocks_to_integrate,
    sampling_frequency, intermediate_frequency,
)
    vals = sats.values
    n = length(vals)
    n == 0 && return nothing
    @batch for i = 1:n
        @inbounds vals[i] = _update_tracked_sat_correlator(
            vals[i],
            dc,
            signal,
            num_samples_signal,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )
    end
    return nothing
end

function downconvert_and_correlate!(
    dc::CPUThreadedDownconvertAndCorrelator,
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
)
    num_samples_signal = get_num_samples(signal)
    _foreach_system!(
        _dc_one_system_threaded!, track_state.tracked_systems,
        dc, signal, num_samples_signal, preferred_num_code_blocks_to_integrate,
        sampling_frequency, intermediate_frequency,
    )
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
