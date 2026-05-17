# Per-thread scratch byte buffers. The code-replica path uses
# `code_replica` for the first signal of a sat; multi-signal sats grow
# `extra_code_replicas` lazily to add slots 2..N. The fused kernel's
# `AbstractVector`-shifts fallback uses `tile_re` / `tile_im`
# concurrently; the multi-signal tile-share kernel shares them too.
# Holding separate buffers per role (instead of bump-arena offset math
# on one) keeps the call sites simple. Both CPU backends share this
# struct: the single-threaded one holds a single `ScratchBuffers`, the
# threaded one holds a `Vector{ScratchBuffers}` (one per thread).
mutable struct ScratchBuffers
    code_replica::Vector{UInt8}
    extra_code_replicas::Vector{Vector{UInt8}}
    tile_re::Vector{UInt8}
    tile_im::Vector{UInt8}
end

ScratchBuffers() = ScratchBuffers(UInt8[], Vector{UInt8}[], UInt8[], UInt8[])

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

# Grow `extra_code_replicas` to hold at least `n_extra` byte vectors.
@inline function _ensure_extra_code_replicas!(bufs::ScratchBuffers, n_extra::Int)
    while length(bufs.extra_code_replicas) < n_extra
        push!(bufs.extra_code_replicas, UInt8[])
    end
end

# Pick the `i`-th code-replica byte buffer for this thread: slot 1 is the
# primary `code_replica`; slots 2..N are extra slots, grown lazily.
@inline function _code_replica_slot(bufs::ScratchBuffers, i::Int)
    i == 1 && return bufs.code_replica
    _ensure_extra_code_replicas!(bufs, i - 1)
    bufs.extra_code_replicas[i - 1]
end

# Multi-signal tile-share scratch handle: yields ScratchViews over the
# shared `tile_re` / `tile_im` SoA buffers, sized for `num_samples`
# Float32s each.
@inline function _with_tile_buffers(
    f,
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    num_samples::Int,
)
    bufs = _scratch_buffers(dc)
    _with_scratch_view(bufs.tile_re, Float32, num_samples) do tile_re
        _with_scratch_view(bufs.tile_im, Float32, num_samples) do tile_im
            f(tile_re, tile_im)
        end
    end
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
    new_signals_data = _correlate_signals(
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

# Single-signal path: gen one code replica + run the in-register fused
# kernel. Returns a one-tuple of `(new_correlator, is_integration_completed)`.
# This is the hot path for N=1 — the fused kernel keeps downconverted
# samples in registers and is ~24% faster than the tile-share kernel
# below at N=1.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal},
    per_signal_completed::Tuple{Bool},
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
    head = signals[1]
    s = head.signal
    correlator = head.correlator
    code_frequency = code_doppler + get_code_frequency(s)
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
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
    ((new_corr, per_signal_completed[1]),)
end

# Multi-signal path (N >= 2): gen N code replicas, then call the tile-
# share fused kernel that does one downconvert + a sample-outer fused
# correlate over all N×NC accumulators. Beats fused-N-times by ~40% at
# N=2 and ~53% at N=3 — see the design doc.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal, TrackedSignal, Vararg{TrackedSignal}},
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
    # Generate each signal's code replica into its per-thread slot, then
    # call the tuple kernel with N replicas + the shared tile buffers.
    code_replicas =
        _gen_all_code_replicas(
            signals, dc, code_doppler, code_phase, sampling_frequency,
            signal_start_sample, samples_to_integrate, prn, num_samples_signal,
        )
    correlators = _signal_correlators(signals)
    sample_shifts_tuple = _signal_sample_shifts(
        signals, code_doppler, sampling_frequency,
    )
    new_correlators = _with_tile_buffers(dc, num_samples_signal) do tile_re, tile_im
        downconvert_and_correlate_fused_tuple!(
            correlators,
            signal,
            code_replicas,
            sample_shifts_tuple,
            carrier_frequency,
            sampling_frequency,
            carrier_phase,
            signal_start_sample,
            samples_to_integrate,
            tile_re,
            tile_im,
        )
    end
    _zip_correlators_with_completed(new_correlators, per_signal_completed)
end

# Generate code replicas for each signal into the per-thread scratch
# slots. Returns a tuple of `ScratchView`s suitable for passing to
# `downconvert_and_correlate_fused_tuple!`. The buffer for signal `i`
# comes from `_code_replica_slot(_, i)`.
#
# Implemented via an enumerated `map` over the heterogeneous signals
# tuple — this stays type-stable and inlines across tuple lengths,
# whereas recursive splatting bails out of the small-N inline path
# around N=3 and introduces per-call boxing.
@inline function _gen_all_code_replicas(
    signals::Tuple,
    dc,
    code_doppler,
    code_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_samples_signal,
)
    bufs = _scratch_buffers(dc)
    idx_tuple = ntuple(identity, length(signals))
    map(signals, idx_tuple) do head, i
        s = head.signal
        code_frequency = code_doppler + get_code_frequency(s)
        sample_shifts = get_correlator_sample_shifts(
            head.correlator, sampling_frequency, code_frequency,
        )
        code_replica_size =
            num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
        per_signal_phase = mod(code_phase, get_code_length(s))
        slot = _code_replica_slot(bufs, i)
        CT = get_code_type(s)
        nbytes = code_replica_size * sizeof(CT)
        length(slot) < nbytes && resize!(slot, nbytes)
        view = GC.@preserve slot ScratchView{CT}(Ptr{CT}(pointer(slot)), code_replica_size)
        gen_code_replica!(
            view, s, code_frequency, sampling_frequency, per_signal_phase,
            signal_start_sample, samples_to_integrate, sample_shifts, prn,
        )
        view
    end
end

@inline _signal_correlators(signals::Tuple) = map(s -> s.correlator, signals)

@inline _signal_sample_shifts(
    signals::Tuple, code_doppler, sampling_frequency,
) = map(signals) do head
    s = head.signal
    code_frequency = code_doppler + get_code_frequency(s)
    get_correlator_sample_shifts(head.correlator, sampling_frequency, code_frequency)
end

@inline _zip_correlators_with_completed(corrs::Tuple, completed::Tuple) =
    map(tuple, corrs, completed)

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
    measurements::Measurements,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
)
    new_track_state = TrackState(
        track_state;
        groups = _copy_groups_slot_vectors(track_state.groups),
    )
    downconvert_and_correlate!(
        dc, measurements, new_track_state,
        preferred_num_code_blocks_to_integrate,
    )
end

"""
$(SIGNATURES)

In-place version: walks each system's `Vector{TrackedSat}` and overwrites
slots in place. Returns the same `track_state`. Allocation-free in steady
state — see [`track!`](@ref).
"""
# Per-group body for the single-threaded backend. Pulled out so
# `_foreach_system!` can call it on each `SignalGroup` in the
# (possibly heterogeneous) `groups` tuple without dynamic dispatch /
# boxing. Routes to this group's band's `Measurement` for the signal
# buffer and front-end metadata.
@inline function _dc_one_system!(
    g::SignalGroup, dc::CPUDownconvertAndCorrelator,
    measurements::Measurements, preferred_num_code_blocks_to_integrate,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    m = measurements[band_key(g.band)]
    signal = m.samples
    num_samples_signal = get_num_samples(m)
    sampling_frequency = m.sampling_frequency
    intermediate_frequency = m.intermediate_frequency
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
    measurements::Measurements,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
)
    _foreach_system!(
        _dc_one_system!, track_state.groups,
        dc, measurements, preferred_num_code_blocks_to_integrate,
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
    measurements::Measurements,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
)
    new_track_state = TrackState(
        track_state;
        groups = _copy_groups_slot_vectors(track_state.groups),
    )
    downconvert_and_correlate!(
        dc, measurements, new_track_state,
        preferred_num_code_blocks_to_integrate,
    )
end

"""
$(SIGNATURES)

In-place version: writes new `TrackedSat` values directly into each system's
existing `Vector{TrackedSat}` backing storage. Different `@batch` iterations
write to disjoint slots, so no synchronization is needed. Returns the same
`track_state`. Allocation-free in steady state — see [`track!`](@ref).
"""
# Per-group body for the threaded backend. Each `@batch` writes to
# disjoint slots in this group's `Vector{TrackedSat}`. Pulled out so
# `_foreach_system!` can call it without boxing on heterogeneous
# group tuples. Routes to this group's band's `Measurement`.
@inline function _dc_one_system_threaded!(
    g::SignalGroup, dc::CPUThreadedDownconvertAndCorrelator,
    measurements::Measurements, preferred_num_code_blocks_to_integrate,
)
    vals = g.satellites.values
    n = length(vals)
    n == 0 && return nothing
    m = measurements[band_key(g.band)]
    signal = m.samples
    num_samples_signal = get_num_samples(m)
    sampling_frequency = m.sampling_frequency
    intermediate_frequency = m.intermediate_frequency
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
    measurements::Measurements,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
)
    _foreach_system!(
        _dc_one_system_threaded!, track_state.groups,
        dc, measurements, preferred_num_code_blocks_to_integrate,
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
