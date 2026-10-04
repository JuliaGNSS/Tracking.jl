# Per-thread scratch byte buffers, one per role: `code_replica` for a sat's first
# signal, `extra_code_replicas` (grown lazily) for signals 2..N, and `tile_re` /
# `tile_im` for the dynamic-shifts fused kernel and the multi-signal tile-share
# kernel. The single-threaded backend holds one, the threaded backend one per
# thread. Buffers are only `resize!`d / `push!`ed in place, hence immutable.
struct ScratchBuffers
    code_replica::Vector{UInt8}
    extra_code_replicas::Vector{Vector{UInt8}}
    tile_re::Vector{UInt8}
    tile_im::Vector{UInt8}
end

ScratchBuffers() = ScratchBuffers(UInt8[], Vector{UInt8}[], UInt8[], UInt8[])

# Isbits pointer-and-length typed view over a `Vector{UInt8}` slot (free to
# build). Implements just enough `AbstractArray` for the kernels: linear
# indexing, `size`, and `pointer` (for SIMD `vstore`).
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

CPU-based implementation of downconversion and correlation. Holds one set of
long-lived scratch buffers (code replica and tile) that grow lazily on first use,
so a hoisted instance has zero allocations per `track!` call in steady state.

For real-time loops, construct the correlator **once outside** the
`track!` loop and pass it via the `downconvert_and_correlator` keyword
argument — the default value rebuilds the buffers on every call.
"""
struct CPUDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::ScratchBuffers
end

CPUDownconvertAndCorrelator() = CPUDownconvertAndCorrelator(ScratchBuffers())

# Residual allocation of the threaded backend: Polyester copies everything a
# `@batch` region references into one argument tuple per launch, elided only if
# every member is `isbits` (bare `Array`s become `PtrArray`s); otherwise it is
# heap-allocated at its `sizeof`, immutable structs by value, mutable objects as
# one pointer. Here the non-`isbits` member is the GNSS signal (a code `Matrix` and
# a `SignalLUT`), needed for code-replica generation: ~96 B per launch. Noise
# despreads ride the same loop through one pointer (`_park_noise_items!`), at no
# extra residual. Measured alternatives: generating from the LUT's bare arrays
# inside the loop reaches 0 B/launch at full throughput but needs a bare-array
# `gen_code!` in GNSSSignals; a serial code-gen pre-pass reaches 0 B but
# serializes ~30 % of the work. Neither is worth ~96 B per code block.
"""
$(SIGNATURES)

Multi-threaded CPU downconvert and correlate, parallelized over the
satellites (PRNs) of each group. Holds one set of scratch buffers per thread,
indexed by `Threads.threadid()` (stable within a `@batch` iteration), so a
hoisted instance's scratch is allocation-free in steady state, apart from a
small residual (~96 B per completed code block, with more than one thread and
more than one work item in a group's loop; a measured noise density counts as
one).

For real-time loops, construct the correlator **once outside** the
`track!` loop and pass it via the `downconvert_and_correlator` keyword
argument — the default value rebuilds the per-thread buffers on every
call.
"""
struct CPUThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{ScratchBuffers}
end

CPUThreadedDownconvertAndCorrelator() =
    CPUThreadedDownconvertAndCorrelator([ScratchBuffers() for _ = 1:Threads.maxthreadid()],)

# The active `ScratchBuffers` for this thread. The index stays bounds-checked:
# foreign threads adopted (`jl_adopt_thread`) after construction can exceed
# `Threads.maxthreadid()` and must fail loudly rather than read out of bounds.
@inline _scratch_buffers(dc::CPUDownconvertAndCorrelator) = dc.buffers
@inline _scratch_buffers(dc::CPUThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

# Grow the role's byte buffer to fit `n` elements of `T` and hand a typed
# `ScratchView` of exactly that size to `f`; allocation-free once the buffer has
# reached its working size. `GC.@preserve` roots the buffer while `f` runs, since
# the view holds a raw `Ptr`.
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
    bufs.extra_code_replicas[i-1]
end

# Yields `ScratchView`s over the `tile_re` / `tile_im` SoA buffers, sized for
# `num_ants * num_samples` Float32s each (antennas back-to-back).
@inline function _with_tile_buffers(
    f,
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    num_samples::Int,
    num_ants::Int = 1,
)
    bufs = _scratch_buffers(dc)
    n = num_samples * num_ants
    _with_scratch_view(bufs.tile_re, Float32, n) do tile_re
        _with_scratch_view(bufs.tile_im, Float32, n) do tile_im
            f(tile_re, tile_im)
        end
    end
end

# THE per-backend correlation primitive: one despread of one signal over one range
# of samples, returning the updated correlator (same type). Every despread goes
# through it — the per-satellite single-signal path (`_correlate_signals`) and the
# open-loop noise reference (`update_noise!`) alike. That makes the noise
# measurement model-free by construction: the reference takes the identical
# quantise → downconvert → despread → accumulate path as the prompt, so `N₀`
# carries the quantisation loss and code amplitude with no per-backend correction
# (and is on the same scale even where the one-/two-bit accumulators are popcount
# counts rather than sample sums).
#
# `code_phase` is modded into the signal's code wrap; `carrier_phase` is the NCO's
# phase on the satellite path and `0.0` for the noise reference.
#
# `code_replica_size` is only used by backends that generate a replica (the CPU
# ones; Int16 and bit-wise pack the code inside the kernel). Size it for the whole
# sample buffer plus the tap spread, not the slice: `gen_code_replica!` writes at
# `start_sample` and the kernel reads offset by `start_sample - 1`. The replica
# comes from the per-thread `ScratchBuffers`, shared by satellites and the noise
# reference; each work item uses it only within its own iteration on its thread.
#
# `use_band_cache = false` (noise reference) keeps the bit-wise backends from
# reading their per-group packed band sign planes, which may belong to another
# band or buffer. Backends without such a cache ignore it.
@inline function _despread_one_signal!(
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    correlator,
    samples,
    signal_type,
    prn,
    sample_shifts,
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    start_sample,
    num_samples,
    code_replica_size,
    use_band_cache::Bool = true,
)
    # GNSSSignals' embedded-LUT `gen_code!` is Int8-only, so the replica is `Int8`
    # (the fused kernel widens it to Float32); CBOC (Galileo E1B) is its Int8
    # integer-amplitude approximation.
    _with_code_replica_buffer(dc, Int8, code_replica_size) do code_replica
        gen_code_replica!(
            code_replica,
            signal_type,
            code_frequency,
            sampling_frequency,
            code_phase,
            start_sample,
            num_samples,
            sample_shifts,
            prn,
        )
        _fused_with_tile_scratch!(
            dc,
            correlator,
            samples,
            code_replica,
            sample_shifts,
            carrier_frequency,
            sampling_frequency,
            carrier_phase,
            start_sample,
            num_samples,
        )::typeof(correlator)
    end
end

# Fused-kernel dispatch: static (`SVector`) shifts use the in-register kernel,
# which needs no tiles; dynamic (`AbstractVector`) shifts get SoA tile buffers
# from the backend's `ScratchBuffers`, keeping that path allocation-free.
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
    downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        start_sample,
        num_samples,
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
    _with_tile_buffers(dc, num_samples, M) do tile_re, tile_im
        downconvert_and_correlate_fused!(
            correlator,
            signal,
            code_replica,
            sample_shifts,
            carrier_frequency,
            sampling_frequency,
            carrier_phase,
            start_sample,
            num_samples,
            tile_re,
            tile_im,
        )::typeof(correlator)
    end
end

# Backend-less variant of `_fused_with_tile_scratch!` for the public
# single-satellite `downconvert_and_correlate!`: dynamic shifts allocate their
# tile buffers per call here (convenience path; the hot paths use pooled scratch).
@inline function _fused_standalone!(
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
    downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        start_sample,
        num_samples,
    )::typeof(correlator)
end

@inline function _fused_standalone!(
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
    tile_re = Vector{Float32}(undef, num_samples * M)
    tile_im = Vector{Float32}(undef, num_samples * M)
    downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        start_sample,
        num_samples,
        tile_re,
        tile_im,
    )::typeof(correlator)
end

# Per-sat downconvert+correlate, shared by every backend; returns the updated
# `TrackedSat`. Kernel choice is dispatched via `_correlate_signals` /
# `_despread_one_signal!`. Each sub-step's window is the MIN samples-to-next-
# boundary across the sat's signals, clamped to the chunk end; carrier/code
# Doppler and start sample are sat-shared, code phase and correlator per signal.
function _update_tracked_sat_correlator(
    sat::TrackedSat,
    dc::AbstractDownconvertAndCorrelator,
    signal,
    num_samples_signal,
    chunk_last_sample,
    sampling_frequency,
    intermediate_frequency,
    stop_before_partial::Bool = false,
)
    # `update` snapshots a `CorrelatorOutput` and resets the accumulator for each
    # signal completing on a sub-step, so a chunk yields 0..n outputs per signal;
    # a partial integration carries in the accumulator. The NCO Doppler is
    # untouched here. For `stop_before_partial`, see `downconvert_and_correlate!`:
    # it makes every completed integration run at a single Doppler.
    while sat.signal_start_sample <= chunk_last_sample
        samples_to_integrate, per_signal_completed = _calc_min_samples_and_completed(
            sat.signals,
            sat.signal_start_sample,
            sampling_frequency,
            sat.code_doppler,
            sat.code_phase,
            num_samples_signal,
            chunk_last_sample,
        )
        samples_to_integrate == 0 && break
        # A sub-step that completes no signal is exactly the chunk-clamped
        # trailing partial — defer it to the post-estimate pass.
        stop_before_partial && !_any_of_tuple(per_signal_completed) && break
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
        sat = update(
            sat,
            samples_to_integrate,
            intermediate_frequency,
            sampling_frequency,
            new_signals_data,
        )
    end
    return sat
end

# `(samples_to_integrate, per_signal_completed)`: the MIN samples-to-next-boundary
# across signals, clamped to the chunk end, and which signals it completes.
@inline function _calc_min_samples_and_completed(
    signals::Tuple,
    signal_start_sample,
    sampling_frequency,
    code_doppler,
    code_phase,
    num_samples_signal,
    chunk_last_sample,
)
    per_signal_to_boundary = _per_signal_samples_to_boundary(
        signals,
        signal_start_sample,
        sampling_frequency,
        code_doppler,
        code_phase,
        num_samples_signal,
    )
    samples_to_integrate = _min_of_tuple(per_signal_to_boundary)
    # Clamp to the chunk end; replica sizing still uses the true buffer length.
    samples_left = chunk_last_sample - signal_start_sample + 1
    samples_to_integrate = min(samples_to_integrate, samples_left)
    per_signal_completed = _flag_completed(per_signal_to_boundary, samples_to_integrate)
    return samples_to_integrate, per_signal_completed
end

# Per-signal samples to the end of its current integration window. Tuple
# recursion keeps the heterogeneous walk inferable.
@inline _per_signal_samples_to_boundary(::Tuple{}, _, _, _, _, _) = ()
@inline function _per_signal_samples_to_boundary(
    signals::Tuple,
    signal_start_sample,
    sampling_frequency,
    code_doppler,
    code_phase,
    num_samples_signal,
)
    head = first(signals)
    # Measured against the secondary-aware wrap (see `_replica_code_wrap`) so a
    # multi-block window ends on a data-symbol boundary.
    p = _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    )
    n = calc_num_samples_left_to_integrate(
        head.signal,
        p.n_blocks,
        sampling_frequency,
        code_doppler,
        p.signal_code_phase,
    )
    (
        n,
        _per_signal_samples_to_boundary(
            Base.tail(signals),
            signal_start_sample,
            sampling_frequency,
            code_doppler,
            code_phase,
            num_samples_signal,
        )...,
    )
end

@inline _min_of_tuple(t::Tuple{Any}) = first(t)
@inline _min_of_tuple(t::Tuple) = min(first(t), _min_of_tuple(Base.tail(t)))

@inline _any_of_tuple(::Tuple{}) = false
@inline _any_of_tuple(t::Tuple) = first(t) || _any_of_tuple(Base.tail(t))

# For each signal, `is_completed = (chosen_samples == samples_to_boundary)`.
@inline _flag_completed(::Tuple{}, _) = ()
@inline _flag_completed(t::Tuple, chosen) =
    (first(t) == chosen, _flag_completed(Base.tail(t), chosen)...)

# Code-phase wrap for replica generation and integration-window sizing. Past sync
# it is the full secondary-/bit-period length, so `gen_code!` bakes the secondary
# code at the right chip and its sign is wiped per block in the replica — at any
# integration length, including one block (issue #125); a multi-block window then
# also ends on a data-symbol boundary instead of summing against a misaligned
# overlay. At one block the boundary is unchanged, as `calc_num_chips_to_integrate`
# re-mods by the primary length. Pre-sync only the primary phase is known, so the
# wrap is the primary length. `num_blocks` is kept for call-site symmetry.
@inline _replica_code_wrap(tsig::TrackedSignal, num_blocks::Integer) =
    has_bit_or_secondary_code_been_found(tsig) ? _post_sync_code_length(tsig) :
    get_code_length(tsig.signal)

# All per-signal replica/kernel parameters, derived in one place so the boundary
# calc, replica sizing/generation and kernel tap offsets cannot diverge
# (issue #133).
@inline function _signal_replica_params(
    tsig::TrackedSignal,
    code_doppler,
    code_phase,
    sampling_frequency,
    num_samples_signal,
)
    s = tsig.signal
    code_frequency = code_doppler + get_code_frequency(s)
    sample_shifts =
        get_correlator_sample_shifts(tsig.correlator, sampling_frequency, code_frequency)
    code_replica_size = num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
    n_blocks = calc_num_code_blocks_to_integrate(
        s,
        tsig.preferred_num_code_blocks_to_integrate,
        has_bit_or_secondary_code_been_found(tsig),
    )
    signal_code_phase = mod(code_phase, _replica_code_wrap(tsig, n_blocks))
    (; code_frequency, sample_shifts, code_replica_size, n_blocks, signal_code_phase)
end

# Single-signal path, shared by every backend (untyped `dc`; the backend
# difference lives in `_despread_one_signal!`). Returns a one-tuple of
# `(new_correlator, is_integration_completed)`. Not routed through the tile-share
# kernel: the in-register fused kernel is ~24% faster at N=1.
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
    p = _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    )
    new_corr = _despread_one_signal!(
        dc,
        head.correlator,
        signal,
        head.signal,
        prn,
        p.sample_shifts,
        p.signal_code_phase,
        carrier_phase,
        p.code_frequency,
        carrier_frequency,
        sampling_frequency,
        signal_start_sample,
        samples_to_integrate,
        p.code_replica_size,
    )
    ((new_corr, per_signal_completed[1]),)
end

# Multi-signal path (N >= 2): generate N code replicas, then one tile-share
# kernel call (see `downconvert_and_correlate_fused_tuple!`).
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal,TrackedSignal,Vararg{TrackedSignal}},
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
    # The replica views hold raw pointers into `dc`'s scratch, so root `dc` from
    # generation through the kernel call.
    new_correlators = GC.@preserve dc begin
        code_replicas = _gen_all_code_replicas(
            signals,
            dc,
            code_doppler,
            code_phase,
            sampling_frequency,
            signal_start_sample,
            samples_to_integrate,
            prn,
            num_samples_signal,
        )
        correlators = _signal_correlators(signals)
        sample_shifts_tuple = _signal_sample_shifts(
            signals,
            code_doppler,
            code_phase,
            sampling_frequency,
            num_samples_signal,
        )
        # All correlators in this sat agree on M (enforced by SignalGroup);
        # read it from the first correlator's type.
        num_ants = get_num_ants(correlators[1])
        _with_tile_buffers(dc, num_samples_signal, num_ants) do tile_re, tile_im
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
    end
    _zip_correlators_with_completed(new_correlators, per_signal_completed)
end

# Generate each signal's code replica into its per-thread slot
# (`_code_replica_slot(_, i)`) and return the `ScratchView`s. An enumerated `map`
# stays type-stable across tuple lengths; recursive splatting boxed around N=3.
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
        p = _signal_replica_params(
            head,
            code_doppler,
            code_phase,
            sampling_frequency,
            num_samples_signal,
        )
        slot = _code_replica_slot(bufs, i)
        CT = Int8  # embedded-LUT `gen_code!` is Int8-only; see `_despread_one_signal!`
        nbytes = p.code_replica_size * sizeof(CT)
        length(slot) < nbytes && resize!(slot, nbytes)
        # The pointer outlives this function; the caller's `GC.@preserve dc` roots it.
        view = ScratchView{CT}(Ptr{CT}(pointer(slot)), p.code_replica_size)
        gen_code_replica!(
            view,
            s,
            p.code_frequency,
            sampling_frequency,
            p.signal_code_phase,
            signal_start_sample,
            samples_to_integrate,
            p.sample_shifts,
            prn,
        )
        view
    end
end

@inline _signal_correlators(signals::Tuple) = map(s -> s.correlator, signals)

# Kernel tap offsets — derived through the same `_signal_replica_params`
# the replica generation uses, so the two can never skew apart.
@inline _signal_sample_shifts(
    signals::Tuple,
    code_doppler,
    code_phase,
    sampling_frequency,
    num_samples_signal,
) = map(signals) do head
    _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    ).sample_shifts
end

@inline _zip_correlators_with_completed(corrs::Tuple, completed::Tuple) =
    map(tuple, corrs, completed)

"""
$(SIGNATURES)

Downconvert and correlate all available satellites. Defined on
`AbstractDownconvertAndCorrelator`, so every backend — the CPU ones, the Int16
one and both bit-wise ones — is served by this method; a backend customises
what happens inside it through `_despread_one_signal!`, `_dc_one_group!`,
`_check_sample_type` and `_threading`. Returns a new
`TrackState` whose slot *values* are detached from the input, but whose key set
(`Indices`) is *shared* with it to avoid copying the hash table every iteration;
the key set is detached once at the `track` boundary
(`reset_start_sample_and_bit_buffer`). Do not
`add_satellite!`/`remove_satellite!` on this function's direct output, or
you will corrupt the input's keys (#123). The copy is otherwise shallow:
per-sat scratch vectors are shared, see [`track`](@ref). The only per-call
allocation is the slot-value copy.

The **noise estimators are shared too**, and this call advances them: a
[`CorrelatorNoiseEstimator`](@ref)'s sliding window and RNG stream are written in
place, so the input's window grows and its draws advance. They are deliberately
not copied: a window holds ~1000 observations per signal, and copying it per call
would cost O(K) and reset the reference's PRN rotation. Consequently:

  - Branching two outputs from one input leaves them sharing one noise reference;
    the same samples enter its window twice.
  - Advancing two such states **concurrently** races on that window and RNG. Give
    each thread its own `TrackState` (built separately, not branched), or its own
    `noise_estimators` entry.
"""
function downconvert_and_correlate(
    dc::AbstractDownconvertAndCorrelator,
    measurements::BandMeasurements,
    track_state::TrackState;
    kwargs...,
)
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    downconvert_and_correlate!(dc, measurements, new_track_state; kwargs...)
end

"""
$(SIGNATURES)

In-place version: writes new `TrackedSat` values directly into each
group's existing `Vector{TrackedSat}` backing storage. On the threaded
backend, different `@batch` iterations write to disjoint slots, so no
synchronization is needed. Returns the same `track_state`.
Allocation-free in steady state — see [`track!`](@ref).

`chunk_duration` and `chunk_index` restrict the pass to one chunk of a fixed
per-band time grid: each satellite integrates only up to sample
`min(round((chunk_index + 1) * chunk_duration * sampling_frequency), num_samples)`.
The boundary is re-anchored to the absolute `chunk_index` each call (not
accumulated), so rounding never drifts and different bands stay time-aligned.
`chunk_duration = nothing` (the default) disables chunking — the whole buffer
is consumed as one chunk. Every completed integration is snapshotted into its
signal's `correlator_outputs`; only the trailing partial stays in the live
correlator.

`stop_before_partial = true` additionally stops each satellite at its last
completed code-block boundary inside the chunk instead of integrating the
chunk-clamped trailing partial. `track!`'s per-chunk pass uses this so the
residue is integrated by the *next* chunk's pass, entirely at the Doppler the
estimator writes in between; leave it `false` (the default) to consume the
whole window.

`samples_unchanged = true` promises that every band's sample buffer holds the
same content as on the previous call with this `dc`; backends may then reuse
sample-derived caches (the bit-wise backends skip re-packing their shared
band sign planes). `track!` passes it on every chunk after the first so the
pack happens once per call; leave it `false` (the default) whenever the
buffers may have been refilled.

`measure_noise = false` skips the per-signal noise measurement, leaving every
[`CorrelatorNoiseEstimator`](@ref)'s window untouched. Pass it on a call that
re-covers already-measured samples, as `track!`'s final drain pass does.

A band **absent** from `measurements` is not measured either, and does not throw:
a multi-band `TrackState` may be advanced one band at a time. Groups holding
satellites on a band still require that band's measurement.

The per-signal noise estimators are mutable state shared with the input on the
out-of-place `downconvert_and_correlate` and [`track`](@ref) — see there.
"""
function downconvert_and_correlate!(
    dc::AbstractDownconvertAndCorrelator,
    measurements::BandMeasurements,
    track_state::TrackState;
    chunk_index::Int = 0,
    chunk_duration = nothing,
    stop_before_partial::Bool = false,
    samples_unchanged::Bool = false,
    measure_noise::Bool = true,
)
    # `measure_noise = false` skips resolving the descriptors, not just running
    # them. Branching here (not in `_noise_items`) keeps the descriptor tuple's type
    # independent of a runtime value; the `()` arm reuses the group-tail
    # specialisation.
    groups = Tuple(track_state.groups)
    if measure_noise
        noise_items =
            _noise_items(dc, track_state, measurements, chunk_index, chunk_duration)
        _dc_groups!(
            groups,
            noise_items,
            _num_noise_items(noise_items),
            dc,
            measurements,
            chunk_index,
            chunk_duration,
            stop_before_partial,
            samples_unchanged,
        )
    else
        _dc_groups!(
            groups,
            (),
            0,
            dc,
            measurements,
            chunk_index,
            chunk_duration,
            stop_before_partial,
            samples_unchanged,
        )
    end
    return track_state
end

# Walk the groups, handing the chunk's noise work to the first. The items do not
# belong to that group (each carries its own band measurement and signal); they
# only ride its per-satellite loop, so a threaded backend runs them alongside the
# satellite correlations. Resolving once and running in one loop keeps each
# signal despread exactly once however satellites are grouped. A `TrackState` with
# no groups has no noise estimators, so nothing is dropped.
@inline _dc_groups!(::Tuple{}, args::Vararg{Any,8}) = nothing
@inline function _dc_groups!(
    groups::Tuple{SignalGroup,Vararg{SignalGroup}},
    noise_items,
    n_noise::Int,
    dc,
    measurements,
    chunk_index,
    chunk_duration,
    stop_before_partial,
    samples_unchanged,
)
    _dc_one_group!(
        first(groups),
        dc,
        measurements,
        chunk_index,
        chunk_duration,
        stop_before_partial,
        samples_unchanged,
        noise_items,
        n_noise,
    )
    _dc_groups!(
        Base.tail(groups),
        (),
        0,
        dc,
        measurements,
        chunk_index,
        chunk_duration,
        stop_before_partial,
        samples_unchanged,
    )
end

# The chunk's noise work: one descriptor per measured signal, keyed off
# `noise_estimators` (so one per signal however satellites are grouped). Resolved
# on the caller's thread, so `_check_sample_type` throws its `ArgumentError` there
# rather than from a worker, and the `@batch` closure captures one pointer (see
# `_park_noise_items!`).
#
# Walked by tuple recursion because estimator types differ per signal; a runtime
# loop would be type-unstable and allocate. Only type-domain decisions are made
# here (signal not tracked, band absent), so the tuple's length depends on types
# alone; runtime ones (`measure_noise`, an empty chunk via `update_noise!`'s
# `num_samples > 0` guard) are left to the caller and the callee.
@inline function _noise_items(
    dc::AbstractDownconvertAndCorrelator,
    track_state::TrackState,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
)
    noise_estimators = track_state.noise_estimators
    _park_noise_items!(
        track_state.noise_descriptor,
        _noise_items_walk(
            Val(keys(noise_estimators)),
            Tuple(noise_estimators),
            dc,
            track_state.groups,
            measurements,
            chunk_index,
            chunk_duration,
        ),
    )
end

# Park the descriptors in the `TrackState`'s reusable cell and return the box
# holding them. Polyester copies immutable `@batch` captures by value into its
# per-launch argument tuple (see `CPUThreadedDownconvertAndCorrelator`): the
# descriptors by value cost 152 B (~96 → 240 B per launch), a mutable box 8 B.
#
# The box's type names the caller's `BandMeasurement`, unknown until the first
# call, so the cell is a type-erased `Base.RefValue{Any}` holding a concretely
# typed box, created once and reused (a new one only if the sample type changes).
# Storing into it is allocation-free and it roots the descriptors; it keeps the
# last chunk's sample buffer reachable until the next call. Nothing to measure
# stays `()`, adding no launch cost.
@inline _park_noise_items!(::Base.RefValue{Any}, ::Tuple{}) = ()
@inline function _park_noise_items!(
    cell::Base.RefValue{Any},
    items::T,
) where {T<:Tuple{Any,Vararg{Any}}}
    box = cell[]
    box isa Base.RefValue{T} || return _new_noise_items_box!(cell, items)
    box[] = items
    box
end

# First call with this descriptor type (or the first after the caller changed one).
# `@noinline` so the hot path above stays a load, a type check and a store.
@noinline function _new_noise_items_box!(
    cell::Base.RefValue{Any},
    items::T,
) where {T<:Tuple}
    box = Base.RefValue{T}(items)
    cell[] = box
    box
end

# How many work items the loop must add to the satellites' range. A property of the
# descriptor set's *type*, so it folds to a literal at every call site.
@inline _num_noise_items(::Tuple{}) = 0
@inline _num_noise_items(::Base.RefValue{<:Tuple{Vararg{Any,N}}}) where {N} = N

# What a group loop is handed for the chunk's noise work: `()` when there is nothing
# to measure — including for every group after the one that carries it — and
# otherwise the single box the descriptors were parked in.
const _NoiseWork = Union{Tuple{},Base.RefValue{<:Tuple}}

@inline _noise_items_walk(
    ::Val{()},
    ::Tuple{},
    ::AbstractDownconvertAndCorrelator,
    ::SignalGroups,
    ::BandMeasurements,
    ::Int,
    _,
) = ()
@inline function _noise_items_walk(
    ::Val{K},
    estimators::Tuple,
    dc::AbstractDownconvertAndCorrelator,
    groups::SignalGroups,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
) where {K}
    (
        _noise_item(
            first(estimators),
            Val(first(K)),
            dc,
            groups,
            measurements,
            chunk_index,
            chunk_duration,
        )...,
        _noise_items_walk(
            Val(Base.tail(K)),
            Base.tail(estimators),
            dc,
            groups,
            measurements,
            chunk_index,
            chunk_duration,
        )...,
    )
end

# One signal's descriptor, as a 0- or 1-tuple so the "nothing to measure" case is
# an empty tuple in the *type* domain rather than a `Nothing` to branch on later.
@inline function _noise_item(
    estimator::AbstractNoiseEstimator,
    ::Val{K},
    dc::AbstractDownconvertAndCorrelator,
    groups::SignalGroups,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
) where {K}
    # A signal no group tracks measures nothing (the correlator-ingest case, filled
    # by `append_noise_observation!` instead).
    found = _first_signal_with_id(Tuple(groups), Val(K))
    isempty(found) && return ()
    signal = first(found)
    # Samples are per band, the despreading code per signal. A band absent from
    # `measurements` measures nothing (bands may be advanced one at a time); a
    # group with satellites on it still throws in `_dc_one_group!`.
    found_m =
        _band_measurement_if_present(measurements, Val(_signal_band_id(typeof(signal))))
    isempty(found_m) && return ()
    m = first(found_m)
    _check_sample_type(dc, m)
    num_samples = get_num_samples(m)
    ((
        estimator,
        m,
        _chunk_first_sample(chunk_duration, chunk_index, m.sampling_frequency, num_samples),
        _chunk_last_sample(chunk_duration, chunk_index, m.sampling_frequency, num_samples),
        signal,
        chunk_index,
    ),)
end

# The `NoiseUpdateContext` is built here, not stored in the descriptor, because it
# carries the backend, which the `@batch` closure already captures.
@inline _apply_noise_item!(item::Tuple, dc) = (
    update_noise!(
        item[1],
        item[2],
        item[3],
        item[4],
        NoiseUpdateContext(item[5], item[6], dc),
    );
    nothing
)

# Run the `j`-th descriptor. `j` is a runtime index into a heterogeneous tuple, so
# the walk unrolls to `n` comparisons (one or two in practice).
@inline _run_noise_item!(::Int, ::Tuple{}, ::Any) = nothing
@inline _run_noise_item!(j::Int, box::Base.RefValue{<:Tuple}, dc) =
    _run_noise_item!(j, box[], dc)
@inline _run_noise_item!(j::Int, items::Tuple, dc) =
    j == 1 ? _apply_noise_item!(first(items), dc) :
    _run_noise_item!(j - 1, Base.tail(items), dc)

# All of them, in order (serial loop).
@inline _run_noise_items!(::Tuple{}, ::Any) = nothing
@inline _run_noise_items!(box::Base.RefValue{<:Tuple}, dc) = _run_noise_items!(box[], dc)
@inline _run_noise_items!(items::Tuple, dc) =
    (_apply_noise_item!(first(items), dc); _run_noise_items!(Base.tail(items), dc))

# Band `B`'s measurement as a 0- or 1-tuple, so absence stays in the type domain;
# folds at compile time.
@inline _band_measurement_if_present(
    measurements::NamedTuple{K,<:Tuple{Vararg{BandMeasurement}}},
    ::Val{B},
) where {K,B} = B in K ? (measurements[B],) : ()

# An instance of signal `K` from any group, as a 0- or 1-tuple so not-found stays
# type-stable. Which group carries it does not matter. Folds at compile time.
@inline _first_signal_with_id(::Tuple{}, ::Val) = ()
@inline function _first_signal_with_id(
    groups::Tuple{SignalGroup,Vararg{SignalGroup}},
    v::Val,
)
    found = _first_signal_with_id(first(groups).signals, v)
    isempty(found) ? _first_signal_with_id(Base.tail(groups), v) : found
end
@inline function _first_signal_with_id(
    signals::Tuple{AbstractGNSSSignal,Vararg{AbstractGNSSSignal}},
    ::Val{K},
) where {K}
    s = first(signals)
    get_signal_id(typeof(s)) === K ? (s,) :
    _first_signal_with_id(Base.tail(signals), Val(K))
end

# Last sample (inclusive) of the current chunk; `nothing` means no chunking.
# Re-anchored to the absolute `chunk_index` so rounding never drifts and bands
# stay time-aligned.
@inline _chunk_last_sample(::Nothing, chunk_index, sampling_frequency, num_samples) =
    num_samples
@inline function _chunk_last_sample(
    chunk_duration,
    chunk_index,
    sampling_frequency,
    num_samples,
)
    grid = uconvert(NoUnits, (chunk_index + 1) * chunk_duration * sampling_frequency)
    min(round(Int, grid), num_samples)
end

# First sample (inclusive) of the chunk, where the noise measurement starts. The
# `nothing` method is needed: the generic one would yield `num_samples + 1`.
@inline _chunk_first_sample(::Nothing, chunk_index, sampling_frequency, num_samples) = 1
@inline _chunk_first_sample(chunk_duration, chunk_index, sampling_frequency, num_samples) =
    _chunk_last_sample(chunk_duration, chunk_index - 1, sampling_frequency, num_samples) + 1

# Optional per-backend sample-type check, run once per group and once per noise
# item. No-op by default; the integer backends reject non-`Complex{Int16}` samples
# with an `ArgumentError`.
@inline _check_sample_type(::AbstractDownconvertAndCorrelator, m) = nothing

# Per-group body shared by every backend: runs the group's loop on its band's
# `BandMeasurement`. `noise_items` / `n_noise` are non-empty for one group only
# (see `_dc_groups!`); an empty group still runs when it carries them.
@inline function _dc_one_group!(
    g::SignalGroup,
    dc::AbstractDownconvertAndCorrelator,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
    stop_before_partial::Bool,
    # The float backends derive nothing from the raw samples worth caching;
    # only the bit-wise backends act on `samples_unchanged`.
    samples_unchanged::Bool,
    noise_items::_NoiseWork,
    n_noise::Int,
)
    vals = g.satellites.values
    isempty(vals) && n_noise == 0 && return nothing
    m = measurements[get_band_id(g.band)]
    _check_sample_type(dc, m)
    num_samples = get_num_samples(m)
    chunk_last_sample =
        _chunk_last_sample(chunk_duration, chunk_index, m.sampling_frequency, num_samples)
    _dc_group_loop!(
        dc,
        vals,
        noise_items,
        n_noise,
        m.samples,
        num_samples,
        chunk_last_sample,
        m.sampling_frequency,
        m.intermediate_frequency,
        stop_before_partial,
    )
end

# Threading strategy for the per-sat group loop, declared per backend via
# `_threading`, so `_dc_group_loop!` has just two bodies shared by all backends.
struct _SerialLoop end
struct _BatchLoop end

# Serial by default; only the multi-threaded backends opt into `@batch`.
@inline _threading(::AbstractDownconvertAndCorrelator) = _SerialLoop()
@inline _threading(::CPUThreadedDownconvertAndCorrelator) = _BatchLoop()

@inline _dc_group_loop!(
    dc::AbstractDownconvertAndCorrelator,
    vals,
    noise_items,
    n_noise::Int,
    args::Vararg{Any,6},
) = _dc_group_loop!(_threading(dc), dc, vals, noise_items, n_noise, args...)

# Serial: noise despreads first, then the satellites.
@inline function _dc_group_loop!(
    ::_SerialLoop,
    dc,
    vals,
    noise_items,
    n_noise::Int,
    args::Vararg{Any,6},
)
    n_noise > 0 && _run_noise_items!(noise_items, dc)
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

# Threaded: the noise despreads are extra work items after the satellites' index
# range, each costing about one satellite's correlation, so `@batch` spreads them
# across idle threads. Thread-safe because each item is the sole writer of its
# estimator's window, totals and RNG (draw order, hence reproducibility, is
# unchanged), and replica scratch is per-thread.
@inline function _dc_group_loop!(
    ::_BatchLoop,
    dc,
    vals,
    noise_items,
    n_noise::Int,
    args::Vararg{Any,6},
)
    n = length(vals)
    @batch for i = 1:(n+n_noise)
        if i <= n
            @inbounds vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
        else
            _run_noise_item!(i - n, noise_items, dc)
        end
    end
    return nothing
end

"""
$(SIGNATURES)

Downconvert and correlate a single satellite on the CPU.
"""
function downconvert_and_correlate!(
    signal_type,
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
        signal_type,
        code_frequency,
        sampling_frequency,
        code_phase,
        signal_start_sample,
        num_samples_left,
        sample_shifts,
        prn,
    )
    _fused_standalone!(
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
