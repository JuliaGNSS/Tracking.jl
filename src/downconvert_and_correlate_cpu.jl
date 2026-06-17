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

CPUThreadedDownconvertAndCorrelator() =
    CPUThreadedDownconvertAndCorrelator([ScratchBuffers() for _ = 1:Threads.maxthreadid()],)

# Look up the active `ScratchBuffers` for this thread. Single-threaded
# backend has just one; threaded backend indexes by `Threads.threadid()`
# (stable under Polyester `@batch`). The index stays bounds-checked:
# `Threads.maxthreadid()` at construction time is not a lifetime bound —
# foreign threads adopted via `@ccallable`/`jl_adopt_thread` after the
# correlator was built get larger ids, which must fail loudly here
# rather than read out of bounds.
@inline _scratch_buffers(dc::CPUDownconvertAndCorrelator) = dc.buffers
@inline _scratch_buffers(dc::CPUThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

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
    bufs.extra_code_replicas[i-1]
end

# Multi-signal tile-share scratch handle: yields ScratchViews over the
# shared `tile_re` / `tile_im` SoA buffers, sized for `M * num_samples`
# Float32s each (M antennas laid out back-to-back).
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

# Per-signal correlation kernel, shared by both CPU backends: generate
# the code replica into the supplied `code_replica` buffer (from the
# calling correlator's `ScratchBuffers`), then run the fused kernel via
# `_fused_with_tile_scratch!` (which supplies SoA tile scratch when the
# shifts are a runtime-sized `AbstractVector`). Returns a value of the
# same correlator type as `correlator`.
#
# Args are per-signal scalars: the caller is responsible for picking the
# right `signal_type`, `correlator`, `code_phase` (modded into the signal's
# primary period), `code_frequency`, and `sample_shifts` for this signal.
@inline function _correlate_one_signal!(
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
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
    # Dynamic-shifts fallback: pull SoA tile buffers from the calling
    # correlator's per-thread `ScratchBuffers` (one re-half + one
    # im-half) so the fused kernel doesn't allocate them per call.
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

# Backend-less standalone fused dispatch for the public single-satellite
# `downconvert_and_correlate!`, which has no `ScratchBuffers` to draw on.
# Static (`SVector`) shifts use the allocation-free in-register kernel;
# dynamic (`AbstractVector`) shifts need SoA tile buffers, allocated per
# call here. This is a convenience entry point — the hot paths run through
# `_fused_with_tile_scratch!` with pooled scratch instead, so the dynamic
# fused kernel itself stays buffer-taking only (no allocating overload).
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
    sampling_frequency,
    intermediate_frequency,
)
    # MIN samples-to-next-boundary across all signals on this sat, clamped
    # to remaining buffer. Each signal's coherent-integration length comes
    # from its own `preferred_num_code_blocks_to_integrate` field; the
    # signal-level boundary calc reads it per signal.
    samples_to_integrate, per_signal_completed = _calc_min_samples_and_completed(
        sat.signals,
        sat.signal_start_sample,
        sampling_frequency,
        sat.code_doppler,
        sat.code_phase,
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
    num_samples_signal,
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
    signal_samples_left = num_samples_signal - signal_start_sample + 1
    samples_to_integrate = min(samples_to_integrate, signal_samples_left)
    per_signal_completed = _flag_completed(per_signal_to_boundary, samples_to_integrate)
    return samples_to_integrate, per_signal_completed
end

# For each signal, return its samples-to-next-primary-code-boundary using
# the signal-specific primary-code-relative phase. Tuple recursion keeps
# the heterogeneous walk inline / inference-friendly.
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
    # Chips-to-next-boundary must be measured against the same wrap the
    # multi-block window aligns to. For an `n_blocks > 1` window past sync it
    # must land on the secondary-/bit-period boundary (one NH10 period = one
    # L5I data symbol), so `_signal_replica_params` uses the secondary-aware
    # wrap; otherwise an N-block window started off the snapped secondary
    # phase straddles the data-symbol boundary and the coherent sum cancels
    # on data transitions. For a single block (or pre-sync) this is just the
    # primary code length.
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

# For each signal, `is_completed = (chosen_samples == samples_to_boundary)`.
@inline _flag_completed(::Tuple{}, _) = ()
@inline _flag_completed(t::Tuple, chosen) =
    (first(t) == chosen, _flag_completed(Base.tail(t), chosen)...)

# Code-phase wrap to use when generating a signal's code replica / sizing its
# integration window. Past sync the wrap is the full secondary-/bit-period
# length so that `gen_code!` bakes the secondary code at the correct chip for
# the current primary-code period: the replica's start phase carries the
# secondary-period offset, so the secondary sign is wiped per block right in
# the replica. This holds at ANY integration length — including the default
# `num_blocks == 1`. The previous code kept the N=1 replica primary-only and
# left the per-block secondary sign for the bit buffer to handle, but the
# post-sync bit decoder never did that wipe, so a `data × Σ(secondary signs)`
# collapse cost ~14 dB of decision margin at the default 1-block integration
# (issue #125). Wiping in the replica matches what `master` did via
# `update_code_phase` and is correct for multi-block windows too (otherwise the
# N blocks sum against a misaligned overlay and the coherent sum cancels).
# The boundary calc is unaffected at N=1: `calc_num_chips_to_integrate` re-mods
# the phase by the primary length, so the wider wrap yields the same per-block
# boundary. Before sync only the primary phase is known, so the wrap stays at
# the primary length. (Signals without a baked secondary code — e.g. GPS L1
# C/A — are unaffected: the primary repeats every `get_code_length` chips,
# secondary length 1.) `num_blocks` is retained for call-site symmetry.
@inline _replica_code_wrap(tsig::TrackedSignal, num_blocks::Integer) =
    has_bit_or_secondary_code_been_found(tsig) ? _post_sync_code_length(tsig) :
    get_code_length(tsig.signal)

# All per-signal replica/kernel parameters, derived in exactly one place
# so the boundary calc, replica sizing/generation, and kernel tap offsets
# can never diverge (issue #133): the code frequency (sat code Doppler +
# this signal's chip rate), the correlator tap `sample_shifts` for that
# frequency, the replica buffer size covering the worst-case tap spread,
# the coherent-integration length `n_blocks`, and the signal-relative
# `code_phase` modded by the secondary-aware `_replica_code_wrap`.
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
    p = _signal_replica_params(
        head,
        code_doppler,
        code_phase,
        sampling_frequency,
        num_samples_signal,
    )
    new_corr =
        _with_code_replica_buffer(dc, get_code_type(s), p.code_replica_size) do code_replica
            _correlate_one_signal!(
                dc,
                code_replica,
                s,
                head.correlator,
                signal,
                p.sample_shifts,
                p.signal_code_phase,
                carrier_phase,
                p.code_frequency,
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
    # Generate each signal's code replica into its per-thread slot, then
    # call the tuple kernel with N replicas + the shared tile buffers.
    # The replica views carry raw pointers into `dc`'s scratch byte
    # vectors, so root `dc` for their entire lifetime — generation
    # through the kernel call.
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
        p = _signal_replica_params(
            head,
            code_doppler,
            code_phase,
            sampling_frequency,
            num_samples_signal,
        )
        slot = _code_replica_slot(bufs, i)
        CT = get_code_type(s)
        nbytes = p.code_replica_size * sizeof(CT)
        length(slot) < nbytes && resize!(slot, nbytes)
        # The raw pointer in this view outlives this function — the caller
        # (`_correlate_signals`) wraps replica generation and the kernel in
        # `GC.@preserve dc`, which roots `slot` (reachable through the
        # correlator's ScratchBuffers) for the views' entire lifetime.
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

Downconvert and correlate all available satellites on the CPU (both the
single-threaded and the multi-threaded backend). Returns a new
`TrackState` whose slot *values* are detached from the input (the input's
per-sat tracking values are left untouched), but whose key set
(`Indices`) is *shared* with the input — this step never changes the key
set, so sharing avoids copying the hash table on every `track` loop
iteration. Detaching the key set happens once at the `track` boundary
(`reset_start_sample_and_bit_buffer`); do not
`add_satellite!`/`remove_satellite!` on this function's direct output, or
you will corrupt the input's keys (#123). The copy is otherwise shallow:
per-sat scratch vectors are shared, see [`track`](@ref). Per-call
code-replica scratch comes from the correlator's (per-thread)
`ScratchBuffers` so the kernel itself stays allocation-free; the only
per-call allocation is the slot-value copy.
"""
function downconvert_and_correlate(
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    measurements::BandMeasurements,
    track_state::TrackState,
)
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    downconvert_and_correlate!(dc, measurements, new_track_state)
end

"""
$(SIGNATURES)

In-place version: writes new `TrackedSat` values directly into each
group's existing `Vector{TrackedSat}` backing storage. On the threaded
backend, different `@batch` iterations write to disjoint slots, so no
synchronization is needed. Returns the same `track_state`.
Allocation-free in steady state — see [`track!`](@ref).
"""
function downconvert_and_correlate!(
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    measurements::BandMeasurements,
    track_state::TrackState,
)
    _foreach_group!(_dc_one_group!, track_state.groups, dc, measurements)
    return track_state
end

# Per-group body shared by both CPU backends. Pulled out so
# `_foreach_group!` can call it on each `SignalGroup` in the (possibly
# heterogeneous) `groups` tuple without dynamic dispatch / boxing.
# Routes to this group's band's `BandMeasurement` for the signal buffer and
# front-end metadata; the serial-vs-`@batch` loop choice is dispatched
# per backend in `_dc_group_loop!`.
@inline function _dc_one_group!(
    g::SignalGroup,
    dc::Union{CPUDownconvertAndCorrelator,CPUThreadedDownconvertAndCorrelator},
    measurements::BandMeasurements,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    m = measurements[band_key(g.band)]
    _dc_group_loop!(
        dc,
        vals,
        m.samples,
        get_num_samples(m),
        m.sampling_frequency,
        m.intermediate_frequency,
    )
end

@inline function _dc_group_loop!(dc::CPUDownconvertAndCorrelator, vals, args::Vararg{Any,4})
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

@inline function _dc_group_loop!(
    dc::CPUThreadedDownconvertAndCorrelator,
    vals,
    args::Vararg{Any,4},
)
    @batch for i = 1:length(vals)
        @inbounds vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
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
