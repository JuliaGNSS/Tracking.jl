# Scratch byte buffers, one per role: `code_replica` for a sat's first signal,
# `extra_code_replicas` (grown lazily) for signals 2..N, and `tile_re` / `tile_im`
# for the `AbstractVector`-shifts fused kernel and the multi-signal tile-share
# kernel. The single-threaded backend holds one, the threaded one one per thread.
# The buffers are only resized in place, so the struct can be immutable.
struct ScratchBuffers
    code_replica::Vector{UInt8}
    extra_code_replicas::Vector{Vector{UInt8}}
    tile_re::Vector{UInt8}
    tile_im::Vector{UInt8}
end

ScratchBuffers() = ScratchBuffers(UInt8[], Vector{UInt8}[], UInt8[], UInt8[])

# Isbits typed view over a `Vector{UInt8}` slot, free to build. Implements just
# enough of `AbstractArray` for the kernels: linear indexing, `size`, and
# `pointer` (for SIMD `vstore`).
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

CPU-based implementation of downconversion and correlation. Its scratch buffers
grow lazily on first use and are reused thereafter, so a hoisted instance is
allocation-free per `track!` call in steady state.

For real-time loops, construct the correlator **once outside** the `track!` loop and
pass it via the `downconvert_and_correlator` keyword argument; the default value
rebuilds the buffers on every call.
"""
struct CPUDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::ScratchBuffers
end

CPUDownconvertAndCorrelator() = CPUDownconvertAndCorrelator(ScratchBuffers())

"""
$(SIGNATURES)

Multi-threaded CPU downconvert and correlate, parallelized over the satellites
(PRNs) of each group. Holds one set of scratch buffers per thread, indexed by
`Threads.threadid()` inside `@batch` (which pins each iteration to a thread), grown
lazily and reused thereafter.

With more than one thread and more than one work item, each `@batch` launch (one per
`downconvert_and_correlate!` call, so one per code block inside `track!`) keeps a
Polyester allocation of ~96 B. It is bounded per launch and does not scale with the
input buffer. Measuring a noise density keeps it at ~96 B (the noise despreads ride
the loop as extra items behind one pointer, see `_park_noise_items!`), but a
single-satellite state then pays it too, as its loop has two items; in exchange the
despread usually costs a fraction of its ~1.75 µs serial time.

Polyester copies everything the `@batch` region references into one argument tuple
per launch. A tuple of only `isbits` members (bare `Array`s qualify, as `PtrArray`s)
is elided; otherwise it is heap-allocated at its `sizeof`, with immutable structs
copied by value and mutable objects as one pointer. The non-`isbits` member here is
the GNSS signal (`GPSL1CA` wraps a `Matrix{Int16}` and a non-`isbits` `SignalLUT`),
which code-replica generation needs. A prototype that generates from the LUT's bare
arrays passed as separate arguments measures 0 B/launch at full throughput, but needs
a GNSSSignals API that does not re-wrap them in a struct; a serial code-gen pre-pass
also reaches 0 but serializes ~30 % of the work. Neither is worth ~96 B per launch.

For real-time loops, construct the correlator **once outside** the `track!` loop and
pass it via the `downconvert_and_correlator` keyword argument; the default value
rebuilds the per-thread buffers on every call.
"""
struct CPUThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{ScratchBuffers}
end

CPUThreadedDownconvertAndCorrelator() =
    CPUThreadedDownconvertAndCorrelator([ScratchBuffers() for _ = 1:Threads.maxthreadid()],)

# This thread's `ScratchBuffers`. Bounds-checked: threads adopted after
# construction (`jl_adopt_thread`) can exceed the `Threads.maxthreadid()` it was
# sized for, and must fail loudly.
@inline _scratch_buffers(dc::CPUDownconvertAndCorrelator) = dc.buffers
@inline _scratch_buffers(dc::CPUThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

# Grow `buf` to fit `n` elements of `T` and hand `f` a `ScratchView` of exactly
# that size, keeping `buf` rooted while `f` runs (the view holds a raw `Ptr`).
@inline function _with_scratch_view(f, buf::Vector{UInt8}, ::Type{T}, n::Int) where {T}
    nbytes = n * sizeof(T)
    length(buf) < nbytes && resize!(buf, nbytes)
    GC.@preserve buf begin
        f(ScratchView{T}(Ptr{T}(pointer(buf)), n))
    end
end

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

# Views over the `tile_re` / `tile_im` SoA buffers, `num_ants * num_samples`
# Float32s each (antennas back-to-back).
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
# of samples, accumulating into (and returning) `correlator`. Every despread goes
# through it, the per-satellite single-signal path (`_correlate_signals`) and the
# open-loop noise reference (`update_noise!`) alike. That makes the noise measurement
# model-free by construction: the reference runs the prompt's own quantise →
# downconvert → despread path, so its `N₀` carries the quantisation loss and code
# amplitude, and is on the same scale even where the one- and two-bit accumulators
# are popcount counts rather than sample sums.
#
# `carrier_phase` is where the caller's carrier stands (`0.0` for the open-loop
# reference). `code_replica_size` is the replica the backend needs if it generates
# one (the CPU kernels do; the Int16 and bit-wise kernels ignore it): the whole sample
# buffer plus the tap spread, since `gen_code_replica!` writes *at* `start_sample`.
# The CPU method draws it from the backend's per-thread `ScratchBuffers`, shared by
# satellites and noise reference: each loop iteration, serial or `@batch`-pinned to a
# thread, uses the slot alone.
#
# `use_band_cache = false` keeps a bit-wise backend from reading the band sign planes
# packed once per group: the noise reference passes it, as the planes belong to the
# group whose loop it rides, possibly on another band. Other backends ignore it.
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
    # GNSSSignals' embedded-LUT `gen_code!` is Int8-only; the fused kernel widens
    # it to Float32. CBOC (Galileo E1B) is the Int8 integer-amplitude approximation.
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

# Fused kernel call with SoA tile buffers from the backend's `ScratchBuffers`, which
# only the `AbstractVector`-shifts fallback needs; `SVector` shifts need none.
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

# The same for the public single-satellite `downconvert_and_correlate!`, which has no
# `ScratchBuffers`: `AbstractVector` shifts allocate their tile buffers per call.
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

# Per-sat downconvert+correlate for every backend; returns the updated TrackedSat.
# The kernel is chosen in `_correlate_signals` / `_despread_one_signal!`.
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
    # Each sub-step runs to the nearest code-block boundary of any signal (or the
    # chunk end), so a chunk yields 0, 1 or several outputs per signal; a partial at
    # the chunk end carries in the accumulator. The NCO Doppler is untouched here.
    #
    # `stop_before_partial = true` leaves the trailing partial to the next chunk's
    # pass, which then integrates boundary to boundary at the Doppler the estimator
    # wrote in between: every completed integration sees a single Doppler. A final
    # pass without the flag drains the buffer's trailing partial.
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
        # A sub-step that completes no signal is the chunk-clamped trailing partial.
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

# The MIN over all signals' samples-to-next-boundary, clamped to the chunk end
# (boundaries and replica sizing still see the whole buffer), and which signals
# that completes.
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
    samples_left = chunk_last_sample - signal_start_sample + 1
    samples_to_integrate = min(samples_to_integrate, samples_left)
    per_signal_completed = _flag_completed(per_signal_to_boundary, samples_to_integrate)
    return samples_to_integrate, per_signal_completed
end

# Each signal's samples to the end of its integration window.
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
    # Measured against the secondary-aware wrap (`_replica_code_wrap`), so a
    # multi-block window past sync ends on a data-symbol boundary instead of
    # straddling one.
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

# Code-phase wrap for a signal's replica and integration window. Past sync it is
# the secondary-/bit-period length, so `gen_code!` bakes the secondary chip of the
# current block into the replica and wipes the overlay at every integration length
# (#125); before sync only the primary phase is known. At one block the boundary is
# unchanged, as `calc_num_chips_to_integrate` re-mods by the primary length.
# `num_blocks` is unused.
@inline _replica_code_wrap(tsig::TrackedSignal, num_blocks::Integer) =
    has_bit_or_secondary_code_been_found(tsig) ? _post_sync_code_length(tsig) :
    get_code_length(tsig.signal)

# All per-signal replica/kernel parameters, derived in one place so the boundary
# calc, replica generation and kernel tap offsets cannot diverge (#133).
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

# Single-signal path for every backend (the per-backend part is
# `_despread_one_signal!`). Returns `((new_correlator, completed),)`.
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

# Multi-signal path (N >= 2): N code replicas, then the tile-share kernel
# `downconvert_and_correlate_fused_tuple!`. Per-backend, as the tile-share kernels
# have no single-signal form.
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
    # The replica views hold raw pointers into `dc`'s scratch, so root `dc`.
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
        # All correlators of a sat share M (enforced by SignalGroup).
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

# Each signal's code replica in its scratch slot, as a tuple of `ScratchView`s. An
# enumerated `map`, as recursive splatting stops inlining around N=3 and boxes.
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
        # Outlives this function; the caller's `GC.@preserve dc` roots `slot`.
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

Downconvert and correlate all available satellites, for every backend (a backend
customises `_despread_one_signal!`, `_dc_one_group!`, `_check_sample_type` and
`_threading`). Returns a new `TrackState` whose slot *values* are detached from the
input but whose key set (`Indices`) is *shared* with it, which saves copying the hash
table on every `track` iteration; `track` detaches the keys once. Do not
`add_satellite!`/`remove_satellite!` on this function's direct output, or you corrupt
the input's keys (#123). The copy is otherwise shallow: per-sat scratch vectors are
shared, see [`track`](@ref). The only per-call allocation is the slot-value copy.

The **noise estimators are shared too**, and this call advances them: a
[`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator)'s window
and RNG stream are written in place, so the input's window grows and its draws
advance. Copying a window (~1000 observations per signal) per call would cost O(K)
and reset the reference's PRN rotation. Hence:

  - Two outputs branched from one input share its noise reference: the same samples
    enter one window twice and every branch divides by the same figure.
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

In-place version: writes the new `TrackedSat` values into each group's existing
storage (threaded `@batch` iterations write disjoint slots) and returns
`track_state`. Allocation-free in steady state, see [`track!`](@ref).

`chunk_duration` and `chunk_index` restrict the pass to one chunk of a fixed per-band
time grid: each satellite integrates up to sample
`min(round((chunk_index + 1) * chunk_duration * sampling_frequency), num_samples)`,
re-anchored to the absolute `chunk_index` so rounding never drifts and bands stay
aligned. `chunk_duration = nothing` (the default) consumes the whole buffer. Every
completed integration is snapshotted into its signal's `correlator_outputs`; only the
trailing partial stays in the live correlator.

`stop_before_partial = true` stops each satellite at its last code-block boundary
inside the chunk, leaving the trailing partial to the next chunk's pass at the
Doppler the estimator writes in between (as `track!` does). Default `false`.

`samples_unchanged = true` promises every band's buffer holds the same content as on
the previous call with this `dc`, so backends may reuse sample-derived caches (the
bit-wise backends skip re-packing their band sign planes). `track!` passes it on every
chunk after the first; keep the default `false` whenever a buffer may have been
refilled.

`measure_noise = false` (default `true`) leaves every
[`CorrelatorNoiseEstimator`](@extref TrackingLoops.CorrelatorNoiseEstimator)'s window
untouched. Pass it on a call that re-covers samples already measured, as `track!`'s
final drain pass does, so they do not enter the window twice.

A band **absent** from `measurements` is not measured and does not throw, so a
multi-band `TrackState` may be advanced one band at a time; groups with satellites on
that band still require its measurement.

The noise estimators are mutable state shared with the input on the out-of-place
`downconvert_and_correlate` and [`track`](@ref), see there.
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
    # `measure_noise = false` skips resolving the noise work, not just running it.
    # Branching here, not in `_noise_items`, keeps each arm's `_NoiseWork` concretely
    # typed; the empty arm is the `()` specialisation the group tail compiles anyway.
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

# Walk the groups, handing the chunk's noise work to the first. The items carry
# their own band measurement and signal; they only ride that group's loop, so on a
# threaded backend the despread runs alongside the satellite correlations. Handing
# them to one group (not every `_dc_one_group!`) despreads each signal exactly once
# however the satellites are grouped. A state with no groups has no noise estimators.
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
# `noise_estimators`, so a signal two groups carry is measured once. Resolved on the
# caller's thread, so `_check_sample_type` still throws its `ArgumentError` before
# any kernel, and parked behind one pointer for the `@batch` closure (see
# `_park_noise_items!`).
#
# Walked by tuple recursion, as estimator types differ per signal and a runtime loop
# would allocate. Only type-domain decisions are made here (a signal no group tracks,
# a band without samples), so the tuple's length depends on types alone; an
# out-of-range chunk is left to `update_noise!`'s `num_samples > 0` guard.
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

# Park the descriptors in the `TrackState`'s reusable cell and return the box. The
# `@batch` argument tuple copies an immutable descriptor by value (152 B for a
# `CorrelatorNoiseEstimator`, a `BandMeasurement` and a `GPSL1CA`, which raised a
# launch from ~96 B to 240 B) but a mutable box as one 8 B pointer.
#
# The box is reused across launches. Its type names the caller's `BandMeasurement`,
# unknown until the first call, so the cell is a `Base.RefValue{Any}` holding a
# concretely typed box, replaced only when the sample type changes. The box roots
# the descriptors, and so keeps the last chunk's sample buffer reachable until the
# next call. Nothing to measure stays `()`, adding nothing to the launch.
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

# First call with this descriptor type; `@noinline` keeps the hot path small.
@noinline function _new_noise_items_box!(
    cell::Base.RefValue{Any},
    items::T,
) where {T<:Tuple}
    box = Base.RefValue{T}(items)
    cell[] = box
    box
end

# Work items to add to the satellites' range; folds to a literal.
@inline _num_noise_items(::Tuple{}) = 0
@inline _num_noise_items(::Base.RefValue{<:Tuple{Vararg{Any,N}}}) where {N} = N

# A group loop's noise work: `()` (nothing to measure, or not the carrying group)
# or the box the descriptors are parked in.
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

# One signal's descriptor as a 0- or 1-tuple, keeping "nothing to measure" in the
# type domain.
@inline function _noise_item(
    estimator::AbstractNoiseEstimator,
    ::Val{K},
    dc::AbstractDownconvertAndCorrelator,
    groups::SignalGroups,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
) where {K}
    # A signal no group tracks measures nothing: the correlator-ingest case, whose
    # window `append_noise_observation!` fills.
    found = _first_signal_with_id(Tuple(groups), Val(K))
    isempty(found) && return ()
    signal = first(found)
    # The samples are per band; only the despreading code is per signal. A band
    # this call brings no samples for measures nothing, as a multi-band state may be
    # advanced one band at a time; a group with satellites there still throws in
    # `_dc_one_group!`.
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

# Builds the `NoiseUpdateContext` here, as the loop already captures the backend.
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

# Run the `j`-th descriptor: an unrolled walk, as the tuple is heterogeneous.
@inline _run_noise_item!(::Int, ::Tuple{}, ::Any) = nothing
@inline _run_noise_item!(j::Int, box::Base.RefValue{<:Tuple}, dc) =
    _run_noise_item!(j, box[], dc)
@inline _run_noise_item!(j::Int, items::Tuple, dc) =
    j == 1 ? _apply_noise_item!(first(items), dc) :
    _run_noise_item!(j - 1, Base.tail(items), dc)

# All of them, in order: the serial backends' path.
@inline _run_noise_items!(::Tuple{}, ::Any) = nothing
@inline _run_noise_items!(box::Base.RefValue{<:Tuple}, dc) = _run_noise_items!(box[], dc)
@inline _run_noise_items!(items::Tuple, dc) =
    (_apply_noise_item!(first(items), dc); _run_noise_items!(Base.tail(items), dc))

# Band `B`'s measurement as a 0- or 1-tuple; folds at compile time.
@inline _band_measurement_if_present(
    measurements::NamedTuple{K,<:Tuple{Vararg{BandMeasurement}}},
    ::Val{B},
) where {K,B} = B in K ? (measurements[B],) : ()

# An instance of signal `K` from any group, as a 0- or 1-tuple; folds at compile
# time.
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

# Last sample (inclusive) of the chunk on the per-band time grid, anchored to the
# absolute `chunk_index` so rounding never drifts; `nothing` is the whole buffer.
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

# First sample (inclusive) of the chunk, where the noise measurement starts. Not
# derived from `_chunk_last_sample(…, chunk_index - 1, …) + 1` without chunking, which
# would be `num_samples + 1` and measure nothing.
@inline _chunk_first_sample(::Nothing, chunk_index, sampling_frequency, num_samples) = 1
@inline _chunk_first_sample(chunk_duration, chunk_index, sampling_frequency, num_samples) =
    _chunk_last_sample(chunk_duration, chunk_index - 1, sampling_frequency, num_samples) + 1

# Per-backend sample-type check before any kernel runs; the integer backends reject
# non-`Complex{Int16}` buffers with an `ArgumentError`.
@inline _check_sample_type(::AbstractDownconvertAndCorrelator, m) = nothing

# Per-group body shared by every backend. `noise_items` / `n_noise` are the chunk's
# noise despreads, non-empty for one group (see `_dc_groups!`); they carry their own
# band measurement, and an empty group still runs to carry them.
@inline function _dc_one_group!(
    g::SignalGroup,
    dc::AbstractDownconvertAndCorrelator,
    measurements::BandMeasurements,
    chunk_index::Int,
    chunk_duration,
    stop_before_partial::Bool,
    # Only the bit-wise backends act on it.
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

# Threading strategy for the group loop, declared per backend via `_threading`.
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

# Serial: the noise despreads run first.
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

# Threaded: the noise despreads are extra items at the end of the satellites' range.
# One costs about what a satellite's correlation does, and nearly nothing in wall
# time where the satellites leave a thread idle. Thread-safe as each item is the only
# writer of its estimator's window and RNG (so seeded runs stay reproducible), and
# replica scratch is per thread.
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
