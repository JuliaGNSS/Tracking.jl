using BenchmarkTools
using GNSSSignals
using GNSSSignals: GalileoE1B
using Unitful: Hz
using Tracking
using Tracking:
    EarlyPromptLateCorrelator,
    get_correlator_sample_shifts,
    get_code_type,
    NumAnts,
    gen_code_replica!,
    TrackState,
    downconvert_and_correlate,
    BitBuffer
using StaticArrays

# GNSSSignals v1 → v2 signal-name shims. AirspeedVelocity uses HEAD's
# `benchmark/benchmarks.jl` to bench *every* rev (via `--bench-on`), so
# this file must work against both v1 (loaded for master's Tracking
# 1.5.0) and v2 (loaded for HEAD's Tracking 2.0.0). Alias the v2 names
# to v1's symbols when v2 isn't loaded.
const GPSL1CA = isdefined(GNSSSignals, :GPSL1CA) ?
    getfield(GNSSSignals, :GPSL1CA) : getfield(GNSSSignals, :GPSL1)
const GPSL5I = isdefined(GNSSSignals, :GPSL5I) ?
    getfield(GNSSSignals, :GPSL5I) : getfield(GNSSSignals, :GPSL5)

# Code-replica buffer element type. GNSSSignals' embedded-LUT `gen_code!` (PR #90)
# is Int8-only; older GNSSSignals emit `get_code_type(signal)` (Int16/Float32).
# Detect the era so this script builds against both (AirspeedVelocity diffs revs).
_code_buf_type(sig) = isdefined(GNSSSignals, :code_engine) ? Int8 : get_code_type(sig)

# Branch-portable per-system storage construction. Three flavours coexist:
#   * Master:           `SystemSatsState(sys, ::Vector{SatState})`
#   * Wrapper branch:   `TrackedSystem(estimator, sys, ::Vector{SatState})`
#   * Multi-signal:     plain `Dictionary{Int, TrackedSat}` (no per-system
#                       container struct at all)
# Detect at load time so the script runs unchanged across all three via
# AirspeedVelocity's `--bench-on=$HEAD_SHA`.
const _HAS_TRACKED_SIGNAL = isdefined(Tracking, :TrackedSignal)
const _HAS_TRACKED_SAT    = isdefined(Tracking, :TrackedSat)
const _HAS_SAT_STATE      = isdefined(Tracking, :SatState)

if !_HAS_TRACKED_SIGNAL
    const _TrackedSystem = isdefined(Tracking, :TrackedSystem) ?
        Tracking.TrackedSystem : Tracking.SystemSatsState
end

# Field / property name on `TrackState` for the per-system tuple.
# Three flavours of branch coexist:
#   * Master:                       real field `multiple_system_sats_state`
#   * Wrapper branch:               real field `satellites`
#   * Multi-band branch (current):  real field `groups`, where each entry is a
#                                   `SignalGroup` bundling a `satellites` dict
#                                   (plus band / signals / num_ants).
# On the multi-band branch the bench needs a NamedTuple-of-dicts shape — derive
# it inline from `groups` and rebuild via fresh `SignalGroup`s after the map.
const _MULTIBAND = _HAS_TRACKED_SIGNAL
const _SYSTEMS_FIELD = _MULTIBAND ? :groups :
    isdefined(Tracking, :TrackedSystem) ? :satellites :
    :multiple_system_sats_state

@inline _get_systems(track_state) =
    _MULTIBAND ? map(g -> g.satellites, track_state.groups) :
                 getproperty(track_state, _SYSTEMS_FIELD)

@inline function _track_state_with_systems(track_state, systems)
    if _MULTIBAND
        new_groups = map(track_state.groups, systems) do g, sats
            Tracking.SignalGroup(g; satellites = sats)
        end
        TrackState(track_state; groups = new_groups)
    else
        TrackState(track_state; NamedTuple{(_SYSTEMS_FIELD,)}((systems,))...)
    end
end

# Branch-portable per-sat construction. Returns an object suitable for
# whatever the loaded Tracking expects in its per-system storage:
#   * Master / wrapper branch: a `SatState`.
#   * Multi-signal branch: a `TrackedSat` (with one TrackedSignal inside,
#     and the doppler_estimator_state already seeded).
function _make_initial_sat(sys, prn, code_phase, carrier_doppler; num_ants = NumAnts(1))
    if _HAS_TRACKED_SIGNAL
        return Tracking.TrackedSat(
            sys, prn, code_phase, carrier_doppler;
            doppler_estimator = Tracking.ConventionalAssistedPLLAndDLL(),
            num_ants,
        )
    else
        return Tracking.SatState(sys, prn, code_phase, carrier_doppler; num_ants)
    end
end

_make_initial_sat_with_num_ants(sys, prn, code_phase, carrier_doppler, num_ants) =
    _make_initial_sat(sys, prn, code_phase, carrier_doppler; num_ants)

const SUITE = BenchmarkGroup()

# ── Helper: set up common benchmark state ──────────────────────────────────

function setup_benchmark(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    gnss_signal = GPSL1CA(),
    num_ants = 1,
)
    code_phase = 10.5
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(gnss_signal)
    code_frequency = code_doppler + get_code_frequency(gnss_signal)

    correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants))
    static_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)

    signal =
        num_ants == 1 ? rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    code_replica = Vector{_code_buf_type(gnss_signal)}(
        undef,
        num_samples + maximum(static_shifts) - minimum(static_shifts),
    )
    _gen_code_replica!(
        code_replica,
        gnss_signal,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples,
        static_shifts,
        1,
    )

    return (;
        correlator,
        signal,
        code_replica,
        static_shifts,
        dynamic_shifts,
        sampling_frequency,
        carrier_doppler,
        code_phase,
        gnss_signal,
        num_samples,
    )
end

# ── Branch-portable constructors ──────────────────────────────────────────
# Tracking ≥ 2.0 dropped the `Val{MESF}` arg from `CPUDownconvertAndCorrelator`,
# `CPUThreadedDownconvertAndCorrelator`, and `gen_code_replica!`. Detect at
# load time and dispatch.

const _NEEDS_VAL = !applicable(Tracking.CPUDownconvertAndCorrelator)

_make_cpu_dc(sampling_frequency) =
    _NEEDS_VAL ? Tracking.CPUDownconvertAndCorrelator(Val(sampling_frequency)) :
    Tracking.CPUDownconvertAndCorrelator()

_make_cpu_threaded_dc(sampling_frequency) =
    _NEEDS_VAL ? Tracking.CPUThreadedDownconvertAndCorrelator(Val(sampling_frequency)) :
    Tracking.CPUThreadedDownconvertAndCorrelator()

function _gen_code_replica!(
    code_replica, gnss_signal, code_frequency, sampling_frequency,
    code_phase, start_sample, num_samples, sample_shifts, prn,
)
    if _NEEDS_VAL
        gen_code_replica!(
            code_replica, gnss_signal, code_frequency, sampling_frequency,
            code_phase, start_sample, num_samples, sample_shifts, prn,
            Val(sampling_frequency),
        )
    else
        gen_code_replica!(
            code_replica, gnss_signal, code_frequency, sampling_frequency,
            code_phase, start_sample, num_samples, sample_shifts, prn,
        )
    end
end

# ── High-level downconvert_and_correlate (full pipeline) ───────────────────

function bench_downconvert_and_correlate(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    gnss_signal = GPSL1CA(),
    num_ants = 1,
)
    downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
    track_state = TrackState(
        gnss_signal,
        [_make_initial_sat_with_num_ants(gnss_signal, 1, 10.5, 1000.0Hz, NumAnts(num_ants))],
    )
    signal =
        num_ants == 1 ? rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    # Multi-band branch uses `BandMeasurement(samples, fs, if)` (named
    # `Measurement` before the issue-#132 rename). The
    # coherent-integration length used to be a trailing argument
    # (`downconvert_and_correlate(dc, measurements, ts, prefer)`) but moved
    # to a per-signal `TrackedSignal` field, leaving a 3-arg call. Pick the
    # name and arity that exist so the script benches both shapes across revs.
    # Master / wrapper branches use the legacy 6-arg form.
    @static if isdefined(Tracking, :BandMeasurement) || isdefined(Tracking, :Measurement)
        measurement_type =
            isdefined(Tracking, :BandMeasurement) ? Tracking.BandMeasurement :
            Tracking.Measurement
        # Measurement NamedTuple key = the running revision's band id. benchpkg
        # runs HEAD's script against every rev, and the key spelling is a
        # breaking change: Tracking ≤ 2.3 keyed by its own `Tracking.band_key`
        # (`:l1`), HEAD keys by `GNSSSignals.get_band_id` (`:L1`). Derive it from
        # whichever accessor the loaded Tracking has so both revs resolve.
        band = get_band(gnss_signal)
        @static if isdefined(Tracking, :band_key)
            band_id = Tracking.band_key(band)
        else
            band_id = get_band_id(band)
        end
        measurements =
            NamedTuple{(band_id,)}((measurement_type(signal, sampling_frequency, 0.0Hz),))
        if hasmethod(
            Tracking.downconvert_and_correlate,
            Tuple{typeof(downconvert_and_correlator),typeof(measurements),typeof(track_state)},
        )
            @benchmarkable Tracking.downconvert_and_correlate(
                $downconvert_and_correlator,
                $measurements,
                $track_state,
            )
        else
            @benchmarkable Tracking.downconvert_and_correlate(
                $downconvert_and_correlator,
                $measurements,
                $track_state,
                1,
            )
        end
    else
        @benchmarkable Tracking.downconvert_and_correlate(
            $downconvert_and_correlator,
            $signal,
            $track_state,
            1,
            $sampling_frequency,
            $(0.0Hz),
        )
    end
end

# ── Fused kernel microbenchmarks ───────────────────────────────────────────

function bench_fused_kernel(;
    signal_type = Float32,
    num_samples = 2000,
    num_ants = 1,
    shifts = :static,
)
    s = setup_benchmark(; signal_type, num_samples, num_ants)
    sample_shifts = shifts == :static ? s.static_shifts : s.dynamic_shifts
    # The dynamic-shifts fallback now requires caller-supplied SoA tile
    # buffers (`tile_re`/`tile_im`). Master's signature is the old 9-arg
    # form (Bumper allocates the tiles inside the kernel); the new
    # branch requires 11 args. Detect via `hasmethod`.
    static_arg_types = Tuple{
        typeof(s.correlator), typeof(s.signal), typeof(s.code_replica),
        typeof(sample_shifts), typeof(s.carrier_doppler),
        typeof(s.sampling_frequency), Float64, Int, Int,
    }
    if shifts == :dynamic && !hasmethod(Tracking.downconvert_and_correlate_fused!, static_arg_types)
        # New branch: 11-arg signature with hoisted tile buffers
        tile_re = Vector{Float32}(undef, num_samples * num_ants)
        tile_im = Vector{Float32}(undef, num_samples * num_ants)
        Tracking.downconvert_and_correlate_fused!(
            s.correlator, s.signal, s.code_replica, sample_shifts,
            s.carrier_doppler, s.sampling_frequency, 0.0,
            1, s.num_samples, tile_re, tile_im,
        )
        return @benchmarkable Tracking.downconvert_and_correlate_fused!(
            $(s.correlator), $(s.signal), $(s.code_replica), $sample_shifts,
            $(s.carrier_doppler), $(s.sampling_frequency), 0.0,
            1, $(s.num_samples), $tile_re, $tile_im,
        )
    end
    # Master path or static-shifts overload: 9-arg signature.
    Tracking.downconvert_and_correlate_fused!(
        s.correlator,
        s.signal,
        s.code_replica,
        sample_shifts,
        s.carrier_doppler,
        s.sampling_frequency,
        0.0,
        1,
        s.num_samples,
    )
    @benchmarkable Tracking.downconvert_and_correlate_fused!(
        $(s.correlator),
        $(s.signal),
        $(s.code_replica),
        $sample_shifts,
        $(s.carrier_doppler),
        $(s.sampling_frequency),
        0.0,
        1,
        $(s.num_samples),
    )
end

# ── Register benchmarks ───────────────────────────────────────────────────

# Full pipeline: CPU, various signal types
foreach((Int16, Int32, Float32, Float64)) do signal_type
    SUITE["downconvert and correlate"]["CPU"][string(signal_type)] =
        bench_downconvert_and_correlate(; signal_type)
end

# Full pipeline: multi-antenna
SUITE["downconvert and correlate"]["CPU"]["Float32 4ant"] =
    bench_downconvert_and_correlate(; num_ants = 4)
SUITE["downconvert and correlate"]["CPU"]["Int16 4ant"] =
    bench_downconvert_and_correlate(; signal_type = Int16, num_ants = 4)

# Full pipeline: track
function bench_track(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
)
    gnss_signal = GPSL1CA()
    downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
    track_state = TrackState(gnss_signal, [_make_initial_sat(gnss_signal, 1, 0.0, 1000Hz)])
    signal = rand(Complex{signal_type}, num_samples)
    @benchmarkable track(
        $signal,
        $track_state,
        $sampling_frequency;
        downconvert_and_correlator = $downconvert_and_correlator,
    )
end
SUITE["track"]["1. Float32/2K – track"] = bench_track()

# In-place track! (only on branches that define it). Mirrors bench_track so
# the comparison report shows them side by side.
function bench_track_inplace(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
)
    gnss_signal = GPSL1CA()
    downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
    track_state = TrackState(gnss_signal, [_make_initial_sat(gnss_signal, 1, 0.0, 1000Hz)])
    signal = rand(Complex{signal_type}, num_samples)
    @benchmarkable Tracking.track!(
        $signal,
        $track_state,
        $sampling_frequency;
        downconvert_and_correlator = $downconvert_and_correlator,
    )
end
if isdefined(Tracking, :track!)
    SUITE["track"]["1. Float32/2K – track!"] = bench_track_inplace()
end

# Fused kernel microbenchmarks (only available on branches with the fused kernel)
if isdefined(Tracking, :downconvert_and_correlate_fused!)
    SUITE["fused kernel"]["1-ant static taps"] = bench_fused_kernel(; shifts = :static)
    SUITE["fused kernel"]["1-ant dynamic taps"] = bench_fused_kernel(; shifts = :dynamic)
    SUITE["fused kernel"]["4-ant static taps"] =
        bench_fused_kernel(; num_ants = 4, shifts = :static)
    SUITE["fused kernel"]["4-ant dynamic taps"] =
        bench_fused_kernel(; num_ants = 4, shifts = :dynamic)
end

# ── Tuple-kernel microbenchmark (multi-signal tile-share path) ─────────
# Exercises `downconvert_and_correlate_fused_tuple!` directly: one
# downconvert into the SoA tile, one sample-outer correlate pass over
# all N×NC accumulators. Covers single-antenna (the existing hot path
# for multi-signal-per-sat) and multi-antenna (the extended path), so a
# PR comparison surfaces regressions on either side.
function bench_fused_tuple_kernel(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    gnss_signal = GPSL1CA(),
    num_ants = 1,
    n_signals = 2,
)
    code_phase = 10.5
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(gnss_signal)
    code_frequency = code_doppler + get_code_frequency(gnss_signal)

    correlator_template = EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants))
    sample_shifts =
        get_correlator_sample_shifts(correlator_template, sampling_frequency, code_frequency)
    code_replica_size =
        num_samples + maximum(sample_shifts) - minimum(sample_shifts)

    signal = num_ants == 1 ?
        rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    correlators = ntuple(
        _ -> EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants)),
        n_signals,
    )
    code_replicas = ntuple(n_signals) do _
        cr = Vector{_code_buf_type(gnss_signal)}(undef, code_replica_size)
        _gen_code_replica!(
            cr, gnss_signal, code_frequency, sampling_frequency, code_phase,
            1, num_samples, sample_shifts, 1,
        )
        cr
    end
    sample_shifts_tuple = ntuple(_ -> sample_shifts, n_signals)
    tile_re = Vector{Float32}(undef, num_samples * num_ants)
    tile_im = Vector{Float32}(undef, num_samples * num_ants)

    @benchmarkable Tracking.downconvert_and_correlate_fused_tuple!(
        $correlators, $signal, $code_replicas, $sample_shifts_tuple,
        $(carrier_doppler + 0.0Hz), $sampling_frequency, 0.0,
        1, $num_samples, $tile_re, $tile_im,
    )
end

# Register only on branches with the tuple kernel (multi-signal branch).
if isdefined(Tracking, :downconvert_and_correlate_fused_tuple!)
    for n_signals in 2:3
        SUITE["fused tuple kernel"]["1-ant N=$n_signals"] =
            bench_fused_tuple_kernel(; num_ants = 1, n_signals)
        SUITE["fused tuple kernel"]["2-ant N=$n_signals"] =
            bench_fused_tuple_kernel(; num_ants = 2, n_signals)
        SUITE["fused tuple kernel"]["4-ant N=$n_signals"] =
            bench_fused_tuple_kernel(; num_ants = 4, n_signals)
    end
end

# ── Per-system multi-satellite track / track! benchmarks ─────────────────
#
# These exercise the full tracking pipeline (downconvert + correlate +
# doppler estimator) on realistic per-system workloads, with
# `bit_buffer.found = true` so the doppler estimator's post-bit-edge
# code path is hit. Four variants are registered per system: the
# {immutable, in-place} × {single-threaded, threaded} cross product.
# The same setup pipeline feeds all four; the only differences between
# entries are the `Tracking.track` vs `Tracking.track!` call and the
# choice of CPU backend.
#
# The `track!` variants are gated on `isdefined(Tracking, :track!)` so
# the script also loads cleanly against master (which has neither
# `track!` nor `CPUThreadedDownconvertAndCorrelator()` zero-arg form).

# Branch-portable per-system container construction. Master accepts a raw
# `Vector{SatState}` via `SystemSatsState(sys, sats)`; the wrapper branch
# needs an estimator to wrap each sat into a `TrackedSat` first via
# `TrackedSystem(estimator, sys, sats)`; the multi-signal branch stores a
# plain `Dictionary{Int, TrackedSat}` directly.
function _build_tracked_system(sys, sats)
    if _HAS_TRACKED_SIGNAL
        return Tracking.to_dictionary(sats)
    elseif _HAS_TRACKED_SAT
        return _TrackedSystem(Tracking.ConventionalAssistedPLLAndDLL(), sys, sats)
    else
        return _TrackedSystem(sys, sats)
    end
end

# Branch-portable "set bit_buffer.found = true" on whatever the dict holds.
# Master holds SatState directly; the wrapper branch holds TrackedSat that
# itself holds SatState; the multi-signal branch holds flat TrackedSat with
# the bit_buffer on `signals[1]`.
#
# BitBuffer has five historical layouts:
#   1. Non-parametric, 7 fields (master + early multi-signal commits).
#   2. Parametric `BitBuffer{B<:Unsigned}`, still 7 fields (Step 1 of the
#      sync-detection redesign).
#   3. Parametric, 9 fields with `secondary_phase::Int` and `polarity::Int8`
#      inserted between `found` and `buffer` (Step 4 of the redesign).
#   4. Parametric, 10 fields with a trailing `soft_bits::Vector{Float32}`.
#   5. Parametric, 11 fields with a further trailing
#      `phase_acc::PhaseAccumulators` (soft-decision L1CA bit-edge
#      detector, issue #124).
# Detect at load time via `hasfield`.
const _HAS_PARAMETRIC_BITBUFFER =
    _HAS_TRACKED_SIGNAL && Tracking.BitBuffer isa UnionAll
const _HAS_BITBUFFER_PHASE_FIELDS =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :secondary_phase)
const _HAS_BITBUFFER_SOFT_BITS =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :soft_bits)
const _HAS_BITBUFFER_PHASE_ACC =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :phase_acc)

# Rebuild a `found = true` BitBuffer matching the existing buffer's
# layout. On parametric branches we pull `B` off the live buffer; on the
# non-parametric branch the type is just `BitBuffer` and the integer
# fields default to `UInt128`.
if _HAS_BITBUFFER_PHASE_ACC
    _bb_int_type(::Tracking.BitBuffer{B}) where {B<:Unsigned} = B
    @inline _make_found_bit_buffer(old_bb) = typeof(old_bb)(
        zero(_bb_int_type(old_bb)),        # code_block_buffer::B
        20,                                # code_block_buffer_lengh
        true,                              # found
        0,                                 # secondary_phase
        Int8(0),                           # polarity
        UInt128(0),                        # buffer
        0,                                 # length
        complex(0.0, 0.0),                 # prompt_accumulator
        0,                                 # prompt_accumulator_integrated_code_blocks
        Float32[],                         # soft_bits
        Tracking.PhaseAccumulators(),      # phase_acc
    )
elseif _HAS_BITBUFFER_SOFT_BITS
    _bb_int_type(::Tracking.BitBuffer{B}) where {B<:Unsigned} = B
    @inline _make_found_bit_buffer(old_bb) = typeof(old_bb)(
        zero(_bb_int_type(old_bb)),        # code_block_buffer::B
        20,                                # code_block_buffer_lengh
        true,                              # found
        0,                                 # secondary_phase
        Int8(0),                           # polarity
        UInt128(0),                        # buffer
        0,                                 # length
        complex(0.0, 0.0),                 # prompt_accumulator
        0,                                 # prompt_accumulator_integrated_code_blocks
        Float32[],                         # soft_bits
    )
elseif _HAS_BITBUFFER_PHASE_FIELDS
    _bb_int_type(::Tracking.BitBuffer{B}) where {B<:Unsigned} = B
    @inline _make_found_bit_buffer(old_bb) = typeof(old_bb)(
        zero(_bb_int_type(old_bb)),        # code_block_buffer::B
        20,                                # code_block_buffer_lengh
        true,                              # found
        0,                                 # secondary_phase
        Int8(0),                           # polarity
        UInt128(0),                        # buffer
        0,                                 # length
        complex(0.0, 0.0),                 # prompt_accumulator
        0,                                 # prompt_accumulator_integrated_code_blocks
    )
elseif _HAS_PARAMETRIC_BITBUFFER
    _bb_int_type(::Tracking.BitBuffer{B}) where {B<:Unsigned} = B
    @inline _make_found_bit_buffer(old_bb) = typeof(old_bb)(
        zero(_bb_int_type(old_bb)),
        20, true, UInt128(0), 0, complex(0.0, 0.0), 0,
    )
elseif _HAS_TRACKED_SIGNAL
    @inline _make_found_bit_buffer(_old_bb) = Tracking.BitBuffer(
        UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0,
    )
end

if _HAS_TRACKED_SIGNAL
    function _with_found_bit_buffer(t, _bb_template)
        sig = only(t.signals)
        new_bb = _make_found_bit_buffer(sig.bit_buffer)
        new_sig = Tracking.TrackedSignal(sig; bit_buffer = new_bb)
        Tracking.TrackedSat(
            t.prn, t.code_phase, t.code_doppler,
            t.carrier_phase, t.carrier_doppler,
            t.signal_start_sample, (new_sig,), t.doppler_estimator_state,
        )
    end
elseif _HAS_TRACKED_SAT
    _with_found_bit_buffer(t, bb) = Tracking.TrackedSat(
        Tracking.SatState(t.sat_state; bit_buffer = bb), t.estimator_state,
    )
else
    _with_found_bit_buffer(s, bb) = Tracking.SatState(s; bit_buffer = bb)
end

# Branch-portable "map a function over the per-sat storage".
if _HAS_TRACKED_SIGNAL
    _map_sats(f, sats_storage) = map(f, sats_storage)
    _rebuild_system_storage(_sss, new_sats) = new_sats
else
    _map_sats(f, sss) = map(f, sss.states)
    _rebuild_system_storage(sss, new_sats) = _TrackedSystem(sss, new_sats)
end

function _make_multi_sat_state(;
    systems,
    nsats_list,
    nsamp,
    prn_max = 32,
    code_dop = 1000.0,
)
    all_sss = []
    total_sats = 0
    for (si, sys) in enumerate(systems)
        ns = nsats_list[min(si, length(nsats_list))]
        pm = sys isa GPSL1CA ? 32 : prn_max
        cd = sys isa GPSL1CA ? 1000.0 : code_dop
        sats = [
            _make_initial_sat(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz)
            for i = 1:ns
        ]
        push!(all_sss, _build_tracked_system(sys, sats))
        total_sats += ns
    end
    _make_track_state(all_sss), rand(ComplexF32, nsamp), total_sats
end

# Branch-portable `TrackState(...)` from a list of per-system / per-group
# storage entries. The multi-band branch's `SignalGroups` type parameter
# requires a `NamedTuple` (not a plain `Tuple`), so build one with
# auto-generated keys. Master / wrapper branches accept either.
if _HAS_TRACKED_SIGNAL
    @inline function _make_track_state(all_sss)
        n = length(all_sss)
        names = ntuple(i -> Symbol("sys", i), n)
        TrackState(NamedTuple{names}(Tuple(all_sss)))
    end
else
    @inline _make_track_state(all_sss) = TrackState(Tuple(all_sss))
end

# Build a `TrackState` for the given system mix with `bit_buffer.found = true`
# on every sat, ready to feed `track` / `track!` for steady-state-style
# benchmarks.
function _make_steady_state_track_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    ts, signal, _ =
        _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    # Pre-multi-signal branches share a single found-true `BitBuffer`
    # across sats; the multi-signal branch builds a per-signal-typed one
    # inside `_with_found_bit_buffer` (parametric `B` per signal). Pass a
    # plain non-parametric template for the older branches and `nothing`
    # for the new one (the new path ignores its second argument).
    found_bb = _HAS_PARAMETRIC_BITBUFFER ? nothing :
        BitBuffer(UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0)
    new_mss = map(_get_systems(ts)) do sss
        new_sats = _map_sats(s -> _with_found_bit_buffer(s, found_bb), sss)
        _rebuild_system_storage(sss, new_sats)
    end
    _track_state_with_systems(ts, new_mss), signal
end

function bench_track_steady_state(
    inplace::Bool, threaded::Bool;
    systems, nsats_list, sfreq, nsamp, prn_max = 32, code_dop = 1000.0,
)
    ts, signal = _make_steady_state_track_state(;
        systems, nsats_list, nsamp, prn_max, code_dop,
    )
    dc = if threaded && isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        _make_cpu_threaded_dc(sfreq)
    else
        _make_cpu_dc(sfreq)
    end
    if inplace
        @benchmarkable Tracking.track!(
            $signal, $ts, $sfreq; downconvert_and_correlator = $dc,
        )
    else
        @benchmarkable Tracking.track(
            $signal, $ts, $sfreq; downconvert_and_correlator = $dc,
        )
    end
end

# Three per-system workloads exercised across four entry-point variants
# (immutable / in-place × single-threaded / threaded). Naming convention:
#
#   track[!][-threaded] / <system> <Nsats>sat/<Nsamp>
#
# Master only registers the immutable variants (no `track!`); the
# threaded suffix is also master-compatible since the underlying
# `CPUThreadedDownconvertAndCorrelator` exists in both.
# Flat keys with numeric prefixes so AirspeedVelocity's alphabetical sort
# clusters all four {track, track!, track-threaded, track!-threaded} variants
# of each case together in the PR comparison table. All entries live under a
# single top-level "track" group; the variant is part of the case name so the
# case prefix dominates the sort order.
const _TRACK_BENCH_CASES = let gpsl1 = GPSL1CA(), gal = GalileoE1B()
    [
        ("2. L1 8sat/5K",   (systems = (gpsl1,),     nsats_list = [8],    sfreq = 5e6Hz,  nsamp = 5000)),
        ("3. E1B 4sat/25K", (systems = (gal,),       nsats_list = [4],    sfreq = 25e6Hz, nsamp = 25000, prn_max = 50, code_dop = 100.0)),
        ("4. 8L1+8E1B/25K", (systems = (gpsl1, gal), nsats_list = [8, 8], sfreq = 25e6Hz, nsamp = 25000, prn_max = 50, code_dop = 100.0)),
    ]
end

for (key, kw) in _TRACK_BENCH_CASES
    SUITE["track"]["$key – track"] = bench_track_steady_state(false, false; kw...)
    SUITE["track"]["$key – track-threaded"] = bench_track_steady_state(false, true; kw...)
    if isdefined(Tracking, :track!)
        SUITE["track"]["$key – track!"] = bench_track_steady_state(true, false; kw...)
        SUITE["track"]["$key – track!-threaded"] = bench_track_steady_state(true, true; kw...)
    end
end

# ── Multi-signal-per-sat track benchmark ─────────────────────────────────
# Measures track cost as N signals stack on a single satellite (the new
# multi-signal hot path). Only registered when TrackedSignal is available.
#
# Step 3 ships the fused-N-times baseline: one independent
# `downconvert_and_correlate_fused!` call per signal in the sat's
# `signals` tuple. Linear scaling with N. Step 4 will add a tile-share
# kernel that does one downconvert + N correlate-from-tile passes,
# expected to save ~20% at N=2 and ~38% at N=3 — see the design doc and
# the microbenchmark probes in `claude_scratch/`.
#
# Stacking N copies of GPSL1CA on one sat is not a real-world scenario
# (a real sat carries one signal of each kind), but it isolates the
# tuple-walk cost on identical N for a clean per-signal-cost comparison.
if _HAS_TRACKED_SIGNAL
    function _make_multi_signal_track_state(; n_signals, nsamp, sfreq)
        gpsl1 = GPSL1CA()
        estimator = Tracking.ConventionalAssistedPLLAndDLL()
        signals = ntuple(
            _ -> Tracking.TrackedSignal(
                gpsl1;
                num_ants = NumAnts(1),
                correlator = Tracking.EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                post_corr_filter = Tracking.DefaultPostCorrFilter(),
            ),
            n_signals,
        )
        carrier_doppler = 1000.0Hz
        # Prefer the public multi-signal `TrackedSat((sigs...), prn, ...)`
        # constructor (added in #133); fall back to the hand-rolled
        # bare-sat → init_estimator_state → rebuild build on revisions that
        # predate it, so this script still benchmarks the baseline.
        sat = if hasmethod(
            Tracking.TrackedSat,
            Tuple{typeof(signals),Int,Float64,typeof(carrier_doppler)},
        )
            Tracking.TrackedSat(
                signals, 1, 10.5, carrier_doppler;
                doppler_estimator = estimator,
            )
        else
            code_doppler =
                carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(gpsl1)
            bare = Tracking.TrackedSat(
                1, 10.5, code_doppler, 0.0, carrier_doppler, 1, signals, nothing,
            )
            de_state = Tracking.init_estimator_state(estimator, bare)
            Tracking.TrackedSat(
                bare.prn, bare.code_phase, bare.code_doppler,
                bare.carrier_phase, bare.carrier_doppler,
                bare.signal_start_sample, bare.signals, de_state,
            )
        end
        ts = TrackState(gpsl1, sat; doppler_estimator = estimator)
        signal = rand(Complex{Float32}, nsamp)
        ts, signal
    end

    for n_signals = 1:3
        prefix = "$(4 + n_signals). multi-signal N=$n_signals/5K"
        ts, signal =
            _make_multi_signal_track_state(; n_signals, nsamp = 5000, sfreq = 5e6Hz)
        dc = _make_cpu_dc(5e6Hz)
        SUITE["track"]["$prefix – track"] = @benchmarkable Tracking.track(
            $signal, $ts, $(5e6Hz); downconvert_and_correlator = $dc,
        )
        if isdefined(Tracking, :track!)
            ts_ip, signal_ip = _make_multi_signal_track_state(;
                n_signals, nsamp = 5000, sfreq = 5e6Hz,
            )
            dc_ip = _make_cpu_dc(5e6Hz)
            SUITE["track"]["$prefix – track!"] = @benchmarkable Tracking.track!(
                $signal_ip, $ts_ip, $(5e6Hz); downconvert_and_correlator = $dc_ip,
            )
        end
    end
end

# ── Int16 vs Float32 backend, full track! ─────────────────────────────────────
# Head-to-head of the Float32 default (CPUThreadedDownconvertAndCorrelator) and
# the integer Complex{Int16} backend (Int16ThreadedDownconvertAndCorrelator)
# through the full `track!` pipeline, on the SAME Complex{Int16} (12-bit ADC)
# capture so the comparison is apples-to-apples. Both accept Complex{Int16}: the
# Float32 path widens to Float32; the integer path requires it (and is the point
# of the backend). Threaded — the real-time default. Each case registers
# "Float32" and "Int16" under `track! Int16 vs Float32/<case>` so the benchmark
# comment can pair them into one Int16-vs-Float32 speedup row. Only registered on
# branches with the Int16 backend (skipped against master so AirspeedVelocity can
# still diff).
if isdefined(Tracking, :Int16ThreadedDownconvertAndCorrelator)
    # Random 12-bit-ADC samples. The kernel's run time is content-independent, but
    # the magnitude must stay within the ±2^11 range the Int16 carrier wipe assumes.
    function _int16_capture(nsamp)
        lim = Int16(2048)
        complex.(rand((-lim):(lim - one(Int16)), nsamp), rand((-lim):(lim - one(Int16)), nsamp))
    end

    # Self-explanatory scenario names (system · sat count · sampling rate). NO
    # "/" in a name — the benchmark-table script keys leaves by "/"-joined path,
    # so a slash would be read as an extra nesting level. Each scenario registers
    # a "Float32" and an "Int16" leaf under the same parent, so they sort
    # adjacently and the table pairs them into one Int16-vs-Float32 speedup row.
    for (name, systems, nsats_list, sfreq, nsamp, prn_max) in (
        ("GPS L1CA, 8 sats @ 5 MHz", (GPSL1CA(),), [8], 5e6Hz, 5000, 32),
        ("GPS L1CA, 8 sats @ 40 MHz", (GPSL1CA(),), [8], 40e6Hz, 40000, 32),
        ("Galileo E1B, 4 sats @ 25 MHz", (GalileoE1B(),), [4], 25e6Hz, 25000, 50),
    )
        sig16 = _int16_capture(nsamp)
        # Two independent track states (track! mutates) so each backend starts
        # from an identical steady-state (bit_buffer.found = true) configuration.
        ts_f, _ = _make_steady_state_track_state(;
            systems, nsats_list, nsamp, prn_max, code_dop = 100.0)
        ts_i, _ = _make_steady_state_track_state(;
            systems, nsats_list, nsamp, prn_max, code_dop = 100.0)
        dc_f = _make_cpu_threaded_dc(sfreq)
        dc_i = Tracking.Int16ThreadedDownconvertAndCorrelator()
        g = SUITE["track! Int16 vs Float32"][name]
        g["Float32"] = @benchmarkable Tracking.track!(
            $sig16, $ts_f, $sfreq; downconvert_and_correlator = $dc_f,
        )
        g["Int16"] = @benchmarkable Tracking.track!(
            $sig16, $ts_i, $sfreq; downconvert_and_correlator = $dc_i,
        )
    end
end
