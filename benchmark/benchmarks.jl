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

# AirspeedVelocity benches every rev with HEAD's copy of this file (`--bench-on`),
# so it must load against every Tracking/GNSSSignals version in use; the shims
# below detect the loaded API. GNSSSignals v1 → v2 signal names:
const GPSL1CA =
    isdefined(GNSSSignals, :GPSL1CA) ? getfield(GNSSSignals, :GPSL1CA) :
    getfield(GNSSSignals, :GPSL1)
const GPSL5I =
    isdefined(GNSSSignals, :GPSL5I) ? getfield(GNSSSignals, :GPSL5I) :
    getfield(GNSSSignals, :GPSL5)

# Code-replica buffer element type. GNSSSignals' embedded-LUT `gen_code!` (PR #90)
# is Int8-only; older GNSSSignals emit `get_code_type(signal)` (Int16/Float32).
# Detect the era so this script builds against both (AirspeedVelocity diffs revs).
_code_buf_type(sig) = isdefined(GNSSSignals, :code_engine) ? Int8 : get_code_type(sig)

# Per-system storage: `SystemSatsState` (master), `TrackedSystem` (wrapper branch),
# or a plain `Dictionary{Int, TrackedSat}` (multi-signal).
const _HAS_TRACKED_SIGNAL = isdefined(Tracking, :TrackedSignal)
const _HAS_TRACKED_SAT = isdefined(Tracking, :TrackedSat)
const _HAS_SAT_STATE = isdefined(Tracking, :SatState)

if !_HAS_TRACKED_SIGNAL
    const _TrackedSystem =
        isdefined(Tracking, :TrackedSystem) ? Tracking.TrackedSystem :
        Tracking.SystemSatsState
end

# `TrackState` field for the per-system tuple: `multiple_system_sats_state`
# (master), `satellites` (wrapper branch), or `groups` of `SignalGroup`s
# (multi-band), which is mapped to and from a NamedTuple of dicts.
const _MULTIBAND = _HAS_TRACKED_SIGNAL
const _SYSTEMS_FIELD =
    _MULTIBAND ? :groups :
    isdefined(Tracking, :TrackedSystem) ? :satellites : :multiple_system_sats_state

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

# Per-sat construction: a `SatState` (master / wrapper branch) or a seeded
# `TrackedSat` (multi-signal). `cn0_estimator = nothing` keeps the package
# default and is splatted away, so revisions without that keyword still work.
function _make_initial_sat(
    sys,
    prn,
    code_phase,
    carrier_doppler;
    num_ants = NumAnts(1),
    cn0_estimator = nothing,
)
    cn0_kw = isnothing(cn0_estimator) ? (;) : (; cn0_estimator)
    if _HAS_TRACKED_SIGNAL
        return Tracking.TrackedSat(
            sys,
            prn,
            code_phase,
            carrier_doppler;
            doppler_estimator = Tracking.ConventionalAssistedPLLAndDLL(),
            num_ants,
            cn0_kw...,
        )
    else
        return Tracking.SatState(sys, prn, code_phase, carrier_doppler; num_ants, cn0_kw...)
    end
end

_make_initial_sat_with_num_ants(
    sys,
    prn,
    code_phase,
    carrier_doppler,
    num_ants;
    cn0_estimator = nothing,
) = _make_initial_sat(sys, prn, code_phase, carrier_doppler; num_ants, cn0_estimator)

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
    code_doppler =
        carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(gnss_signal)
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
# `CPUThreadedDownconvertAndCorrelator`, and `gen_code_replica!`.

const _NEEDS_VAL = !applicable(Tracking.CPUDownconvertAndCorrelator)

_make_cpu_dc(sampling_frequency) =
    _NEEDS_VAL ? Tracking.CPUDownconvertAndCorrelator(Val(sampling_frequency)) :
    Tracking.CPUDownconvertAndCorrelator()

_make_cpu_threaded_dc(sampling_frequency) =
    _NEEDS_VAL ? Tracking.CPUThreadedDownconvertAndCorrelator(Val(sampling_frequency)) :
    Tracking.CPUThreadedDownconvertAndCorrelator()

function _gen_code_replica!(
    code_replica,
    gnss_signal,
    code_frequency,
    sampling_frequency,
    code_phase,
    start_sample,
    num_samples,
    sample_shifts,
    prn,
)
    if _NEEDS_VAL
        gen_code_replica!(
            code_replica,
            gnss_signal,
            code_frequency,
            sampling_frequency,
            code_phase,
            start_sample,
            num_samples,
            sample_shifts,
            prn,
            Val(sampling_frequency),
        )
    else
        gen_code_replica!(
            code_replica,
            gnss_signal,
            code_frequency,
            sampling_frequency,
            code_phase,
            start_sample,
            num_samples,
            sample_shifts,
            prn,
        )
    end
end

# ── High-level downconvert_and_correlate (one PRN's pipeline) ──────────────
#
# The per-satellite pipeline (quantise, downconvert, despread, accumulate) for one
# PRN, between the `fused kernel` and `track` rows. The sat uses
# `MomentsCN0Estimator` (same signature on every rev) so no noise reference is
# provisioned: the default `NoiseRefCN0Estimator` would despread a second
# correlator channel here, roughly doubling these rows. That cost is measured by
# the `noise estimation` group and the `track` rows.
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
        [
            _make_initial_sat_with_num_ants(
                gnss_signal,
                1,
                10.5,
                1000.0Hz,
                NumAnts(num_ants);
                cn0_estimator = Tracking.MomentsCN0Estimator(100),
            ),
        ],
    )
    signal =
        num_ants == 1 ? rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    # Measurement type (`BandMeasurement`, `Measurement` before #132) and call
    # arity (3-arg, or 4-arg with a trailing integration length) differ by rev;
    # master / wrapper branches use the legacy 6-arg form.
    @static if isdefined(Tracking, :BandMeasurement) || isdefined(Tracking, :Measurement)
        measurement_type =
            isdefined(Tracking, :BandMeasurement) ? Tracking.BandMeasurement :
            Tracking.Measurement
        # Band key: `Tracking.band_key` (`:l1`, Tracking ≤ 2.3) or
        # `GNSSSignals.get_band_id` (`:L1`).
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
            Tuple{
                typeof(downconvert_and_correlator),
                typeof(measurements),
                typeof(track_state),
            },
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
    # Dynamic shifts: newer revs take caller-supplied SoA tile buffers (11 args);
    # master uses the 9-arg form.
    static_arg_types = Tuple{
        typeof(s.correlator),
        typeof(s.signal),
        typeof(s.code_replica),
        typeof(sample_shifts),
        typeof(s.carrier_doppler),
        typeof(s.sampling_frequency),
        Float64,
        Int,
        Int,
    }
    if shifts == :dynamic &&
       !hasmethod(Tracking.downconvert_and_correlate_fused!, static_arg_types)
        # New branch: 11-arg signature with hoisted tile buffers
        tile_re = Vector{Float32}(undef, num_samples * num_ants)
        tile_im = Vector{Float32}(undef, num_samples * num_ants)
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
            tile_re,
            tile_im,
        )
        return @benchmarkable Tracking.downconvert_and_correlate_fused!(
            $(s.correlator),
            $(s.signal),
            $(s.code_replica),
            $sample_shifts,
            $(s.carrier_doppler),
            $(s.sampling_frequency),
            0.0,
            1,
            $(s.num_samples),
            $tile_re,
            $tile_im,
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
SUITE["track"]["GPS L1CA, 1 sat, 2K @ 5 MHz – out-of-place, Float32"] = bench_track()

# In-place track! (only on branches that define it).
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
    SUITE["track"]["GPS L1CA, 1 sat, 2K @ 5 MHz – in-place, Float32"] =
        bench_track_inplace()
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
# `downconvert_and_correlate_fused_tuple!` directly, single- and multi-antenna.
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
    code_doppler =
        carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(gnss_signal)
    code_frequency = code_doppler + get_code_frequency(gnss_signal)

    correlator_template = EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants))
    sample_shifts = get_correlator_sample_shifts(
        correlator_template,
        sampling_frequency,
        code_frequency,
    )
    code_replica_size = num_samples + maximum(sample_shifts) - minimum(sample_shifts)

    signal =
        num_ants == 1 ? rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    correlators =
        ntuple(_ -> EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants)), n_signals)
    code_replicas = ntuple(n_signals) do _
        cr = Vector{_code_buf_type(gnss_signal)}(undef, code_replica_size)
        _gen_code_replica!(
            cr,
            gnss_signal,
            code_frequency,
            sampling_frequency,
            code_phase,
            1,
            num_samples,
            sample_shifts,
            1,
        )
        cr
    end
    sample_shifts_tuple = ntuple(_ -> sample_shifts, n_signals)
    tile_re = Vector{Float32}(undef, num_samples * num_ants)
    tile_im = Vector{Float32}(undef, num_samples * num_ants)

    @benchmarkable Tracking.downconvert_and_correlate_fused_tuple!(
        $correlators,
        $signal,
        $code_replicas,
        $sample_shifts_tuple,
        $(carrier_doppler + 0.0Hz),
        $sampling_frequency,
        0.0,
        1,
        $num_samples,
        $tile_re,
        $tile_im,
    )
end

# Register only on branches with the tuple kernel (multi-signal branch).
if isdefined(Tracking, :downconvert_and_correlate_fused_tuple!)
    for n_signals = 2:3
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
# Full pipeline on realistic per-system workloads, with `bit_buffer.found = true`
# so the post-bit-edge path is hit. Variants: {out-of-place, in-place} ×
# {single-threaded, threaded}; `track!` rows are gated on `isdefined`.

# Per-system container; see the storage-flavour note at the top.
function _build_tracked_system(sys, sats)
    if _HAS_TRACKED_SIGNAL
        return Tracking.to_dictionary(sats)
    elseif _HAS_TRACKED_SAT
        return _TrackedSystem(Tracking.ConventionalAssistedPLLAndDLL(), sys, sats)
    else
        return _TrackedSystem(sys, sats)
    end
end

# Branch-portable "set bit_buffer.found = true" (on `signals[1]` for the
# multi-signal branch). BitBuffer layouts across revs, detected via `hasfield`:
#   1. Non-parametric, 7 fields.
#   2. Parametric `BitBuffer{B<:Unsigned}`, 7 fields.
#   3. + `secondary_phase`, `polarity` after `found` (9 fields).
#   4. + trailing `soft_bits` (10).
#   5. + trailing `phase_acc` (11, issue #124).
#   6. Hard-bit `buffer` / `length` pair removed (9).
const _HAS_PARAMETRIC_BITBUFFER = _HAS_TRACKED_SIGNAL && Tracking.BitBuffer isa UnionAll
const _HAS_BITBUFFER_PHASE_FIELDS =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :secondary_phase)
const _HAS_BITBUFFER_SOFT_BITS =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :soft_bits)
const _HAS_BITBUFFER_PHASE_ACC =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :phase_acc)
const _HAS_BITBUFFER_HARD_BITS =
    _HAS_PARAMETRIC_BITBUFFER && hasfield(Tracking.BitBuffer, :buffer)

# Rebuild a `found = true` BitBuffer in the live buffer's layout (and `B`).
if _HAS_BITBUFFER_PHASE_ACC && !_HAS_BITBUFFER_HARD_BITS
    _bb_int_type(::Tracking.BitBuffer{B}) where {B<:Unsigned} = B
    @inline _make_found_bit_buffer(old_bb) = typeof(old_bb)(
        zero(_bb_int_type(old_bb)),        # code_block_buffer::B
        20,                                # code_block_buffer_length
        true,                              # found
        0,                                 # secondary_phase
        Int8(0),                           # polarity
        complex(0.0, 0.0),                 # prompt_accumulator
        0,                                 # prompt_accumulator_integrated_code_blocks
        Float32[],                         # soft_bits
        Tracking.PhaseAccumulators(),      # phase_acc
    )
elseif _HAS_BITBUFFER_PHASE_ACC
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
        20,
        true,
        UInt128(0),
        0,
        complex(0.0, 0.0),
        0,
    )
elseif _HAS_TRACKED_SIGNAL
    @inline _make_found_bit_buffer(_old_bb) =
        Tracking.BitBuffer(UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0)
end

if _HAS_TRACKED_SIGNAL
    function _with_found_bit_buffer(t, _bb_template)
        sig = only(t.signals)
        new_bb = _make_found_bit_buffer(sig.bit_buffer)
        new_sig = Tracking.TrackedSignal(sig; bit_buffer = new_bb)
        Tracking.TrackedSat(
            t.prn,
            t.code_phase,
            t.code_doppler,
            t.carrier_phase,
            t.carrier_doppler,
            t.signal_start_sample,
            (new_sig,),
            t.doppler_estimator_state,
        )
    end
elseif _HAS_TRACKED_SAT
    _with_found_bit_buffer(t, bb) = Tracking.TrackedSat(
        Tracking.SatState(t.sat_state; bit_buffer = bb),
        t.estimator_state,
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
            _make_initial_sat(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for
            i = 1:ns
        ]
        push!(all_sss, _build_tracked_system(sys, sats))
        total_sats += ns
    end
    _make_track_state(all_sss), rand(ComplexF32, nsamp), total_sats
end

# `TrackState` from per-system storage; the multi-band branch needs a NamedTuple.
if _HAS_TRACKED_SIGNAL
    @inline function _make_track_state(all_sss)
        n = length(all_sss)
        names = ntuple(i -> Symbol("sys", i), n)
        TrackState(NamedTuple{names}(Tuple(all_sss)))
    end
else
    @inline _make_track_state(all_sss) = TrackState(Tuple(all_sss))
end

# `TrackState` for the given system mix with `bit_buffer.found = true` on every sat.
function _make_steady_state_track_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    ts, signal, _ = _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    # Older branches share one found-true template; parametric ones build their
    # own in `_with_found_bit_buffer` and ignore it.
    found_bb =
        _HAS_PARAMETRIC_BITBUFFER ? nothing :
        BitBuffer(UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0)
    new_mss = map(_get_systems(ts)) do sss
        new_sats = _map_sats(s -> _with_found_bit_buffer(s, found_bb), sss)
        _rebuild_system_storage(sss, new_sats)
    end
    _track_state_with_systems(ts, new_mss), signal
end

function bench_track_steady_state(
    inplace::Bool,
    threaded::Bool;
    systems,
    nsats_list,
    sfreq,
    nsamp,
    prn_max = 32,
    code_dop = 1000.0,
)
    ts, signal =
        _make_steady_state_track_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    dc = if threaded && isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        _make_cpu_threaded_dc(sfreq)
    else
        _make_cpu_dc(sfreq)
    end
    if inplace
        @benchmarkable Tracking.track!(
            $signal,
            $ts,
            $sfreq;
            downconvert_and_correlator = $dc,
        )
    else
        @benchmarkable Tracking.track(
            $signal,
            $ts,
            $sfreq;
            downconvert_and_correlator = $dc,
        )
    end
end

# Naming convention for every leaf under the "track" group:
#
#   <constellation>, <Nsats> sats[, <extra>], <Nsamp> @ <rate> – <variant>
#
# where `<variant>` is `out-of-place` (`track`) or `in-place` (`track!`), then the
# backend (`Float32` / `Int16` / `OneBit` / `TwoBit`, always named), then
# `, threaded` on the rows that are. The scenario carries sat count, sample count
# and rate because the comparison table is read without this file. No numeric
# prefixes: the table sorts leaves alphabetically on the "/"-joined path, and this
# text keeps each scenario's variants (and a backend's Float32 sibling) adjacent.
# ("500K" sorts before "5K"; left as is.)
const _TRACK_BENCH_CASES = let gpsl1 = GPSL1CA(), gal = GalileoE1B()
    [
        (
            "GPS L1CA, 8 sats, 5K @ 5 MHz",
            (systems = (gpsl1,), nsats_list = [8], sfreq = 5e6Hz, nsamp = 5000),
        ),
        (
            "Galileo E1B, 4 sats, 25K @ 25 MHz",
            (
                systems = (gal,),
                nsats_list = [4],
                sfreq = 25e6Hz,
                nsamp = 25000,
                prn_max = 50,
                code_dop = 100.0,
            ),
        ),
        (
            "GPS L1CA + Galileo E1B, 8+8 sats, 25K @ 25 MHz",
            (
                systems = (gpsl1, gal),
                nsats_list = [8, 8],
                sfreq = 25e6Hz,
                nsamp = 25000,
                prn_max = 50,
                code_dop = 100.0,
            ),
        ),
    ]
end

for (key, kw) in _TRACK_BENCH_CASES
    SUITE["track"]["$key – out-of-place, Float32"] =
        bench_track_steady_state(false, false; kw...)
    SUITE["track"]["$key – out-of-place, Float32, threaded"] =
        bench_track_steady_state(false, true; kw...)
    if isdefined(Tracking, :track!)
        SUITE["track"]["$key – in-place, Float32"] =
            bench_track_steady_state(true, false; kw...)
        SUITE["track"]["$key – in-place, Float32, threaded"] =
            bench_track_steady_state(true, true; kw...)
    end
end

# ── Multi-signal-per-sat track benchmark ─────────────────────────────────
# Track cost as N signals stack on one satellite. N identical GPSL1CA copies are
# not realistic but isolate the per-signal cost of the tuple walk.
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
        # Public tuple constructor (#133), or a hand-rolled build on older revs.
        sat =
            if hasmethod(
                Tracking.TrackedSat,
                Tuple{typeof(signals),Int,Float64,typeof(carrier_doppler)},
            )
                Tracking.TrackedSat(
                    signals,
                    1,
                    10.5,
                    carrier_doppler;
                    doppler_estimator = estimator,
                )
            else
                code_doppler =
                    carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(gpsl1)
                bare = Tracking.TrackedSat(
                    1,
                    10.5,
                    code_doppler,
                    0.0,
                    carrier_doppler,
                    1,
                    signals,
                    nothing,
                )
                de_state = Tracking.init_estimator_state(estimator, bare)
                Tracking.TrackedSat(
                    bare.prn,
                    bare.code_phase,
                    bare.code_doppler,
                    bare.carrier_phase,
                    bare.carrier_doppler,
                    bare.signal_start_sample,
                    bare.signals,
                    de_state,
                )
            end
        ts = TrackState(gpsl1, sat; doppler_estimator = estimator)
        signal = rand(Complex{Float32}, nsamp)
        ts, signal
    end

    for n_signals = 1:3
        # Sweeps signals per sat, so the label carries both counts.
        prefix = "GPS L1CA, 1 sat, $n_signals signal$(n_signals == 1 ? "" : "s"), 5K @ 5 MHz"
        ts, signal =
            _make_multi_signal_track_state(; n_signals, nsamp = 5000, sfreq = 5e6Hz)
        dc = _make_cpu_dc(5e6Hz)
        SUITE["track"]["$prefix – out-of-place, Float32"] = @benchmarkable Tracking.track(
            $signal,
            $ts,
            $(5e6Hz);
            downconvert_and_correlator = $dc,
        )
        if isdefined(Tracking, :track!)
            ts_ip, signal_ip =
                _make_multi_signal_track_state(; n_signals, nsamp = 5000, sfreq = 5e6Hz)
            dc_ip = _make_cpu_dc(5e6Hz)
            SUITE["track"]["$prefix – in-place, Float32"] = @benchmarkable Tracking.track!(
                $signal_ip,
                $ts_ip,
                $(5e6Hz);
                downconvert_and_correlator = $dc_ip,
            )
        end
    end
end

# ── Multi-code-period allocation guard (track!) ───────────────────────────────
# A warm `track!` must be allocation-free per iteration, so `memory` must stay flat
# as the block count grows; a per-iteration leak shows as a jump between lengths.
# Pre-sync runs to 200 blocks because the CFAR detector only scores (and evaluates
# its Student-t threshold) from `2 * blocks_per_bit` = 40 blocks; the hard
# assertion is the pre-sync guard in test/track_in_place.jl. `setup` builds a fresh
# state and warms it with one call, so the measured call is the second; with
# `evals = 1` every measured call gets a fresh state, which keeps the loop filter
# from drifting to NaN (see the Int16 rows).
if _HAS_TRACKED_SIGNAL && isdefined(Tracking, :track!)
    let sfreq = 5e6Hz, gpsl1 = GPSL1CA()
        samples_per_block = 5000               # 1 ms GPS L1CA code period @ 5 MHz
        for n_blocks in (2, 20, 200)
            nsamp = samples_per_block * n_blocks
            sig = rand(ComplexF32, nsamp)
            dc = _make_cpu_dc(sfreq)
            SUITE["track"]["GPS L1CA, 1 sat, bit sync pending, $(n_blocks) blk @ 5 MHz – in-place, Float32"] = @benchmarkable(
                Tracking.track!($sig, ts, $sfreq; downconvert_and_correlator = $dc),
                setup = (
                    ts = first(
                        _make_multi_sat_state(;
                            systems = ($gpsl1,),
                            nsats_list = [1],
                            nsamp = $nsamp,
                        ),
                    );
                    Tracking.track!($sig, ts, $sfreq; downconvert_and_correlator = $dc)
                ),
                evals = 1,
            )
            # Post-sync there is no per-block search to leak from.
            n_blocks == 200 && continue
            SUITE["track"]["GPS L1CA, 1 sat, bit sync found, $(n_blocks) blk @ 5 MHz – in-place, Float32"] = @benchmarkable(
                Tracking.track!($sig, ts, $sfreq; downconvert_and_correlator = $dc),
                setup = (
                    ts = first(
                        _make_steady_state_track_state(;
                            systems = ($gpsl1,),
                            nsats_list = [1],
                            nsamp = $nsamp,
                            prn_max = 32,
                            code_dop = 1000.0,
                        ),
                    );
                    Tracking.track!($sig, ts, $sfreq; downconvert_and_correlator = $dc)
                ),
                evals = 1,
            )
        end
    end
end

# ── Int16 vs Float32 backend, full track! ─────────────────────────────────────
# Threaded backends through the full `track!` on the same Complex{Int16} (12-bit
# ADC) capture. Each case registers sibling leaves under
# `track! Int16 vs Float32/<case>`, which the benchmark comment pairs into
# speedup rows.
if isdefined(Tracking, :Int16ThreadedDownconvertAndCorrelator)
    # Random 12-bit-ADC samples. The kernel's run time is content-independent, but
    # the magnitude must stay within the ±2^11 range the Int16 carrier wipe assumes.
    function _int16_capture(nsamp)
        lim = Int16(2048)
        complex.(rand((-lim):(lim-one(Int16)), nsamp), rand((-lim):(lim-one(Int16)), nsamp))
    end

    # `max_meas` is positional now; older revs took it as a keyword.
    _make_int16_threaded_dc(max_meas) =
        applicable(Tracking.Int16ThreadedDownconvertAndCorrelator, max_meas) ?
        Tracking.Int16ThreadedDownconvertAndCorrelator(max_meas) :
        Tracking.Int16ThreadedDownconvertAndCorrelator(; max_meas)
    _make_int16_dc(max_meas) =
        applicable(Tracking.Int16DownconvertAndCorrelator, max_meas) ?
        Tracking.Int16DownconvertAndCorrelator(max_meas) :
        Tracking.Int16DownconvertAndCorrelator(; max_meas)

    # Scenario names must not contain "/" (the table script nests on it).
    # 8-PRN GPS L1CA capture matching `_make_multi_sat_state`'s per-sat pattern, for
    # the 100 ms case: its ~100 loop updates in one call would drift to a NaN
    # Doppler over pure noise. Amplitudes stay well inside ±2048.
    function _l1ca_composite_capture(nsamp, sfreq)
        sys = GPSL1CA()
        fsn = sfreq / Hz
        acc = zeros(ComplexF64, nsamp)
        for i = 1:8
            d = 1000.0 + i * 10
            code = gen_code(
                nsamp,
                sys,
                i,
                sfreq,
                get_code_frequency(sys) + (d / 1540)Hz,
                10.5 + i * 0.1,
            )
            acc .+= 64.0 .* code .* cis.(2π .* d .* (0:(nsamp-1)) ./ fsn)
        end
        acc .+= 32.0 .* randn(nsamp) .+ 32.0im .* randn(nsamp)
        complex.(round.(Int16, real.(acc)), round.(Int16, imag.(acc)))
    end

    # The 100 ms case catches a per-step (instead of per-call) band repack in the
    # bit backends, invisible at 1 ms.
    for (name, systems, nsats_list, sfreq, nsamp, prn_max, capture) in (
        ("GPS L1CA, 8 sats @ 5 MHz", (GPSL1CA(),), [8], 5e6Hz, 5000, 32, _int16_capture),
        (
            "GPS L1CA, 8 sats @ 5 MHz, 100 ms buffer",
            (GPSL1CA(),),
            [8],
            5e6Hz,
            500_000,
            32,
            n -> _l1ca_composite_capture(n, 5e6Hz),
        ),
        ("GPS L1CA, 8 sats @ 40 MHz", (GPSL1CA(),), [8], 40e6Hz, 40000, 32, _int16_capture),
        (
            "Galileo E1B, 4 sats @ 25 MHz",
            (GalileoE1B(),),
            [4],
            25e6Hz,
            25000,
            50,
            _int16_capture,
        ),
    )
        sig16 = capture(nsamp)
        dc_f = _make_cpu_threaded_dc(sfreq)
        # max_meas = 2^11 matches the ±2048 full-scale of `_int16_capture` above.
        dc_i = _make_int16_threaded_dc(2^11)
        g = SUITE["track! Int16 vs Float32"][name]
        # Tracking random samples drifts the loop filter to a NaN Doppler
        # (`InexactError`) over many evals, so `setup` rebuilds a fresh synced state
        # per sample and `evals = 1` runs each measured call once.
        g["Float32"] = @benchmarkable(
            Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_f),
            setup = (
                ts = first(
                    _make_steady_state_track_state(;
                        systems = $systems,
                        nsats_list = $nsats_list,
                        nsamp = $nsamp,
                        prn_max = $prn_max,
                        code_dop = 100.0,
                    ),
                )
            ),
            evals = 1,
        )
        g["Int16"] = @benchmarkable(
            Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_i),
            setup = (
                ts = first(
                    _make_steady_state_track_state(;
                        systems = $systems,
                        nsats_list = $nsats_list,
                        nsamp = $nsamp,
                        prn_max = $prn_max,
                        code_dop = 100.0,
                    ),
                )
            ),
            evals = 1,
        )
        # One-bit (bit-wise) backend, same capture. BPSK-only (it errors on CBOC), so
        # register it only for the non-Galileo cases. Guarded for base revs without it.
        if isdefined(Tracking, :OneBitThreadedDownconvertAndCorrelator) &&
           !(first(systems) isa GalileoE1B)
            dc_b = Tracking.OneBitThreadedDownconvertAndCorrelator()
            g["OneBit"] = @benchmarkable(
                Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_b),
                setup = (
                    ts = first(
                        _make_steady_state_track_state(;
                            systems = $systems,
                            nsats_list = $nsats_list,
                            nsamp = $nsamp,
                            prn_max = $prn_max,
                            code_dop = 100.0,
                        ),
                    )
                ),
                evals = 1,
            )
        end
        # Two-bit (sign+magnitude bit-wise) backend, same capture and BPSK-only gate.
        if isdefined(Tracking, :TwoBitThreadedDownconvertAndCorrelator) &&
           !(first(systems) isa GalileoE1B)
            dc_t = Tracking.TwoBitThreadedDownconvertAndCorrelator()
            g["TwoBit"] = @benchmarkable(
                Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_t),
                setup = (
                    ts = first(
                        _make_steady_state_track_state(;
                            systems = $systems,
                            nsats_list = $nsats_list,
                            nsamp = $nsamp,
                            prn_max = $prn_max,
                            code_dop = 100.0,
                        ),
                    )
                ),
                evals = 1,
            )
        end
    end
end

# ── Backend axes: multi-signal / multi-antenna / dynamic taps ─────────────────
# Head-to-head rows over multi-signal, multi-antenna and dynamic tap counts, under
# the existing `INT16_GROUP` because the workflow runs the base ref's
# `bench_table.jl`, which would not render a new group. Multi-signal / multi-antenna
# time `downconvert_and_correlate`; dynamic taps are unreachable through the pipeline
# (correlators pass `SVector` shifts), so that row times the kernel fallbacks.
if isdefined(Tracking, :OneBitThreadedDownconvertAndCorrelator) &&
   isdefined(Tracking, :Int16ThreadedDownconvertAndCorrelator)
    const _AXES_SIG = GPSL1CA()
    const _AXES_FS = 5e6Hz
    const _AXES_NSAMP = 5000
    _int16_capture_mat(nsamp, M) = repeat(_int16_capture(nsamp); outer = (1, M))

    # A multi-signal sat: N GPS L1CA signals sharing one carrier + measurement.
    function _axes_multisignal_state(n_signals)
        est = Tracking.ConventionalAssistedPLLAndDLL()
        signals = ntuple(
            _ -> Tracking.TrackedSignal(
                _AXES_SIG;
                num_ants = NumAnts(1),
                correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                post_corr_filter = Tracking.DefaultPostCorrFilter(),
            ),
            n_signals,
        )
        sat = Tracking.TrackedSat(signals, 1, 10.5, 1000.0Hz; doppler_estimator = est)
        TrackState(_AXES_SIG, sat; doppler_estimator = est)
    end

    # A single-signal sat with M antenna channels.
    _axes_multiantenna_state(num_ants) = TrackState(
        _AXES_SIG,
        [_make_initial_sat(_AXES_SIG, 1, 10.5, 1000.0Hz; num_ants = NumAnts(num_ants))],
    )

    # Revs with the one-bit backend key bands by `GNSSSignals.get_band_id` (`:L1`).
    _axes_dc(dc, cap, ts) =
        let meas = NamedTuple{(get_band_id(_AXES_SIG),)}((
                Tracking.BandMeasurement(cap, _AXES_FS, 0.0Hz),
            ))
            @benchmarkable Tracking.downconvert_and_correlate($dc, $meas, $ts)
        end

    # multi-signal (M=1) and multi-antenna (N=1) via the uniform pipeline.
    for n_signals in (3,)
        cap = _int16_capture(_AXES_NSAMP)
        g = SUITE["track! Int16 vs Float32"]["GPS L1CA, 1 sat, $n_signals signals @ 5 MHz"]
        g["Float32"] = _axes_dc(
            _make_cpu_threaded_dc(_AXES_FS),
            cap,
            _axes_multisignal_state(n_signals),
        )
        g["Int16"] =
            _axes_dc(_make_int16_threaded_dc(2^11), cap, _axes_multisignal_state(n_signals))
        g["OneBit"] = _axes_dc(
            Tracking.OneBitThreadedDownconvertAndCorrelator(),
            cap,
            _axes_multisignal_state(n_signals),
        )
        if isdefined(Tracking, :TwoBitThreadedDownconvertAndCorrelator)
            g["TwoBit"] = _axes_dc(
                Tracking.TwoBitThreadedDownconvertAndCorrelator(),
                cap,
                _axes_multisignal_state(n_signals),
            )
        end
    end
    for num_ants in (4,)
        cap = _int16_capture_mat(_AXES_NSAMP, num_ants)
        g = SUITE["track! Int16 vs Float32"]["GPS L1CA, 1 sat, $num_ants antennas @ 5 MHz"]
        g["Float32"] = _axes_dc(
            _make_cpu_threaded_dc(_AXES_FS),
            cap,
            _axes_multiantenna_state(num_ants),
        )
        g["Int16"] =
            _axes_dc(_make_int16_threaded_dc(2^11), cap, _axes_multiantenna_state(num_ants))
        g["OneBit"] = _axes_dc(
            Tracking.OneBitThreadedDownconvertAndCorrelator(),
            cap,
            _axes_multiantenna_state(num_ants),
        )
        if isdefined(Tracking, :TwoBitThreadedDownconvertAndCorrelator)
            g["TwoBit"] = _axes_dc(
                Tracking.TwoBitThreadedDownconvertAndCorrelator(),
                cap,
                _axes_multiantenna_state(num_ants),
            )
        end
    end

    # dynamic (runtime AbstractVector) tap count — kernel level, all three backends.
    let g = SUITE["track! Int16 vs Float32"]["dynamic taps @ 5 MHz (kernel)"]
        code_doppler = 1000.0Hz * GNSSSignals.get_code_center_frequency_ratio(_AXES_SIG)
        code_frequency = code_doppler + get_code_frequency(_AXES_SIG)
        correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
        dyn_shifts =
            collect(get_correlator_sample_shifts(correlator, _AXES_FS, code_frequency))
        cap = _int16_capture(_AXES_NSAMP)
        span = maximum(dyn_shifts) - minimum(dyn_shifts)
        # Float32 fused kernel, dynamic-shifts (11-arg) path with hoisted SoA tiles.
        code_replica_f = Vector{Int8}(undef, _AXES_NSAMP + span)
        _gen_code_replica!(
            code_replica_f,
            _AXES_SIG,
            code_frequency,
            _AXES_FS,
            10.5,
            1,
            _AXES_NSAMP,
            dyn_shifts,
            1,
        )
        tile_re = Vector{Float32}(undef, _AXES_NSAMP)
        tile_im = Vector{Float32}(undef, _AXES_NSAMP)
        g["Float32"] = @benchmarkable Tracking.downconvert_and_correlate_fused!(
            $correlator,
            $cap,
            $code_replica_f,
            $dyn_shifts,
            $(1000.0Hz + 0.0Hz),
            $_AXES_FS,
            0.0,
            1,
            $_AXES_NSAMP,
            $tile_re,
            $tile_im,
        )
        # Int16 hybrid-blocked kernel, dynamic-shifts (AbstractVector) fallback.
        dc_i = _make_int16_dc(2^11)
        g["Int16"] = @benchmarkable Tracking._int16_hybrid_blocked!(
            $dc_i,
            $cap,
            NumAnts{1}(),
            $_AXES_SIG,
            1,
            $dyn_shifts,
            10.5,
            0.0,
            $code_frequency,
            $(1000.0Hz + 0.0Hz),
            $_AXES_FS,
            1,
            $_AXES_NSAMP,
        )
        # One-bit bit-wise kernel, dynamic-shifts (AbstractVector) fallback.
        dc_b = Tracking.OneBitDownconvertAndCorrelator()
        g["OneBit"] = @benchmarkable Tracking._onebit_hybrid_blocked!(
            $dc_b,
            $cap,
            NumAnts{1}(),
            $_AXES_SIG,
            1,
            $dyn_shifts,
            10.5,
            0.0,
            $code_frequency,
            $(1000.0Hz + 0.0Hz),
            $_AXES_FS,
            1,
            $_AXES_NSAMP,
        )
        # Two-bit bit-wise kernel, dynamic-shifts (AbstractVector) fallback.
        if isdefined(Tracking, :TwoBitDownconvertAndCorrelator)
            dc_t = Tracking.TwoBitDownconvertAndCorrelator()
            g["TwoBit"] = @benchmarkable Tracking._twobit_hybrid_blocked!(
                $dc_t,
                $cap,
                NumAnts{1}(),
                $_AXES_SIG,
                1,
                $dyn_shifts,
                10.5,
                0.0,
                $code_frequency,
                $(1000.0Hz + 0.0Hz),
                $_AXES_FS,
                1,
                $_AXES_NSAMP,
            )
        end
    end
end

# ── Int16 / OneBit `track!` in the base-vs-head regression table ───────────────
# The head-to-head group is rendered head-only, so register a few Int16/OneBit
# `track!` rows under `track`, which `bench_table.jl` diffs base-vs-head. Same
# rebuild-state-per-eval pattern as above.
if isdefined(Tracking, :Int16DownconvertAndCorrelator) && isdefined(Tracking, :track!)
    # The 100 ms case (per-call band pack, see above) is repeated here so the
    # base-vs-head table tracks it; it uses the locked composite capture and adds
    # a TwoBit leaf.
    for (key, systems, nsats_list, sfreq, nsamp, prn_max, onebit, twobit, capture) in (
        (
            "GPS L1CA, 8 sats, 5K @ 5 MHz",
            (GPSL1CA(),),
            [8],
            5e6Hz,
            5000,
            32,
            true,
            false,
            _int16_capture,
        ),
        (
            "Galileo E1B, 4 sats, 25K @ 25 MHz",
            (GalileoE1B(),),
            [4],
            25e6Hz,
            25000,
            50,
            false,
            false,
            _int16_capture,
        ),
        (
            "GPS L1CA, 8 sats, 500K @ 5 MHz",
            (GPSL1CA(),),
            [8],
            5e6Hz,
            500_000,
            32,
            true,
            true,
            n -> _l1ca_composite_capture(n, 5e6Hz),
        ),
    )
        sig16 = capture(nsamp)
        dc_i = _make_int16_dc(2^11)   # max_meas = 2^11 matches _int16_capture's ±2048
        SUITE["track"]["$key – in-place, Int16"] = @benchmarkable(
            Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_i),
            setup = (
                ts = first(
                    _make_steady_state_track_state(;
                        systems = $systems,
                        nsats_list = $nsats_list,
                        nsamp = $nsamp,
                        prn_max = $prn_max,
                        code_dop = 100.0,
                    ),
                )
            ),
            evals = 1,
        )
        if onebit && isdefined(Tracking, :OneBitDownconvertAndCorrelator)
            dc_b = Tracking.OneBitDownconvertAndCorrelator()
            SUITE["track"]["$key – in-place, OneBit"] = @benchmarkable(
                Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_b),
                setup = (
                    ts = first(
                        _make_steady_state_track_state(;
                            systems = $systems,
                            nsats_list = $nsats_list,
                            nsamp = $nsamp,
                            prn_max = $prn_max,
                            code_dop = 100.0,
                        ),
                    )
                ),
                evals = 1,
            )
        end
        if twobit && isdefined(Tracking, :TwoBitDownconvertAndCorrelator)
            dc_t = Tracking.TwoBitDownconvertAndCorrelator()
            SUITE["track"]["$key – in-place, TwoBit"] = @benchmarkable(
                Tracking.track!($sig16, ts, $sfreq; downconvert_and_correlator = $dc_t),
                setup = (
                    ts = first(
                        _make_steady_state_track_state(;
                            systems = $systems,
                            nsats_list = $nsats_list,
                            nsamp = $nsamp,
                            prn_max = $prn_max,
                            code_dop = 100.0,
                        ),
                    )
                ),
                evals = 1,
            )
        end
    end
end

# ── Per-signal noise estimation ──────────────────────────────────────────────
# `update_noise!` is the per-signal, per-chunk O(N) despread; it should sit near a
# per-satellite correlate, not a whole `track!`. `append_noise_observation!` is the
# hardware ingest path and must be allocation-free. The default `track` rows
# include the noise reference on purpose (end-to-end cost).
if isdefined(Tracking, :CorrelatorNoiseEstimator)
    function bench_update_noise(; num_samples = 4000, sampling_frequency = 4e6Hz)
        gnss_signal = GPSL1CA()
        estimator = Tracking.CorrelatorNoiseEstimator()
        measurement =
            Tracking.BandMeasurement(rand(ComplexF32, num_samples), sampling_frequency)
        context = Tracking.NoiseUpdateContext(
            gnss_signal,
            0,
            Tracking.CPUDownconvertAndCorrelator(),
        )
        @benchmarkable Tracking.update_noise!(
            $estimator,
            $measurement,
            1,
            $num_samples,
            $context,
        )
    end
    SUITE["noise estimation"]["update_noise! – 1 ms @ 4 MHz"] = bench_update_noise()
    SUITE["noise estimation"]["update_noise! – 1 ms @ 20 MHz"] =
        bench_update_noise(; num_samples = 20_000, sampling_frequency = 20e6Hz)

    # Antenna-array reference: despreads all `M` columns. Sub-linear in `M`
    # (`M = 4` ≈2× the 1-antenna row; bandwidth-bound and vectorised across
    # antennas); this row catches that ratio degrading, e.g. if the `M²` per-tap
    # `Σ b·bᴴ` pooling stops being dominated by the despread. `M = 4` is the
    # largest allocation-free size (an 8×8 `SMatrix` leaves the stack).
    if hasmethod(Tracking.CorrelatorNoiseEstimator, Tuple{}, (:num_ants,))
        function bench_update_noise_array(;
            num_samples = 4000,
            sampling_frequency = 4e6Hz,
            num_ants = 4,
        )
            gnss_signal = GPSL1CA()
            estimator = Tracking.CorrelatorNoiseEstimator(; num_ants = NumAnts(num_ants))
            measurement = Tracking.BandMeasurement(
                rand(ComplexF32, num_samples, num_ants),
                sampling_frequency,
            )
            context = Tracking.NoiseUpdateContext(
                gnss_signal,
                0,
                Tracking.CPUDownconvertAndCorrelator(),
            )
            @benchmarkable Tracking.update_noise!(
                $estimator,
                $measurement,
                1,
                $num_samples,
                $context,
            )
        end
        SUITE["noise estimation"]["update_noise! – 1 ms @ 4 MHz, 4 antennas"] =
            bench_update_noise_array()
    end

    function bench_append_noise_observation()
        estimator = Tracking.CorrelatorNoiseEstimator()
        observation = Tracking.noise_observation_from_samples(4000.0, 4000, 4e6Hz)
        # Fill the window first: the interesting cost is the steady state, where
        # every push is paired with a `popfirst!`.
        for _ = 1:2000
            Tracking.append_noise_observation!(estimator, observation)
        end
        @benchmarkable Tracking.append_noise_observation!($estimator, $observation)
    end
    SUITE["noise estimation"]["append_noise_observation! – steady state"] =
        bench_append_noise_observation()

    function bench_track_noise_ref(;
        signal_type = Float32,
        num_samples = 2000,
        sampling_frequency = 5e6Hz,
    )
        gnss_signal = GPSL1CA()
        downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
        track_state = TrackState(
            gnss_signal,
            [
                Tracking.TrackedSat(
                    gnss_signal,
                    1,
                    0.0,
                    1000Hz;
                    cn0_estimator = Tracking.NoiseRefCN0Estimator(),
                ),
            ],
        )
        signal = rand(Complex{signal_type}, num_samples)
        @benchmarkable Tracking.track!(
            $signal,
            $track_state,
            $sampling_frequency;
            downconvert_and_correlator = $downconvert_and_correlator,
        )
    end
    SUITE["noise estimation"]["track! – GPS L1CA, 1 sat, 2K @ 5 MHz with a noise reference"] =
        bench_track_noise_ref()
end
