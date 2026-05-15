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
    SatState,
    TrackState,
    downconvert_and_correlate,
    BitBuffer
using StaticArrays

# Branch-portable alias for the per-system tracking container. On
# Tracking 2.0+ this is `TrackedSystem`; on master it was
# `SystemSatsState`. The benchmark script runs against both via
# AirspeedVelocity's `--bench-on=$HEAD_SHA`, so detect at load time.
const _TrackedSystem = isdefined(Tracking, :TrackedSystem) ?
    Tracking.TrackedSystem : Tracking.SystemSatsState

# Field name on `TrackState` for the per-system tuple. Master named it
# `multiple_system_sats_state`; on the new branch it's `tracked_systems`.
const _SYSTEMS_FIELD = isdefined(Tracking, :TrackedSystem) ?
    :tracked_systems : :multiple_system_sats_state

@inline _get_systems(track_state) = getfield(track_state, _SYSTEMS_FIELD)
@inline _track_state_with_systems(track_state, systems) =
    TrackState(track_state; NamedTuple{(_SYSTEMS_FIELD,)}((systems,))...)

const SUITE = BenchmarkGroup()

# ── Helper: set up common benchmark state ──────────────────────────────────

function setup_benchmark(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1CA(),
    num_ants = 1,
)
    code_phase = 10.5
    carrier_doppler = 1000.0Hz
    code_doppler = carrier_doppler * GNSSSignals.get_code_center_frequency_ratio(system)
    code_frequency = code_doppler + get_code_frequency(system)

    correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants))
    static_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)

    signal =
        num_ants == 1 ? rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    code_replica = Vector{get_code_type(system)}(
        undef,
        num_samples + maximum(static_shifts) - minimum(static_shifts),
    )
    _gen_code_replica!(
        code_replica,
        system,
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
        system,
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
    code_replica, system, code_frequency, sampling_frequency,
    code_phase, start_sample, num_samples, sample_shifts, prn,
)
    if _NEEDS_VAL
        gen_code_replica!(
            code_replica, system, code_frequency, sampling_frequency,
            code_phase, start_sample, num_samples, sample_shifts, prn,
            Val(sampling_frequency),
        )
    else
        gen_code_replica!(
            code_replica, system, code_frequency, sampling_frequency,
            code_phase, start_sample, num_samples, sample_shifts, prn,
        )
    end
end

# ── High-level downconvert_and_correlate (full pipeline) ───────────────────

function bench_downconvert_and_correlate(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1CA(),
    num_ants = 1,
)
    downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
    track_state = TrackState(
        system,
        [SatState(system, 1, 10.5, 1000.0Hz; num_ants = NumAnts(num_ants))],
    )
    signal =
        num_ants == 1 ? rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    @benchmarkable Tracking.downconvert_and_correlate(
        $downconvert_and_correlator,
        $signal,
        $track_state,
        1,
        $sampling_frequency,
        $(0.0Hz),
    )
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
    system = GPSL1CA()
    downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
    track_state = TrackState(system, [SatState(system, 1, 0.0, 1000Hz)])
    signal = rand(Complex{signal_type}, num_samples)
    @benchmarkable track(
        $signal,
        $track_state,
        $sampling_frequency;
        downconvert_and_correlator = $downconvert_and_correlator,
    )
end
SUITE["track"]["Float32"] = bench_track()

# In-place track! (only on branches that define it). Mirrors bench_track so
# the comparison report shows them side by side.
function bench_track_inplace(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
)
    system = GPSL1CA()
    downconvert_and_correlator = _make_cpu_dc(sampling_frequency)
    track_state = TrackState(system, [SatState(system, 1, 0.0, 1000Hz)])
    signal = rand(Complex{signal_type}, num_samples)
    @benchmarkable Tracking.track!(
        $signal,
        $track_state,
        $sampling_frequency;
        downconvert_and_correlator = $downconvert_and_correlator,
    )
end
if isdefined(Tracking, :track!)
    SUITE["track!"]["Float32"] = bench_track_inplace()
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

# Branch-portable per-system container construction. Master accepts a
# raw Vector{SatState} via `SystemSatsState(sys, sats)`; the wrapper
# branch needs an estimator to wrap each sat into a TrackedSat first
# via `TrackedSystem(estimator, sys, sats)`.
function _build_tracked_system(sys, sats::Vector{<:SatState})
    if isdefined(Tracking, :TrackedSat)
        return _TrackedSystem(Tracking.ConventionalAssistedPLLAndDLL(), sys, sats)
    else
        return _TrackedSystem(sys, sats)
    end
end

# Branch-portable "set bit_buffer.found = true" on whatever the dict holds.
# Master holds SatState directly; the wrapper branch holds TrackedSat.
_with_found_bit_buffer(s::SatState, bb) = SatState(s; bit_buffer = bb)
if isdefined(Tracking, :TrackedSat)
    _with_found_bit_buffer(t::Tracking.TrackedSat, bb) =
        Tracking.TrackedSat(SatState(t.sat_state; bit_buffer = bb), t.estimator_state)
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
        sats = [SatState(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for i = 1:ns]
        push!(all_sss, _build_tracked_system(sys, sats))
        total_sats += ns
    end
    TrackState(Tuple(all_sss)), rand(ComplexF32, nsamp), total_sats
end

# Build a `TrackState` for the given system mix with `bit_buffer.found = true`
# on every sat, ready to feed `track` / `track!` for steady-state-style
# benchmarks.
function _make_steady_state_track_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    ts, signal, _ =
        _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    found_bb = BitBuffer(UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0)
    new_mss = map(_get_systems(ts)) do sss
        new_sats = map(s -> _with_found_bit_buffer(s, found_bb), sss.states)
        _TrackedSystem(sss, new_sats)
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
const _TRACK_BENCH_CASES = let gpsl1 = GPSL1CA(), gal = GalileoE1B()
    [
        ("L1 8sat/5K",       (systems = (gpsl1,),     nsats_list = [8],    sfreq = 5e6Hz,  nsamp = 5000)),
        ("E1B 4sat/25K",     (systems = (gal,),       nsats_list = [4],    sfreq = 25e6Hz, nsamp = 25000, prn_max = 50, code_dop = 100.0)),
        ("8L1+8E1B/25K",     (systems = (gpsl1, gal), nsats_list = [8, 8], sfreq = 25e6Hz, nsamp = 25000, prn_max = 50, code_dop = 100.0)),
    ]
end

for (key, kw) in _TRACK_BENCH_CASES
    SUITE["track"][key] = bench_track_steady_state(false, false; kw...)
    SUITE["track-threaded"][key] = bench_track_steady_state(false, true; kw...)
    if isdefined(Tracking, :track!)
        SUITE["track!"][key] = bench_track_steady_state(true, false; kw...)
        SUITE["track!-threaded"][key] = bench_track_steady_state(true, true; kw...)
    end
end
