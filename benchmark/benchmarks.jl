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
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate,
    BitBuffer
using StaticArrays

const SUITE = BenchmarkGroup()

# ── Helper: set up common benchmark state ──────────────────────────────────

function setup_benchmark(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1(),
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
    system = GPSL1(),
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
    system = GPSL1()
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
    system = GPSL1()
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

# ── Multi-satellite benchmarks (threaded if available, CPU fallback) ──────

# Branch-portable SystemSatsState construction. Master accepts a raw
# Vector{SatState}; the wrapper branch needs an estimator to wrap each sat
# into a TrackedSat first.
function _build_system_sats_state(sys, sats::Vector{<:SatState})
    if isdefined(Tracking, :TrackedSat)
        return SystemSatsState(Tracking.ConventionalAssistedPLLAndDLL(), sys, sats)
    else
        return SystemSatsState(sys, sats)
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
        pm = sys isa GPSL1 ? 32 : prn_max
        cd = sys isa GPSL1 ? 1000.0 : code_dop
        sats = [SatState(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for i = 1:ns]
        push!(all_sss, _build_system_sats_state(sys, sats))
        total_sats += ns
    end
    TrackState(Tuple(all_sss)), rand(ComplexF32, nsamp), total_sats
end

function bench_multi_sat(
    threaded::Bool;
    systems,
    nsats_list,
    sfreq,
    nsamp,
    prn_max = 32,
    code_dop = 1000.0,
)
    ts, signal, _ =
        _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    if threaded && isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        dc = _make_cpu_threaded_dc(sfreq)
    else
        dc = _make_cpu_dc(sfreq)
    end
    @benchmarkable downconvert_and_correlate($dc, $signal, $ts, 1, $sfreq, $(0.0Hz))
end

# ── Full track loop with multiple satellites in steady-state ──────────────
# Exercises estimate_dopplers_and_filter_prompt with bit_buffer.found == true

function bench_track_steady_state(;
    systems,
    nsats_list,
    sfreq,
    nsamp,
    prn_max = 32,
    code_dop = 1000.0,
)
    ts, signal, total_sats =
        _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    dc = _make_cpu_dc(sfreq)
    # Set bit_buffer.found = true to simulate steady-state tracking
    # (random signal data never triggers bit detection on its own)
    found_bb = BitBuffer(UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0)
    new_mss = map(ts.multiple_system_sats_state) do sss
        new_sats = map(s -> _with_found_bit_buffer(s, found_bb), sss.states)
        SystemSatsState(sss, new_sats)
    end
    ts = TrackState(ts; multiple_system_sats_state = new_mss)
    @benchmarkable track($signal, $ts, $sfreq; downconvert_and_correlator = $dc)
end

SUITE["track steady-state"]["L1 8sat/5K"] = bench_track_steady_state(;
    systems = (GPSL1(),), nsats_list = [8], sfreq = 5e6Hz, nsamp = 5000,
)

# In-place track! steady-state (only on branches that define it). Reuses the
# steady-state preparation done by `bench_track_steady_state`'s setup
# (BitBuffer.found = true) by simply re-running the same setup pipeline,
# except it returns a `track!` benchmarkable instead of `track`.
#
# Pass `threaded = true` for the threaded backend variant.
function bench_track_inplace_steady_state(
    threaded::Bool = false;
    systems,
    nsats_list,
    sfreq,
    nsamp,
    prn_max = 32,
    code_dop = 1000.0,
)
    ts, signal, _ =
        _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    dc = if threaded && isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        _make_cpu_threaded_dc(sfreq)
    else
        _make_cpu_dc(sfreq)
    end
    found_bb = BitBuffer(UInt128(0), 20, true, UInt128(0), 0, complex(0.0, 0.0), 0)
    new_mss = map(ts.multiple_system_sats_state) do sss
        new_sats = map(s -> _with_found_bit_buffer(s, found_bb), sss.states)
        SystemSatsState(sss, new_sats)
    end
    ts = TrackState(ts; multiple_system_sats_state = new_mss)
    @benchmarkable Tracking.track!(
        $signal, $ts, $sfreq; downconvert_and_correlator = $dc,
    )
end

if isdefined(Tracking, :track!)
    SUITE["track! steady-state"]["L1 8sat/5K"] = bench_track_inplace_steady_state(
        false;
        systems = (GPSL1(),), nsats_list = [8], sfreq = 5e6Hz, nsamp = 5000,
    )
    if isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        SUITE["track! steady-state"]["L1 8sat/5K threaded"] =
            bench_track_inplace_steady_state(
                true;
                systems = (GPSL1(),), nsats_list = [8], sfreq = 5e6Hz, nsamp = 5000,
            )
    end
end

let gpsl1 = GPSL1(), gal = GalileoE1B()
    SUITE["multi-sat"]["L1 8sat/5K"] = bench_multi_sat(
        false;
        systems = (gpsl1,),
        nsats_list = [8],
        sfreq = 5e6Hz,
        nsamp = 5000,
    )
    SUITE["multi-sat"]["E1B 4sat/25K"] = bench_multi_sat(
        false;
        systems = (gal,),
        nsats_list = [4],
        sfreq = 25e6Hz,
        nsamp = 25000,
        prn_max = 50,
        code_dop = 100.0,
    )
    SUITE["multi-sat"]["8L1+8E1B/25K"] = bench_multi_sat(
        false;
        systems = (gpsl1, gal),
        nsats_list = [8, 8],
        sfreq = 25e6Hz,
        nsamp = 25000,
        prn_max = 50,
        code_dop = 100.0,
    )
    SUITE["multi-sat-threaded"]["L1 8sat/5K"] = bench_multi_sat(
        true;
        systems = (gpsl1,),
        nsats_list = [8],
        sfreq = 5e6Hz,
        nsamp = 5000,
    )
    SUITE["multi-sat-threaded"]["E1B 4sat/25K"] = bench_multi_sat(
        true;
        systems = (gal,),
        nsats_list = [4],
        sfreq = 25e6Hz,
        nsamp = 25000,
        prn_max = 50,
        code_dop = 100.0,
    )
    SUITE["multi-sat-threaded"]["8L1+8E1B/25K"] = bench_multi_sat(
        true;
        systems = (gpsl1, gal),
        nsats_list = [8, 8],
        sfreq = 25e6Hz,
        nsamp = 25000,
        prn_max = 50,
        code_dop = 100.0,
    )
end
