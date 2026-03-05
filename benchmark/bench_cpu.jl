# CPU downconvert-and-correlate + track benchmarks

using BenchmarkTools
using GNSSSignals
using Unitful: Hz
using Tracking
using Tracking:
    CPUDownconvertAndCorrelator,
    EarlyPromptLateCorrelator,
    NumAnts,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate,
    get_correlator_sample_shifts,
    get_code_type,
    gen_code_replica!
using StaticArrays

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
    static_shifts = get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)

    signal = num_ants == 1 ?
        rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    code_replica = Vector{get_code_type(system)}(
        undef, num_samples + maximum(static_shifts) - minimum(static_shifts),
    )
    gen_code_replica!(
        code_replica, system, code_frequency, sampling_frequency,
        code_phase, 1, num_samples, static_shifts, 1, Val(sampling_frequency),
    )

    return (;
        correlator, signal, code_replica,
        static_shifts, dynamic_shifts,
        sampling_frequency, carrier_doppler, code_phase,
        system, num_samples,
    )
end

# ── High-level downconvert_and_correlate (full pipeline) ───────────────────

function bench_downconvert_and_correlate(;
    signal_type = Float32,
    num_samples = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1(),
    num_ants = 1,
)
    downconvert_and_correlator = CPUDownconvertAndCorrelator(Val(sampling_frequency))
    system_sats_state = SystemSatsState(
        system,
        [SatState(system, 1, 10.5, 1000.0Hz; num_ants = NumAnts(num_ants))],
    )
    track_state = TrackState((system_sats_state,))
    signal = num_ants == 1 ?
        rand(Complex{signal_type}, num_samples) :
        rand(Complex{signal_type}, num_samples, num_ants)

    @benchmarkable Tracking.downconvert_and_correlate(
        $downconvert_and_correlator, $signal, $track_state, 1,
        $sampling_frequency, $(0.0Hz),
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
    # Warmup to trigger compilation
    Tracking.downconvert_and_correlate_fused!(
        s.correlator, s.signal, s.code_replica, sample_shifts,
        s.carrier_doppler, s.sampling_frequency, 0.0, 1, s.num_samples,
    )
    @benchmarkable Tracking.downconvert_and_correlate_fused!(
        $(s.correlator), $(s.signal), $(s.code_replica), $sample_shifts,
        $(s.carrier_doppler), $(s.sampling_frequency), 0.0, 1, $(s.num_samples),
    )
end

function cpu_suite()
    suite = BenchmarkGroup()

    # Full pipeline: CPU, various signal types
    foreach((Int16, Int32, Float32, Float64)) do signal_type
        suite["downconvert and correlate"]["CPU"][string(signal_type)] =
            bench_downconvert_and_correlate(; signal_type)
    end

    # Full pipeline: multi-antenna
    suite["downconvert and correlate"]["CPU"]["Float32 4ant"] =
        bench_downconvert_and_correlate(; num_ants = 4)
    suite["downconvert and correlate"]["CPU"]["Int16 4ant"] =
        bench_downconvert_and_correlate(; signal_type = Int16, num_ants = 4)

    # Full pipeline: track
    system = GPSL1()
    sampling_frequency = 5e6Hz
    downconvert_and_correlator = CPUDownconvertAndCorrelator(Val(sampling_frequency))
    track_state = TrackState(system, [SatState(system, 1, 0.0, 1000Hz)])
    signal = rand(ComplexF32, 2000)
    suite["track"]["Float32"] = @benchmarkable track(
        $signal, $track_state, $sampling_frequency;
        downconvert_and_correlator = $downconvert_and_correlator,
    )

    # Fused kernel microbenchmarks (only available on branches with the fused kernel)
    if isdefined(Tracking, :downconvert_and_correlate_fused!)
        suite["fused kernel"]["1-ant static taps"]  = bench_fused_kernel(; shifts = :static)
        suite["fused kernel"]["1-ant dynamic taps"] = bench_fused_kernel(; shifts = :dynamic)
        suite["fused kernel"]["4-ant static taps"]  = bench_fused_kernel(; num_ants = 4, shifts = :static)
        suite["fused kernel"]["4-ant dynamic taps"] = bench_fused_kernel(; num_ants = 4, shifts = :dynamic)
    end

    suite
end
