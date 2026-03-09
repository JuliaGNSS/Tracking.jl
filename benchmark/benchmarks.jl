using BenchmarkTools
using GNSSSignals
using GNSSSignals: GalileoE1B
using Unitful: Hz
using CUDA
using Tracking
using Tracking: EarlyPromptLateCorrelator, get_correlator_sample_shifts, get_code_type,
    NumAnts, gen_code_replica!, SystemSatsState, SatState, TrackState,
    downconvert_and_correlate
using StaticArrays

function bench_downconvert_and_correlate(
    type;
    signal_type = Float32,
    num_samples_signal = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1(),
)
    type == :GPU && !CUDA.functional() && return

    code_phase = 10.5
    intermediate_frequency = 0.0Hz

    maximum_expected_sampling_frequency = Val(sampling_frequency)

    system_sats_state =
        PACKAGE_VERSION <= v"0.15.6" ?
        SystemSatsState(
            system,
            [SatState(system, 1, sampling_frequency, code_phase, 1000.0Hz)],
        ) : SystemSatsState(system, [SatState(system, 1, code_phase, 1000.0Hz)])

    multiple_system_sats_state = (system_sats_state,)

    downconvert_and_correlator =
        PACKAGE_VERSION <= v"0.15.4" ?
        (type == :CPU ? CPUDownconvertAndCorrelator() : GPUDownconvertAndCorrelator()) :
        (
            PACKAGE_VERSION <= v"0.15.5" ?
            (
                type == :CPU ?
                CPUDownconvertAndCorrelator(
                    maximum_expected_sampling_frequency,
                    multiple_system_sats_state,
                    num_samples_signal,
                ) :
                GPUDownconvertAndCorrelator(multiple_system_sats_state, num_samples_signal)
            ) :
            (
                type == :CPU ?
                CPUDownconvertAndCorrelator(maximum_expected_sampling_frequency) :
                GPUDownconvertAndCorrelator()
            )
        )

    array_transform = type == :CPU ? Array : cu

    track_state =
        PACKAGE_VERSION <= v"0.15.4" ?
        TrackState(
            multiple_system_sats_state;
            num_samples = num_samples_signal,
            downconvert_and_correlator,
        ) :
        (
            PACKAGE_VERSION <= v"0.15.5" ?
            TrackState(
                multiple_system_sats_state;
                num_samples = num_samples_signal,
                downconvert_and_correlator,
                maximum_expected_sampling_frequency,
            ) : TrackState(multiple_system_sats_state)
        )

    signal = array_transform(rand(Complex{signal_type}, num_samples_signal))

    type == :GPU && !CUDA.functional() && CUDA.versioninfo()
    if PACKAGE_VERSION <= v"0.15.3"
        system_sats_sample_params =
            Tracking.init_sample_params(multiple_system_sats_state, 1)
        next_system_sats_sample_params = Tracking.calc_sample_params(
            multiple_system_sats_state,
            system_sats_sample_params,
            num_samples_signal,
            sampling_frequency,
            1,
        )

        return if PACKAGE_VERSION <= v"0.15.2"
            @benchmarkable Tracking.downconvert_and_correlate(
                $signal,
                $track_state,
                $next_system_sats_sample_params,
                $sampling_frequency,
                $intermediate_frequency,
                $num_samples_signal,
            )
        else
            @benchmarkable Tracking.downconvert_and_correlate(
                $signal,
                $track_state,
                $next_system_sats_sample_params,
                $sampling_frequency,
                $intermediate_frequency,
                $num_samples_signal,
                $(Val(sampling_frequency)),
            )
        end
    else
        preferred_num_code_blocks_to_integrate = 1
        return if PACKAGE_VERSION <= v"0.15.4"
            @benchmarkable Tracking.downconvert_and_correlate(
                $signal,
                $track_state,
                $preferred_num_code_blocks_to_integrate,
                $sampling_frequency,
                $intermediate_frequency,
                $num_samples_signal,
                $(Val(sampling_frequency)),
            )
        else
            return if PACKAGE_VERSION <= v"0.15.5"
                @benchmarkable Tracking.downconvert_and_correlate(
                    $signal,
                    $track_state,
                    $preferred_num_code_blocks_to_integrate,
                    $sampling_frequency,
                    $intermediate_frequency,
                    $num_samples_signal,
                )
            else
                @benchmarkable Tracking.downconvert_and_correlate(
                    $downconvert_and_correlator,
                    $signal,
                    $track_state,
                    $preferred_num_code_blocks_to_integrate,
                    $sampling_frequency,
                    $intermediate_frequency,
                )
            end
        end
    end
end

function bench_track(;
    signal_type,
    num_samples_signal = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1(),
)
    track_state =
        PACKAGE_VERSION <= v"0.15.4" ?
        TrackState(
            system,
            [SatState(system, 1, sampling_frequency, 0.0, 1000Hz)];
            num_samples = num_samples_signal,
        ) :
        (
            PACKAGE_VERSION <= v"0.15.5" ?
            TrackState(
                system,
                [SatState(system, 1, sampling_frequency, 0.0, 1000Hz)];
                num_samples = num_samples_signal,
                maximum_expected_sampling_frequency = Val(sampling_frequency),
            ) :
            PACKAGE_VERSION <= v"0.15.6" ?
            TrackState(system, [SatState(system, 1, sampling_frequency, 0.0, 1000Hz)]) :
            TrackState(system, [SatState(system, 1, 0.0, 1000Hz)])
        )
    signal = rand(Complex{signal_type}, num_samples_signal)
    return if PACKAGE_VERSION <= v"0.15.2"
        @benchmarkable track($signal, $track_state, $sampling_frequency)
    else
        return if PACKAGE_VERSION <= v"0.15.4"
            @benchmarkable track(
                $signal,
                $track_state,
                $sampling_frequency,
                maximum_expected_sampling_frequency = $(Val(sampling_frequency)),
            )
        else
            return if PACKAGE_VERSION <= v"0.15.5"
                @benchmarkable track($signal, $track_state, $sampling_frequency)
            else
                downconvert_and_correlator =
                    CPUDownconvertAndCorrelator(Val(sampling_frequency))
                @benchmarkable track(
                    $signal,
                    $track_state,
                    $sampling_frequency;
                    downconvert_and_correlator = $downconvert_and_correlator,
                )
            end
        end
    end
end

const SUITE = BenchmarkGroup()

foreach((Int16, Int32, Float32, Float64)) do signal_type
    SUITE["downconvert and correlate"]["CPU"][string(signal_type)] =
        bench_downconvert_and_correlate(:CPU; signal_type)
end
if (CUDA.functional())
    SUITE["downconvert and correlate"]["GPU"] = bench_downconvert_and_correlate(:GPU)
end
SUITE["track"]["Float32"] = bench_track(; signal_type = Float32)

# ── Multi-satellite benchmarks (threaded if available, CPU fallback) ──────

function _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max=32, code_dop=1000.0)
    all_sss = []
    total_sats = 0
    for (si, sys) in enumerate(systems)
        ns = nsats_list[min(si, length(nsats_list))]
        pm = sys isa GPSL1 ? 32 : prn_max
        cd = sys isa GPSL1 ? 1000.0 : code_dop
        sats = [SatState(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for i = 1:ns]
        push!(all_sss, SystemSatsState(sys, sats))
        total_sats += ns
    end
    TrackState(Tuple(all_sss)), rand(ComplexF32, nsamp), total_sats
end

function bench_multi_sat(; systems, nsats_list, sfreq, nsamp, prn_max=32, code_dop=1000.0)
    ts, signal, total_sats = _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    if isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        dc = Tracking.CPUThreadedDownconvertAndCorrelator(
            systems, Val(sfreq); max_sats=max(total_sats, 4), max_num_samples=nsamp,
        )
    else
        dc = CPUDownconvertAndCorrelator(Val(sfreq))
    end
    @benchmarkable downconvert_and_correlate($dc, $signal, $ts, 1, $sfreq, $(0.0Hz))
end

let gpsl1 = GPSL1(), gal = GalileoE1B()
    SUITE["multi-sat"]["L1 8sat/5K"] =
        bench_multi_sat(; systems=(gpsl1,), nsats_list=[8], sfreq=5e6Hz, nsamp=5000)
    SUITE["multi-sat"]["E1B 4sat/25K"] =
        bench_multi_sat(; systems=(gal,), nsats_list=[4], sfreq=25e6Hz, nsamp=25000, prn_max=50, code_dop=100.0)
    SUITE["multi-sat"]["8L1+8E1B/25K"] =
        bench_multi_sat(; systems=(gpsl1, gal), nsats_list=[8, 8], sfreq=25e6Hz, nsamp=25000, prn_max=50, code_dop=100.0)
end
