using BenchmarkTools
using GNSSSignals
using Unitful: Hz
using CUDA
using Tracking

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

    arch_function = (
        CPU = (
            downconvert_and_correlator = CPUDownconvertAndCorrelator(),
            array_transform = Array,
        ),
        GPU = (
            downconvert_and_correlator = GPUDownconvertAndCorrelator(),
            array_transform = cu,
        ),
    )

    system_sats_state = SystemSatsState(
        system,
        [SatState(system, 1, sampling_frequency, code_phase, 1000.0Hz)];
    )
    multiple_system_sats_state = (system_sats_state,)

    track_state = TrackState(
        multiple_system_sats_state;
        num_samples = num_samples_signal,
        downconvert_and_correlator = arch_function[type].downconvert_and_correlator,
    )

    signal =
        arch_function[type].array_transform(rand(Complex{signal_type}, num_samples_signal))

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
        return @benchmarkable Tracking.downconvert_and_correlate(
            $signal,
            $track_state,
            $preferred_num_code_blocks_to_integrate,
            $sampling_frequency,
            $intermediate_frequency,
            $num_samples_signal,
            $(Val(sampling_frequency)),
        )
    end
end

function bench_track(;
    signal_type,
    num_samples_signal = 2000,
    sampling_frequency = 5e6Hz,
    system = GPSL1(),
)
    track_state = TrackState(
        system,
        [SatState(system, 1, sampling_frequency, 0.0, 1000Hz)];
        num_samples = num_samples_signal,
    )
    signal = rand(Complex{signal_type}, num_samples_signal)
    return if PACKAGE_VERSION <= v"0.15.2"
        @benchmarkable track($signal, $track_state, $sampling_frequency)
    else
        @benchmarkable track(
            $signal,
            $track_state,
            $sampling_frequency,
            maximum_expected_sampling_frequency = $(Val(sampling_frequency)),
        )
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