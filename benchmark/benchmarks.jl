using BenchmarkTools
using GNSSSignals
using Unitful: Hz
using CUDA
using Tracking

function bench_downconvert_and_correlate(
    type;
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

    system_sats_sample_params = Tracking.init_sample_params(multiple_system_sats_state, 1)
    next_system_sats_sample_params = Tracking.calc_sample_params(
        multiple_system_sats_state,
        system_sats_sample_params,
        num_samples_signal,
        sampling_frequency,
        1,
    )

    signal = arch_function[type].array_transform(randn(ComplexF32, num_samples_signal))

    type == :GPU && !CUDA.functional() && CUDA.versioninfo()
    @benchmarkable Tracking.downconvert_and_correlate(
        $signal,
        $track_state,
        $next_system_sats_sample_params,
        $sampling_frequency,
        $intermediate_frequency,
        $num_samples_signal,
    ) evals = 10 samples = 10000
end

const SUITE = BenchmarkGroup()

SUITE["downconvert and correlate"]["CPU"] = bench_downconvert_and_correlate(:CPU)
if (CUDA.functional())
    SUITE["downconvert and correlate"]["GPU"] = bench_downconvert_and_correlate(:GPU)
end