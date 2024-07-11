using BenchmarkTools, InteractiveUtils
@testset "Benchmark with $type" for type in (:CPU, :GPU)
    type == :GPU && !CUDA.functional() && return

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 1000
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
        gpsl1,
        [SatState(gpsl1, 1, sampling_frequency, code_phase, 1000.0Hz)];
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

    println("Benchmark downconvert and correlate with $type:")
    versioninfo()
    if (type == :GPU && CUDA.functional())
        display(CUDA.versioninfo())
        display(CUDA.device())
    end
    @btime Tracking.downconvert_and_correlate(
        $signal,
        $track_state,
        $next_system_sats_sample_params,
        $sampling_frequency,
        $intermediate_frequency,
        $num_samples_signal,
    )
end