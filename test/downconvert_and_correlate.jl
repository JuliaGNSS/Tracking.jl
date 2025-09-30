module DownconvertAndCorrelateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using CUDA: CUDA, cu
using GNSSSignals: GPSL1, gen_code, get_code_frequency, get_code_center_frequency_ratio
using Pkg
using Bumper: SlabBuffer
using Tracking:
    CPUDownconvertAndCorrelator,
    GPUDownconvertAndCorrelator,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate,
    get_last_fully_integrated_correlator

# Helper function to check CUDA.jl version
function is_cuda_below_version(version_string::String)
    cuda_version = Pkg.installed()["CUDA"]
    return cuda_version < VersionNumber(version_string)
end

@testset "Downconvert and Correlator" begin
    downconvert_and_correlator = CPUDownconvertAndCorrelator(Val(5e6Hz))
    @test downconvert_and_correlator.buffer isa SlabBuffer
end

@testset "Downconvert and correlate with $type" for type in (:CPU, :GPU)
    type == :GPU && !CUDA.functional() && return

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    intermediate_frequency = 0.0Hz

    system_sats_state = SystemSatsState(
        gpsl1,
        [SatState(gpsl1, 1, code_phase, 1000.0Hz), SatState(gpsl1, 2, 11.0, 500.0Hz)];
    )
    multiple_system_sats_state = (system_sats_state,)

    maximum_expected_sampling_frequency = Val(sampling_frequency)

    downconvert_and_correlator =
        type == :CPU ? CPUDownconvertAndCorrelator(maximum_expected_sampling_frequency) :
        GPUDownconvertAndCorrelator(multiple_system_sats_state, num_samples_signal)

    array_transform = type == :CPU ? Array : cu

    track_state = TrackState(multiple_system_sats_state)

    preferred_num_code_blocks_to_integrate = 1

    signal = array_transform(
        gen_code(
            num_samples_signal,
            gpsl1,
            1,
            sampling_frequency,
            get_code_frequency(gpsl1) + 1000Hz * get_code_center_frequency_ratio(gpsl1),
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 1000.0Hz / sampling_frequency),
    )

    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    # GPU uses floating point arithmetic and might differ a little with the fixed point arithmetic
    if type == :CPU
        @test real.(
            get_last_fully_integrated_correlator(next_track_state, 1).accumulators
        ) ≈ [2921, 4949, 2917]
    else
        if is_cuda_below_version("5.8.0")
            @test real.(
                get_last_fully_integrated_correlator(next_track_state, 1).accumulators
            ) ≈ [2921, 4949, 2921]
        else
            # Unfortuntely, the texture index rounding is off with CUDA 5.9 TODO: report
            accumulator =
                real.(
                    get_last_fully_integrated_correlator(next_track_state, 1).accumulators
                )
            @test 2900 < accumulator[1] < 3000
            @test 4900 < accumulator[2] < 5000
            @test 2900 < accumulator[3] < 3000
        end
    end

    signal = array_transform(
        gen_code(
            num_samples_signal,
            gpsl1,
            2,
            sampling_frequency,
            get_code_frequency(gpsl1) + 500Hz * get_code_center_frequency_ratio(gpsl1),
            11.0,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 500.0Hz / sampling_frequency),
    )

    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    # GPU uses floating point arithmetic and might differ a little with the fixed point arithmetic
    if type == :CPU
        @test real.(
            get_last_fully_integrated_correlator(next_track_state, 2).accumulators
        ) ≈ [2919, 4947, 2915]
    else
        if is_cuda_below_version("5.8.0")
            @test real.(
                get_last_fully_integrated_correlator(next_track_state, 2).accumulators
            ) ≈ [2919, 4947, 2919]
        else
            # Unfortuntely, the texture index rounding is off with CUDA 5.9 TODO: report
            accumulator =
                real.(
                    get_last_fully_integrated_correlator(next_track_state, 2).accumulators
                )
            @test 2900 < accumulator[1] < 3000
            @test 4900 < accumulator[2] < 5000
            @test 2900 < accumulator[3] < 3000
        end
    end
end

end
