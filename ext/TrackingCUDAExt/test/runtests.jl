module TrackingCUDAExtTests

using Test: @test, @testset, @inferred
using Unitful: Hz
using Tracking
using CUDA: CUDA, cu
using GNSSSignals: GPSL1, gen_code, get_code_frequency, get_code_center_frequency_ratio, get_code
using Pkg
using Bumper: SlabBuffer
using Tracking:
    CPUDownconvertAndCorrelator,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate,
    get_last_fully_integrated_correlator,
    track,
    get_code_phase,
    get_carrier_phase,
    get_carrier_doppler,
    get_code_doppler,
    get_last_fully_integrated_filtered_prompt,
    NumAnts,
    get_default_correlator

# Access extension types via Base.get_extension
const TrackingCUDAExt = Base.get_extension(Tracking, :TrackingCUDAExt)
const GPUDownconvertAndCorrelator = TrackingCUDAExt.GPUDownconvertAndCorrelator

# Helper function to check CUDA.jl version
function is_cuda_below_version(version_string::String)
    cuda_version = Pkg.installed()["CUDA"]
    return cuda_version < VersionNumber(version_string)
end

@testset "GPU Downconvert and Correlate" begin
    !CUDA.functional() && return

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

    downconvert_and_correlator =
        GPUDownconvertAndCorrelator(multiple_system_sats_state, num_samples_signal)

    track_state = TrackState(multiple_system_sats_state)

    preferred_num_code_blocks_to_integrate = 1

    signal = cu(
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

    signal = cu(
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

@testset "Track multiple signals with GPU" begin
    !CUDA.functional() && return
    gpsl1 = GPSL1()
    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + get_code_frequency(gpsl1)
    sampling_frequency = 4e6Hz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2

    num_samples = 4000
    num_ants = NumAnts(3)

    correlator = get_default_correlator(gpsl1, num_ants)
    sat_states =
        [SatState(gpsl1, 1, start_code_phase, carrier_doppler - 20Hz; num_ants, correlator)]

    system_sats_state = SystemSatsState(gpsl1, sat_states)

    # TODO: Why doesn't @inferred work here?
    downconvert_and_correlator =
        GPUDownconvertAndCorrelator((system_sats_state,), num_samples)

    track_state = @inferred TrackState(system_sats_state)

    signal =
        cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ start_carrier_phase) .*
        get_code.(
            gpsl1,
            code_frequency .* range ./ sampling_frequency .+ start_code_phase,
            prn,
        )
    signal_mat = cu(repeat(signal; outer = (1, 3)))

    track_state = @inferred track(
        signal_mat,
        track_state,
        sampling_frequency;
        downconvert_and_correlator,
    )

    iterations = 2000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    tracked_prompts = zeros(ComplexF64, iterations)
    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * i / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        code_phase =
            mod(code_frequency * 4000 * i / sampling_frequency + start_code_phase, 1023)
        signal =
            cis.(2π .* carrier_doppler .* range ./ sampling_frequency .+ carrier_phase) .*
            get_code.(
                gpsl1,
                code_frequency .* range ./ sampling_frequency .+ code_phase,
                prn,
            )
        signal_mat = cu(repeat(signal; outer = (1, 3)))
        track_state = @inferred track(
            signal_mat,
            track_state,
            sampling_frequency;
            downconvert_and_correlator,
        )
        comp_carrier_phase =
            mod2pi(
                2π * carrier_doppler * 4000 * (i + 1) / sampling_frequency +
                start_carrier_phase +
                π,
            ) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sampling_frequency + start_code_phase,
            1023,
        )
        tracked_code_phases[i] = get_code_phase(track_state)
        tracked_carrier_phases[i] = get_carrier_phase(track_state)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_state) / Hz
        tracked_code_dopplers[i] = get_code_doppler(track_state) / Hz
        tracked_prompts[i] = get_last_fully_integrated_filtered_prompt(track_state)
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 1e-2
    @test tracked_carrier_phases[end] + π ≈ carrier_phases[end] atol = 5e-5
end

end
