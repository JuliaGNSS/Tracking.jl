module DownconvertAndCorrelateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, gen_code, get_code_frequency, get_code_center_frequency_ratio
using Bumper: SlabBuffer
import Tracking
using Tracking:
    CPUDownconvertAndCorrelator,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate,
    get_correlator

@testset "Downconvert and Correlator" begin
    downconvert_and_correlator = CPUDownconvertAndCorrelator(Val(5e6Hz))
    @test downconvert_and_correlator.buffer isa SlabBuffer
end

@testset "Downconvert and correlate with CPU" begin
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
        CPUDownconvertAndCorrelator(maximum_expected_sampling_frequency)

    track_state = TrackState(multiple_system_sats_state)

    preferred_num_code_blocks_to_integrate = 1

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            1,
            sampling_frequency,
            get_code_frequency(gpsl1) + 1000Hz * get_code_center_frequency_ratio(gpsl1),
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 1000.0Hz / sampling_frequency)

    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    @test real.(get_correlator(next_track_state, 1).accumulators) ≈ [2921, 4949, 2917] rtol=1e-3

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            2,
            sampling_frequency,
            get_code_frequency(gpsl1) + 500Hz * get_code_center_frequency_ratio(gpsl1),
            11.0,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 500.0Hz / sampling_frequency)

    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    @test real.(get_correlator(next_track_state, 2).accumulators) ≈ [2919, 4947, 2915] rtol=1e-3
end

@testset "Fused downconvert and correlate" begin
    using StructArrays: StructVector
    using GNSSSignals: get_code_center_frequency_ratio, get_code_frequency, get_code_type
    using Tracking:
        EarlyPromptLateCorrelator,
        NumAnts,
        get_accumulators,
        get_correlator_sample_shifts

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    carrier_doppler = 1000.0Hz
    intermediate_frequency = 0.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + intermediate_frequency
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    correlator = EarlyPromptLateCorrelator()
    sample_shifts = get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)

    # Generate a realistic signal
    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal - 1)) * carrier_doppler / sampling_frequency)

    # Generate code replica
    code_replica_length = num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples_signal,
        sample_shifts,
        prn,
        Val(sampling_frequency),
    )

    # Downconvert buffer for split path
    downconvert_buffer = StructVector{ComplexF32}(undef, num_samples_signal)

    carrier_phase = 0.0
    start_sample = 1
    num_samples = num_samples_signal

    # Split path: downconvert then correlate
    Tracking.downconvert!(
        downconvert_buffer,
        signal,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        start_sample,
        num_samples,
        NumAnts{1}(),
    )
    split_correlator = Tracking.correlate(
        correlator,
        downconvert_buffer,
        sample_shifts,
        code_replica,
        start_sample,
        num_samples,
    )

    # Fused path
    fused_correlator = Tracking.downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        start_sample,
        num_samples,
    )

    @test get_accumulators(fused_correlator) ≈ get_accumulators(split_correlator) rtol = 1e-3
end

end
