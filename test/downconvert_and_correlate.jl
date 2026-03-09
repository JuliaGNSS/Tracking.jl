module DownconvertAndCorrelateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1, gen_code, get_code_frequency, get_code_center_frequency_ratio, get_code_type
using Bumper: SlabBuffer
import Tracking
using Tracking:
    AbstractCorrelator,
    CPUDownconvertAndCorrelator,
    EarlyPromptLateCorrelator,
    NumAnts,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate,
    get_accumulators,
    get_correlator,
    get_correlator_sample_shifts,
    update_accumulator

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
    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    correlator = EarlyPromptLateCorrelator()
    sample_shifts = get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal - 1)) * carrier_doppler / sampling_frequency)

    code_replica_length = num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica, gpsl1, code_frequency, sampling_frequency,
        code_phase, 1, num_samples_signal, sample_shifts, prn,
        Val(sampling_frequency),
    )

    # Scalar reference: plain Julia loop, no SIMD, no tricks
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i in 1:num_samples_signal
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = signal[i] * carrier
        for tap in 1:NC
            code_idx = i + sample_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Static (@generated) path
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator, signal, code_replica, sample_shifts,
        carrier_frequency, sampling_frequency, 0.0, 1, num_samples_signal,
    )

    @test get_accumulators(fused_static) ≈ ref_accumulators rtol = 1e-4

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator, signal, code_replica, dynamic_shifts,
        carrier_frequency, sampling_frequency, 0.0, 1, num_samples_signal,
    )

    @test get_accumulators(fused_dynamic) ≈ ref_accumulators rtol = 1e-4

    # Static and dynamic paths should produce nearly identical results
    @test get_accumulators(fused_static) ≈ get_accumulators(fused_dynamic) rtol = 1e-6
end

@testset "Fused downconvert and correlate with Vector accumulators" begin
    # Minimal correlator with Vector (not SVector) accumulators
    struct DynamicCorrelator <: AbstractCorrelator{1}
        accumulators::Vector{ComplexF64}
        shifts::Vector{Int}
    end
    Tracking.get_accumulators(c::DynamicCorrelator) = c.accumulators
    Tracking.update_accumulator(c::DynamicCorrelator, acc) =
        DynamicCorrelator(collect(acc), c.shifts)

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    # Use the same shifts as EarlyPromptLateCorrelator but as a Vector
    epl = EarlyPromptLateCorrelator()
    static_shifts = get_correlator_sample_shifts(epl, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)

    dynamic_correlator = DynamicCorrelator(zeros(ComplexF64, 3), dynamic_shifts)

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal - 1)) * carrier_doppler / sampling_frequency)

    code_replica_length = num_samples_signal + maximum(dynamic_shifts) - minimum(dynamic_shifts)
    code_replica = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica, gpsl1, code_frequency, sampling_frequency,
        code_phase, 1, num_samples_signal, static_shifts, prn,
        Val(sampling_frequency),
    )

    # Scalar reference
    NC = length(dynamic_shifts)
    min_shift = minimum(dynamic_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i in 1:num_samples_signal
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = signal[i] * carrier
        for tap in 1:NC
            code_idx = i + dynamic_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Test: use the dynamic Vector accumulator path
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        dynamic_correlator, signal, code_replica, dynamic_shifts,
        carrier_frequency, sampling_frequency, 0.0, 1, num_samples_signal,
    )

    @test get_accumulators(fused_dynamic) ≈ ref_accumulators rtol = 1e-4
end

end
