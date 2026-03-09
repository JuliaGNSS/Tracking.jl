module DownconvertAndCorrelateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals:
    GPSL1, GalileoE1B, gen_code, get_code_frequency, get_code_center_frequency_ratio, get_code_type
using Bumper: SlabBuffer
import Tracking
using Tracking:
    AbstractCorrelator,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    EarlyPromptLateCorrelator,
    KADownconvertAndCorrelator,
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
    @test get_accumulators(fused_static) ≈ get_accumulators(fused_dynamic) rtol = 1e-4
end

@testset "Fused downconvert and correlate multi-antenna" begin
    using StaticArrays: SVector

    gpsl1 = GPSL1()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    carrier_doppler = 1000.0Hz
    prn = 1
    num_ants = 2

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants))
    sample_shifts = get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)

    # Multi-antenna signal: each antenna gets slightly different data
    signal_ant1 =
        gen_code(
            num_samples_signal, gpsl1, prn, sampling_frequency,
            code_frequency, code_phase,
        ) .* cis.(2π * (0:(num_samples_signal - 1)) * carrier_doppler / sampling_frequency)
    signal_ant2 = signal_ant1 .* cis(0.3)  # phase-shifted copy
    signal = hcat(signal_ant1, signal_ant2)

    code_replica_length = num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica, gpsl1, code_frequency, sampling_frequency,
        code_phase, 1, num_samples_signal, sample_shifts, prn,
        Val(sampling_frequency),
    )

    # Scalar reference for multi-antenna
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = [zeros(ComplexF64, num_ants) for _ in 1:NC]
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i in 1:num_samples_signal
        carrier = cis(-carrier_step * (i - 1))
        for j in 1:num_ants
            dc_sample = signal[i, j] * carrier
            for tap in 1:NC
                code_idx = i + sample_shifts[tap] - min_shift
                ref_accumulators[tap][j] += dc_sample * code_replica[code_idx]
            end
        end
    end

    # Static (@generated) path
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator, signal, code_replica, sample_shifts,
        carrier_frequency, sampling_frequency, 0.0, 1, num_samples_signal,
    )
    for tap in 1:NC
        @test collect(get_accumulators(fused_static)[tap]) ≈ ref_accumulators[tap] rtol = 1e-4
    end

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator, signal, code_replica, dynamic_shifts,
        carrier_frequency, sampling_frequency, 0.0, 1, num_samples_signal,
    )
    for tap in 1:NC
        @test collect(get_accumulators(fused_dynamic)[tap]) ≈ ref_accumulators[tap] rtol = 1e-4
    end

    # Static and dynamic should match
    for k in 1:NC
        @test collect(get_accumulators(fused_static)[k]) ≈ collect(get_accumulators(fused_dynamic)[k]) rtol = 1e-4
    end

    # Dynamic path with non-SIMD-aligned signal length (covers scalar remainder loop)
    num_samples_odd = 5003  # not a multiple of SIMD width
    signal_odd_ant1 =
        gen_code(
            num_samples_odd, gpsl1, prn, sampling_frequency,
            code_frequency, code_phase,
        ) .* cis.(2π * (0:(num_samples_odd - 1)) * carrier_doppler / sampling_frequency)
    signal_odd = hcat(signal_odd_ant1, signal_odd_ant1 .* cis(0.3))
    code_replica_odd = Vector{get_code_type(gpsl1)}(undef, num_samples_odd + maximum(sample_shifts) - minimum(sample_shifts))
    Tracking.gen_code_replica!(
        code_replica_odd, gpsl1, code_frequency, sampling_frequency,
        code_phase, 1, num_samples_odd, sample_shifts, prn,
        Val(sampling_frequency),
    )
    ref_odd = [zeros(ComplexF64, num_ants) for _ in 1:NC]
    for i in 1:num_samples_odd
        carrier = cis(-carrier_step * (i - 1))
        for j in 1:num_ants
            dc_sample = signal_odd[i, j] * carrier
            for tap in 1:NC
                ref_odd[tap][j] += dc_sample * code_replica_odd[i + sample_shifts[tap] - min_shift]
            end
        end
    end
    fused_odd = Tracking.downconvert_and_correlate_fused!(
        correlator, signal_odd, code_replica_odd, dynamic_shifts,
        carrier_frequency, sampling_frequency, 0.0, 1, num_samples_odd,
    )
    for tap in 1:NC
        @test collect(get_accumulators(fused_odd)[tap]) ≈ ref_odd[tap] rtol = 1e-4
    end
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

@testset "Downconvert and correlate with CPUThreaded" begin
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
        CPUThreadedDownconvertAndCorrelator((gpsl1,), Val(sampling_frequency))

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

    next_track_state = downconvert_and_correlate(
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

    next_track_state = downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    @test real.(get_correlator(next_track_state, 2).accumulators) ≈ [2919, 4947, 2915] rtol=1e-3

    # Test early return when signal_start_sample is past the signal end
    sat_past_end = SatState(
        SatState(gpsl1, 1, code_phase, 1000.0Hz);
        signal_start_sample = num_samples_signal + 1,
    )
    sss_skip = SystemSatsState(gpsl1, [sat_past_end])
    ts_skip = TrackState((sss_skip,))
    result_skip = downconvert_and_correlate(
        downconvert_and_correlator, signal, ts_skip, 1, sampling_frequency, intermediate_frequency,
    )
    @test get_correlator(result_skip, 1).accumulators == sat_past_end.correlator.accumulators
end

@testset "Downconvert and correlate with KA (CPU backend)" begin
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

    downconvert_and_correlator = KADownconvertAndCorrelator((gpsl1,), Array)

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

    next_track_state = downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    @test real.(get_correlator(next_track_state, 1).accumulators) ≈ [2921, 4949, 2917] atol = 25

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            2,
            sampling_frequency,
            get_code_frequency(gpsl1) + 500Hz * get_code_center_frequency_ratio(gpsl1),
            11.0,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 500.0Hz / sampling_frequency)

    next_track_state = downconvert_and_correlate(
        downconvert_and_correlator,
        signal,
        track_state,
        preferred_num_code_blocks_to_integrate,
        sampling_frequency,
        intermediate_frequency,
    )

    @test real.(get_correlator(next_track_state, 2).accumulators) ≈ [2919, 4947, 2915] atol = 25
end

@testset "Downconvert and correlate GalileoE1B: KA vs CPU" begin
    gal = GalileoE1B()
    sampling_frequency = 25e6Hz
    code_phase = 10.5
    num_samples_signal = 100000
    intermediate_frequency = 0.0Hz

    system_sats_state = SystemSatsState(
        gal,
        [SatState(gal, 1, code_phase, 100.0Hz), SatState(gal, 2, 11.0, 200.0Hz)];
    )
    track_state = TrackState((system_sats_state,))

    preferred_num_code_blocks_to_integrate = 1

    signal =
        gen_code(
            num_samples_signal,
            gal,
            1,
            sampling_frequency,
            get_code_frequency(gal) + 100Hz * get_code_center_frequency_ratio(gal),
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 100.0Hz / sampling_frequency)

    cpu_dc = CPUDownconvertAndCorrelator(Val(sampling_frequency))
    cpu_result = downconvert_and_correlate(
        cpu_dc, signal, track_state, preferred_num_code_blocks_to_integrate,
        sampling_frequency, intermediate_frequency,
    )

    ka_dc = KADownconvertAndCorrelator((gal,), Array)
    ka_result = downconvert_and_correlate(
        ka_dc, signal, track_state, preferred_num_code_blocks_to_integrate,
        sampling_frequency, intermediate_frequency,
    )

    cpu_accum = real.(get_correlator(cpu_result, 1).accumulators)
    ka_accum = real.(get_correlator(ka_result, 1).accumulators)
    @test ka_accum ≈ cpu_accum rtol = 0.01

    # Prompt correlator (index 3 for VeryEarlyPromptLate) should be peak
    prompt_idx = length(cpu_accum) ÷ 2 + 1
    @test abs(cpu_accum[prompt_idx]) == maximum(abs.(cpu_accum))
end

@testset "Downconvert and correlate multi-system GPSL1+GalileoE1B (KA CPU)" begin
    gpsl1 = GPSL1()
    gal = GalileoE1B()
    sampling_frequency = 25e6Hz
    num_samples_signal = 25000
    intermediate_frequency = 0.0Hz

    sss_l1 = SystemSatsState(
        gpsl1,
        [SatState(gpsl1, 1, 10.5, 1000.0Hz)];
    )
    sss_gal = SystemSatsState(
        gal,
        [SatState(gal, 1, 10.5, 100.0Hz)];
    )
    track_state = TrackState((sss_l1, sss_gal))

    # Generate GPSL1 signal for PRN 1
    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            1,
            sampling_frequency,
            get_code_frequency(gpsl1) + 1000Hz * get_code_center_frequency_ratio(gpsl1),
            10.5,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 1000.0Hz / sampling_frequency)

    ka_dc = KADownconvertAndCorrelator((gpsl1, gal), Array)
    result = downconvert_and_correlate(
        ka_dc, signal, track_state, 1, sampling_frequency, intermediate_frequency,
    )

    # GPSL1 sat should correlate well (prompt peak at index 2 for EarlyPromptLate)
    l1_accum = real.(get_correlator(result, 1, 1).accumulators)
    l1_prompt_idx = length(l1_accum) ÷ 2 + 1
    @test l1_accum[l1_prompt_idx] > l1_accum[1]
    @test l1_accum[l1_prompt_idx] > l1_accum[end]
    @test l1_accum[l1_prompt_idx] > 4000  # strong correlation

    # GalileoE1B sat should NOT correlate with GPSL1 signal (noise-level)
    gal_accum = real.(get_correlator(result, 2, 1).accumulators)
    gal_prompt_idx = length(gal_accum) ÷ 2 + 1
    @test abs(gal_accum[gal_prompt_idx]) < abs(l1_accum[l1_prompt_idx]) / 10
end

const CUDA_AVAILABLE = try
    @eval using CUDA: CUDA, CuArray, cu
    CUDA.functional()
catch
    false
end

if CUDA_AVAILABLE
    @testset "Downconvert and correlate with KA (CUDA backend)" begin
        gpsl1 = GPSL1()
        sampling_frequency = 5e6Hz
        code_phase = 10.5
        num_samples_signal = 5000
        intermediate_frequency = 0.0Hz

        system_sats_state = SystemSatsState(
            gpsl1,
            [
                SatState(gpsl1, 1, code_phase, 1000.0Hz),
                SatState(gpsl1, 2, 11.0, 500.0Hz),
            ];
        )
        multiple_system_sats_state = (system_sats_state,)

        downconvert_and_correlator = KADownconvertAndCorrelator((gpsl1,), CuArray)

        track_state = TrackState(multiple_system_sats_state)

        preferred_num_code_blocks_to_integrate = 1

        signal_cpu =
            gen_code(
                num_samples_signal,
                gpsl1,
                1,
                sampling_frequency,
                get_code_frequency(gpsl1) +
                1000Hz * get_code_center_frequency_ratio(gpsl1),
                code_phase,
            ) .* cis.(2π * (0:(num_samples_signal-1)) * 1000.0Hz / sampling_frequency)
        signal = CuArray(ComplexF32.(signal_cpu))

        next_track_state = downconvert_and_correlate(
            downconvert_and_correlator,
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )

        @test real.(get_correlator(next_track_state, 1).accumulators) ≈
              [2921, 4949, 2917] atol = 25

        signal_cpu =
            gen_code(
                num_samples_signal,
                gpsl1,
                2,
                sampling_frequency,
                get_code_frequency(gpsl1) +
                500Hz * get_code_center_frequency_ratio(gpsl1),
                11.0,
            ) .* cis.(2π * (0:(num_samples_signal-1)) * 500.0Hz / sampling_frequency)
        signal = CuArray(ComplexF32.(signal_cpu))

        next_track_state = downconvert_and_correlate(
            downconvert_and_correlator,
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )

        @test real.(get_correlator(next_track_state, 2).accumulators) ≈
              [2919, 4947, 2915] atol = 25
    end
end

const AMDGPU_AVAILABLE = try
    @eval using AMDGPU
    AMDGPU.functional()
catch
    false
end

if AMDGPU_AVAILABLE
    @testset "Downconvert and correlate with KA (AMDGPU backend)" begin
        gpsl1 = GPSL1()
        sampling_frequency = 5e6Hz
        code_phase = 10.5
        num_samples_signal = 5000
        intermediate_frequency = 0.0Hz

        system_sats_state = SystemSatsState(
            gpsl1,
            [
                SatState(gpsl1, 1, code_phase, 1000.0Hz),
                SatState(gpsl1, 2, 11.0, 500.0Hz),
            ];
        )
        multiple_system_sats_state = (system_sats_state,)

        downconvert_and_correlator = KADownconvertAndCorrelator((gpsl1,), ROCArray)

        track_state = TrackState(multiple_system_sats_state)

        preferred_num_code_blocks_to_integrate = 1

        signal_cpu =
            gen_code(
                num_samples_signal,
                gpsl1,
                1,
                sampling_frequency,
                get_code_frequency(gpsl1) +
                1000Hz * get_code_center_frequency_ratio(gpsl1),
                code_phase,
            ) .* cis.(2π * (0:(num_samples_signal-1)) * 1000.0Hz / sampling_frequency)
        signal = ROCArray(ComplexF32.(signal_cpu))

        next_track_state = downconvert_and_correlate(
            downconvert_and_correlator,
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )

        @test real.(get_correlator(next_track_state, 1).accumulators) ≈
              [2921, 4949, 2917] atol = 25

        signal_cpu =
            gen_code(
                num_samples_signal,
                gpsl1,
                2,
                sampling_frequency,
                get_code_frequency(gpsl1) +
                500Hz * get_code_center_frequency_ratio(gpsl1),
                11.0,
            ) .* cis.(2π * (0:(num_samples_signal-1)) * 500.0Hz / sampling_frequency)
        signal = ROCArray(ComplexF32.(signal_cpu))

        next_track_state = downconvert_and_correlate(
            downconvert_and_correlator,
            signal,
            track_state,
            preferred_num_code_blocks_to_integrate,
            sampling_frequency,
            intermediate_frequency,
        )

        @test real.(get_correlator(next_track_state, 2).accumulators) ≈
              [2919, 4947, 2915] atol = 25
    end
end

end
