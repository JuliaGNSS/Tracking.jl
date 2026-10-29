module DownconvertAndCorrelateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA, gen_code, get_chip_rate, get_chip_rate_to_center_frequency_ratio
import Tracking
using Tracking:
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    TrackedSat,
    TrackState,
    BandMeasurement,
    downconvert_and_correlate,
    get_correlator_bank,
    get_correlator_bank_outputs,
    get_sat_state
import TrackingLoops
using TrackingLoops:
    AbstractCorrelatorBank,
    EarlyPromptLateCorrelatorBank,
    NumAntennas,
    get_correlator_outputs,
    calc_correlator_sample_offsets

@testset "Downconvert and Correlator" begin
    # Scratch grows lazily, so default construction needs no setup.
    @test CPUDownconvertAndCorrelator() isa CPUDownconvertAndCorrelator
    @test CPUThreadedDownconvertAndCorrelator() isa CPUThreadedDownconvertAndCorrelator
end

@testset "Downconvert and correlate with $DC" for DC in [
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
]
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples = 5000
    intermediate_frequency = 0.0Hz

    downconvert_and_correlator = DC()

    # One code period (~4949 samples) completes in the 5000-sample buffer, so the
    # result is read from `correlator_bank_outputs` (see `_advance_sat` in
    # src/downconvert_and_correlate.jl). Each check uses its own TrackState so the
    # reused `correlator_bank_outputs` buffer is not carried between calls.
    only_output(track_state, prn) =
        only(get_correlator_bank_outputs(only(get_sat_state(track_state, prn).signals)))

    samples =
        gen_code(
            num_samples,
            gpsl1,
            1,
            sampling_frequency,
            get_chip_rate(gpsl1) + 1000Hz * get_chip_rate_to_center_frequency_ratio(gpsl1),
            code_phase,
        ) .* cis.(2π * (0:(num_samples-1)) * 1000.0Hz / sampling_frequency)

    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 1, code_phase, 1000.0Hz)])
    measurements =
        (L1 = BandMeasurement(samples, sampling_frequency, intermediate_frequency),)
    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        measurements,
        track_state,
    )

    @test real.(only_output(next_track_state, 1).correlator_bank.correlator_outputs) ≈
          [2921, 4949, 2917] rtol = 1e-3

    samples =
        gen_code(
            num_samples,
            gpsl1,
            2,
            sampling_frequency,
            get_chip_rate(gpsl1) + 500Hz * get_chip_rate_to_center_frequency_ratio(gpsl1),
            11.0,
        ) .* cis.(2π * (0:(num_samples-1)) * 500.0Hz / sampling_frequency)

    track_state = TrackState(gpsl1, [TrackedSat(gpsl1, 2, 11.0, 500.0Hz)])
    measurements =
        (L1 = BandMeasurement(samples, sampling_frequency, intermediate_frequency),)
    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        measurements,
        track_state,
    )

    @test real.(only_output(next_track_state, 2).correlator_bank.correlator_outputs) ≈
          [2919, 4947, 2915] rtol = 1e-3
end

@testset "Early return when signal_start_sample past end with $DC" for DC in [
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
]
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples = 5000
    intermediate_frequency = 0.0Hz

    samples = ones(ComplexF64, num_samples)
    sat_past_end = TrackedSat(
        TrackedSat(gpsl1, 1, code_phase, 1000.0Hz);
        signal_start_sample = num_samples + 1,
    )
    ts_skip = TrackState(gpsl1, [sat_past_end])
    downconvert_and_correlator = DC()

    measurements =
        (L1 = BandMeasurement(samples, sampling_frequency, intermediate_frequency),)

    # Immutable form
    result_skip =
        downconvert_and_correlate(downconvert_and_correlator, measurements, ts_skip)
    @test get_correlator_bank(result_skip, 1).correlator_outputs ==
          get_correlator_bank(sat_past_end).correlator_outputs

    # In-place form takes the same early return; the sat stays unchanged.
    Tracking.downconvert_and_correlate!(downconvert_and_correlator, measurements, ts_skip)
    @test get_correlator_bank(ts_skip, 1).correlator_outputs ==
          get_correlator_bank(sat_past_end).correlator_outputs
end

@testset "Fused downconvert and correlate" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)

    correlator_bank = EarlyPromptLateCorrelatorBank()
    sample_shifts =
        calc_correlator_sample_offsets(correlator_bank, sampling_frequency, code_frequency)

    samples =
        gen_code(num_samples, gpsl1, prn, sampling_frequency, code_frequency, code_phase) .*
        cis.(2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency)

    code_replica_length = num_samples + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = Vector{Int8}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples,
        sample_shifts,
        prn,
    )

    # Scalar reference: plain Julia loop, no SIMD, no tricks
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = samples[i] * carrier
        for tap = 1:NC
            code_idx = i + sample_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Static (@generated) path
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        samples,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples,
    )

    @test get_correlator_outputs(fused_static) ≈ ref_accumulators rtol = 1e-4

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    tile_re = Vector{Float32}(undef, num_samples)
    tile_im = Vector{Float32}(undef, num_samples)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        samples,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples,
        tile_re,
        tile_im,
    )

    @test get_correlator_outputs(fused_dynamic) ≈ ref_accumulators rtol = 1e-4

    # Static and dynamic paths should produce nearly identical results
    @test get_correlator_outputs(fused_static) ≈ get_correlator_outputs(fused_dynamic) rtol =
        1e-4
end

@testset "Fused downconvert and correlate multi-antenna" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples = 5000
    carrier_doppler = 1000.0Hz
    prn = 1
    num_antennas = 2

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)

    correlator_bank =
        EarlyPromptLateCorrelatorBank(; num_antennas = NumAntennas(num_antennas))
    sample_shifts =
        calc_correlator_sample_offsets(correlator_bank, sampling_frequency, code_frequency)

    # Multi-antenna signal: each antenna gets slightly different data
    samples_ant1 =
        gen_code(num_samples, gpsl1, prn, sampling_frequency, code_frequency, code_phase) .*
        cis.(2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency)
    samples_ant2 = samples_ant1 .* cis(0.3)  # phase-shifted copy
    samples = hcat(samples_ant1, samples_ant2)

    code_replica_length = num_samples + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = Vector{Int8}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples,
        sample_shifts,
        prn,
    )

    # Scalar reference for multi-antenna
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = [zeros(ComplexF64, num_antennas) for _ = 1:NC]
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples
        carrier = cis(-carrier_step * (i - 1))
        for j = 1:num_antennas
            dc_sample = samples[i, j] * carrier
            for tap = 1:NC
                code_idx = i + sample_shifts[tap] - min_shift
                ref_accumulators[tap][j] += dc_sample * code_replica[code_idx]
            end
        end
    end

    # Static (@generated) path
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        samples,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples,
    )
    for tap = 1:NC
        @test collect(get_correlator_outputs(fused_static)[tap]) ≈ ref_accumulators[tap] rtol =
            1e-4
    end

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    tile_re = Vector{Float32}(undef, num_samples * num_antennas)
    tile_im = Vector{Float32}(undef, num_samples * num_antennas)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        samples,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples,
        tile_re,
        tile_im,
    )
    for tap = 1:NC
        @test collect(get_correlator_outputs(fused_dynamic)[tap]) ≈ ref_accumulators[tap] rtol =
            1e-4
    end

    # Static and dynamic should match
    for k = 1:NC
        @test collect(get_correlator_outputs(fused_static)[k]) ≈
              collect(get_correlator_outputs(fused_dynamic)[k]) rtol = 1e-4
    end

    # Dynamic path with non-SIMD-aligned sample count (covers scalar remainder loop)
    num_samples_odd = 5003  # not a multiple of SIMD width
    signal_odd_ant1 =
        gen_code(
            num_samples_odd,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_odd-1)) * carrier_doppler / sampling_frequency)
    signal_odd = hcat(signal_odd_ant1, signal_odd_ant1 .* cis(0.3))
    code_replica_odd = Vector{Int8}(
        undef,
        num_samples_odd + maximum(sample_shifts) - minimum(sample_shifts),
    )
    Tracking.gen_code_replica!(
        code_replica_odd,
        gpsl1,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples_odd,
        sample_shifts,
        prn,
    )
    ref_odd = [zeros(ComplexF64, num_antennas) for _ = 1:NC]
    for i = 1:num_samples_odd
        carrier = cis(-carrier_step * (i - 1))
        for j = 1:num_antennas
            dc_sample = signal_odd[i, j] * carrier
            for tap = 1:NC
                ref_odd[tap][j] +=
                    dc_sample * code_replica_odd[i+sample_shifts[tap]-min_shift]
            end
        end
    end
    tile_re_odd = Vector{Float32}(undef, num_samples_odd * num_antennas)
    tile_im_odd = Vector{Float32}(undef, num_samples_odd * num_antennas)
    fused_odd = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        signal_odd,
        code_replica_odd,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples_odd,
        tile_re_odd,
        tile_im_odd,
    )
    for tap = 1:NC
        @test collect(get_correlator_outputs(fused_odd)[tap]) ≈ ref_odd[tap] rtol = 1e-4
    end
end

@testset "Fused downconvert and correlate with Vector accumulators" begin
    # Minimal correlator with Vector (not SVector) accumulators
    struct DynamicCorrelator <: AbstractCorrelatorBank{1}
        correlator_outputs::Vector{ComplexF64}
        shifts::Vector{Int}
    end
    TrackingLoops.get_correlator_outputs(c::DynamicCorrelator) = c.correlator_outputs
    TrackingLoops.with_correlator_outputs(c::DynamicCorrelator, acc) =
        DynamicCorrelator(collect(acc), c.shifts)

    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)

    # Use the same shifts as EarlyPromptLateCorrelatorBank but as a Vector
    epl = EarlyPromptLateCorrelatorBank()
    static_shifts = calc_correlator_sample_offsets(epl, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)

    dynamic_correlator = DynamicCorrelator(zeros(ComplexF64, 3), dynamic_shifts)

    samples =
        gen_code(num_samples, gpsl1, prn, sampling_frequency, code_frequency, code_phase) .*
        cis.(2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency)

    code_replica_length = num_samples + maximum(dynamic_shifts) - minimum(dynamic_shifts)
    code_replica = Vector{Int8}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples,
        static_shifts,
        prn,
    )

    # Scalar reference
    NC = length(dynamic_shifts)
    min_shift = minimum(dynamic_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = samples[i] * carrier
        for tap = 1:NC
            code_idx = i + dynamic_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Test: use the dynamic Vector accumulator path
    tile_re = Vector{Float32}(undef, num_samples)
    tile_im = Vector{Float32}(undef, num_samples)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        dynamic_correlator,
        samples,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples,
        tile_re,
        tile_im,
    )

    @test get_correlator_outputs(fused_dynamic) ≈ ref_accumulators rtol = 1e-4
end

@testset "Fused dynamic-shifts kernel with start_sample > 1" begin
    # Issue #126 (a): the AbstractVector-shifts kernel must read the code replica at
    # the absolute window offset `start_sample`, like the static kernels.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_buffer_samples = 5000
    start_sample = 3001
    num_samples = 2000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)

    correlator_bank = EarlyPromptLateCorrelatorBank()
    sample_shifts =
        calc_correlator_sample_offsets(correlator_bank, sampling_frequency, code_frequency)

    samples =
        gen_code(
            num_buffer_samples,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_buffer_samples-1)) * carrier_doppler / sampling_frequency)

    # Code phase of the first sample of the window at start_sample
    window_code_phase =
        code_phase + (start_sample - 1) * Float64(code_frequency / sampling_frequency)
    code_replica_length =
        start_sample - 1 + num_samples + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = zeros(Int8, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        window_code_phase,
        start_sample,
        num_samples,
        sample_shifts,
        prn,
    )

    # Scalar reference over the window only
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = samples[start_sample-1+i] * carrier
        for tap = 1:NC
            code_idx = start_sample - 1 + i + sample_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    tile_re = Vector{Float32}(undef, num_samples)
    tile_im = Vector{Float32}(undef, num_samples)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        samples,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        start_sample,
        num_samples,
        tile_re,
        tile_im,
    )
    @test get_correlator_outputs(fused_dynamic) ≈ ref_accumulators rtol = 1e-4

    # Must match the static in-register kernel on the same window
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator_bank,
        samples,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        start_sample,
        num_samples,
    )
    @test get_correlator_outputs(fused_dynamic) ≈ get_correlator_outputs(fused_static) rtol =
        1e-4
end

# Custom correlator whose `calc_correlator_sample_offsets` returns a plain
# `Vector` — a supported, exported extension point (issue #126 (b)).
struct VectorShiftsCorrelator <: AbstractCorrelatorBank{1}
    correlator_outputs::Vector{ComplexF64}
    shifts::Vector{Int}
end
TrackingLoops.get_correlator_outputs(c::VectorShiftsCorrelator) = c.correlator_outputs
TrackingLoops.with_correlator_outputs(c::VectorShiftsCorrelator, acc) =
    VectorShiftsCorrelator(collect(acc), c.shifts)
TrackingLoops.calc_correlator_sample_offsets(
    c::VectorShiftsCorrelator,
    sampling_frequency,
    code_frequency,
) = c.shifts

@testset "Vector-shifts correlator through CPU backend paths" begin
    # Issue #126 (b): both CPU backends and the public single-satellite
    # `downconvert_and_correlate!` must handle `Vector` shifts on a window not
    # starting at sample 1.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_buffer_samples = 5000
    start_sample = 3001
    num_samples = 2000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)

    epl = EarlyPromptLateCorrelatorBank()
    static_shifts = calc_correlator_sample_offsets(epl, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)
    correlator_bank = VectorShiftsCorrelator(zeros(ComplexF64, 3), dynamic_shifts)

    samples =
        gen_code(
            num_buffer_samples,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_buffer_samples-1)) * carrier_doppler / sampling_frequency)

    window_code_phase =
        code_phase + (start_sample - 1) * Float64(code_frequency / sampling_frequency)
    code_replica_length =
        start_sample - 1 + num_samples + maximum(dynamic_shifts) - minimum(dynamic_shifts)

    # Scalar reference over the window
    ref_code_replica = zeros(Int8, code_replica_length)
    Tracking.gen_code_replica!(
        ref_code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        window_code_phase,
        start_sample,
        num_samples,
        dynamic_shifts,
        prn,
    )
    NC = length(dynamic_shifts)
    min_shift = minimum(dynamic_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = samples[start_sample-1+i] * carrier
        for tap = 1:NC
            code_idx = start_sample - 1 + i + dynamic_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * ref_code_replica[code_idx]
        end
    end

    # Public single-satellite entry point (12-arg downconvert_and_correlate!)
    code_replica = zeros(Int8, code_replica_length)
    public_result = Tracking.downconvert_and_correlate!(
        gpsl1,
        samples,
        correlator_bank,
        code_replica,
        window_code_phase,
        0.0,
        code_frequency,
        carrier_frequency,
        sampling_frequency,
        start_sample,
        num_samples,
        prn,
    )
    @test get_correlator_outputs(public_result) ≈ ref_accumulators rtol = 1e-4

    # Both backends' despread primitive `_despread_one_signal!`.
    for dc in (CPUDownconvertAndCorrelator(), CPUThreadedDownconvertAndCorrelator())
        backend_result = Tracking._despread_one_signal!(
            dc,
            correlator_bank,
            samples,
            gpsl1,
            prn,
            dynamic_shifts,
            window_code_phase,
            0.0,
            code_frequency,
            carrier_frequency,
            sampling_frequency,
            start_sample,
            num_samples,
            code_replica_length,
        )
        @test get_correlator_outputs(backend_result) ≈ ref_accumulators rtol = 1e-4
    end

    # Same entry point with SVector shifts: covers the static (in-register) branch.
    epl_static = EarlyPromptLateCorrelatorBank()
    static_code_replica = zeros(Int8, code_replica_length)
    static_result = Tracking.downconvert_and_correlate!(
        gpsl1,
        samples,
        epl_static,
        static_code_replica,
        window_code_phase,
        0.0,
        code_frequency,
        carrier_frequency,
        sampling_frequency,
        start_sample,
        num_samples,
        prn,
    )
    @test get_correlator_outputs(static_result) ≈ ref_accumulators rtol = 1e-4
end

@testset "Fused downconvert and correlate tuple kernel" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_chip_rate_to_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_chip_rate(gpsl1)

    @testset "Single-antenna, N=$N signals" for N = 2:3
        correlator_banks = ntuple(_ -> EarlyPromptLateCorrelatorBank(), N)
        sample_shifts_each = calc_correlator_sample_offsets(
            correlator_banks[1],
            sampling_frequency,
            code_frequency,
        )
        all_sample_shifts = ntuple(_ -> sample_shifts_each, N)

        samples =
            gen_code(
                num_samples,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase,
            ) .* cis.(2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency)

        # Per-signal code replicas (same code in this scenario; just reuse).
        code_replica_length =
            num_samples + maximum(sample_shifts_each) - minimum(sample_shifts_each)
        code_replicas = ntuple(N) do _
            cr = Vector{Int8}(undef, code_replica_length)
            Tracking.gen_code_replica!(
                cr,
                gpsl1,
                code_frequency,
                sampling_frequency,
                code_phase,
                1,
                num_samples,
                sample_shifts_each,
                prn,
            )
            cr
        end

        # Scalar reference (same for every signal, since they share code+replica).
        NC = length(sample_shifts_each)
        min_shift = minimum(sample_shifts_each)
        ref_accumulators = zeros(ComplexF64, NC)
        carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
        for i = 1:num_samples
            carrier = cis(-carrier_step * (i - 1))
            dc_sample = samples[i] * carrier
            for tap = 1:NC
                code_idx = i + sample_shifts_each[tap] - min_shift
                ref_accumulators[tap] += dc_sample * code_replicas[1][code_idx]
            end
        end

        tile_re = Vector{Float32}(undef, num_samples)
        tile_im = Vector{Float32}(undef, num_samples)
        new_correlator_banks = Tracking.downconvert_and_correlate_fused_tuple!(
            correlator_banks,
            samples,
            code_replicas,
            all_sample_shifts,
            carrier_frequency,
            sampling_frequency,
            0.0,
            1,
            num_samples,
            tile_re,
            tile_im,
        )

        @test length(new_correlator_banks) == N
        for i = 1:N
            @test get_correlator_outputs(new_correlator_banks[i]) ≈ ref_accumulators rtol =
                1e-4
        end
    end

    @testset "Multi-antenna ($num_antennas ants), N=$N signals" for num_antennas in (2, 4),
        N = 2:3

        correlator_banks = ntuple(
            _ -> EarlyPromptLateCorrelatorBank(; num_antennas = NumAntennas(num_antennas)),
            N,
        )
        sample_shifts_each = calc_correlator_sample_offsets(
            correlator_banks[1],
            sampling_frequency,
            code_frequency,
        )
        all_sample_shifts = ntuple(_ -> sample_shifts_each, N)

        # Multi-antenna signal: each antenna gets a phase-shifted copy of the
        # base signal so accumulator values per antenna differ in a known way.
        samples_ant1 =
            gen_code(
                num_samples,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase,
            ) .* cis.(2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency)
        ant_phase_offsets = collect(0.0:0.3:(0.3*(num_antennas-1)))
        samples = hcat([samples_ant1 .* cis(ϕ) for ϕ in ant_phase_offsets]...)

        code_replica_length =
            num_samples + maximum(sample_shifts_each) - minimum(sample_shifts_each)
        code_replicas = ntuple(N) do _
            cr = Vector{Int8}(undef, code_replica_length)
            Tracking.gen_code_replica!(
                cr,
                gpsl1,
                code_frequency,
                sampling_frequency,
                code_phase,
                1,
                num_samples,
                sample_shifts_each,
                prn,
            )
            cr
        end

        # Scalar reference per antenna.
        NC = length(sample_shifts_each)
        min_shift = minimum(sample_shifts_each)
        ref_accumulators = [zeros(ComplexF64, num_antennas) for _ = 1:NC]
        carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
        for i = 1:num_samples
            carrier = cis(-carrier_step * (i - 1))
            for j = 1:num_antennas
                dc_sample = samples[i, j] * carrier
                for tap = 1:NC
                    code_idx = i + sample_shifts_each[tap] - min_shift
                    ref_accumulators[tap][j] += dc_sample * code_replicas[1][code_idx]
                end
            end
        end

        tile_re = Vector{Float32}(undef, num_samples * num_antennas)
        tile_im = Vector{Float32}(undef, num_samples * num_antennas)
        new_correlator_banks = Tracking.downconvert_and_correlate_fused_tuple!(
            correlator_banks,
            samples,
            code_replicas,
            all_sample_shifts,
            carrier_frequency,
            sampling_frequency,
            0.0,
            1,
            num_samples,
            tile_re,
            tile_im,
        )

        @test length(new_correlator_banks) == N
        for i = 1:N
            for tap = 1:NC
                @test collect(get_correlator_outputs(new_correlator_banks[i])[tap]) ≈
                      ref_accumulators[tap] rtol = 1e-4
            end
        end
    end
end

end
