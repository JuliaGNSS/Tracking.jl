module DownconvertAndCorrelateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals:
    GPSL1CA, gen_code, get_code_frequency, get_code_center_frequency_ratio, get_code_type
import Tracking
using Tracking:
    AbstractCorrelator,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    EarlyPromptLateCorrelator,
    NumAnts,
    TrackedSat,
    TrackState,
    BandMeasurement,
    downconvert_and_correlate,
    get_accumulators,
    get_correlator,
    get_correlator_sample_shifts,
    update_accumulator

@testset "Downconvert and Correlator" begin
    # Both backends own long-lived `ScratchBuffers` (one per thread for
    # the threaded backend) that grow lazily on first use, so default
    # construction is cheap and valid without any further setup.
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
    num_samples_signal = 5000
    intermediate_frequency = 0.0Hz

    sats = [TrackedSat(gpsl1, 1, code_phase, 1000.0Hz), TrackedSat(gpsl1, 2, 11.0, 500.0Hz)]
    downconvert_and_correlator = DC()
    track_state = TrackState(gpsl1, sats)

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            1,
            sampling_frequency,
            get_code_frequency(gpsl1) + 1000Hz * get_code_center_frequency_ratio(gpsl1),
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * 1000.0Hz / sampling_frequency)

    measurements =
        (l1 = BandMeasurement(signal, sampling_frequency, intermediate_frequency),)
    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        measurements,
        track_state,
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

    measurements =
        (l1 = BandMeasurement(signal, sampling_frequency, intermediate_frequency),)
    next_track_state = @inferred downconvert_and_correlate(
        downconvert_and_correlator,
        measurements,
        track_state,
    )

    @test real.(get_correlator(next_track_state, 2).accumulators) ≈ [2919, 4947, 2915] rtol=1e-3
end

@testset "Early return when signal_start_sample past end with $DC" for DC in [
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
]
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    intermediate_frequency = 0.0Hz

    signal = ones(ComplexF64, num_samples_signal)
    sat_past_end = TrackedSat(
        TrackedSat(gpsl1, 1, code_phase, 1000.0Hz);
        signal_start_sample = num_samples_signal + 1,
    )
    ts_skip = TrackState(gpsl1, [sat_past_end])
    downconvert_and_correlator = DC()

    measurements =
        (l1 = BandMeasurement(signal, sampling_frequency, intermediate_frequency),)

    # Immutable form
    result_skip =
        downconvert_and_correlate(downconvert_and_correlator, measurements, ts_skip)
    @test get_correlator(result_skip, 1).accumulators ==
          get_correlator(sat_past_end).accumulators

    # In-place form takes the same `signal_samples_to_integrate == 0` early
    # return per sat. The TrackedSat is reassigned to itself unchanged.
    Tracking.downconvert_and_correlate!(downconvert_and_correlator, measurements, ts_skip)
    @test get_correlator(ts_skip, 1).accumulators ==
          get_correlator(sat_past_end).accumulators
end

@testset "Fused downconvert and correlate" begin
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    correlator = EarlyPromptLateCorrelator()
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)

    code_replica_length =
        num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
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
    )

    # Scalar reference: plain Julia loop, no SIMD, no tricks
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples_signal
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = signal[i] * carrier
        for tap = 1:NC
            code_idx = i + sample_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Static (@generated) path
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples_signal,
    )

    @test get_accumulators(fused_static) ≈ ref_accumulators rtol = 1e-4

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    tile_re = Vector{Float32}(undef, num_samples_signal)
    tile_im = Vector{Float32}(undef, num_samples_signal)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples_signal,
        tile_re,
        tile_im,
    )

    @test get_accumulators(fused_dynamic) ≈ ref_accumulators rtol = 1e-4

    # Static and dynamic paths should produce nearly identical results
    @test get_accumulators(fused_static) ≈ get_accumulators(fused_dynamic) rtol = 1e-4
end

@testset "Fused downconvert and correlate multi-antenna" begin
    using StaticArrays: SVector

    gpsl1 = GPSL1CA()
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
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)

    # Multi-antenna signal: each antenna gets slightly different data
    signal_ant1 =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)
    signal_ant2 = signal_ant1 .* cis(0.3)  # phase-shifted copy
    signal = hcat(signal_ant1, signal_ant2)

    code_replica_length =
        num_samples_signal + maximum(sample_shifts) - minimum(sample_shifts)
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
    )

    # Scalar reference for multi-antenna
    NC = length(sample_shifts)
    min_shift = minimum(sample_shifts)
    ref_accumulators = [zeros(ComplexF64, num_ants) for _ = 1:NC]
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples_signal
        carrier = cis(-carrier_step * (i - 1))
        for j = 1:num_ants
            dc_sample = signal[i, j] * carrier
            for tap = 1:NC
                code_idx = i + sample_shifts[tap] - min_shift
                ref_accumulators[tap][j] += dc_sample * code_replica[code_idx]
            end
        end
    end

    # Static (@generated) path
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples_signal,
    )
    for tap = 1:NC
        @test collect(get_accumulators(fused_static)[tap]) ≈ ref_accumulators[tap] rtol =
            1e-4
    end

    # Dynamic (AbstractVector) path
    dynamic_shifts = collect(sample_shifts)
    tile_re = Vector{Float32}(undef, num_samples_signal * num_ants)
    tile_im = Vector{Float32}(undef, num_samples_signal * num_ants)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples_signal,
        tile_re,
        tile_im,
    )
    for tap = 1:NC
        @test collect(get_accumulators(fused_dynamic)[tap]) ≈ ref_accumulators[tap] rtol =
            1e-4
    end

    # Static and dynamic should match
    for k = 1:NC
        @test collect(get_accumulators(fused_static)[k]) ≈
              collect(get_accumulators(fused_dynamic)[k]) rtol = 1e-4
    end

    # Dynamic path with non-SIMD-aligned signal length (covers scalar remainder loop)
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
    code_replica_odd = Vector{get_code_type(gpsl1)}(
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
    ref_odd = [zeros(ComplexF64, num_ants) for _ = 1:NC]
    for i = 1:num_samples_odd
        carrier = cis(-carrier_step * (i - 1))
        for j = 1:num_ants
            dc_sample = signal_odd[i, j] * carrier
            for tap = 1:NC
                ref_odd[tap][j] +=
                    dc_sample * code_replica_odd[i+sample_shifts[tap]-min_shift]
            end
        end
    end
    tile_re_odd = Vector{Float32}(undef, num_samples_odd * num_ants)
    tile_im_odd = Vector{Float32}(undef, num_samples_odd * num_ants)
    fused_odd = Tracking.downconvert_and_correlate_fused!(
        correlator,
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

    gpsl1 = GPSL1CA()
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
        ) .* cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)

    code_replica_length =
        num_samples_signal + maximum(dynamic_shifts) - minimum(dynamic_shifts)
    code_replica = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
    Tracking.gen_code_replica!(
        code_replica,
        gpsl1,
        code_frequency,
        sampling_frequency,
        code_phase,
        1,
        num_samples_signal,
        static_shifts,
        prn,
    )

    # Scalar reference
    NC = length(dynamic_shifts)
    min_shift = minimum(dynamic_shifts)
    ref_accumulators = zeros(ComplexF64, NC)
    carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
    for i = 1:num_samples_signal
        carrier = cis(-carrier_step * (i - 1))
        dc_sample = signal[i] * carrier
        for tap = 1:NC
            code_idx = i + dynamic_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * code_replica[code_idx]
        end
    end

    # Test: use the dynamic Vector accumulator path
    tile_re = Vector{Float32}(undef, num_samples_signal)
    tile_im = Vector{Float32}(undef, num_samples_signal)
    fused_dynamic = Tracking.downconvert_and_correlate_fused!(
        dynamic_correlator,
        signal,
        code_replica,
        dynamic_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        1,
        num_samples_signal,
        tile_re,
        tile_im,
    )

    @test get_accumulators(fused_dynamic) ≈ ref_accumulators rtol = 1e-4
end

@testset "Fused dynamic-shifts kernel with start_sample > 1" begin
    # Issue #126 (a): the AbstractVector-shifts fused kernel must read the
    # code replica at the absolute window offset, matching where
    # gen_code_replica! writes it (index start_sample), like the static
    # in-register kernel and the tuple tile-share kernel do.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    start_sample = 3001
    num_samples = 2000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    correlator = EarlyPromptLateCorrelator()
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)

    # Code phase of the first sample of the window at start_sample
    window_code_phase =
        code_phase + (start_sample - 1) * Float64(code_frequency / sampling_frequency)
    code_replica_length =
        start_sample - 1 + num_samples + maximum(sample_shifts) - minimum(sample_shifts)
    code_replica = zeros(get_code_type(gpsl1), code_replica_length)
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
        dc_sample = signal[start_sample-1+i] * carrier
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
        correlator,
        signal,
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
    @test get_accumulators(fused_dynamic) ≈ ref_accumulators rtol = 1e-4

    # Must match the static in-register kernel on the same window
    fused_static = Tracking.downconvert_and_correlate_fused!(
        correlator,
        signal,
        code_replica,
        sample_shifts,
        carrier_frequency,
        sampling_frequency,
        0.0,
        start_sample,
        num_samples,
    )
    @test get_accumulators(fused_dynamic) ≈ get_accumulators(fused_static) rtol = 1e-4
end

# Custom correlator whose `get_correlator_sample_shifts` returns a plain
# `Vector` — a supported, exported extension point (issue #126 (b)).
struct VectorShiftsCorrelator <: AbstractCorrelator{1}
    accumulators::Vector{ComplexF64}
    shifts::Vector{Int}
end
Tracking.get_accumulators(c::VectorShiftsCorrelator) = c.accumulators
Tracking.update_accumulator(c::VectorShiftsCorrelator, acc) =
    VectorShiftsCorrelator(collect(acc), c.shifts)
Tracking.get_correlator_sample_shifts(
    c::VectorShiftsCorrelator,
    sampling_frequency,
    code_frequency,
) = c.shifts

@testset "Vector-shifts correlator through CPU backend paths" begin
    # Issue #126 (b): both CPU backends and the public single-satellite
    # downconvert_and_correlate! must handle correlators whose sample
    # shifts are a runtime-sized Vector, and produce correct results for
    # windows that do not start at sample 1.
    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    start_sample = 3001
    num_samples = 2000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    epl = EarlyPromptLateCorrelator()
    static_shifts = get_correlator_sample_shifts(epl, sampling_frequency, code_frequency)
    dynamic_shifts = collect(static_shifts)
    correlator = VectorShiftsCorrelator(zeros(ComplexF64, 3), dynamic_shifts)

    signal =
        gen_code(
            num_samples_signal,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            code_phase,
        ) .* cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)

    window_code_phase =
        code_phase + (start_sample - 1) * Float64(code_frequency / sampling_frequency)
    code_replica_length =
        start_sample - 1 + num_samples + maximum(dynamic_shifts) - minimum(dynamic_shifts)

    # Scalar reference over the window
    ref_code_replica = zeros(get_code_type(gpsl1), code_replica_length)
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
        dc_sample = signal[start_sample-1+i] * carrier
        for tap = 1:NC
            code_idx = start_sample - 1 + i + dynamic_shifts[tap] - min_shift
            ref_accumulators[tap] += dc_sample * ref_code_replica[code_idx]
        end
    end

    # Public single-satellite entry point (12-arg downconvert_and_correlate!)
    code_replica = zeros(get_code_type(gpsl1), code_replica_length)
    public_result = Tracking.downconvert_and_correlate!(
        gpsl1,
        signal,
        correlator,
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
    @test get_accumulators(public_result) ≈ ref_accumulators rtol = 1e-4

    # Both backends' per-signal kernel
    for dc in (CPUDownconvertAndCorrelator(), CPUThreadedDownconvertAndCorrelator())
        backend_code_replica = zeros(get_code_type(gpsl1), code_replica_length)
        backend_result = Tracking._correlate_one_signal!(
            dc,
            backend_code_replica,
            gpsl1,
            correlator,
            signal,
            dynamic_shifts,
            window_code_phase,
            0.0,
            code_frequency,
            carrier_frequency,
            sampling_frequency,
            start_sample,
            num_samples,
            prn,
        )
        @test get_accumulators(backend_result) ≈ ref_accumulators rtol = 1e-4
    end

    # The same public entry point with a standard SVector-shifts correlator
    # must route through the in-register kernel and agree with the scalar
    # reference (covers the static branch of the standalone fused dispatch).
    epl_static = EarlyPromptLateCorrelator()
    static_code_replica = zeros(get_code_type(gpsl1), code_replica_length)
    static_result = Tracking.downconvert_and_correlate!(
        gpsl1,
        signal,
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
    @test get_accumulators(static_result) ≈ ref_accumulators rtol = 1e-4
end

@testset "Fused downconvert and correlate tuple kernel" begin
    using StaticArrays: SVector

    gpsl1 = GPSL1CA()
    sampling_frequency = 5e6Hz
    code_phase = 10.5
    num_samples_signal = 5000
    carrier_doppler = 1000.0Hz
    prn = 1

    carrier_frequency = carrier_doppler + 0.0Hz
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
    code_frequency = code_doppler + get_code_frequency(gpsl1)

    @testset "Single-antenna, N=$N signals" for N = 2:3
        correlators = ntuple(_ -> EarlyPromptLateCorrelator(), N)
        sample_shifts_each =
            get_correlator_sample_shifts(correlators[1], sampling_frequency, code_frequency)
        all_sample_shifts = ntuple(_ -> sample_shifts_each, N)

        signal =
            gen_code(
                num_samples_signal,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase,
            ) .*
            cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)

        # Per-signal code replicas (same code in this scenario; just reuse).
        code_replica_length =
            num_samples_signal + maximum(sample_shifts_each) - minimum(sample_shifts_each)
        code_replicas = ntuple(N) do _
            cr = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
            Tracking.gen_code_replica!(
                cr,
                gpsl1,
                code_frequency,
                sampling_frequency,
                code_phase,
                1,
                num_samples_signal,
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
        for i = 1:num_samples_signal
            carrier = cis(-carrier_step * (i - 1))
            dc_sample = signal[i] * carrier
            for tap = 1:NC
                code_idx = i + sample_shifts_each[tap] - min_shift
                ref_accumulators[tap] += dc_sample * code_replicas[1][code_idx]
            end
        end

        tile_re = Vector{Float32}(undef, num_samples_signal)
        tile_im = Vector{Float32}(undef, num_samples_signal)
        new_correlators = Tracking.downconvert_and_correlate_fused_tuple!(
            correlators,
            signal,
            code_replicas,
            all_sample_shifts,
            carrier_frequency,
            sampling_frequency,
            0.0,
            1,
            num_samples_signal,
            tile_re,
            tile_im,
        )

        @test length(new_correlators) == N
        for i = 1:N
            @test get_accumulators(new_correlators[i]) ≈ ref_accumulators rtol = 1e-4
        end
    end

    @testset "Multi-antenna ($num_ants ants), N=$N signals" for num_ants in (2, 4), N = 2:3
        correlators =
            ntuple(_ -> EarlyPromptLateCorrelator(; num_ants = NumAnts(num_ants)), N)
        sample_shifts_each =
            get_correlator_sample_shifts(correlators[1], sampling_frequency, code_frequency)
        all_sample_shifts = ntuple(_ -> sample_shifts_each, N)

        # Multi-antenna signal: each antenna gets a phase-shifted copy of the
        # base signal so accumulator values per antenna differ in a known way.
        signal_ant1 =
            gen_code(
                num_samples_signal,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase,
            ) .*
            cis.(2π * (0:(num_samples_signal-1)) * carrier_doppler / sampling_frequency)
        ant_phase_offsets = collect(0.0:0.3:(0.3*(num_ants-1)))
        signal = hcat([signal_ant1 .* cis(ϕ) for ϕ in ant_phase_offsets]...)

        code_replica_length =
            num_samples_signal + maximum(sample_shifts_each) - minimum(sample_shifts_each)
        code_replicas = ntuple(N) do _
            cr = Vector{get_code_type(gpsl1)}(undef, code_replica_length)
            Tracking.gen_code_replica!(
                cr,
                gpsl1,
                code_frequency,
                sampling_frequency,
                code_phase,
                1,
                num_samples_signal,
                sample_shifts_each,
                prn,
            )
            cr
        end

        # Scalar reference per antenna.
        NC = length(sample_shifts_each)
        min_shift = minimum(sample_shifts_each)
        ref_accumulators = [zeros(ComplexF64, num_ants) for _ = 1:NC]
        carrier_step = 2π * Float64(carrier_frequency / sampling_frequency)
        for i = 1:num_samples_signal
            carrier = cis(-carrier_step * (i - 1))
            for j = 1:num_ants
                dc_sample = signal[i, j] * carrier
                for tap = 1:NC
                    code_idx = i + sample_shifts_each[tap] - min_shift
                    ref_accumulators[tap][j] += dc_sample * code_replicas[1][code_idx]
                end
            end
        end

        tile_re = Vector{Float32}(undef, num_samples_signal * num_ants)
        tile_im = Vector{Float32}(undef, num_samples_signal * num_ants)
        new_correlators = Tracking.downconvert_and_correlate_fused_tuple!(
            correlators,
            signal,
            code_replicas,
            all_sample_shifts,
            carrier_frequency,
            sampling_frequency,
            0.0,
            1,
            num_samples_signal,
            tile_re,
            tile_im,
        )

        @test length(new_correlators) == N
        for i = 1:N
            for tap = 1:NC
                @test collect(get_accumulators(new_correlators[i])[tap]) ≈
                      ref_accumulators[tap] rtol = 1e-4
            end
        end
    end
end

end
