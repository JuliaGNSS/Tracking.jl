module VectorMultiSignalTest

using Test: @test, @testset, @test_throws
using Unitful: Hz, MHz
using StaticArrays: SVector
using Logging: with_logger, NullLogger
using Dictionaries: dictionary
using GNSSSignals: GPSL1CA, GPSL5I, GPSL5Q, get_carrier_phase_offset
using Tracking:
    Tracking,
    TrackedSignal,
    TrackState,
    VectorPLLAndDLL,
    CorrelatorOutput,
    EarlyPromptLateCorrelator,
    add_satellite!,
    append_correlator_output!,
    enable_vt!,
    estimate_dopplers_and_filter_prompt!,
    get_correlator,
    get_doppler_estimator_state,
    get_sat_state,
    mean_carrier_discr,
    mean_code_discr,
    reset_code_discr_acc!,
    reset_carrier_discr_acc!,
    set_carrier_freq_updates!,
    set_code_freq_updates!

const FS = 5e6Hz
const SAMPLES = 5000

function vector_state(signals; vt_on = true)
    ts = TrackState(; signals = (g = signals,), doppler_estimator = VectorPLLAndDLL())
    add_satellite!(ts; prn = 1, group = :g, code_phase = 0.0, carrier_doppler = 0.0Hz)
    if vt_on
        enable_vt!(ts, (1,))
        set_carrier_freq_updates!(ts, dictionary((1 => 3.0Hz,)))
        set_code_freq_updates!(ts, dictionary((1 => -0.5Hz,)))
    end
    ts
end

# One record for signal `i` with carrier phase error `phase_error`, carrying the
# signal's own nominal carrier phase offset as the correlate phase would see it;
# `late` skews the late tap to give the DLL something to read.
function append_record!(ts, i, phase_error; sample_index = SAMPLES, late = 1.0)
    signal = Tracking.get_signal(get_sat_state(ts, :g, 1), i)
    rotation = SAMPLES * cis(phase_error + get_carrier_phase_offset(signal))
    correlator = get_correlator(ts, :g, 1, i)
    taps = SVector(0.5, 1.0, 0.5 * late) .* rotation
    output = CorrelatorOutput(
        EarlyPromptLateCorrelator(
            taps,
            correlator.preferred_early_late_to_prompt_code_shift,
        ),
        SAMPLES,
        sample_index,
    )
    append_correlator_output!(ts, output, :g, 1, i)
end

# The records carry no noise observations; the C/N₀ warning that follows is
# beside the point here.
estimate!(ts) = with_logger(NullLogger()) do
    estimate_dopplers_and_filter_prompt!(ts, (L5 = FS, L1 = FS))
end

# One chunk of two records per signal, so the second FLL reading has a previous
# prompt.
function feed!(ts, errors, lates; passenger_offset = 0)
    for (i, (phase_error, late)) in enumerate(zip(errors, lates))
        offset = i == 1 ? 0 : passenger_offset
        append_record!(ts, i, phase_error; sample_index = SAMPLES + offset, late)
        append_record!(ts, i, phase_error + 0.05; sample_index = 2SAMPLES + offset, late)
    end
    estimate!(ts)
    ts
end

@testset "every signal hands the navigation filter its own measurements" begin
    driver_alone = feed!(vector_state((GPSL5Q(),)), (0.1,), (1.2,))
    passenger_alone = feed!(vector_state((GPSL5I(),)), (0.3,), (0.8,))
    # The passenger's records ending on the driver's samples or not.
    for passenger_offset in (0, -1)
        ts = feed!(
            vector_state((GPSL5Q(), GPSL5I())),
            (0.1, 0.3),
            (1.2, 0.8);
            passenger_offset,
        )
        sat = get_sat_state(ts, :g, 1)
        state = get_doppler_estimator_state(sat)
        @test map(first, state.code_discr_acc) == (2, 2)
        @test map(first, state.carrier_discr_acc) == (2, 2)
        # Each slot holds what that signal reads tracked alone.
        @test mean_code_discr(ts, :g, 1, GPSL5Q) == mean_code_discr(driver_alone, :g, 1)
        @test mean_carrier_discr(sat, 1) == mean_carrier_discr(driver_alone, :g, 1)
        if passenger_offset == 0
            @test mean_code_discr(ts, :g, 1, GPSL5I) ==
                  mean_code_discr(passenger_alone, :g, 1)
            @test mean_carrier_discr(sat, 2) == mean_carrier_discr(passenger_alone, :g, 1)
        end
        # The navigation filter's reads stay per signal.
        @test_throws ArgumentError mean_code_discr(state)
        reset_code_discr_acc!(ts)
        reset_carrier_discr_acc!(ts)
        state = get_doppler_estimator_state(get_sat_state(ts, :g, 1))
        @test state.code_discr_acc == ((0, 0.0), (0, 0.0))
        @test isnothing(mean_carrier_discr(ts, :g, 1, 2))
    end
end

@testset "a sampling frequency in MHz accumulates the same measurements" begin
    function measurements(sampling_frequency)
        ts = vector_state((GPSL5Q(), GPSL5I()))
        for (i, (phase_error, late)) in enumerate(zip((0.1, 0.3), (1.2, 0.8)))
            append_record!(ts, i, phase_error; sample_index = SAMPLES, late)
            append_record!(ts, i, phase_error + 0.05; sample_index = 2SAMPLES, late)
        end
        with_logger(NullLogger()) do
            estimate_dopplers_and_filter_prompt!(ts, (L5 = sampling_frequency,))
        end
        sat = get_sat_state(ts, :g, 1)
        map(i -> (mean_code_discr(sat, i), mean_carrier_discr(sat, i)), (1, 2))
    end
    in_hz = measurements(FS)
    in_mhz = measurements(5.0MHz)
    for i = 1:2
        @test in_mhz[i][1] ≈ in_hz[i][1]
        @test in_mhz[i][2] ≈ in_hz[i][2]
    end
end

@testset "nothing is accumulated in the scalar fallback" begin
    ts = feed!(vector_state((GPSL5Q(), GPSL5I()); vt_on = false), (0.1, 0.3), (1.2, 0.8))
    state = get_doppler_estimator_state(get_sat_state(ts, :g, 1))
    @test state.code_discr_acc == ((0, 0.0), (0, 0.0))
    @test state.carrier_discr_acc == ((0, 0.0Hz), (0, 0.0Hz))
end

@testset "the driver's loops are unchanged by a passenger" begin
    alone = feed!(vector_state((GPSL5Q(),)), (0.1,), (1.2,))
    pair = feed!(vector_state((GPSL5Q(), GPSL5I())), (0.1, 0.3), (1.2, 0.8))
    @test Tracking.get_carrier_doppler(pair, :g, 1) ==
          Tracking.get_carrier_doppler(alone, :g, 1)
    @test Tracking.get_code_doppler(pair, :g, 1) == Tracking.get_code_doppler(alone, :g, 1)
end

# Fold `n` record pairs through the per-satellite update, as the estimate phase
# does. A function, not module scope, so `@allocated` measures the fold rather
# than global lookups.
function fold_records!(ts, outputs, n)
    sats = Tracking.get_sat_states(ts, :g)
    for _ = 1:n
        sat = sats[1]
        foreach(sat.signals, outputs) do signal, output
            push!(Tracking.get_correlator_outputs(signal), output)
            empty!(Tracking.get_filtered_prompts(signal))
        end
        sats[1] = Tracking._update_tracked_sat_doppler(
            sat,
            FS,
            ((nothing, false), (nothing, false)),
        )
    end
end

@static if VERSION >= v"1.11"
    @testset "per-signal accumulation is allocation-free" begin
        ts = vector_state((GPSL5Q(), GPSL5I()))
        sat = get_sat_state(ts, :g, 1)
        outputs = map(sat.signals) do signal
            rotation = SAMPLES * cis(0.1 + get_carrier_phase_offset(signal.signal))
            CorrelatorOutput(
                EarlyPromptLateCorrelator(SVector(0.5, 1.0, 0.6) .* rotation, 0.5),
                SAMPLES,
                SAMPLES,
            )
        end
        fold_records!(ts, outputs, 10)
        @test (@allocated fold_records!(ts, outputs, 100)) == 0
    end
end

end
