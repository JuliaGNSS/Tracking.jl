module DiscriminatorCombiningTest

using Test: @test, @testset, @test_throws
using Unitful: Hz, s, ns, m
using StaticArrays: SVector
using Logging: with_logger, NullLogger
using Dictionaries: dictionary
using GNSSSignals:
    GPSL1CA,
    GalileoE1B,
    GalileoE1C,
    GPSL5I,
    GPSL5Q,
    GPSL1C_D,
    GPSL1C_P,
    get_carrier_phase_offset
using Tracking:
    Tracking,
    TrackedSignal,
    TrackState,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    CorrelatorOutput,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    add_satellite!,
    append_correlator_output!,
    estimate_dopplers_and_filter_prompt!,
    get_carrier_doppler,
    get_code_doppler,
    get_correlator,
    get_correlator_outputs,
    get_filtered_prompts,
    get_group_delay,
    set_group_delay!,
    VectorPLLAndDLL,
    enable_vt!,
    set_code_freq_updates!,
    set_carrier_freq_updates!,
    mean_code_discr,
    mean_carrier_discr

@testset "group delay" begin
    @test isnothing(get_group_delay(TrackedSignal(GalileoE1B())))
    @test get_group_delay(TrackedSignal(GalileoE1B(); group_delay = 1.5ns)) ≈ 1.5e-9s

    ts = TrackState(; signals = (e1 = (GalileoE1C(), GalileoE1B()),))
    add_satellite!(ts; prn = 11, group = :e1, code_phase = 0.0, carrier_doppler = 0.0Hz)
    @test isnothing(get_group_delay(ts, :e1, 11, GalileoE1B))

    set_group_delay!(ts, :e1, 11, GalileoE1B, -0.5ns)
    @test get_group_delay(ts, :e1, 11, GalileoE1B) === -5.0e-10s
    @test isnothing(get_group_delay(ts, :e1, 11, GalileoE1C))

    set_group_delay!(ts, 11, 1, 0.0s)
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # `nothing` clears again: unknown, not zero.
    set_group_delay!(ts, 11, GalileoE1B, nothing)
    @test isnothing(get_group_delay(ts, :e1, 11, 2))
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # A group delay is a time.
    @test_throws Exception set_group_delay!(ts, 11, 2, 1.0m)
    @test_throws Exception set_group_delay!(ts, 11, 2, 1.0)
end

const FS = 5e6Hz
const SAMPLES = 5000

function combining_state(
    signals;
    combining = true,
    assisted = false,
    estimator = assisted ?
                ConventionalAssistedPLLAndDLL(; discriminator_combining = combining) :
                ConventionalPLLAndDLL(; discriminator_combining = combining),
)
    ts = TrackState(; signals = (g = signals,), doppler_estimator = estimator)
    add_satellite!(ts; prn = 1, group = :g, code_phase = 0.0, carrier_doppler = 0.0Hz)
end

with_taps(c::EarlyPromptLateCorrelator, taps) = EarlyPromptLateCorrelator(
    SVector(taps[2], taps[3], taps[4]),
    c.preferred_early_late_to_prompt_code_shift,
)
with_taps(c::VeryEarlyPromptLateCorrelator, taps) = VeryEarlyPromptLateCorrelator(
    SVector(taps...),
    c.preferred_early_late_to_prompt_code_shift,
    c.preferred_very_early_late_to_prompt_code_shift,
)

# One record for signal `i` with carrier phase error `phase_error` against the
# driver's frame: the signal's own nominal carrier phase offset is on its taps,
# exactly as the correlate phase would see it. `late` skews the late taps to
# give the DLL something to read.
function append_record!(
    ts,
    i,
    phase_error;
    sample_index = SAMPLES,
    late = 1.0,
    samples = SAMPLES,
)
    signal = Tracking.get_signal(Tracking.get_sat_state(ts, :g, 1), i)
    rotation = SAMPLES * cis(phase_error + get_carrier_phase_offset(signal))
    taps = (0.25, 0.5, 1.0, 0.5 * late, 0.25 * late) .* rotation
    output = CorrelatorOutput(
        with_taps(get_correlator(ts, :g, 1, i), taps),
        samples,
        sample_index,
    )
    append_correlator_output!(ts, output, :g, 1, i)
end

# The records carry no noise observations; the C/N₀ warning that follows is beside
# the point here.
estimate!(ts) = with_logger(NullLogger()) do
    estimate_dopplers_and_filter_prompt!(ts, (L5 = FS, L1 = FS))
end

@testset "discriminator combining" begin
    @testset "off, or single-signal: bit-identical to the driver alone" begin
        alone = combining_state((GPSL5Q(),); combining = false)
        single = combining_state((GPSL5Q(),))
        off = combining_state((GPSL5Q(), GPSL5I()); combining = false)
        for ts in (alone, single, off)
            append_record!(ts, 1, 0.1; late = 1.2)
        end
        append_record!(off, 2, -0.3; late = 0.8)
        foreach(estimate!, (alone, single, off))
        for ts in (single, off)
            @test get_carrier_doppler(ts, :g, 1) == get_carrier_doppler(alone, :g, 1)
            @test get_code_doppler(ts, :g, 1) == get_code_doppler(alone, :g, 1)
        end
        # The passenger's record was still applied to the passenger.
        @test length(get_filtered_prompts(off, :g, 1, 2)) == 1
    end

    @testset "carrier: power-weighted mean in the driver's phase frame" begin
        # L5Q (−π/2) drives, L5I (0) passes: equal power, so the PLL reads the
        # plain mean of the two phase errors once L5I is de-rotated.
        both = combining_state((GPSL5Q(), GPSL5I()))
        append_record!(both, 1, 0.1)
        append_record!(both, 2, 0.3)
        mean_alone = combining_state((GPSL5Q(),))
        append_record!(mean_alone, 1, 0.2)
        foreach(estimate!, (both, mean_alone))
        @test get_carrier_doppler(both, :g, 1) ≈ get_carrier_doppler(mean_alone, :g, 1)
        @test isempty(get_correlator_outputs(both, :g, 1, 2))
        @test length(get_filtered_prompts(both, :g, 1, 2)) == 1

        # L1C-P 0.75, L1C-D 0.25.
        l1c = combining_state((GPSL1C_P(), GPSL1C_D()))
        append_record!(l1c, 1, 0.1)
        append_record!(l1c, 2, 0.5)
        l1c_alone = combining_state((GPSL1C_P(),))
        append_record!(l1c_alone, 1, 0.75 * 0.1 + 0.25 * 0.5)
        foreach(estimate!, (l1c, l1c_alone))
        @test get_carrier_doppler(l1c, :g, 1) ≈ get_carrier_doppler(l1c_alone, :g, 1)
    end

    @testset "records that do not coincide are not combined" begin
        for (index, samples) in ((SAMPLES - 1, SAMPLES), (SAMPLES, SAMPLES - 1))
            ts = combining_state((GPSL5Q(), GPSL5I()))
            append_record!(ts, 1, 0.1)
            append_record!(ts, 2, 0.3; sample_index = index, samples)
            alone = combining_state((GPSL5Q(),))
            append_record!(alone, 1, 0.1)
            foreach(estimate!, (ts, alone))
            @test get_carrier_doppler(ts, :g, 1) == get_carrier_doppler(alone, :g, 1)
            @test length(get_filtered_prompts(ts, :g, 1, 2)) == 1
            @test isempty(get_correlator_outputs(ts, :g, 1, 2))
        end
    end

    @testset "code: only with group delays on both ends, referred to the driver" begin
        function code_doppler(delays; passenger_late = 1.2)
            ts = combining_state((GPSL5Q(), GPSL5I()))
            for (i, delay) in enumerate(delays)
                set_group_delay!(ts, :g, 1, i, delay)
            end
            append_record!(ts, 1, 0.0; late = 1.2)
            append_record!(ts, 2, 0.0; late = passenger_late)
            estimate!(ts)
            get_code_doppler(ts, :g, 1)
        end
        alone = combining_state((GPSL5Q(),))
        append_record!(alone, 1, 0.0; late = 1.2)
        estimate!(alone)
        driver_only = get_code_doppler(alone, :g, 1)

        # Unknown on either end: the passenger stays out of the code loop.
        @test code_doppler((nothing, nothing); passenger_late = 0.8) == driver_only
        @test code_doppler((0.0s, nothing); passenger_late = 0.8) == driver_only
        @test code_doppler((nothing, 0.0s); passenger_late = 0.8) == driver_only
        # Known and equal, identical taps: the mean is the driver's value.
        @test code_doppler((0.0s, 0.0s)) ≈ driver_only
        # Known and different, and a different passenger reading: it takes part.
        @test !(code_doppler((0.0s, 0.0s); passenger_late = 0.8) ≈ driver_only)
        # A passenger *less* delayed than the driver sits at a larger code phase,
        # so its reading is lowered before averaging: with identical taps the
        # combined code error, hence the code Doppler, drops.
        @test code_doppler((1.0ns, 0.0s)) < driver_only
        @test code_doppler((0.0s, 1.0ns)) > driver_only
    end

    @testset "assisted carrier loop combines the FLL too" begin
        # Two records per chunk so the second has a previous prompt.
        function carrier_doppler(signals, errors)
            ts = combining_state(signals; assisted = true)
            for (i, (first_error, second_error)) in enumerate(errors)
                append_record!(ts, i, first_error; sample_index = SAMPLES)
                append_record!(ts, i, second_error; sample_index = 2SAMPLES)
            end
            estimate!(ts)
            get_carrier_doppler(ts, :g, 1)
        end
        both = carrier_doppler((GPSL5Q(), GPSL5I()), ((0.0, 0.1), (0.0, 0.3)))
        alone = carrier_doppler((GPSL5Q(),), ((0.0, 0.2),))
        @test both ≈ alone
    end
end

# Fold `n` coincident record pairs through the per-satellite update, as the
# estimate phase does. A function, not module scope, so `@allocated` measures the
# fold rather than global lookups.
function fold_records!(ts, outputs, n)
    sats = Tracking.get_sat_states(ts, :g)
    for _ = 1:n
        sat = sats[1]
        foreach(sat.signals, outputs) do signal, output
            push!(Tracking.get_correlator_outputs(signal), output)
            empty!(get_filtered_prompts(signal))
        end
        sats[1] = Tracking._update_tracked_sat_doppler(
            sat,
            FS,
            ((nothing, false), (nothing, false)),
        )
    end
end

@static if VERSION >= v"1.11"
    @testset "combining is allocation-free (assisted = $assisted)" for assisted in
                                                                       (false, true)
        ts = combining_state((GPSL5Q(), GPSL5I()); assisted)
        set_group_delay!(ts, :g, 1, 1, 0.0s)
        set_group_delay!(ts, :g, 1, 2, 0.5ns)
        sat = Tracking.get_sat_state(ts, :g, 1)
        outputs = map(sat.signals) do signal
            rotation = SAMPLES * cis(0.1 + get_carrier_phase_offset(signal.signal))
            taps = (0.25, 0.5, 1.0, 0.6, 0.25) .* rotation
            CorrelatorOutput(with_taps(signal.correlator, taps), SAMPLES, SAMPLES)
        end
        fold_records!(ts, outputs, 10)
        @test (@allocated fold_records!(ts, outputs, 100)) == 0
    end
end

@testset "vector tracking" begin
    signals = (GPSL5Q(), GPSL5I())
    # Two chunks of two records each, so the FLL has a previous prompt and the
    # passenger is applied across chunk boundaries.
    function feed!(ts, errors; late = (1.2, 0.8), passenger_index = SAMPLES, chunks = 2)
        for chunk = 1:chunks
            for (i, phase_error) in enumerate(errors)
                index = i == 1 ? SAMPLES : passenger_index
                append_record!(ts, i, phase_error; sample_index = index, late = late[i])
                append_record!(
                    ts,
                    i,
                    phase_error + 0.05;
                    sample_index = index + SAMPLES,
                    late = late[i],
                )
            end
            estimate!(ts)
        end
        ts
    end
    vector_state(signals; combining = true) = combining_state(
        signals;
        estimator = VectorPLLAndDLL(; discriminator_combining = combining),
    )
    function with_delays!(ts)
        foreach(i -> set_group_delay!(ts, :g, 1, i, (i - 1) * 1.0ns), 1:2)
        ts
    end
    function vt_on!(ts)
        enable_vt!(ts, (1,))
        set_carrier_freq_updates!(ts, dictionary((1 => 3.0Hz,)))
        set_code_freq_updates!(ts, dictionary((1 => -0.5Hz,)))
        ts
    end

    @testset "vt_on = false combines exactly like the conventional estimator" begin
        conventional =
            feed!(with_delays!(combining_state(signals; assisted = true)), (0.1, 0.3))
        vector = feed!(with_delays!(vector_state(signals)), (0.1, 0.3))
        @test get_carrier_doppler(vector, :g, 1) == get_carrier_doppler(conventional, :g, 1)
        @test get_code_doppler(vector, :g, 1) == get_code_doppler(conventional, :g, 1)
        # Nothing is accumulated while the scalar fallback runs.
        @test isnothing(mean_code_discr(vector, :g, 1, 2))
    end

    @testset "vt_on = true combines only the carrier phase" begin
        both = feed!(vt_on!(with_delays!(vector_state(signals))), (0.1, 0.3))
        alone = feed!(vt_on!(vector_state((GPSL5Q(),))), (0.2,))
        # PLL: the mean phase error; FLL and code: the navigation filter's.
        @test get_carrier_doppler(both, :g, 1) ≈ get_carrier_doppler(alone, :g, 1)
        @test get_code_doppler(both, :g, 1) ≈ get_code_doppler(alone, :g, 1)
    end

    @testset "every signal accumulates its own raw discriminators" begin
        # The passenger's slot must hold what the passenger reads tracked alone,
        # unshifted by the group delay, whether it was combined or not. One
        # chunk, so both satellites read with the same code Doppler.
        passenger_alone =
            feed!(vt_on!(vector_state((GPSL5I(),))), (0.3,); late = (0.8,), chunks = 1)
        driver_alone = feed!(vt_on!(vector_state((GPSL5Q(),))), (0.1,); chunks = 1)
        for (combining, passenger_index) in
            ((true, SAMPLES), (false, SAMPLES), (true, SAMPLES - 1))
            ts = feed!(
                vt_on!(with_delays!(vector_state(signals; combining))),
                (0.1, 0.3);
                passenger_index,
                chunks = 1,
            )
            sat = Tracking.get_sat_state(ts, :g, 1)
            state = Tracking.get_doppler_estimator_state(sat)
            @test state.code_discr_acc[1][1] == 2
            @test state.code_discr_acc[2][1] == 2
            @test state.carrier_discr_acc[2][1] == 2
            @test mean_code_discr(ts, :g, 1, GPSL5Q) == mean_code_discr(driver_alone, :g, 1)
            if passenger_index == SAMPLES
                @test mean_code_discr(ts, :g, 1, GPSL5I) ==
                      mean_code_discr(passenger_alone, :g, 1)
                @test mean_carrier_discr(sat, 2) ==
                      mean_carrier_discr(passenger_alone, :g, 1)
            end
        end
        # A multi-signal satellite needs a signal selector.
        multi = feed!(vt_on!(vector_state(signals)), (0.1, 0.3))
        @test_throws Exception mean_code_discr(
            Tracking.get_doppler_estimator_state(Tracking.get_sat_state(multi, :g, 1)),
        )
    end
end

@static if VERSION >= v"1.11"
    @testset "vector combining is allocation-free (vt_on = $vt_on)" for vt_on in
                                                                        (false, true)
        ts = combining_state(
            (GPSL5Q(), GPSL5I());
            estimator = VectorPLLAndDLL(; discriminator_combining = true),
        )
        set_group_delay!(ts, :g, 1, 1, 0.0s)
        set_group_delay!(ts, :g, 1, 2, 0.5ns)
        vt_on && enable_vt!(ts, (1,))
        sat = Tracking.get_sat_state(ts, :g, 1)
        outputs = map(sat.signals) do signal
            rotation = SAMPLES * cis(0.1 + get_carrier_phase_offset(signal.signal))
            taps = (0.25, 0.5, 1.0, 0.6, 0.25) .* rotation
            CorrelatorOutput(with_taps(signal.correlator, taps), SAMPLES, SAMPLES)
        end
        fold_records!(ts, outputs, 10)
        @test (@allocated fold_records!(ts, outputs, 100)) == 0
    end
end

@testset "a record correlated before a same-chunk sync is excluded" begin
    # Synced by a record earlier in this chunk: the bit buffer reports `found`
    # although the chunk started unsynced (`found_before_fold = false`).
    synced(t) = TrackedSignal(
        t;
        bit_buffer = typeof(t.bit_buffer)(
            zero(typeof(t.bit_buffer.code_block_buffer)),
            0,
            true,
            0,
            Int8(+1),
            complex(0.0, 0.0),
            0,
            Float32[],
            Tracking.PhaseAccumulators(),
        ),
    )
    excluded(signal; found_before_fold) = last(
        Tracking._apply_passenger_record(
            synced(TrackedSignal(signal)),
            CorrelatorOutput(get_correlator(TrackedSignal(signal)), SAMPLES, SAMPLES),
            found_before_fold,
            1,
            FS,
            nothing,
            false,
            0.0,
        ),
    )
    # GPS L5I carries a secondary code, whose wipe-off the record lacked.
    @test excluded(GPSL5I(); found_before_fold = false)
    @test !excluded(GPSL5I(); found_before_fold = true)
    # GPS L1 C/A has none, so its record correlated the same either side.
    @test !excluded(GPSL1CA(); found_before_fold = false)
end

end
