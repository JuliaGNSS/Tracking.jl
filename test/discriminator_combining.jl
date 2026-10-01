module DiscriminatorCombiningTest

using Test: @test, @testset, @test_throws
using Unitful: Hz, MHz, s, ns, m, DimensionError
using StaticArrays: SVector
using Logging: with_logger, NullLogger
using GNSSSignals:
    GalileoE1B, GalileoE1C, GPSL5I, GPSL5Q, GPSL1C_D, GPSL1C_P, get_carrier_phase_offset
using Tracking:
    Tracking,
    TrackedSignal,
    TrackState,
    ConventionalAssistedPLLAndDLL,
    ConventionalPLLAndDLL,
    CorrelatorOutput,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    ThirdOrderBilinearLF,
    ThirdOrderAssistedBilinearLF,
    SecondOrderBilinearLF,
    SatConventionalPLLAndDLL,
    add_satellite!,
    append_correlator_output!,
    estimate_dopplers_and_filter_prompt!,
    get_carrier_doppler,
    get_code_doppler,
    get_correlator,
    get_correlator_outputs,
    get_filtered_prompts,
    get_group_delay,
    set_group_delay!

@testset "combine_discriminators is a type parameter" begin
    @test ConventionalPLLAndDLL() isa
          ConventionalPLLAndDLL{ThirdOrderBilinearLF,SecondOrderBilinearLF,false}
    combining = ConventionalPLLAndDLL(; combine_discriminators = true)
    @test combining isa
          ConventionalPLLAndDLL{ThirdOrderBilinearLF,SecondOrderBilinearLF,true}
    @test ConventionalAssistedPLLAndDLL(; combine_discriminators = true) ===
          ConventionalPLLAndDLL(ThirdOrderAssistedBilinearLF; combine_discriminators = true)
    # The bandwidth-update constructor keeps it.
    @test ConventionalPLLAndDLL(combining; code_loop_filter_bandwidth = 2.0Hz) isa
          ConventionalPLLAndDLL{ThirdOrderBilinearLF,SecondOrderBilinearLF,true}

    # The per-sat state carries it, with one unknown group delay per signal.
    ts = TrackState(;
        signals = (e1 = (GalileoE1C(), GalileoE1B()),),
        doppler_estimator = combining,
    )
    add_satellite!(ts; prn = 11, group = :e1, code_phase = 0.0, carrier_doppler = 0.0Hz)
    state = Tracking.get_doppler_estimator_state(Tracking.get_sat_state(ts, :e1, 11))
    @test state isa
          SatConventionalPLLAndDLL{<:ThirdOrderBilinearLF,<:SecondOrderBilinearLF,2,true}
    @test isconcretetype(typeof(state))
    @test all(isnan, state.group_delays)

    # By keyword, the delays are given in their public form.
    by_keyword = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        group_delays = (0.0s, nothing),
        combine_discriminators = true,
    )
    @test by_keyword isa SatConventionalPLLAndDLL{<:Any,<:Any,2,true}
    @test get_group_delay(by_keyword, 1) === 0.0s
    @test isnothing(get_group_delay(by_keyword, 2))
    @test_throws ArgumentError SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        group_delays = (NaN * s,),
    )
end

@testset "group delay" begin
    ts = TrackState(;
        signals = (e1 = (GalileoE1C(), GalileoE1B()),),
        doppler_estimator = ConventionalPLLAndDLL(; combine_discriminators = true),
    )
    add_satellite!(ts; prn = 11, group = :e1, code_phase = 0.0, carrier_doppler = 0.0Hz)
    @test isnothing(get_group_delay(ts, :e1, 11, GalileoE1B))

    set_group_delay!(ts, :e1, 11, GalileoE1B, -0.5ns)
    @test get_group_delay(ts, :e1, 11, GalileoE1B) === -5.0e-10s
    @test isnothing(get_group_delay(ts, :e1, 11, GalileoE1C))

    set_group_delay!(ts, :e1, 11, 1, 0.0s)
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # `nothing` clears again: unknown, not zero.
    set_group_delay!(ts, :e1, 11, GalileoE1B, nothing)
    @test isnothing(get_group_delay(ts, :e1, 11, 2))
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # A group delay is a time.
    @test_throws DimensionError set_group_delay!(ts, :e1, 11, 2, 1.0m)
    @test_throws DimensionError set_group_delay!(ts, :e1, 11, 2, 1.0)
    # NaN is the stored form of "unknown", so it is not a value one can set.
    @test_throws ArgumentError set_group_delay!(ts, :e1, 11, 2, NaN * s)

    # Any time unit is accepted, and a reset keeps the values.
    set_group_delay!(ts, :e1, 11, 2, 1.5ns)
    @test get_group_delay(ts, :e1, 11, 2) ≈ 1.5e-9s
    Tracking.reset_loop_filters!(ts)
    @test get_group_delay(ts, :e1, 11, 2) ≈ 1.5e-9s
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # The per-satellite form takes the signal too.
    sat = Tracking.get_sat_state(ts, :e1, 11)
    @test get_group_delay(sat, GalileoE1B) ≈ 1.5e-9s
    @test get_group_delay(sat, 1) === 0.0s

    # A satellite that does not combine stores its group delays all the same.
    conventional = TrackState(; signals = (e1 = (GalileoE1C(), GalileoE1B()),))
    add_satellite!(
        conventional;
        prn = 11,
        group = :e1,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    set_group_delay!(conventional, :e1, 11, 2, 0.5ns)
    @test get_group_delay(conventional, :e1, 11, 2) ≈ 5.0e-10s
    @test isnothing(get_group_delay(conventional, :e1, 11, 1))
end

const FS = 5e6Hz
const SAMPLES = 5000

# `combining = false` is the conventional estimator, with the same loop filters.
function combining_state(signals; combining = true, assisted = false)
    carrier_filter = assisted ? ThirdOrderAssistedBilinearLF : ThirdOrderBilinearLF
    estimator = ConventionalPLLAndDLL(carrier_filter; combine_discriminators = combining)
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
    @testset "nothing to combine: bit-identical to the conventional estimator" begin
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
function fold_records!(ts, outputs, n, sampling_frequency = FS)
    sats = Tracking.get_sat_states(ts, :g)
    for _ = 1:n
        sat = sats[1]
        foreach(sat.signals, outputs) do signal, output
            push!(Tracking.get_correlator_outputs(signal), output)
            empty!(get_filtered_prompts(signal))
        end
        sats[1] = Tracking._update_tracked_sat_doppler(
            sat,
            sampling_frequency,
            ((nothing, false), (nothing, false)),
        )
    end
end

@static if VERSION >= v"1.11"
    # A sampling frequency in `MHz` still gives FLL readings in `Hz` (`fll_disc`),
    # so the passenger sums keep their concrete type.
    @testset "combining is allocation-free (assisted = $assisted, fs = $fs)" for (
        assisted,
        fs,
    ) in Iterators.product(
        (false, true),
        (FS, 5.0MHz),
    )
        ts = combining_state((GPSL5Q(), GPSL5I()); assisted)
        set_group_delay!(ts, :g, 1, 1, 0.0s)
        set_group_delay!(ts, :g, 1, 2, 0.5ns)
        sat = Tracking.get_sat_state(ts, :g, 1)
        outputs = map(sat.signals) do signal
            rotation = SAMPLES * cis(0.1 + get_carrier_phase_offset(signal.signal))
            taps = (0.25, 0.5, 1.0, 0.6, 0.25) .* rotation
            CorrelatorOutput(with_taps(signal.correlator, taps), SAMPLES, SAMPLES)
        end
        fold_records!(ts, outputs, 10, fs)
        @test (@allocated fold_records!(ts, outputs, 100, fs)) == 0
    end
end

end
