module SatStateTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms, upreferred
using Dictionaries: dictionary
using GNSSSignals:
    GNSSSignals,
    AbstractGNSSSignal,
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    get_code_center_frequency_ratio,
    get_code_length,
    get_code_frequency
using Acquisition: Acquisition, AcquisitionResults
import Tracking
using Tracking:
    TrackedSat,
    get_prn,
    get_code_phase,
    get_code_doppler,
    get_carrier_phase,
    get_carrier_doppler,
    get_integrated_samples,
    get_signal_start_sample,
    get_correlator,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_filtered_prompt,
    get_last_fully_integrated_num_code_blocks,
    get_last_fully_integrated_integration_time,
    TrackState,
    get_signals,
    get_sat_state,
    has_bit_or_secondary_code_been_found,
    get_bit_buffer,
    to_dictionary,
    max_code_length

# `_make_acq`: shared Acquisition-version shim for building
# `AcquisitionResults` — see test/acquisition_test_helpers.jl.
include("acquisition_test_helpers.jl")

@testset "Satellite state" begin
    gpsl1 = GPSL1CA()
    sat_state = @inferred TrackedSat(gpsl1, 1, 10.0, 500.0Hz)
    @test get_prn(sat_state) == 1
    @test get_code_phase(sat_state) == 10.0
    @test get_code_doppler(sat_state) == 500.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 500.0Hz
    @test get_integrated_samples(sat_state) == 0
    @test get_signal_start_sample(sat_state) == 1
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test has_bit_or_secondary_code_been_found(sat_state) == false
    @test length(get_bit_buffer(sat_state)) == 0

    # `_make_acq` is the shared Acquisition-version shim — see
    # test/acquisition_test_helpers.jl.
    acq = _make_acq(gpsl1, 5, 524.6, 100.0Hz)
    sat_state = @inferred TrackedSat(acq)
    @test get_prn(sat_state) == 5
    @test get_code_phase(sat_state) == 524.6
    @test get_code_doppler(sat_state) == 100.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 100.0Hz
    @test get_integrated_samples(sat_state) == 0.0
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_filtered_prompt(sat_state) == complex(0.0, 0.0)
    @test has_bit_or_secondary_code_been_found(sat_state) == false
    @test length(get_bit_buffer(sat_state)) == 0
end

@testset "TrackedSat signal/dict helpers" begin
    gpsl1 = GPSL1CA()
    sat = TrackedSat(gpsl1, 11, 10.5, 100.0Hz)

    # `get_signals` returns the per-signal tuple.
    sigs = @inferred get_signals(sat)
    @test length(sigs) == 1
    @test sigs[1].signal isa GPSL1CA

    # `to_dictionary` on an already-dictionary input is a no-op.
    d = dictionary([11 => sat])
    @test to_dictionary(d) === d

    # `get_sat_state(::Dictionary)` (no identifier) returns the only sat.
    @test get_sat_state(d).prn == 11

    # `max_code_length` is the upper-bound wrap for L1 C/A: 1023 chips
    # × 20 blocks per data bit = 20460. The runtime wrap returned by
    # `current_code_wrap` shrinks to 1023 until bit-edge sync, then
    # widens to 20460 — exercised below.
    @test @inferred(max_code_length(sat.signals)) == 20460

    # Before sync `bit_buffer.found == false`, so `current_code_wrap`
    # returns just the primary code length.
    @test @inferred(Tracking.current_code_wrap(sat.signals)) == 1023

    # After sync the runtime wrap widens to the full data-bit period.
    synced_sig = Tracking.TrackedSignal(
        only(sat.signals);
        bit_buffer = Tracking.BitBuffer{UInt64}(
            zero(UInt64),
            0,
            true,                # found
            0,
            Int8(+1),
            complex(0.0, 0.0),
            0,
            Float32[],
            Tracking.PhaseAccumulators(),
        ),
    )
    synced_signals = (synced_sig,)
    @test @inferred(Tracking.current_code_wrap(synced_signals)) == 20460

    # `_max_code_length` recursion terminator on an empty tuple — 1, the
    # lcm identity. Not reachable via TrackedSat itself (signals tuple is
    # always non-empty), but exercised here so the base case counts as
    # covered.
    @test Tracking._max_code_length(()) == 1
    @test Tracking._current_code_wrap(()) == 1

    @testset "_post_sync_code_length per signal" begin
        # Worst-case wrap contributions across all supported signals.
        import GNSSSignals: GalileoE1B, GalileoE1B_BOC11, GPSL5I, GPSL1C_D, GPSL1C_P
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1CA())) ==
              1023 * 20
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B())) ==
              4092 * 1   # 1 block per symbol
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B_BOC11())) ==
              4092 * 1   # BOC(1,1) approximation — same 1-block-per-symbol shape
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL5I())) ==
              10230 * 10
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_D())) ==
              10230 * 1
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_P())) ==
              10230 * 1800
    end
end

# Pilot-style signal (data frequency 0) with arbitrary primary / secondary
# code lengths. All real signal pairings have wrap periods that divide each
# other, so the lcm-vs-max distinction in the shared code wrap (issue #129)
# is only observable with fabricated lengths. Pre-sync the signal
# contributes `code_length`; post-sync `code_length × secondary_code_length`.
struct FakeWrapSignal <: AbstractGNSSSignal{Matrix{Int16}}
    code_length::Int
    secondary_code_length::Int
end
GNSSSignals.get_code_length(s::FakeWrapSignal) = s.code_length
GNSSSignals.get_secondary_code_length(s::FakeWrapSignal) = s.secondary_code_length
GNSSSignals.get_data_frequency(::FakeWrapSignal) = 0Hz

@testset "shared code wrap is a common multiple, not a max (issue #129)" begin
    # The per-signal code phase is re-derived as
    # `mod(code_phase, _replica_code_wrap(signal))`, which is only correct
    # when the shared wrap is an integer multiple of *every* signal's
    # replica wrap — `max` is not a common multiple in general.
    base = Tracking.TrackedSignal(GPSL1CA())
    synced_buffer(::Tracking.BitBuffer{B}) where {B} = Tracking.BitBuffer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(+1),
        complex(0.0, 0.0),
        0,
        Float32[],
        Tracking.PhaseAccumulators(),
    )
    fake_tracked_signal(code_length, secondary_length, found) = Tracking.TrackedSignal(
        FakeWrapSignal(code_length, secondary_length),
        base.integrated_samples,
        base.correlator,
        base.last_fully_integrated_correlator,
        base.last_fully_integrated_filtered_prompt,
        base.cn0_estimator,
        found ? synced_buffer(base.bit_buffer) : base.bit_buffer,
        base.post_corr_filter,
        base.filtered_prompts,
        base.correlator_outputs,
        base.preferred_num_code_blocks_to_integrate,
        base.last_fully_integrated_num_code_blocks,
        base.group_delay,
    )

    # Pre-sync wraps 4 and 6: the shared wrap must be 12 (max would give 6,
    # under which a phase of e.g. 5 re-derives as 1 for the 4-chip signal —
    # but so does a phase of 9, silently corrupting the 4-chip signal).
    s4 = fake_tracked_signal(4, 5, false)
    s6 = fake_tracked_signal(6, 1, false)
    @test @inferred(Tracking.current_code_wrap((s4, s6))) == 12

    # Once the 4-chip signal syncs, its wrap widens to 4 × 5 = 20 → lcm 60.
    s4_synced = fake_tracked_signal(4, 5, true)
    @test @inferred(Tracking.current_code_wrap((s4_synced, s6))) == 60

    # The compile-time bound covers the all-synced worst case: lcm(20, 6).
    @test @inferred(max_code_length((s4, s6))) == 60
end

@testset "Public multi-signal TrackedSat constructor" begin
    # Issue #133 item 2: a tuple of signals must be constructible directly,
    # without hand-rolling the bare-sat → init_estimator_state → rebuild
    # two-stage build.
    sigs = (GPSL1C_P(), GPSL1C_D(), GPSL1CA())
    sat = @inferred TrackedSat(sigs, 11, 1234.5, 500.0Hz)
    @test get_prn(sat) == 11
    @test get_code_phase(sat) == 1234.5
    @test get_carrier_doppler(sat) == 500.0Hz
    # The default code doppler scales by the estimator-driver signal's
    # (signals[1]) code/center frequency ratio.
    @test get_code_doppler(sat) == 500.0Hz * get_code_center_frequency_ratio(GPSL1C_P())
    @test get_signal_start_sample(sat) == 1
    sigs_on_sat = get_signals(sat)
    @test length(sigs_on_sat) == 3
    @test sigs_on_sat[1].signal isa GPSL1C_P
    @test sigs_on_sat[2].signal isa GPSL1C_D
    @test sigs_on_sat[3].signal isa GPSL1CA
    @test sat.doppler_estimator_state isa Tracking.SatConventionalPLLAndDLL

    # A one-tuple builds the same concrete type as the scalar-signal form.
    sat_tuple = @inferred TrackedSat((GPSL1CA(),), 1, 10.0, 500.0Hz)
    sat_scalar = @inferred TrackedSat(GPSL1CA(), 1, 10.0, 500.0Hz)
    @test typeof(sat_tuple) === typeof(sat_scalar)

    # Kwargs flow through: explicit carrier phase, code doppler, estimator.
    # The default (auto-bandwidth) estimator sizes the sat's loop from its
    # own driver signal (signals[1] = GPS L1C-P → 1.8 Hz), not from a fixed
    # value on the estimator (which is `nothing` = auto).
    estimator = Tracking.ConventionalAssistedPLLAndDLL()
    @test estimator.carrier_loop_filter_bandwidth === nothing
    sat_kw = TrackedSat(
        sigs,
        7,
        0.25,
        -250.0Hz;
        carrier_phase = 0.5,
        code_doppler = -0.3Hz,
        doppler_estimator = estimator,
    )
    @test get_carrier_phase(sat_kw) ≈ 0.5
    @test get_code_doppler(sat_kw) == -0.3Hz
    @test sat_kw.doppler_estimator_state.carrier_loop_filter_bandwidth ==
          Tracking.default_carrier_loop_filter_bandwidth(GPSL1C_P())
end

@testset "Last fully integrated integration time" begin
    # C/N₀ is processing-independent, so a consumer asking a *detectability*
    # question of the last record — is the peak still above the noise? — needs
    # the record's own integration time, because post-integration SNR is
    # C/N₀ · T. Before these accessors that number was only reachable as a
    # private field, and the obvious public stand-in is a trap: see the
    # `get_integrated_samples` check below.
    gpsl1 = GPSL1CA()
    code_period = get_code_length(gpsl1) / get_code_frequency(gpsl1)
    @test upreferred(code_period) ≈ 1ms

    tsig = Tracking.TrackedSignal(gpsl1)
    # A fresh signal has completed nothing, so it reports one block — the
    # divisor `estimate_cn0` used before it knew any better.
    @test get_last_fully_integrated_num_code_blocks(tsig) == 1
    @test get_last_fully_integrated_integration_time(tsig) ≈ code_period

    for num_blocks in (1, 2, 20)
        t = Tracking.TrackedSignal(tsig; last_fully_integrated_num_code_blocks = num_blocks)
        @test get_last_fully_integrated_num_code_blocks(t) == num_blocks
        @test get_last_fully_integrated_integration_time(t) ≈ num_blocks * code_period
    end

    # A full GPS L1 C/A bit is 20 blocks = 20 ms.
    twenty = Tracking.TrackedSignal(tsig; last_fully_integrated_num_code_blocks = 20)
    @test upreferred(get_last_fully_integrated_integration_time(twenty)) ≈ 20ms

    # Not the same thing as `get_integrated_samples`, which counts the record
    # *currently* accumulating and is reset to zero whenever one completes — so
    # at the moment a consumer reads a completed record it is not the length of
    # anything.
    @test get_integrated_samples(twenty) == 0

    # Forwarded at the satellite and track-state levels like its siblings.
    sat = TrackedSat(gpsl1, 1, 10.0, 500.0Hz)
    @test get_last_fully_integrated_num_code_blocks(sat) == 1
    @test get_last_fully_integrated_integration_time(sat) ≈ code_period
    @test get_last_fully_integrated_num_code_blocks(sat, 1) == 1
    @test get_last_fully_integrated_integration_time(sat, 1) ≈ code_period

    track_state = TrackState(gpsl1, [sat])
    @test get_last_fully_integrated_num_code_blocks(track_state, 1) == 1
    @test get_last_fully_integrated_integration_time(track_state, 1) ≈ code_period
end

# --------------------------------------------------------------------------
# Per-signal group delay
# --------------------------------------------------------------------------

using Test: @test_logs, @test_throws
using Unitful: @u_str, DimensionError, s
using GNSSSignals: GPSL5I, GPSL5Q, GalileoE1B, GalileoE1C, GalileoE5aI, GalileoE5aQ
using Tracking: TrackedSignal, add_satellite!, get_group_delay, set_group_delay!

@testset "group delay: default, plumbing, and addressing" begin
    # Every signal starts unknown, whatever its constellation: this package holds
    # no table of per-signal values, because deciding one needs constellation
    # knowledge and usually a decoded navigation message. A signal defaulting to
    # `0.0` here would be this package guessing on the caller's behalf.
    for signal in (
        GPSL1CA(),
        GPSL1C_D(),
        GPSL1C_P(),
        GPSL5I(),
        GPSL5Q(),
        GalileoE1B(),
        GalileoE1C(),
        GalileoE5aI(),
        GalileoE5aQ(),
    )
        @test get_group_delay(TrackedSignal(signal)) === nothing
    end

    l1cp, l1cd, l1ca = GPSL1C_P(), GPSL1C_D(), GPSL1CA()
    ts = TrackState(; signals = (gps_l1 = (l1cp, l1cd, l1ca),))
    ts = add_satellite!(
        ts;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )

    # Addressed by signal type …
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_D, -1.5e-9s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_D) === -1.5e-9s
    # … and readable straight off the `TrackState`, like every other per-signal
    # setting.
    @test get_group_delay(ts, :gps_l1, 7, GPSL1C_D) === -1.5e-9s
    # … and by index, and the other signals are untouched.
    set_group_delay!(ts, :gps_l1, 7, 3, 3.0e-9s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), 3) === 3.0e-9s
    # Slot 1 is an ordinary slot: untouched by the writes above, it is still the
    # `nothing` every slot starts at, and it is settable like any other.
    @test isnothing(get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_P))
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_P, 0.0s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_P) === 0.0s
    # The group by position, and the single-group form without a group.
    set_group_delay!(ts, 1, 7, GPSL1C_P, 0.5e-9s)
    @test get_group_delay(ts, 1, 7, 1) === 0.5e-9s
    set_group_delay!(ts, 7, GPSL1C_P, 0.0s)
    @test get_group_delay(ts, :gps_l1, 7, 1) === 0.0s

    # `nothing` marks it unknown again — the one case a plain `Maybe` kwarg could
    # not express, hence the `Some` wrapper inside.
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_D, nothing)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_D) === nothing

    # Any time unit is accepted and normalized to seconds as a Float64 quantity,
    # so the field type stays concrete however the caller spells the value.
    set_group_delay!(ts, :gps_l1, 7, GPSL1CA, -1200u"ps")
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1CA) === -1.2e-9s
    set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 0s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1CA) === 0.0s

    # A bare number, and a quantity that is not a time, are both refused: the
    # field carries its unit like every dimensioned quantity here, and at
    # sub-nanosecond scale an assumed unit is a metre-scale mistake.
    @test_throws DimensionError set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 1.2e-9)
    @test_throws DimensionError TrackedSignal(GPSL1CA(); group_delay = 1.2e-9)
    @test_throws DimensionError set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 1.2Hz)
    # …except a length, which gets a sentence naming the conversion and the datum
    # trap that comes with it: an SSR/PPP "code bias" in metres is the likely
    # wrong input here, not an exotic one.
    @test_throws ArgumentError set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 0.45u"m")
    @test_throws ArgumentError TrackedSignal(GPSL1CA(); group_delay = 0.45u"m")
    @test occursin("Unitful.c0", sprint(showerror, try
        Tracking._as_group_delay(0.45u"m")
    catch e
        e
    end))
    # A signal selector is required on a multi-signal satellite.
    @test_throws MethodError set_group_delay!(ts, :gps_l1, 7, 1.0e-9s)
end

@testset "the datum is whichever slot you state, and only differences are read" begin
    # `_driver_relative_group_delay` is the one place the differencing rule is
    # written: datum minus delay, so a passenger delayed relative to the driver
    # (larger value) comes out negative, and either end unknown is `nothing`.
    @test Tracking._driver_relative_group_delay(3.0e-9s, 1.0e-9s) ≈ -2.0e-9s
    @test Tracking._driver_relative_group_delay(1.0e-9s, 3.0e-9s) ≈ 2.0e-9s
    @test isnothing(Tracking._driver_relative_group_delay(nothing, 1.0e-9s))
    @test isnothing(Tracking._driver_relative_group_delay(1.0e-9s, nothing))
    @test isnothing(Tracking._driver_relative_group_delay(nothing, nothing))
    # Shifting both ends by one constant changes nothing — the datum cancels.
    offset = 5.0e-9s
    @test Tracking._driver_relative_group_delay(1.2e-9s + offset, offset) ≈
          Tracking._driver_relative_group_delay(1.2e-9s, 0.0s)
    # Converted to chips at the code frequency in effect.
    @test Tracking._group_delay_to_chips(1.0e-6s, 1.023e6Hz) ≈ 1.023
end

@testset "a pre-built signals tuple states its own datum" begin
    # The `TrackedSignal`-tuple constructor is the full-control form and seeds
    # nothing, so a caller using it owns the datum as it owns the correlators.
    l1cp, l1cd = GPSL1C_P(), GPSL1C_D()
    sat = TrackedSat((TrackedSignal(l1cp), TrackedSignal(l1cd)), 1, 0.0, 1000.0Hz)
    @test isnothing(get_group_delay(sat, 1))
    @test isnothing(get_group_delay(sat, 2))

    sat = Tracking._set_sat_group_delay(sat, 2.0e-9s, 2)
    @test get_group_delay(sat, 2) == 2.0e-9s
    rebuilt = TrackedSat(
        sat;
        # `Some` because `nothing` is a legal value of this field, so the
        # copy-update constructor cannot read a bare `nothing` as "clear it".
        signals = (
            TrackedSignal(
                sat.signals[1];
                group_delay = Some{Union{Nothing,typeof(1.0s)}}(1.0e-9s),
            ),
            sat.signals[2],
        ),
    )
    @test get_group_delay(rebuilt, 1) == 1.0e-9s
    @test get_group_delay(rebuilt, 2) == 2.0e-9s

    # A value on slot 1 at construction is kept as given: no normalisation and
    # nothing refused.
    sat = TrackedSat(
        (
            TrackedSignal(l1cp; group_delay = 1.0e-9s),
            TrackedSignal(l1cd; group_delay = -1.0e-9s),
        ),
        7,
        0.0,
        1000.0Hz,
    )
    @test get_group_delay(sat, GPSL1C_P) == 1.0e-9s
    @test get_group_delay(sat, GPSL1C_D) == -1.0e-9s
    fresh = TrackedSat((l1cp, l1cd), 7, 0.0, 1000.0Hz)
    @test isnothing(get_group_delay(fresh, GPSL1C_P))
end

@testset "setting a group delay is estimator-agnostic and silent" begin
    l1cp, l1cd = GPSL1C_P(), GPSL1C_D()
    for estimator in (Tracking.ConventionalAssistedPLLAndDLL(), Tracking.VectorPLLAndDLL())
        ts = TrackState(; signals = (gps_l1 = (l1cp, l1cd),), doppler_estimator = estimator)
        ts = add_satellite!(
            ts;
            group = :gps_l1,
            prn = 7,
            code_phase = 0.0,
            carrier_doppler = 1000.0Hz,
        )
        @test_logs set_group_delay!(ts, :gps_l1, 7, GPSL1C_D, 1.0e-9s)
        @test get_group_delay(ts, :gps_l1, 7, GPSL1C_D) == 1.0e-9s
    end
end

end
