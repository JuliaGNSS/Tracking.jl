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
    get_chip_rate_to_center_frequency_ratio,
    get_code_length,
    get_chip_rate
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
    get_bit_buffer,
    to_dictionary,
    max_code_length
import TrackingLoops
using TrackingLoops: has_bit_or_secondary_code_been_found

# `_make_acq`, an Acquisition-version shim.
include("acquisition_test_helpers.jl")

@testset "Satellite state" begin
    gpsl1 = GPSL1CA()
    sat_state = @inferred TrackedSat(gpsl1, 1, 10.0, 500.0Hz)
    @test get_prn(sat_state) == 1
    @test get_code_phase(sat_state) == 10.0
    @test get_code_doppler(sat_state) ==
          500.0Hz * get_chip_rate_to_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 500.0Hz
    @test get_integrated_samples(sat_state) == 0
    @test get_signal_start_sample(sat_state) == 1
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test has_bit_or_secondary_code_been_found(sat_state) == false
    @test length(get_bit_buffer(sat_state)) == 0

    acq = _make_acq(gpsl1, 5, 524.6, 100.0Hz)
    sat_state = @inferred TrackedSat(acq)
    @test get_prn(sat_state) == 5
    @test get_code_phase(sat_state) == 524.6
    @test get_code_doppler(sat_state) ==
          100.0Hz * get_chip_rate_to_center_frequency_ratio(gpsl1)
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

    sigs = @inferred get_signals(sat)
    @test length(sigs) == 1
    @test sigs[1].signal isa GPSL1CA

    d = dictionary([11 => sat])
    @test to_dictionary(d) === d

    @test get_sat_state(d).prn == 11

    # The upper bound, 1023 chips × 20 blocks per bit; the runtime wrap is
    # 1023 until sync, then 20460.
    @test @inferred(max_code_length(sat.signals)) == 20460

    @test @inferred(Tracking.current_code_wrap(sat.signals)) == 1023

    synced_sig = Tracking.TrackedSignal(
        only(sat.signals);
        bit_buffer = TrackingLoops.BitBuffer{UInt64}(
            zero(UInt64),
            0,
            true,                # found
            0,
            Int8(+1),
            complex(0.0, 0.0),
            0,
            Float32[],
            TrackingLoops.PhaseAccumulators(),
        ),
    )
    synced_signals = (synced_sig,)
    @test @inferred(Tracking.current_code_wrap(synced_signals)) == 20460

    # The recursion's base case, the lcm identity; a TrackedSat never has an
    # empty signal tuple.
    @test Tracking._max_code_length(()) == 1
    @test Tracking._current_code_wrap(()) == 1

    @testset "_post_sync_code_length per signal" begin
        import GNSSSignals: GalileoE1B, GalileoE1B_BOC11, GPSL5I, GPSL1C_D, GPSL1C_P
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1CA())) ==
              1023 * 20
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B())) ==
              4092 * 1   # 1 block per symbol
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B_BOC11())) ==
              4092 * 1
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL5I())) ==
              10230 * 10
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_D())) ==
              10230 * 1
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_P())) ==
              10230 * 1800
    end
end

# A pilot with arbitrary code lengths: real signal pairings have wraps that
# divide each other, so lcm and max agree on them (#129).
struct FakeWrapSignal <: AbstractGNSSSignal{Matrix{Int16}}
    code_length::Int
    secondary_code_length::Int
end
GNSSSignals.get_code_length(s::FakeWrapSignal) = s.code_length
GNSSSignals.get_secondary_code_length(s::FakeWrapSignal) = s.secondary_code_length
GNSSSignals.get_bit_rate(::Type{FakeWrapSignal}) = 0Hz

@testset "shared code wrap is a common multiple, not a max (issue #129)" begin
    # Each signal re-derives its phase as `mod(code_phase, wrap)`, correct only
    # if the shared wrap is a multiple of every signal's.
    base = Tracking.TrackedSignal(GPSL1CA())
    synced_buffer(::TrackingLoops.BitBuffer{B}) where {B} = TrackingLoops.BitBuffer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(+1),
        complex(0.0, 0.0),
        0,
        Float32[],
        TrackingLoops.PhaseAccumulators(),
    )
    fake_tracked_signal(code_length, secondary_length, found) = Tracking.TrackedSignal(
        FakeWrapSignal(code_length, secondary_length),
        base.integrated_samples,
        base.correlator,
        base.last_fully_integrated_correlator,
        Tracking.TrackedSignal(
            base;
            bit_buffer = found ? synced_buffer(base.loop_state.bit_buffer) :
                         base.loop_state.bit_buffer,
        ).loop_state,
        base.filtered_prompts,
        base.correlator_outputs,
        base.preferred_num_code_blocks_to_integrate,
    )

    # Wraps 4 and 6: 12, where max would give 6, no multiple of 4.
    s4 = fake_tracked_signal(4, 5, false)
    s6 = fake_tracked_signal(6, 1, false)
    @test @inferred(Tracking.current_code_wrap((s4, s6))) == 12

    # Synced, the 4-chip signal wraps at 4 × 5 = 20.
    s4_synced = fake_tracked_signal(4, 5, true)
    @test @inferred(Tracking.current_code_wrap((s4_synced, s6))) == 60

    # The bound is the all-synced case.
    @test @inferred(max_code_length((s4, s6))) == 60
end

@testset "Public multi-signal TrackedSat constructor" begin
    # #133.
    sigs = (GPSL1C_P(), GPSL1C_D(), GPSL1CA())
    sat = @inferred TrackedSat(sigs, 11, 1234.5, 500.0Hz)
    @test get_prn(sat) == 11
    @test get_code_phase(sat) == 1234.5
    @test get_carrier_doppler(sat) == 500.0Hz
    # Code Doppler from the driver, `signals[1]`.
    @test get_code_doppler(sat) ==
          500.0Hz * get_chip_rate_to_center_frequency_ratio(GPSL1C_P())
    @test get_signal_start_sample(sat) == 1
    sigs_on_sat = get_signals(sat)
    @test length(sigs_on_sat) == 3
    @test sigs_on_sat[1].signal isa GPSL1C_P
    @test sigs_on_sat[2].signal isa GPSL1C_D
    @test sigs_on_sat[3].signal isa GPSL1CA
    @test sat.doppler_estimator_state isa TrackingLoops.SatConventionalPLLAndDLL

    sat_tuple = @inferred TrackedSat((GPSL1CA(),), 1, 10.0, 500.0Hz)
    sat_scalar = @inferred TrackedSat(GPSL1CA(), 1, 10.0, 500.0Hz)
    @test typeof(sat_tuple) === typeof(sat_scalar)

    # The auto-bandwidth estimator sizes the loop from the driver, GPS L1C-P.
    estimator = TrackingLoops.ConventionalAssistedPLLAndDLL()
    @test estimator.wide_carrier_loop_filter_bandwidth === nothing
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
    @test sat_kw.doppler_estimator_state.bandwidths.wide_carrier ==
          TrackingLoops.default_wide_carrier_loop_filter_bandwidth(GPSL1C_P())
end

@testset "Last fully integrated integration time" begin
    # Post-integration SNR is C/N₀ · T, so judging the last record needs its
    # own integration time.
    gpsl1 = GPSL1CA()
    code_period = get_code_length(gpsl1) / get_chip_rate(gpsl1)
    @test upreferred(code_period) ≈ 1ms

    tsig = Tracking.TrackedSignal(gpsl1)
    # Nothing completed yet: one block.
    @test get_last_fully_integrated_num_code_blocks(tsig) == 1
    @test get_last_fully_integrated_integration_time(tsig) ≈ code_period

    for num_blocks in (1, 2, 20)
        t = Tracking.TrackedSignal(tsig; last_fully_integrated_num_code_blocks = num_blocks)
        @test get_last_fully_integrated_num_code_blocks(t) == num_blocks
        @test get_last_fully_integrated_integration_time(t) ≈ num_blocks * code_period
    end

    twenty = Tracking.TrackedSignal(tsig; last_fully_integrated_num_code_blocks = 20)
    @test upreferred(get_last_fully_integrated_integration_time(twenty)) ≈ 20ms

    # Not `get_integrated_samples`, which counts the record in progress.
    @test get_integrated_samples(twenty) == 0

    sat = TrackedSat(gpsl1, 1, 10.0, 500.0Hz)
    @test get_last_fully_integrated_num_code_blocks(sat) == 1
    @test get_last_fully_integrated_integration_time(sat) ≈ code_period
    @test get_last_fully_integrated_num_code_blocks(sat, 1) == 1
    @test get_last_fully_integrated_integration_time(sat, 1) ≈ code_period

    track_state = TrackState(gpsl1, [sat])
    @test get_last_fully_integrated_num_code_blocks(track_state, 1) == 1
    @test get_last_fully_integrated_integration_time(track_state, 1) ≈ code_period
end

end
