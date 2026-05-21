module SatStateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using Dictionaries: dictionary
using GNSSSignals: GPSL1CA, get_code_center_frequency_ratio
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
    get_signals,
    get_sat_state,
    has_bit_or_secondary_code_been_found,
    get_bit_buffer,
    to_dictionary,
    max_code_length

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

    acq_args = (
        gpsl1,
        5,
        5e6Hz,
        100.0Hz,
        524.6,
        45.0,
        1.0,
        10.0,
        1,
        randn(100, 100),
        -500:100.0:500,
    )
    acq = pkgversion(Acquisition) >= v"2" ?
        AcquisitionResults(acq_args..., 1, length(-500:100.0:500)) :
        AcquisitionResults(acq_args...)
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
            zero(UInt128),
            0,
            complex(0.0, 0.0),
            0,
        ),
    )
    synced_signals = (synced_sig,)
    @test @inferred(Tracking.current_code_wrap(synced_signals)) == 20460

    # `_max_code_length` recursion terminator on an empty tuple. Not
    # reachable via TrackedSat itself (signals tuple is always non-empty),
    # but exercised here so the base case counts as covered.
    @test Tracking._max_code_length(()) == 0
    @test Tracking._current_code_wrap(()) == 0

    @testset "_post_sync_code_length per signal" begin
        # Worst-case wrap contributions across all supported signals.
        import GNSSSignals: GalileoE1B, GalileoE1B_BOC11, GPSL5I, GPSL1C_D, GPSL1C_P
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1CA())) == 1023 * 20
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B())) == 4092 * 1   # 1 block per symbol
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B_BOC11())) == 4092 * 1   # BOC(1,1) approximation — same 1-block-per-symbol shape
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL5I())) == 10230 * 10
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_D())) == 10230 * 1
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_P())) == 10230 * 1800
    end
end

end
