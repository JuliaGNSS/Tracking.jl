module SatStateTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using Dictionaries: dictionary
using GNSSSignals:
    GNSSSignals, AbstractGNSSSignal, GPSL1CA, get_code_center_frequency_ratio
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
    # Construct across Acquisition's compat range: v1 has 11 fields; v2.0–v2.4
    # appended num_blocks/block_size (13 fields); v2.5+ inserted
    # secondary_code_phase (after code_phase) and appended
    # num_secondary_rotations (15 fields).
    acq = if pkgversion(Acquisition) >= v"2.5"
        AcquisitionResults(
            acq_args[1:5]..., nothing, acq_args[6:end]...,
            1, length(-500:100.0:500), 1,
        )
    elseif pkgversion(Acquisition) >= v"2"
        AcquisitionResults(acq_args..., 1, length(-500:100.0:500))
    else
        AcquisitionResults(acq_args...)
    end
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
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1CA())) == 1023 * 20
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B())) == 4092 * 1   # 1 block per symbol
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GalileoE1B_BOC11())) == 4092 * 1   # BOC(1,1) approximation — same 1-block-per-symbol shape
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL5I())) == 10230 * 10
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_D())) == 10230 * 1
        @test Tracking._post_sync_code_length(Tracking.TrackedSignal(GPSL1C_P())) == 10230 * 1800
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
        zero(B), 0, true, 0, Int8(+1), zero(UInt128), 0, complex(0.0, 0.0), 0, Float32[],
        Tracking.PhaseAccumulators(),
    )
    fake_tracked_signal(code_length, secondary_length, found) = Tracking.TrackedSignal(
        FakeWrapSignal(code_length, secondary_length),
        base.integrated_samples,
        base.is_integration_completed,
        base.correlator,
        base.last_fully_integrated_correlator,
        base.last_fully_integrated_filtered_prompt,
        base.cn0_estimator,
        found ? synced_buffer(base.bit_buffer) : base.bit_buffer,
        base.post_corr_filter,
        base.filtered_prompts,
        base.preferred_num_code_blocks_to_integrate,
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

end
