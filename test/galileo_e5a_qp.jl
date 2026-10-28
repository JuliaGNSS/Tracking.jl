module GalileoE5aQPTest

# Galileo E5a-QP — the E5a quick-acquisition aid (OS SIS ICD v2.2 §2.3.1.4).
#
# Dataless BPSK(5) at 5.115 Mcps, a 330-chip code repeated 31 times in 2 ms with
# no overlay. One block is 64.5 µs, so it integrates whole 31-block cycles, or
# the loop would run at 15.5 kHz (#236).

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz, dBHz, ms, upreferred, ustrip
using GNSSSignals:
    GalileoE5aI,
    GalileoE5aQ,
    GalileoE5aQP,
    GPSL5I,
    gen_code,
    get_band_id,
    get_chip_rate_to_center_frequency_ratio,
    get_chip_rate,
    get_code_length,
    get_bit_rate,
    get_secondary_code_length
import Tracking
using Tracking:
    SignalGroup,
    TrackState,
    TrackedSat,
    TrackedSignal,
    get_carrier_doppler,
    get_code_phase,
    get_last_fully_integrated_num_code_blocks,
    get_num_bits,
    get_preferred_num_code_blocks_to_integrate,
    get_sat_state,
    track
import TrackingLoops
using TrackingLoops:
    default_narrow_carrier_loop_filter_bandwidth,
    effective_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    estimate_cn0,
    default_num_code_blocks_to_integrate,
    max_num_code_blocks_to_integrate,
    has_bit_or_secondary_code_been_found

# 31 blocks = 10230 chips = 2 ms, the ICD's "repeated 31 times within 2 ms".
const BLOCKS_PER_CYCLE = 31

@testset "Galileo E5a-QP" begin
    e5a_qp = GalileoE5aQP()
    prn = 6

    # The reproducer from #236.
    @test TrackedSat(e5a_qp, prn, 0.0, 0.0Hz) isa TrackedSat

    @testset "Signal shape the policy rests on" begin
        @test get_code_length(e5a_qp) == 330
        @test get_chip_rate(e5a_qp) == 5_115_000Hz
        @test get_bit_rate(e5a_qp) == 0Hz
        @test get_secondary_code_length(e5a_qp) == 1
        @test BLOCKS_PER_CYCLE * get_code_length(e5a_qp) == 10230
        @test upreferred(
            BLOCKS_PER_CYCLE * get_code_length(e5a_qp) / get_chip_rate(e5a_qp),
        ) ≈ upreferred(2ms)
        @test get_band_id(e5a_qp) === :L5
        @test get_band_id(GalileoE5aI()) === :L5
    end

    @testset "Short-period integration policy" begin
        # One 2 ms cycle is ceiling and default; longer is a coherence
        # question for the caller.
        @test @inferred(max_num_code_blocks_to_integrate(e5a_qp)) == BLOCKS_PER_CYCLE
        @test @inferred(default_num_code_blocks_to_integrate(e5a_qp)) == BLOCKS_PER_CYCLE
        @test get_preferred_num_code_blocks_to_integrate(TrackedSignal(e5a_qp)) ==
              BLOCKS_PER_CYCLE

        for signal in (GalileoE5aI(), GalileoE5aQ(), GPSL5I())
            @test @inferred(default_num_code_blocks_to_integrate(signal)) == 1
        end
        @test max_num_code_blocks_to_integrate(GalileoE5aI()) == 20
        @test max_num_code_blocks_to_integrate(GalileoE5aQ()) == 100

        # The flat 18 Hz carrier default, uncapped at the 2 ms cycle.
        BL = @inferred default_narrow_carrier_loop_filter_bandwidth(e5a_qp)
        @test BL ≈ 18.0Hz
        cycle_time =
            upreferred(BLOCKS_PER_CYCLE * get_code_length(e5a_qp) / get_chip_rate(e5a_qp))
        @test effective_carrier_loop_filter_bandwidth(BL, cycle_time) ≈ 18.0Hz
        @test @inferred(default_code_loop_filter_bandwidth(e5a_qp)) ≈ 1.0Hz
    end

    @testset "Not a co-tracked E5a component — chip rate rules it out" begin
        # Half the E5a-I/Q chip rate, so it cannot share their group (#151).
        @test get_chip_rate(e5a_qp) != get_chip_rate(GalileoE5aI())
        @test_throws ArgumentError SignalGroup((GalileoE5aI(), e5a_qp))
        @test SignalGroup((e5a_qp,)) isa SignalGroup
    end

    @testset "Tracks a clean replica in whole 2 ms cycles" begin
        start_code_phase = 100
        # 5 samples per chip, so one primary block is exactly 1650 samples.
        sampling_frequency = 5 * get_chip_rate(e5a_qp)
        samples_per_block = round(Int, upreferred(get_code_length(e5a_qp) * 5))
        # One pre-sync block, then three whole 31-block cycles.
        num_blocks = 1 + 3 * BLOCKS_PER_CYCLE
        num_samples = num_blocks * samples_per_block

        code = gen_code(
            num_samples,
            e5a_qp,
            prn,
            sampling_frequency,
            get_chip_rate(e5a_qp),
            start_code_phase,
        )
        samples = ComplexF32.(code)

        track_state = TrackState(e5a_qp, [TrackedSat(e5a_qp, prn, start_code_phase, 0.0Hz)])
        track_state = track(samples, track_state, sampling_frequency)
        sat_state = get_sat_state(track_state, prn)

        @test get_code_phase(sat_state) ≈ start_code_phase atol = 0.1
        @test has_bit_or_secondary_code_been_found(sat_state) == true
        @test get_last_fully_integrated_num_code_blocks(sat_state) == BLOCKS_PER_CYCLE

        @test get_num_bits(sat_state) == 0
    end

    @testset "Pulls a Doppler offset in over 2 ms integrations" begin
        # A 5 Hz error closed over 240 ms of 2 ms integrations at 18 Hz.
        doppler = 200.0Hz
        samples_per_chip = 4
        sampling_frequency = samples_per_chip * get_chip_rate(e5a_qp)
        code_frequency =
            get_chip_rate(e5a_qp) +
            doppler * get_chip_rate_to_center_frequency_ratio(e5a_qp)
        num_cycles = 120
        samples_per_cycle = BLOCKS_PER_CYCLE * samples_per_chip * get_code_length(e5a_qp)
        num_samples = num_cycles * samples_per_cycle

        code = gen_code(num_samples, e5a_qp, prn, sampling_frequency, code_frequency, 0.0)
        samples = ComplexF32.(
            cis.(
                2π .* ustrip(Hz, doppler) .* (0:(num_samples-1)) ./
                ustrip(Hz, sampling_frequency),
            ) .* code,
        )

        track_state = TrackState(e5a_qp, [TrackedSat(e5a_qp, prn, 0.0, 195.0Hz)])
        for cycle = 1:num_cycles
            track_state = track(
                view(samples, ((cycle-1)*samples_per_cycle+1):(cycle*samples_per_cycle)),
                track_state,
                sampling_frequency,
            )
        end
        sat_state = get_sat_state(track_state, prn)

        @test get_carrier_doppler(sat_state) ≈ doppler atol = 2.0Hz
        # On the real axis, with no energy lost to the 2 ms integration.
        prompt = Tracking.get_last_fully_integrated_filtered_prompt(sat_state)
        @test real(prompt) / abs(prompt) > 0.9
        @test abs(prompt) > 0.9
        @test get_last_fully_integrated_num_code_blocks(sat_state) == BLOCKS_PER_CYCLE
        @test get_num_bits(sat_state) == 0
        @test estimate_cn0(track_state, prn) > 40.0dBHz
    end
end

end
