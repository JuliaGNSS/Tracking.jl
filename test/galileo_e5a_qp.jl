module GalileoE5aQPTest

# Galileo E5a-QP — the E5a quick-acquisition aid (OS SIS ICD v2.2 §2.3.1.4).
#
# A dataless BPSK(5) component at 5.115 Mcps whose primary code is only 330
# chips, repeated 31 times within 2 ms with no overlay. Tracking it is a
# *policy* decision, not just a dispatch method: one primary code block is
# 64.5 µs, so the tracker must group whole 31-block (2 ms) code cycles into one
# coherent integration, or the loop would run at 15.5 kHz off a 279 Hz
# reference bandwidth. These tests pin that policy, its arithmetic and the
# resulting end-to-end pass (issue #236).

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz, dBHz, ms, upreferred, ustrip
using GNSSSignals:
    GalileoE5aI,
    GalileoE5aQ,
    GalileoE5aQP,
    GPSL5I,
    gen_code,
    get_band_id,
    get_code_center_frequency_ratio,
    get_code_frequency,
    get_code_length,
    get_data_frequency,
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
    EarlyPromptLateCorrelator,
    NumAnts,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    detect_bit_or_secondary_code_sync,
    estimate_cn0,
    get_code_block_buffer_type,
    get_default_correlator,
    default_num_code_blocks_to_integrate,
    max_num_code_blocks_to_integrate,
    has_bit_or_secondary_code_been_found

# 31 primary code blocks = 10230 chips = 2 ms — the ICD's own "repeated 31
# times within 2 ms", and the unit this package integrates E5a-QP in.
const BLOCKS_PER_CYCLE = 31

@testset "Galileo E5a-QP" begin
    e5a_qp = GalileoE5aQP()
    prn = 6

    # The reproducer from issue #236: default construction used to throw a
    # MethodError from `get_default_correlator`.
    @test TrackedSat(e5a_qp, prn, 0.0, 0.0Hz) isa TrackedSat

    @testset "Signal shape the policy rests on" begin
        @test get_code_length(e5a_qp) == 330
        @test get_code_frequency(e5a_qp) == 5_115_000Hz
        # Dataless and overlay-free: no bit and no secondary code to sync on.
        @test get_data_frequency(e5a_qp) == 0Hz
        @test get_secondary_code_length(e5a_qp) == 1
        # 31 primary blocks are exactly 2 ms, and exactly 10230 chips.
        @test BLOCKS_PER_CYCLE * get_code_length(e5a_qp) == 10230
        @test upreferred(
            BLOCKS_PER_CYCLE * get_code_length(e5a_qp) / get_code_frequency(e5a_qp),
        ) ≈ upreferred(2ms)
        # E5a-QP shares the E5a carrier with E5a-I / E5a-Q.
        @test get_band_id(e5a_qp) === :L5
        @test get_band_id(GalileoE5aI()) === :L5
    end

    @testset "Correlator, sync and search-buffer defaults" begin
        # Plain BPSK(5) (`LOC`) → the C/A-style EarlyPromptLate default.
        for n in (1, 3)
            @test @inferred(get_default_correlator(e5a_qp, NumAnts(n))) ==
                  EarlyPromptLateCorrelator(; num_ants = NumAnts(n))
        end

        # No data bit and no overlay, so there is no boundary to find: the
        # detector reports sync immediately and unconditionally, which is what
        # unlocks the multi-block integration below.
        for (bits, n) in ((UInt8(0x0), 0), (UInt8(0x1), 1), (UInt8(0xff), 97))
            res = @inferred detect_bit_or_secondary_code_sync(e5a_qp, prn, bits, n)
            @test res.found == true
            @test res.phase == 0
            @test res.polarity == +1
        end

        # Neither soft detector applies (one needs a multi-block data bit, the
        # other an overlay), so the hard path above is the one that runs.
        @test TrackingLoops.uses_soft_bit_edge_detection(e5a_qp) == false
        @test TrackingLoops.uses_soft_secondary_code_detection(e5a_qp) == false

        # Nothing to search for, so the packed sign window is dead state.
        @test @inferred(get_code_block_buffer_type(e5a_qp)) === UInt8
    end

    @testset "Short-period integration policy" begin
        # One 2 ms code cycle is both the ceiling and the default: the code
        # repeats with no structure to straddle, so integrating longer is a
        # coherence question for the caller, not a correctness one.
        @test @inferred(max_num_code_blocks_to_integrate(e5a_qp)) == BLOCKS_PER_CYCLE
        @test @inferred(default_num_code_blocks_to_integrate(e5a_qp)) == BLOCKS_PER_CYCLE
        @test get_preferred_num_code_blocks_to_integrate(TrackedSignal(e5a_qp)) ==
              BLOCKS_PER_CYCLE

        # Every other signal keeps the historical single-block default.
        for signal in (GalileoE5aI(), GalileoE5aQ(), GPSL5I())
            @test @inferred(default_num_code_blocks_to_integrate(signal)) == 1
        end
        # …and their ceiling is still the data bit / overlay period.
        @test max_num_code_blocks_to_integrate(GalileoE5aI()) == 20
        @test max_num_code_blocks_to_integrate(GalileoE5aQ()) == 100

        # The carrier default is a per-primary-block *reference* bandwidth, so
        # for a 64.5 µs block it is a large 279 Hz — and the loop's automatic
        # 1/N scaling at the 31-block integration brings it to a sane 9 Hz.
        BL = @inferred default_carrier_loop_filter_bandwidth(e5a_qp)
        @test BL ≈ 0.018 / upreferred(get_code_length(e5a_qp) / get_code_frequency(e5a_qp))
        @test BL / BLOCKS_PER_CYCLE ≈ 9.0Hz rtol = 1e-3
        @test @inferred(default_code_loop_filter_bandwidth(e5a_qp)) ≈ 1.0Hz
    end

    @testset "Not a co-tracked E5a component — chip rate rules it out" begin
        # E5a-QP runs at half the E5a-I/Q chip rate, so it cannot share a
        # satellite's SignalGroup with them (see issue #151); it is acquired
        # and tracked on its own and handed over, not paired.
        @test get_code_frequency(e5a_qp) != get_code_frequency(GalileoE5aI())
        @test_throws ArgumentError SignalGroup((GalileoE5aI(), e5a_qp))
        @test SignalGroup((e5a_qp,)) isa SignalGroup
    end

    @testset "Tracks a clean replica in whole 2 ms cycles" begin
        start_code_phase = 100
        # 5 samples per chip, so one primary block is exactly 1650 samples.
        sampling_frequency = 5 * get_code_frequency(e5a_qp)
        samples_per_block = round(Int, upreferred(get_code_length(e5a_qp) * 5))
        # One pre-sync block, then three whole 31-block cycles.
        num_blocks = 1 + 3 * BLOCKS_PER_CYCLE
        num_samples = num_blocks * samples_per_block

        code = gen_code(
            num_samples,
            e5a_qp,
            prn,
            sampling_frequency,
            get_code_frequency(e5a_qp),
            start_code_phase,
        )
        samples = ComplexF32.(code)

        track_state = TrackState(e5a_qp, [TrackedSat(e5a_qp, prn, start_code_phase, 0.0Hz)])
        track_state = track(samples, track_state, sampling_frequency)
        sat_state = get_sat_state(track_state, prn)

        # Zero Doppler, so the code phase after a whole number of blocks is the
        # one that was seeded.
        @test get_code_phase(sat_state) ≈ start_code_phase atol = 0.1
        # The detector fires on the first block, so every later integration is
        # a full 2 ms cycle.
        @test has_bit_or_secondary_code_been_found(sat_state) == true
        @test get_last_fully_integrated_num_code_blocks(sat_state) == BLOCKS_PER_CYCLE

        # A pilot carries no navigation data: correlation and tracking work,
        # but no bit is ever emitted and none may be required of it.
        @test get_num_bits(sat_state) == 0
    end

    @testset "Pulls a Doppler offset in over 2 ms integrations" begin
        # The software reference check the policy is really about: a 20 Hz
        # carrier-Doppler error on a clean replica, closed by the loop over
        # 120 ms of 2 ms integrations. At the 279 Hz per-block reference
        # bandwidth the estimator's automatic 1/N scaling turns into a 9 Hz
        # effective one, which is a bandwidth a PLL can actually hold.
        doppler = 200.0Hz
        samples_per_chip = 4
        sampling_frequency = samples_per_chip * get_code_frequency(e5a_qp)
        code_frequency =
            get_code_frequency(e5a_qp) + doppler * get_code_center_frequency_ratio(e5a_qp)
        num_cycles = 60
        samples_per_cycle = BLOCKS_PER_CYCLE * samples_per_chip * get_code_length(e5a_qp)
        num_samples = num_cycles * samples_per_cycle

        code = gen_code(num_samples, e5a_qp, prn, sampling_frequency, code_frequency, 0.0)
        samples = ComplexF32.(
            cis.(
                2π .* ustrip(Hz, doppler) .* (0:(num_samples-1)) ./
                ustrip(Hz, sampling_frequency),
            ) .* code,
        )

        # Start 20 Hz off and feed one 2 ms cycle per `track` call.
        track_state = TrackState(e5a_qp, [TrackedSat(e5a_qp, prn, 0.0, 180.0Hz)])
        for cycle = 1:num_cycles
            track_state = track(
                view(samples, ((cycle-1)*samples_per_cycle+1):(cycle*samples_per_cycle)),
                track_state,
                sampling_frequency,
            )
        end
        sat_state = get_sat_state(track_state, prn)

        @test get_carrier_doppler(sat_state) ≈ doppler atol = 2.0Hz
        # The loop has rotated the prompt onto the real axis, and no energy was
        # lost to the 2 ms coherent integration.
        prompt = Tracking.get_last_fully_integrated_filtered_prompt(sat_state)
        @test real(prompt) / abs(prompt) > 0.9
        @test abs(prompt) > 0.9
        @test get_last_fully_integrated_num_code_blocks(sat_state) == BLOCKS_PER_CYCLE
        # Still a pilot: tracking observables, never a navigation bit.
        @test get_num_bits(sat_state) == 0
        @test estimate_cn0(track_state, prn) > 40.0dBHz
    end
end

end
