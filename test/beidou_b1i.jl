module BeiDouB1ITest

using Test: @test, @testset, @inferred
using Unitful: Hz
using Random: MersenneTwister
using GNSSSignals: BeiDouB1I, get_secondary_code, get_secondary_code_length, secondary_value
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

rotl(x::T, r, N) where {T} =
    r == 0 ? x : ((x << r) | (x >> (N - r))) & ((one(T) << N) - one(T))

# PRN 6 is a MEO/IGSO satellite (D1 message) and carries the NH20 overlay;
# PRN 1 is GEO (D2) and carries none.
const MEO_PRN = 6
const GEO_PRN = 1

# Block index of the first symbol boundary in the synthetic streams, i.e. the
# rotation the detector should recover.
const SYMBOL_OFFSET = 7

# Drive the live soft path (the pair `_buffer_find_bit` calls) with a synthetic
# prompt stream: one symbol every `blocks_per_symbol` blocks times the overlay
# chip, plus unit-variance complex noise. Returns the 1-based lock block and
# rotation, or `nothing`.
function soft_sync(prn, blocks_per_symbol; nblocks, amplitude = 3.0, seed = 1)
    signal = BeiDouB1I()
    N = get_secondary_code_length(signal)
    overlay = get_secondary_code(signal)
    rng = MersenneTwister(seed)
    accumulators = Tracking.PhaseAccumulators()
    Tracking._seed_phase_accumulators!(accumulators, N)
    symbol = 1.0
    for i = 0:(nblocks-1)
        (i - SYMBOL_OFFSET) % blocks_per_symbol == 0 && (symbol = rand(rng, (-1.0, 1.0)))
        chip = secondary_value(overlay, prn, mod(i - SYMBOL_OFFSET, N))
        prompt = amplitude * symbol * chip + randn(rng, ComplexF64)
        Tracking._update_secondary_accumulators!(
            accumulators,
            ComplexF64(prompt),
            i,
            N,
            signal,
            prn,
        )
        result = Tracking._detect_secondary_code_cfar(accumulators, N, 0.999, i + 1)
        result.found && return (block = i + 1, rotation = mod(i + 1, N))
    end
    nothing
end

@testset "BeiDou B1I" begin
    b1i = BeiDouB1I()
    N = get_secondary_code_length(b1i)  # 20 (NH20)
    @test N == 20

    # One NH20 period is one D1 symbol. No lock below one full period.
    @test @inferred(detect_bit_or_secondary_code_sync(b1i, MEO_PRN, UInt32(0x0), N - 1)).found ==
          false

    @testset "NH20 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt32, b1i, MEO_PRN)
        # NH20 per BDS-SIS-ICD-B1I-3.0 §5.2.1, packed newest-first.
        @test reference == UInt32(0b00000100110101001110)
        for r in (0, 9, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(b1i, MEO_PRN, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(b1i, MEO_PRN, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    @testset "GEO satellites carry no NH20 overlay" begin
        # GEO PRNs get an all-ones column (see docs/src/bit_sync.md).
        mask = (one(UInt32) << N) - one(UInt32)
        for prn in (1, 5, 59, 63)
            @test Tracking._packed_secondary_code(UInt32, b1i, prn) == mask
        end
        # Rotation-invariant, so the hard sweep would lock at an arbitrary
        # phase. Pinned so a move to the hard path cannot pass silently.
        reference = Tracking._packed_secondary_code(UInt32, b1i, GEO_PRN)
        @test all(rotl(reference, r, N) == reference for r = 0:(N-1))
    end

    @testset "Soft sync — the live path, and what a GEO satellite does to it" begin
        # MEO/IGSO: NH20 at the D1 rate locks at the true period boundary.
        meo = soft_sync(MEO_PRN, 20; nblocks = 2000)
        @test meo !== nothing
        @test meo.rotation == SYMBOL_OFFSET

        # Hypothetical all-ones column at the D1 rate: reduces to a bit-edge
        # search and still finds the symbol boundary.
        geo_at_d1_rate = soft_sync(GEO_PRN, 20; nblocks = 4000)
        @test geo_at_d1_rate !== nothing
        @test geo_at_d1_rate.rotation == SYMBOL_OFFSET

        # Real GEO: no overlay and D2 at 2 blocks/symbol never syncs. Pins the
        # claim in docs/src/bit_sync.md.
        for seed = 1:3
            @test soft_sync(GEO_PRN, 2; nblocks = 20000, seed) === nothing
        end
    end

    # Soft secondary-code detection for every PRN; never the bit-edge detector.
    @test Tracking.uses_soft_secondary_code_detection(b1i) == true
    @test Tracking.uses_soft_bit_edge_detection(b1i) == false

    # BPSK → EarlyPromptLate.
    @test @inferred(get_default_correlator(b1i, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(b1i, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Flat defaults, uncapped at 1 ms.
    @test @inferred(default_carrier_loop_filter_bandwidth(b1i)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(b1i)) ≈ 1.0Hz

    # 20-block NH20 window fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(b1i)) === UInt32
end

end
