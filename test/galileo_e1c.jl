module GalileoE1CTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GalileoE1C, GalileoE1C_BOC11, get_secondary_code_length
import Tracking
using Tracking:
    detect_bit_or_secondary_code_sync,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    VeryEarlyPromptLateCorrelator,
    NumAnts

rotl(x::T, r, N) where {T} =
    r == 0 ? x : ((x << r) | (x >> (N - r))) & ((one(T) << N) - one(T))

# `GalileoE1C` (full CBOC) and `GalileoE1C_BOC11` (BOC(1,1) approximation)
# share primary code, code length, code rate, CS25 secondary code, and band
# — only the modulation differs — so all tracking-side traits behave
# identically. Parameterise over both.
@testset "Galileo E1C ($(nameof(typeof(galileo_e1c))))" for galileo_e1c in (
    GalileoE1C(),
    GalileoE1C_BOC11(),
)
    prn = 1
    N = get_secondary_code_length(galileo_e1c)  # 25 (CS25)
    @test N == 25

    # E1C is the E1 pilot: no data, a 25-chip CS25 secondary code is the sync
    # feature. Below one full period the detector returns `found = false`.
    @test @inferred(
        detect_bit_or_secondary_code_sync(galileo_e1c, prn, UInt32(0x0), N - 1)
    ).found == false

    @testset "CS25 search — clean lock at known phase / polarity" begin
        reference = Tracking._packed_secondary_code(UInt32, galileo_e1c, prn)
        for r in (0, 11, N - 1)
            received = rotl(reference, r, N)
            res = @inferred detect_bit_or_secondary_code_sync(galileo_e1c, prn, received, N)
            @test res.found == true
            @test res.phase == r
            @test res.polarity == +1
        end
        negated = reference ⊻ ((one(UInt32) << N) - one(UInt32))
        res = @inferred detect_bit_or_secondary_code_sync(galileo_e1c, prn, negated, N)
        @test res.found == true
        @test res.phase == 0
        @test res.polarity == -1
    end

    # E1C shares the E1 modulation family with E1B → VeryEarlyPromptLate default.
    @test @inferred(get_default_correlator(galileo_e1c, NumAnts(1))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(galileo_e1c, NumAnts(3))) ==
          VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # 4 ms primary period (4092 chips at 1.023 Mcps), same as E1B → 4.5 Hz / 0.25 Hz.
    @test @inferred(default_carrier_loop_filter_bandwidth(galileo_e1c)) ≈ 4.5Hz
    @test @inferred(default_code_loop_filter_bandwidth(galileo_e1c)) ≈ 0.25Hz

    # 25-block CS25 window fits in a UInt32.
    @test @inferred(get_code_block_buffer_type(galileo_e1c)) === UInt32
end

end
