module GPSL1CATest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA
using Tracking:
    is_upcoming_integration_new_bit,
    get_default_correlator,
    get_code_block_buffer_type,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    EarlyPromptLateCorrelator,
    NumAnts

@testset "GPS L1" begin
    gpsl1 = GPSL1CA()
    # PRN argument is ignored by the L1 C/A detector but required by the
    # uniform `is_upcoming_integration_new_bit(signal, prn, ...)` API.
    prn = 1

    # Matched at negative polarity (20 ones followed by 20 zeros).
    res = @inferred(is_upcoming_integration_new_bit(gpsl1, prn, UInt64(0xfffff00000), 50))
    @test res.found == true
    @test res.polarity == -1

    # Not enough integrations yet — buffer hasn't filled to 40 blocks.
    @test @inferred(is_upcoming_integration_new_bit(gpsl1, prn, UInt64(0xfffff), 10)).found == false

    # Matched at positive polarity (20 zeros followed by 20 ones).
    res = @inferred(is_upcoming_integration_new_bit(gpsl1, prn, UInt64(0xfffff), 40))
    @test res.found == true
    @test res.polarity == +1

    sampling_frequency = 5e6Hz

    @test @inferred(get_default_correlator(gpsl1, NumAnts(1))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    @test @inferred(get_default_correlator(gpsl1, NumAnts(3))) ==
          EarlyPromptLateCorrelator(; num_ants = NumAnts(3))

    # Per-signal default loop bandwidths: sized for the 1 ms primary code
    # period at BL·T ≈ 0.018. Historically the package shipped 18 Hz / 1 Hz
    # as the universal default — that value falls out of the formula here.
    @test @inferred(default_carrier_loop_filter_bandwidth(gpsl1)) ≈ 18.0Hz
    @test @inferred(default_code_loop_filter_bandwidth(gpsl1)) ≈ 1.0Hz

    # 40-block sync window fits in a UInt64.
    @test @inferred(get_code_block_buffer_type(gpsl1)) === UInt64
end

end
