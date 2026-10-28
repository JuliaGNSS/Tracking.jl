module BitBufferTest

# The warning Tracking logs when a correlator record carries the bit accumulator past
# the navigation-bit boundary. TrackingLoops' `apply_record` reports the overshoot
# and its bit buffer drops sync (tested there); the record path through `track!` is in
# `test/external_correlator_producer.jl`.

using Test: @test, @testset, @test_logs
using GNSSSignals: GPSL1CA, get_signal_id
using Tracking: _warn_bit_boundary_overshoot

@testset "The overshoot is reported to the caller, which warns" begin
    # What a user sees when external records cross a bit boundary. A PRN of its own:
    # the warning is `maxlog = 1` per signal and PRN.
    @test_logs (:warn, r"past the navigation-bit boundary") match_mode = :any _warn_bit_boundary_overshoot(
        get_signal_id(GPSL1CA()),
        99,
        21,
        20,
    )
end

end
