module DiscriminatorCombiningTest

using Test: @test, @testset, @test_throws
using Unitful: Hz, s, ns, m
using GNSSSignals: GalileoE1B, GalileoE1C
using Tracking:
    Tracking, TrackedSignal, TrackState, add_satellite!, get_group_delay, set_group_delay!

@testset "group delay" begin
    @test isnothing(get_group_delay(TrackedSignal(GalileoE1B())))
    @test get_group_delay(TrackedSignal(GalileoE1B(); group_delay = 1.5ns)) ≈ 1.5e-9s

    ts = TrackState(; signals = (e1 = (GalileoE1C(), GalileoE1B()),))
    add_satellite!(ts; prn = 11, group = :e1, code_phase = 0.0, carrier_doppler = 0.0Hz)
    @test isnothing(get_group_delay(ts, :e1, 11, GalileoE1B))

    set_group_delay!(ts, :e1, 11, GalileoE1B, -0.5ns)
    @test get_group_delay(ts, :e1, 11, GalileoE1B) === -5.0e-10s
    @test isnothing(get_group_delay(ts, :e1, 11, GalileoE1C))

    set_group_delay!(ts, 11, 1, 0.0s)
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # `nothing` clears again: unknown, not zero.
    set_group_delay!(ts, 11, GalileoE1B, nothing)
    @test isnothing(get_group_delay(ts, :e1, 11, 2))
    @test get_group_delay(ts, :e1, 11, 1) === 0.0s

    # A group delay is a time.
    @test_throws Exception set_group_delay!(ts, 11, 2, 1.0m)
    @test_throws Exception set_group_delay!(ts, 11, 2, 1.0)
end

end
