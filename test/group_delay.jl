module GroupDelayTest

using Test: @test, @testset
using Unitful: Hz, ns, s
using GNSSSignals: GPSL1CA, GalileoE1B, GalileoE1C
import Tracking
using Tracking:
    SignalGroup,
    TrackState,
    TrackedSignal,
    add_satellite!,
    get_group_delay,
    set_group_delay!

_e1_state() = TrackState(; signals = (e1 = (GalileoE1C(), GalileoE1B()),))

@testset "Group delays" begin
    @test isnothing(TrackedSignal(GalileoE1B()).group_delay)
    delayed = TrackedSignal(GalileoE1B(); group_delay = 2.5ns)
    @test get_group_delay(delayed) === 2.5e-9 * 1.0s
    # Any time unit and number type is converted to seconds in `Float64`.
    @test get_group_delay(TrackedSignal(GalileoE1B(); group_delay = 10ns)) ===
          convert(typeof(1.0s), 10ns)
    # The kwarg-update constructor keeps it on `nothing` and takes a new value
    # (`nothing` included) as `Some`.
    @test get_group_delay(TrackedSignal(delayed; integrated_samples = 3)) === 2.5e-9s
    @test isnothing(get_group_delay(TrackedSignal(delayed; group_delay = Some(nothing))))

    track_state = _e1_state()
    add_satellite!(
        track_state;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    set_group_delay!(track_state, :e1, 11, GalileoE1B, 3.0ns)
    @test get_group_delay(track_state, :e1, 11, GalileoE1B) ≈ 3.0e-9s
    @test isnothing(get_group_delay(track_state, :e1, 11, GalileoE1C))
    set_group_delay!(track_state, :e1, 11, 1, 1.0ns)
    @test get_group_delay(track_state, :e1, 11, GalileoE1C) ≈ 1.0e-9s
    set_group_delay!(track_state, :e1, 11, GalileoE1B, nothing)
    @test isnothing(get_group_delay(track_state, :e1, 11, 2))

    single = TrackState(; signal = GPSL1CA())
    add_satellite!(single; prn = 1, code_phase = 0.0, carrier_doppler = 0.0Hz)
    set_group_delay!(single, 1, 4.0ns)
    @test get_group_delay(single, 1) ≈ 4.0e-9s
    set_group_delay!(single, 5.0ns)
    @test get_group_delay(single) ≈ 5.0e-9s
end

end
