module CarrierPhaseTest

# `_carrier_phase_derotation` turns each component's bit prompt into the frame
# of the driver (`signals[1]`), which the loop locks onto the real axis. A
# component in quadrature with it (GPS L5 I/Q, Galileo E5a I/Q) would otherwise
# land on the imaginary axis and lose its `real(prompt)` bit decision.

using Test: @test, @testset
using TrackingLoops: _carrier_phase_derotation
using GNSSSignals:
    get_carrier_phase_offset,
    GPSL5I,
    GPSL5Q,
    GalileoE5aI,
    GalileoE5aQ,
    GalileoE1B,
    GalileoE1C,
    GPSL1CA

@testset "_carrier_phase_derotation" begin
    # GPS L5: the L5Q pilot (−π/2) on the real axis puts L5I data at +j·A.
    rot = _carrier_phase_derotation(get_carrier_phase_offset(GPSL5Q()), GPSL5I())
    @test rot ≈ -im
    data_prompt = complex(0.0, 5.0)   # +j·5
    @test real(data_prompt * rot) ≈ 5.0
    @test imag(data_prompt * rot) ≈ 0.0 atol = 1e-12

    # Galileo E5a: same quadrature relationship.
    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE5aQ()),
        GalileoE5aI(),
    ) ≈ im

    # Exactly 1 for the driver itself, a co-phased pair or a single signal.
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL5Q()), GPSL5Q()) ===
          complex(1.0, 0.0)
    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE5aQ()),
        GalileoE5aQ(),
    ) === complex(1.0, 0.0)

    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE1C()),
        GalileoE1B(),
    ) === complex(1.0, 0.0)

    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL5I()), GPSL5I()) ===
          complex(1.0, 0.0)
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL1CA()), GPSL1CA()) ===
          complex(1.0, 0.0)

    p = complex(1.25, -0.5)
    @test p * _carrier_phase_derotation(get_carrier_phase_offset(GPSL1CA()), GPSL1CA()) ===
          p

    # Relative to the driver, not to the pilot: with the data driving, its own
    # prompt is untouched and the pilot (at −j·A) is rotated onto the real axis.
    driver_offset_data = get_carrier_phase_offset(GPSL5I())         # data drives
    @test _carrier_phase_derotation(driver_offset_data, GPSL5I()) === complex(1.0, 0.0)  # driver no-op
    @test _carrier_phase_derotation(driver_offset_data, GPSL5Q()) ≈ im                   # pilot passenger
    pilot_prompt = complex(0.0, -5.0)   # −j·5
    @test real(pilot_prompt * _carrier_phase_derotation(driver_offset_data, GPSL5Q())) ≈ 5.0
    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE5aI()),
        GalileoE5aI(),
    ) === complex(1.0, 0.0)
    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE5aI()),
        GalileoE5aQ(),
    ) ≈ -im
end

end
