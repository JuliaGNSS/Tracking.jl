module CarrierPhaseTest

# `_carrier_phase_derotation` rotates a component's bit prompt onto the driver's real
# axis, so quadrature data/pilot pairs still decode; see its comment in
# conventional_pll_and_dll.jl.

using Test: @test, @testset
using Tracking: _carrier_phase_derotation
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
    # GPS L5 pilot driver: the data sits at +j·A and must come back to the real axis
    # with its magnitude preserved.
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

    # The driver de-rotates against itself ⇒ exact no-op (cis(0) === 1 + 0im).
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL5Q()), GPSL5Q()) ===
          complex(1.0, 0.0)
    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE5aQ()),
        GalileoE5aQ(),
    ) === complex(1.0, 0.0)

    # Co-phased CBOC pair (Galileo E1C driver, E1B data): both in-phase ⇒ no-op.
    @test _carrier_phase_derotation(
        get_carrier_phase_offset(GalileoE1C()),
        GalileoE1B(),
    ) === complex(1.0, 0.0)

    # Single-signal (data-only) sat: driver == data ⇒ no-op.
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL5I()), GPSL5I()) ===
          complex(1.0, 0.0)
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL1CA()), GPSL1CA()) ===
          complex(1.0, 0.0)

    # A no-op rotation leaves any prompt untouched (bit-identical).
    p = complex(1.25, -0.5)
    @test p * _carrier_phase_derotation(get_carrier_phase_offset(GPSL1CA()), GPSL1CA()) ===
          p

    # The rotation is relative to the driver (`signals[1]`), not to "the pilot": with a
    # data driver, the data is a no-op and the pilot passenger is rotated (harmless).
    driver_offset_data = get_carrier_phase_offset(GPSL5I())         # data drives
    @test _carrier_phase_derotation(driver_offset_data, GPSL5I()) === complex(1.0, 0.0)  # driver no-op
    @test _carrier_phase_derotation(driver_offset_data, GPSL5Q()) ≈ im                   # pilot passenger
    pilot_prompt = complex(0.0, -5.0)   # −j·5
    @test real(pilot_prompt * _carrier_phase_derotation(driver_offset_data, GPSL5Q())) ≈ 5.0
    # Same for Galileo E5a with the data component as driver.
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
