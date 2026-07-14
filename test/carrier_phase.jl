module CarrierPhaseTest

# Tests for the combined pilot+data carrier-quadrature fix: `_carrier_phase_derotation`
# de-rotates each component's bit-buffer prompt onto the driver's (real) phase
# frame, using the per-signal carrier phase from GNSSSignals' `get_carrier_phase_offset`.
#
# In a combined group the loops lock the driver (`signals[1]`) onto the real
# axis. A component in carrier quadrature with the driver (QPSK data/pilot pairs
# — GPS L5 I/Q, Galileo E5a I/Q) then lands on the imaginary axis, so its
# `real(prompt)` bit decision collapses. The de-rotation restores it.

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
    # GPS L5: driver = L5Q pilot (phase −π/2), data = L5I (phase 0). A pilot
    # locked to the real axis puts the data on the imaginary axis (+j·A for
    # amplitude A); the de-rotation must bring it back to the real axis
    # (magnitude preserved, sign resolved downstream by the decoder).
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

    # Single-signal (data-only) sat: driver == data ⇒ no-op, so the pre-fix
    # behaviour is preserved bit-for-bit.
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL5I()), GPSL5I()) ===
          complex(1.0, 0.0)
    @test _carrier_phase_derotation(get_carrier_phase_offset(GPSL1CA()), GPSL1CA()) ===
          complex(1.0, 0.0)

    # A no-op rotation leaves any prompt untouched (bit-identical).
    p = complex(1.25, -0.5)
    @test p * _carrier_phase_derotation(get_carrier_phase_offset(GPSL1CA()), GPSL1CA()) ===
          p

    # Symmetry: the de-rotation is relative to whichever signal drives the loop
    # (signals[1]), not to "the pilot". If the DATA component is the driver and
    # the pilot is the passenger (the reverse of GNSSReceiver's convention), the
    # loop locks the data onto the real axis, so:
    #   - the data driver's own bit prompt is a no-op (already real, decodes
    #     directly, exactly like data-only tracking), and
    #   - the pilot passenger is rotated onto real (harmless — a pilot carries no
    #     data, so its bit prompt never feeds a decoder).
    driver_offset_data = get_carrier_phase_offset(GPSL5I())         # data drives
    @test _carrier_phase_derotation(driver_offset_data, GPSL5I()) === complex(1.0, 0.0)  # driver no-op
    @test _carrier_phase_derotation(driver_offset_data, GPSL5Q()) ≈ im                   # pilot passenger
    # A pilot passenger sits at −90° when the data driver is locked to real (−jA);
    # the +im rotation brings it back to the real axis.
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
