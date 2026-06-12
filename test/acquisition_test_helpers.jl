# Shared helper for building an `AcquisitionResults` across the Acquisition
# versions in test/Project.toml's compat range:
#   * v1            — 11 positional fields.
#   * v2.0 – v2.4   — added num_blocks/block_size (FM-DBZP backend) → 13 fields.
#   * v2.5+         — added secondary_code_phase (after code_phase) and
#                     num_secondary_rotations (last) → 15 fields.
#
# `include`d by each test module that hands acquisition results to the
# tracking API (test/add_satellite.jl, test/sat_state.jl), so the version
# shim is maintained in exactly one place. Expects `Acquisition` and
# `AcquisitionResults` to be imported by the including module.
function _make_acq(signal, prn, code_phase, carrier_doppler; noise_power = 1.0)
    args = (
        signal, prn, 5e6Hz, carrier_doppler, code_phase,
        45.0, noise_power, 10.0, 1, randn(10, 10), -500:100.0:500,
    )
    if pkgversion(Acquisition) >= v"2.5"
        AcquisitionResults(
            args[1:5]..., nothing, args[6:end]...,
            1, length(-500:100.0:500), 1,
        )
    elseif pkgversion(Acquisition) >= v"2"
        AcquisitionResults(args..., 1, length(-500:100.0:500))
    else
        AcquisitionResults(args...)
    end
end
