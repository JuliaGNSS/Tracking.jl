# Shared helper for building an `AcquisitionResults` (Acquisition 3: 15
# positional fields, with `secondary_code_phase` after `code_phase` and
# `num_secondary_rotations` last).
#
# The including module imports `Acquisition` and `AcquisitionResults`.
function _make_acq(signal, prn, code_phase, carrier_doppler; noise_power = 1.0)
    AcquisitionResults(
        signal,
        prn,
        5e6Hz,
        carrier_doppler,
        code_phase,
        nothing,
        45.0,
        noise_power,
        10.0,
        1,
        randn(10, 10),
        -500:100.0:500,
        1,
        length(-500:100.0:500),
        1,
    )
end
