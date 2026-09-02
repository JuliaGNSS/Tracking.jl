"""
$(SIGNATURES)

Calculates the code phase error in chips using the noncoherent early minus late
envelope normalized discriminator.

Uses the generalized normalization for arbitrary early-late spacing `d` (in chips):
`(2 - d) / 2 * (E - L) / (E + L)`
which reduces to `1/2 * (E - L) / (E + L)` for the standard 1-chip spacing.

See: Kaplan & Hegarty, "Understanding GPS: Principles and Applications", 2nd ed.,
Table 5.5; GNSS-SDR tracking_discriminators.cc.
"""
function dll_disc(
    signal::AbstractGNSSSignal,
    correlator::EarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    code_frequency = code_doppler + get_code_frequency(signal)
    code_phase_delta = code_frequency / sampling_frequency
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    distance_between_early_and_late =
        get_early_late_sample_spacing(correlator, sampling_frequency, code_frequency) *
        code_phase_delta
    (2 - distance_between_early_and_late) / 2 * (E - L) / (E + L)
end

# Piecewise-linear sine-BOC(1,1) autocorrelation envelope |R(τ)| (infinite bandwidth,
# long-code approximation): the main peak falls from 1 with slope −3 to the zero crossing
# at 1/3 chips, the (negative) side lobe's envelope rises with slope +3 to 1/2 at
# 0.5 chips and falls with slope −1 to zero at 1 chip.
_boc11_envelope(offset) =
    offset < 1 / 3 ? 1 - 3 * offset :
    offset < 1 / 2 ? 3 * offset - 1 : offset < 1 ? 1 - offset : 0.0

# Envelope derivative as the mean of the one-sided derivatives, so a tap sitting exactly
# on an envelope knot (e.g. the side-lobe vertex at 0.5 chips under coarse sampling, where
# the ± excursions of the tap pair run along different segments) gets the slope those
# excursions actually produce.
function _boc11_envelope_slope(offset)
    right = offset < 1 / 3 ? -3.0 : offset < 1 / 2 ? 3.0 : offset < 1 ? -1.0 : 0.0
    left = offset <= 1 / 3 ? -3.0 : offset <= 1 / 2 ? 3.0 : offset <= 1 ? -1.0 : 0.0
    (left + right) / 2
end

# S-curve slope of the very-early-minus-late discriminator at the origin, for the
# inner (early/late) and outer (very-early/very-late) tap offsets in chips: each tap
# pair's envelope difference contributes −2·slope·τ to the numerator and its envelope
# sum 2·|R| to the denominator. For the default ±0.15/±0.6 taps this is
# 4 / (2 − 3·0.15 − 0.6) ≈ 4.2.
_veml_discriminator_slope(inner_offset, outer_offset) =
    -(_boc11_envelope_slope(inner_offset) + _boc11_envelope_slope(outer_offset)) /
    (_boc11_envelope(inner_offset) + _boc11_envelope(outer_offset))

"""
$(SIGNATURES)

Calculates the code phase error in chips using the noncoherent very early minus late
envelope normalized discriminator `(VE + E - VL - L) / (VE + E + VL + L)` for
BOC(1,1)-dominant signals (Galileo E1, GPS L1C), divided by its S-curve slope so the
output is calibrated in chips — the contract the early-prompt-late method above keeps
with its `(2 - d) / 2` factor. The slope (≈ 4.3 for the default ±0.15/±0.6 chip taps)
is evaluated on the piecewise-linear sine-BOC(1,1) autocorrelation envelope at the
actual sample-quantized tap offsets, the calibration GNSS-SDR applies where it wants
chips from a BOC discriminator (`CalculateSlopeAbs` on `SinBocCorrelationFunction`).
Against the full CBOC/TMBOC modulations a residual gain error of the modulation
mismatch remains, and the discriminator is linear only within the inner tap offset.

Raw discriminator form from:
https://gnss-sdr.org/docs/sp-blocks/tracking/#implementation-galileo_e1_dll_pll_veml_tracking
"""
function dll_disc(
    signal::AbstractGNSSSignal,
    correlator::VeryEarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    code_frequency = code_doppler + get_code_frequency(signal)
    code_phase_delta = upreferred(code_frequency / sampling_frequency)
    inner_offset =
        calc_preferred_code_shift_to_sample_shift(
            correlator.preferred_early_late_to_prompt_code_shift,
            sampling_frequency,
            code_frequency,
        ) * code_phase_delta
    outer_offset =
        calc_preferred_code_shift_to_sample_shift(
            correlator.preferred_very_early_late_to_prompt_code_shift,
            sampling_frequency,
            code_frequency,
        ) * code_phase_delta
    slope = _veml_discriminator_slope(inner_offset, outer_offset)
    VE = abs(get_very_early(correlator))
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    VL = abs(get_very_late(correlator))
    raw = (VE + E - VL - L) / (VE + E + VL + L)
    # A tap layout whose S-curve is locally flat cannot be calibrated — return the raw
    # discriminator rather than dividing by zero. Every functional layout (inner taps on
    # the main peak, outer taps on the side lobe) has a positive slope.
    slope == 0 ? raw : raw / slope
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(signal::AbstractGNSSSignal, correlator)
    p = get_prompt(correlator)
    atan(imag(p) / real(p))
end

"""
$(SIGNATURES)

Calculates the carrier frequency error in Hz.
"""
function fll_disc(signal::AbstractGNSSSignal, correlator, previous_prompt, integration_time)
    if previous_prompt == 0
        # return 0 when there is no previous prompt
        return 0.0/integration_time
    end

    current_prompt = get_prompt(correlator)

    result = conj(previous_prompt) * current_prompt
    cross = imag(result)
    dot = real(result)

    # atan(+-Int) produces valid outputs (+-π / 2)
    return atan(cross / dot) / (2 * pi * integration_time)
end
