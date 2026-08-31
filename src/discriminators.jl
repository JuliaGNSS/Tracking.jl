"""
$(SIGNATURES)

Normalized envelope ratio `numerator / denominator`, or zero when the correlator
carries no energy at all.

Every discriminator here normalizes by a sum of envelopes, so a correlator whose
accumulators are all zero makes it compute `0 / 0`. That happens in practice: a
satellite that has just been added has not integrated anything yet, and a
hardware correlator (GNSSReceiver.jl#107) can deliver a record with
`integrated_samples == 0` after a dropped dump or a channel reassignment.

The resulting `NaN` does not stay local. It flows into the loop filter and out
as `carrier_doppler`/`code_doppler`, and the next iteration of the correlate
loop converts a Doppler-derived sample count to an integer — so the receiver
dies with `InexactError: Int64(NaN)` from inside the tracking task, far from the
zero correlator that caused it. Zero is the honest discriminator output for no
energy: no measurement, no correction.
"""
@inline function normalized_envelope_ratio(numerator, denominator)
    iszero(denominator) ? zero(numerator) / oneunit(denominator) : numerator / denominator
end

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
    # Note the numerator carries the scale factor, so the arithmetic keeps the
    # original `(2 - d) / 2 * (E - L) / (E + L)` left-to-right association and
    # stays bit-identical to it whenever there is energy to normalise by.
    normalized_envelope_ratio((2 - distance_between_early_and_late) / 2 * (E - L), E + L)
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
    raw = normalized_envelope_ratio(VE + E - VL - L, VE + E + VL + L)
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
    # A zero prompt has no phase to measure, and `atan(0 / 0)` is NaN — see
    # `normalized_envelope_ratio` for why that NaN takes the whole tracking task
    # down rather than staying in this satellite's loop. The guard is on the
    # prompt as a whole, not on `real(p)`: a purely imaginary prompt is a
    # legitimate ±π/2, which `atan(±Inf)` already returns.
    iszero(p) && return atan(zero(real(p)))
    atan(imag(p) / real(p))
end

"""
$(SIGNATURES)

Calculates the carrier frequency error in Hz.
"""
function fll_disc(signal::AbstractGNSSSignal, correlator, previous_prompt, integration_time)
    current_prompt = get_prompt(correlator)

    # No previous prompt means no frequency to difference. A zero *current*
    # prompt is the same situation one epoch later and needs the same answer:
    # `result` would be 0, so `atan(0 / 0)` is NaN, and that NaN reaches
    # `carrier_doppler` and kills the tracking task (see
    # `normalized_envelope_ratio`).
    if iszero(previous_prompt) || iszero(current_prompt)
        return 0.0 / integration_time
    end

    result = conj(previous_prompt) * current_prompt
    cross = imag(result)
    dot = real(result)

    # atan(+-Int) produces valid outputs (+-π / 2)
    return atan(cross / dot) / (2 * pi * integration_time)
end
