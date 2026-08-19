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
    code_phase_delta = upreferred(code_frequency / sampling_frequency)
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

# Tap offsets of a VEML correlator's inner (early/late) and outer
# (very-early/very-late) pairs, in chips, at the code frequency in effect — where the
# taps *actually* sit once the preferred chip shifts have been rounded to whole samples.
#
# Shared by `dll_disc` and `dll_disc_noise_gain` rather than computed in each, because
# the gain is the noise gain *of that discriminator* only if both are evaluated at the
# same offsets: the S-curve slope, and with it the chip calibration of the one and the
# variance of the other, moves with them.
@inline function _veml_tap_offsets(
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
    (inner_offset, outer_offset)
end

# S-curve slope of the very-early-minus-late discriminator at the origin, for the
# inner (early/late) and outer (very-early/very-late) tap offsets in chips: each tap
# pair's envelope difference contributes −2·slope·τ to the numerator and its envelope
# sum 2·|R| to the denominator. For the default ±0.15/±0.6 taps this is
# 4 / (2 − 3·0.15 − 0.6) ≈ 4.2.
#
# A tap layout that puts *both* pairs past the correlation support (≥ 1 chip) makes the
# envelope sum zero, and the quotient `-0.0 / 0.0` a `NaN` — which would then defeat the
# `iszero(slope)` guard in `dll_disc` and the `iszero(denominator)` one in
# `dll_disc_noise_gain`, since `NaN` compares equal to nothing. Report `0.0` instead: such a
# layout carries no delay information, which is what those two guards are there to say, and
# a `NaN` weight would otherwise propagate out of the combining accumulator into the
# *driver's* loop filter and wreck a satellite that was tracking fine.
function _veml_discriminator_slope(inner_offset, outer_offset)
    envelope_sum = _boc11_envelope(inner_offset) + _boc11_envelope(outer_offset)
    iszero(envelope_sum) && return 0.0
    -(_boc11_envelope_slope(inner_offset) + _boc11_envelope_slope(outer_offset)) /
    envelope_sum
end

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
    inner_offset, outer_offset =
        _veml_tap_offsets(signal, correlator, code_doppler, sampling_frequency)
    slope = _veml_discriminator_slope(inner_offset, outer_offset)
    VE = abs(get_very_early(correlator))
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    VL = abs(get_very_late(correlator))
    raw = (VE + E - VL - L) / (VE + E + VL + L)
    # A tap layout whose S-curve is locally flat cannot be calibrated — return the raw
    # discriminator rather than dividing by zero. Every functional layout (inner taps on
    # the main peak, outer taps on the side lobe) has a positive slope.
    iszero(slope) ? raw : raw / slope
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

Carrier phase error in radians, measured after rotating the prompt by
`rotation` — the form [multi-signal combining](@ref Multi-signal-discriminator-combining)
uses for a signal that is not the one the loops lock onto.

The PLL locks the estimator-driver signal onto the real axis, so a component
transmitted in carrier quadrature with the driver (GPS L5I vs L5Q, Galileo
E5aI vs E5aQ) sits on the imaginary axis and `atan(imag/real)` would read
`±π/2` instead of that component's share of the common phase error.
`rotation = `[`_carrier_phase_derotation`](@ref)`(driver_phase, signal)` brings
it back onto the driver's frame, where its prompt measures the *same* carrier
phase error the driver's does — the property that makes the two discriminators
combinable at all.

For a co-phased component (and for the driver itself) `rotation` is
`1 + 0im` and this is bit-identical to the two-argument method.

Being a Costas (π-ambiguous) discriminator is also why combining fuses
discriminators rather than prompts: it is blind to the unknown ±1
navigation-data sign, so a data-bearing component's discriminator can be averaged
with a pilot's where their *prompts* would partially cancel. See
[Multi-signal discriminator combining](@ref Multi-signal-discriminator-combining)
for the code loop's separate reason.
"""
function pll_disc(signal::AbstractGNSSSignal, correlator, rotation::Complex)
    p = get_prompt(correlator) * rotation
    atan(imag(p) / real(p))
end

# ---------------------------------------------------------------------------
# Discriminator noise gains — the `G` in `σ² ≈ G / SNR`
# ---------------------------------------------------------------------------
#
# Used to weight one signal's discriminator against another's when several
# signals of one satellite drive a common loop (see
# `_accumulate_passenger_discriminators`). Only *ratios* between the signals of a group
# matter, so every gain below is expressed against the same SNR definition —
# `SNR ≡ (received power) · (integration time) / N₀` — and every per-signal
# factor that is common across a group (the front-end noise density; the code
# amplitude, which `normalize` divides out of signal and noise alike) cancels and
# is never formed.
#
# Getting a gain wrong costs *combining efficiency*, never bias: each
# discriminator is individually calibrated (in chips / radians / Hz), so any
# positive weights produce a consistent combined estimate — the weights only
# decide how close to minimum-variance it lands. That is why the models below
# are allowed to be first-order.
#
# Only the *code* gain is a hook, because only it varies across a group: the
# correlator layout is per signal, and a BOC signal's VEML discriminator is ~8x
# more precise than a BPSK signal's early-late one at equal SNR. The carrier
# gains are hard-coded in `_record_contribution` instead — every shipped
# `pll_disc` is the same Costas `atan(Q/I)` with `sigma² ~ 1 / 2·SNR`, which
# cancels in the normalized mean, and every `fll_disc` is that phase difference
# over `2π·T`, whose only signal-dependent factor is `T` itself. A caller that
# overrides `pll_disc` or `fll_disc` for its own signal type with a
# *differently* precise discriminator — a four-quadrant `atan2` on a pilot, say,
# which is 3 dB better than the Costas form — therefore gets weighted as if it
# had not: correct, but not minimum-variance. Say so here rather than adding two
# more hooks nothing shipped would specialize.

"""
$(SIGNATURES)

Post-integration SNR of one correlator record, up to factors common to a
satellite's signals: the component's **nominal** power share
(`GNSSSignals.get_relative_power`) times the record's `integrated_samples`.

The power is the ICD split rather than the record's own `|P|²`, and it therefore
does *not* notice a component that is absent or unlocked — see
[Multi-signal discriminator combining](@ref Multi-signal-discriminator-combining)
for why that trade is the right way round, and what it obliges the caller to do.

`integrated_samples` rather than a time, because the sampling frequency is one of
the factors shared across a group.
"""
@inline _nominal_record_snr(signal::AbstractGNSSSignal, integrated_samples::Integer) =
    get_relative_power(signal) * integrated_samples

"""
$(SIGNATURES)

Fallback noise gain for a correlator type Tracking.jl ships no model for — a
custom multi-tap or multipath-mitigating correlator with its own
`dll_disc` method.

Returns the 1-chip early-late value `1/4`, i.e. such a discriminator is
weighted as if it had the accuracy of a wide-correlator BPSK early-late one.
That under-weights a sharper discriminator; define this method for your
correlator type to weight it on its merits. It cannot bias the result — see the
note above.
"""
@inline dll_disc_noise_gain(
    ::AbstractGNSSSignal,
    ::AbstractCorrelator,
    code_doppler,
    sampling_frequency,
) = 1 / 4

"""
$(SIGNATURES)

Noise gain `G` of the code discriminator for `correlator`, i.e. the constant in
`σ²_chips ≈ G / SNR`. Multi-signal combining weights each signal's
DLL discriminator by `SNR / G`, so this is what lets a BOC signal's VEML
discriminator outvote a BPSK signal's early-late one at equal C/N₀.

This is a statement about *variance*, and it is not made redundant by the
chip calibration every `dll_disc` method carries (the `(2 - d) / 2` factor here,
the S-curve slope division for VEML). That calibration equalizes the
discriminators' **mean** response — every method answers `1.0 · τ` in chips,
which is what makes two signals' discriminators averageable at all — by dividing
signal and noise alike by the slope. A raw discriminator that was steeper
therefore comes out of it with *less* noise per chip, not the same: `G` is
exactly what is left, and `G` is measured after the calibration, never a
substitute for it. Hence `d / 4 ≈ 0.25` for a 1-chip early-late layout against
`≈ 0.03` for the default VEML taps on a BOC(1,1) peak — an 8× precision
difference that survives both discriminators reading unit slope.

For the noncoherent early-minus-late envelope discriminator with early-late
spacing `d` chips this is `d / 4`, the Kaplan & Hegarty (2nd ed., Table 5.6)
tracking-jitter constant. The `d` dependence comes entirely from the *noise
correlation* between the early and late taps (`ρ ≈ 1 - d` for a triangular
autocorrelation): narrowing the taps shrinks the S-curve slope and the
differenced noise together, and the latter wins.
"""
@inline function dll_disc_noise_gain(
    signal::AbstractGNSSSignal,
    correlator::EarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    code_frequency = code_doppler + get_code_frequency(signal)
    code_phase_delta = upreferred(code_frequency / sampling_frequency)
    d =
        get_early_late_sample_spacing(correlator, sampling_frequency, code_frequency) *
        code_phase_delta
    d / 4
end

"""
$(SIGNATURES)

Noise gain `G` of the very-early-minus-late code discriminator, in the same
`σ²_chips ≈ G / SNR` convention as the early-late method above.

Evaluated on the same piecewise-linear sine-BOC(1,1) autocorrelation envelope
the discriminator's own chip calibration uses, at the same sample-quantized tap
offsets — literally the same, both taken from `_veml_tap_offsets`: with S-curve slope `k` (`_veml_discriminator_slope`) and
envelope values `e_i`, `e_o` at the inner and outer taps, the four taps
contribute `4 · σ²/2` of noise to a numerator whose denominator is
`2 · (e_i + e_o)`, giving `G = 1 / (2 · k² · (e_i + e_o)²)`.

Unlike the early-late gain this does **not** model the noise correlation
between the tap pairs, so it over-states VEML noise and under-weights a VEML
signal — the conservative direction, and only a matter of combining efficiency
(see the note above [`dll_disc_noise_gain`](@ref Tracking.dll_disc_noise_gain)). For the default ±0.15/±0.6
chip taps it lands near `0.03`, about 8× smaller than the 1-chip early-late
gain, matching the ~3× steeper BOC(1,1) main peak.
"""
@inline function dll_disc_noise_gain(
    signal::AbstractGNSSSignal,
    correlator::VeryEarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    inner_offset, outer_offset =
        _veml_tap_offsets(signal, correlator, code_doppler, sampling_frequency)
    slope = _veml_discriminator_slope(inner_offset, outer_offset)
    envelope_sum = _boc11_envelope(inner_offset) + _boc11_envelope(outer_offset)
    # A degenerate tap layout (flat S-curve, or both taps past the correlation
    # support) carries no delay information. Returning `Inf` gives it weight
    # zero instead of a division by zero — the discriminator is uninformative,
    # not infinitely precise.
    denominator = 2 * slope^2 * envelope_sum^2
    iszero(denominator) ? Inf : 1 / denominator
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
