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
    system::AbstractGNSS,
    correlator::EarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    code_frequency = code_doppler + get_code_frequency(system)
    code_phase_delta = code_frequency / sampling_frequency
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    distance_between_early_and_late =
        get_early_late_sample_spacing(correlator, sampling_frequency, code_frequency) *
        code_phase_delta
    (2 - distance_between_early_and_late) / 2 * (E - L) / (E + L)
end

"""
$(SIGNATURES)

Calculates the code phase error in chips using the noncoherent very early minus late
envelope normalized discriminator for BOC signals.

Implementation from:
https://gnss-sdr.org/docs/sp-blocks/tracking/#implementation-galileo_e1_dll_pll_veml_tracking
"""
function dll_disc(
    system::AbstractGNSS,
    correlator::VeryEarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    VE = abs(get_very_early(correlator))
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    VL = abs(get_very_late(correlator))
    (VE + E - VL - L) / (VE + E + VL + L)
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(system::AbstractGNSS, correlator)
    p = get_prompt(correlator)
    atan(imag(p) / real(p))
end

"""
$(SIGNATURES)

Calculates the carrier frequency error in Hz.
"""
function fll_disc(system::AbstractGNSS, correlator, previous_prompt, integration_time)
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
