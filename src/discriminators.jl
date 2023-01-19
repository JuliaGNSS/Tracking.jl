"""
$(SIGNATURES)

Calculates the code phase error in chips.
"""
function dll_disc(
    system::AbstractGNSS,
    correlator::EarlyPromptLateCorrelator,
    code_phase_delta
)
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    distance_between_early_and_late =
        get_early_late_sample_spacing(correlator) *
        code_phase_delta
    (E - L) / (E + L) / (2 * (2 - distance_between_early_and_late))
end

"""
$(SIGNATURES)
DLL discriminator for very early prompt late correlator.
Implementation from:
https://gnss-sdr.org/docs/sp-blocks/tracking/#implementation-galileo_e1_dll_pll_veml_tracking
"""
function dll_disc(
    system::AbstractGNSS,
    correlator::VeryEarlyPromptLateCorrelator,
    code_phase_delta
)
    VE = abs(get_very_early(correlator))
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    VL = abs(get_very_late(correlator))
    distance_between_early_and_late =
        get_early_late_sample_spacing(correlator) *
        code_phase_delta
    (VE + E - VL - L) / (VE + E + VL + L) / (2 * (2 - distance_between_early_and_late))
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(system::AbstractGNSS, correlator)
    p = get_prompt(correlator)
    atan(imag(p) / real(p))
end
