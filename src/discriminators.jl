"""
$(SIGNATURES)

Calculates the code phase error in chips.
"""
function dll_disc(
    system::AbstractGNSS,
    correlator::EarlyPromptLateCorrelator,
    code_doppler,
    sampling_frequency,
)
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    distance_between_early_and_late = get_early_late_code_spacing(correlator)
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
    code_doppler,
    sampling_frequency,
)
    VE = abs(get_very_early(correlator))
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    VL = abs(get_very_late(correlator))
    distance_between_early_and_late = get_early_late_code_spacing(correlator)
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
