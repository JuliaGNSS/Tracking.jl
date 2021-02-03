"""
$(SIGNATURES)

Calculates the code phase error in chips.
"""
function dll_disc(
    system::AbstractGNSS,
    correlator,
    correlator_sample_shifts,
    code_phase_delta
)
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    distance_between_early_and_late = get_early_late_sample_spacing(correlator, correlator_sample_shifts) * code_phase_delta
    (E - L) / (E + L) / (2 * (2 - distance_between_early_and_late))
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(system::AbstractGNSS, correlator)
    p = get_prompt(correlator)
    atan(imag(p) / real(p))
end
