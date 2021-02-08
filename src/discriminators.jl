"""
$(SIGNATURES)

Calculates the code phase error in chips.
"""
function dll_disc(
    system::AbstractGNSS,
    correlator,
    correlator_sample_shifts,
    early_late_index_shift,
    code_phase_delta
)
    E = abs(get_early(correlator, correlator_sample_shifts, early_late_index_shift))
    L = abs(get_late(correlator, correlator_sample_shifts, early_late_index_shift))
    distance_between_early_and_late =
        get_early_late_sample_spacing(correlator_sample_shifts, early_late_index_shift) *
        code_phase_delta
    (E - L) / (E + L) / (2 * (2 - distance_between_early_and_late))
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(system::AbstractGNSS, correlator, correlator_sample_shifts)
    p = get_prompt(correlator, correlator_sample_shifts)
    atan(imag(p) / real(p))
end
