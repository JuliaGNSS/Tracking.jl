"""
$(SIGNATURES)

Calculates the code phase error in chips.
"""
function dll_disc(
    ::Type{S},
    correlator,
    early_late_sample_shift,
    code_phase_delta
) where S <: AbstractGNSSSystem
    E = abs(get_early(correlator))
    L = abs(get_late(correlator))
    distance_between_early_and_late = 2 * early_late_sample_shift * code_phase_delta
    (E - L) / (E + L) / (2 * (2 - distance_between_early_and_late))
end

"""
$(SIGNATURES)

Calculates the carrier phase error in radians.
"""
function pll_disc(::Type{S}, correlator) where S <: AbstractGNSSSystem
    p = get_prompt(correlator)
    atan(imag(p) / real(p))
end
