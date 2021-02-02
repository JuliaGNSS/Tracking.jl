"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for generic BOC signal
"""
function is_upcoming_integration_new_bit(
    boc::BOCcos{T},
    prns,
    num_prns
) where T <: AbstractGNSS
    is_upcoming_integration_new_bit(boc, prns, num_prns)
end

function get_default_correlator(
    boc::BOCcos{T}, 
    numAnts::NumAnts{N}
) where {N, T <: AbstractGNSS}
    EarlyPromptLateCorrelator(numAnts)
end
