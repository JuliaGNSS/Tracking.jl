"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for generic BOC signal
"""
function is_upcoming_integration_new_bit(
    ::Type{<:GenericBOC{T}}, 
    prns,
    num_prns
) where T<:AbstractGNSSSystem
    is_upcoming_integration_new_bit(T,prns,num_prns)
end

function get_default_correlator(
    ::Type{<:GenericBOC{T}}, 
    numAnts::NumAnts{N}
) where {N, T<:AbstractGNSSSystem}
    EarlyPromptLateCorrelator(numAnts)
end
