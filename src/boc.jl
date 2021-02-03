"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for generic BOC signal
"""
function is_upcoming_integration_new_bit(
    boc::BOCcos,
    prns,
    num_prns
)
    is_upcoming_integration_new_bit(boc.system, prns, num_prns)
end

function get_default_correlator(
    ::BOCcos,
    numAnts::NumAnts
)
    EarlyPromptLateCorrelator(numAnts)
end
