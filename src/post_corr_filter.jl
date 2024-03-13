"""
$(SIGNATURES)

Abstract post correlation filter for the prompt signal.
"""
abstract type AbstractPostCorrFilter end

"""
$(SIGNATURES)

This is the default post correlation filter. For a single antenna
channel it will just return 
"""
struct DefaultPostCorrFilter <: AbstractPostCorrFilter end

update(filter::DefaultPostCorrFilter, prompt) = filter

(filter::DefaultPostCorrFilter)(x) = x
(filter::DefaultPostCorrFilter)(x::AbstractVector) = last(x)

"""
$(SIGNATURES)

Default post corr filter
"""
function get_default_post_corr_filter(correlator::AbstractCorrelator{T}) where {T<:SVector}
    x -> x[1]
end

function get_default_post_corr_filter(correlator)
    x -> x
end
