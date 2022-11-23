abstract type AbstractPostCorrFilter end

struct DefaultPostCorrFilter <: AbstractPostCorrFilter end

update(filter::DefaultPostCorrFilter, prompt) = filter

(filter::DefaultPostCorrFilter)(x) = x
(filter::DefaultPostCorrFilter)(x::AbstractVector) = last(x)