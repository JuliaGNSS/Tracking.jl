module Tracking
    using DocStringExtensions

    export dll_disc, pll_disc

    include("discriminators.jl")
    include("loop_filters.jl")
end