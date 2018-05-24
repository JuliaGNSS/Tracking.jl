module Tracking
    using DocStringExtensions

    export dll_disc, pll_disc

    include("discriminators.jl")
    include("loop_filters.jl")
    include("dll_and_pll.jl")
end