module Tracking
    using DocStringExtensions, GNSSSignals

    export dll_disc, pll_disc

    include("discriminators.jl")
    include("loop_filters.jl")
    include("dll_and_pll.jl")
end