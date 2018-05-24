module Tracking
    using DocStringExtensions

    export dll_disc, pll_disc, init_1st_order_loop_filter, init_2nd_order_loop_filter, init_3rd_order_loop_filter
    include("discriminators.jl")
    include("loop_filters.jl")
end