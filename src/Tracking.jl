module Tracking
    using DocStringExtensions, GNSSSignals

    export dll_disc, pll_disc,init_1st_order_loop_filter,init_2nd_order_loop_filter,init_3rd_order_loop_filter, init_PLL, init_DLL, prompt

    include("discriminators.jl")
    include("loop_filters.jl")
    include("dll_and_pll.jl")
    include("tracking_loop.jl")
end