module Tracking
    using DocStringExtensions, GNSSSignals

    export dll_disc,
        pll_disc,
        prompt,
        init_1st_order_loop_filter,
        init_2nd_order_loop_filter,
        init_3rd_order_loop_filter,
        init_code_replica,
        init_carrier_replica,
        init_tracking

    include("discriminators.jl")
    include("loop_filters.jl")
    include("replicas.jl")
    include("tracking_loop.jl")
    include("joined_tracking.jl")
end