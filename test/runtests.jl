using
    Test,
    Tracking,
    GNSSSignals,
    Random,
    StaticArrays,
    TrackingLoopFilters,
    StructArrays,
    Statistics

import Unitful: MHz, kHz, Hz, s, ms, dBHz

include("agc.jl")
include("code_replica.jl")
include("carrier_replica.jl")
include("downconvert.jl")
include("discriminators.jl")
include("gps_l1.jl")
include("gps_l5.jl")
include("galileo_e1b.jl")
include("secondary_code_or_bit_detector.jl")
include("bit_buffer.jl")
include("correlator.jl")
include("tracking_state.jl")
include("tracking_results.jl")
include("tracking_loop.jl")
include("cn0_estimation.jl")
