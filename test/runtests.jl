using
    Test,
    Tracking,
    GNSSSignals,
    Random,
    StaticArrays,
    TrackingLoopFilters,
    StructArrays,
    Statistics,
    CUDA

using Acquisition: AcquisitionResults
using Unitful: MHz, kHz, Hz, s, ms, dBHz

include("code_replica.jl")
include("carrier_replica.jl")
include("downconvert.jl")
include("discriminators.jl")
include("gps_l1.jl")
include("gps_l5.jl")
include("galileo_e1b.jl")
include("boc.jl")
include("secondary_code_or_bit_detector.jl")
include("bit_buffer.jl")
include("correlator.jl")
include("tracking_state.jl")
include("tracking_results.jl")
include("tracking_loop.jl")
include("cn0_estimation.jl")

if CUDA.functional()
    include("cuda/bit_buffer.jl")
    include("cuda/boc.jl")
    include("cuda/cn0_estimation.jl")
    include("cuda/discriminators.jl")
    include("cuda/galileo_e1b.jl")
    include("cuda/gps_l1.jl")
    include("cuda/gps_l5.jl")
    include("cuda/secondary_code_or_bit_detector.jl")
    include("cuda/tracking_loop.jl")
    include("cuda/tracking_results.jl")
    include("cuda/tracking_state.jl")
end
