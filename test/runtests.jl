using CUDA
using GNSSSignals
using Random
using StaticArrays
using Statistics
using StructArrays
using Test
using Tracking
using TrackingLoopFilters

using Acquisition: AcquisitionResults
using Unitful: MHz, kHz, Hz, s, ms, dBHz

include("conventional_pll_and_dll.jl")
include("sat_state.jl")
include("sample_parameters.jl")
include("downconvert_and_correlate.jl")
include("post_corr_filter.jl")
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
include("cn0_estimation.jl")
include("track.jl")
include("tracking_state.jl")

if CUDA.functional()
    include("cuda/bit_buffer.jl")
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
