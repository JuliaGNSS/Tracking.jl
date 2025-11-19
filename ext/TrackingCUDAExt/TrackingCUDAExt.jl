module TrackingCUDAExt

using Tracking
using CUDA

# Import dependencies needed by GPU code
using GNSSSignals
using StaticArrays
using Dictionaries
using Accessors
using DocStringExtensions
using Unitful: Hz
import Adapt

# Include the GPU implementation
include("downconvert_and_correlate_gpu.jl")

# Export GPU types - accessible as Tracking.TrackingCUDAExt.GPUDownconvertAndCorrelator
export GPUSatDownconvertAndCorrelator,
       GPUSystemDownconvertAndCorrelator,
       GPUDownconvertAndCorrelator,
       convert_code_to_texture_memory

end
