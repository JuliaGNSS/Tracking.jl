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

# Re-export GPU functionality to make it available when CUDA is loaded
export GPUSatDownconvertAndCorrelator,
       GPUSystemDownconvertAndCorrelator,
       GPUDownconvertAndCorrelator,
       convert_code_to_texture_memory

# Include the GPU implementation
include("downconvert_and_correlate_gpu.jl")

end
