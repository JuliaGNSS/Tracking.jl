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

# Import functions and types from Tracking that we need to extend or use
import Tracking: downconvert_and_correlate, get_num_samples, update
using Tracking: TrackedSat

# Include the GPU implementation
include("downconvert_and_correlate_gpu.jl")

# Export GPU types from the extension
# Users access these via: Base.get_extension(Tracking, :TrackingCUDAExt).GPUDownconvertAndCorrelator
export GPUSatDownconvertAndCorrelator,
       GPUSystemDownconvertAndCorrelator,
       GPUDownconvertAndCorrelator,
       convert_code_to_texture_memory

end
