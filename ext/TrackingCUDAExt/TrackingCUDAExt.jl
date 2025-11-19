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

# Make GPU types available in the parent Tracking module
# In package extensions, export doesn't propagate to the parent module,
# so we explicitly assign them to the parent module's namespace in __init__
function __init__()
    @eval Tracking begin
        GPUSatDownconvertAndCorrelator = $GPUSatDownconvertAndCorrelator
        GPUSystemDownconvertAndCorrelator = $GPUSystemDownconvertAndCorrelator
        GPUDownconvertAndCorrelator = $GPUDownconvertAndCorrelator
        convert_code_to_texture_memory = $convert_code_to_texture_memory
    end
end

end
