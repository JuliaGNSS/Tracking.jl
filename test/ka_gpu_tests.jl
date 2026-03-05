# KernelAbstractions GPU tests only — run on Buildkite with CUDA/AMDGPU agents.
# Usage: julia --project -e 'using Pkg; Pkg.test()'
# Or directly: julia --project test/ka_gpu_tests.jl

using Test

@testset "KA GPU Tests" begin
    include("downconvert_and_correlate.jl")
    include("track.jl")
end
