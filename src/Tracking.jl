module Tracking
    using DocStringExtensions, GNSSSignals, DataStructures, LinearAlgebra
    import Unitful: Hz, s

    export prompt,
        init_tracking,
        Initials

    struct TrackingResults
        carrier_doppler::typeof(1.0Hz)
        carrier_phase::Float64
        code_doppler::typeof(1.0Hz)
        code_phase::Float64
        prompt::Array{Complex{Float64}, 1}
    end

    struct Initials
        carrier_doppler::typeof(1.0Hz)
        carrier_phase::Float64
        code_doppler::typeof(1.0Hz)
        code_phase::Float64
    end

    function Initials(res::TrackingResults)
        Initials(res.carrier_doppler, res.carrier_phase, res.code_doppler, res.code_phase)
    end

    include("discriminators.jl")
    include("loop_filters.jl")
    include("replicas.jl")
    include("tracking_loop.jl")
    include("joined_tracking.jl")
end
