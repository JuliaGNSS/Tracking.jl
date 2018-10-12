module Tracking
    using DocStringExtensions, GNSSSignals, DataStructures, LinearAlgebra, Statistics, StructArrays, StaticArrays, Destruct
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

    struct TrackingPhases
        carrier::Float64
        code::Float64
    end

    struct TrackingDopplers
        carrier::typeof(1.0Hz)
        code::typeof(1.0Hz)
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

    include("downconvert_and_correlate.jl")
    include("discriminators.jl")
    include("loop_filters.jl")
    include("tracking_loop.jl")
    include("joined_tracking.jl")
end
