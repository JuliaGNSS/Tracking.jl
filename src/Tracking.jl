module Tracking
    using DocStringExtensions, GNSSSignals

    export prompt,
        init_tracking,
        Initials

    struct TrackingResults
        carrier_doppler::Float64
        carrier_phase::Float64
        code_doppler::Float64
        code_phase::Float64
        prompts_correlated_signals::Array{Complex{Float64}, 1}
    end

    struct Initials
        carrier_doppler::Float64
        carrier_phase::Float64
        code_doppler::Float64
        code_phase::Float64
    end
    
    struct JoinedTrackingResults
        l1_results::TrackingResults
        l5_results::TrackingResults
    end    

    include("discriminators.jl")
    include("loop_filters.jl")
    include("replicas.jl")
    include("tracking_loop.jl")
    include("joined_tracking.jl")
end