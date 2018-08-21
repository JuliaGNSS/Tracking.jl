module Tracking
    using DocStringExtensions, GNSSSignals

    export prompt,
        init_tracking,
        GNSSSystem,
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

    struct GNSSSystem
        f₀::Float64 # Carrier center frequency
        code_f₀::Float64 # Code center frequency
        sampling_freq::Float64
        gen_sampled_code::Function
        calc_next_code_phase::Function
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