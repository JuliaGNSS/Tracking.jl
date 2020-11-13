module Tracking
    using
        DocStringExtensions,
        GNSSSignals,
        StaticArrays,
        TrackingLoopFilters,
        StructArrays,
        LoopVectorization,
        CUDA
    using Unitful: upreferred, Hz, dBHz, ms
    import Base.zero, Base.length, Base.resize!, LinearAlgebra.dot

    export
        get_early,
        get_prompt,
        get_late,
        get_correlator,
        get_carrier_doppler,
        get_carrier_phase,
        get_code_doppler,
        get_code_phase,
        get_early_late_sample_shift,
        get_secondary_code_or_bit_found,
        get_state,
        get_cn0,
        track,
        TrackingState,
        NumAnts,
        MomentsCN0Estimator,
        AbstractCN0Estimator,
        get_bits,
        get_num_bits,
        EarlyPromptLateCorrelator,
        VeryEarlyPromptLateCorrelator,
        SecondaryCodeOrBitDetector,
        GainControlledSignal

    struct NumAnts{x}
    end

    SOA = NamedTuple{(:re, :im),Tuple{Array{Float32,1},Array{Float32,1}}}
    SOC = NamedTuple{(:re, :im),Tuple{CuArray{Float32,1},CuArray{Float32,1}}}

    NumAnts(x) = NumAnts{x}()

    include("agc.jl")
    include("code_replica.jl")
    include("carrier_replica.jl")
    include("downconvert.jl")
    include("cn0_estimation.jl")
    include("discriminators.jl")
    include("bit_buffer.jl")
    include("correlator.jl")
    include("secondary_code_or_bit_detector.jl")
    include("tracking_state.jl")
    include("tracking_results.jl")
    include("tracking_loop.jl")
    include("gpsl1.jl")
    include("gpsl5.jl")
    include("galileo_e1b.jl")
end
