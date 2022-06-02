module Tracking

    using CUDA
    using DocStringExtensions
    using GNSSSignals
    using LoopVectorization
    using StaticArrays
    using StructArrays
    using TrackingLoopFilters

    using Acquisition: AcquisitionResults
    using Unitful: upreferred, Hz, dBHz, ms
    import Base.zero, Base.length, Base.resize!, LinearAlgebra.dot

    export get_early, get_prompt, get_late
    export get_early_index,get_prompt_index, get_late_index
    export get_correlator
    export get_accumulator, get_accumulators
    export get_num_accumulators
    export get_carrier_doppler, get_carrier_phase
    export get_code_doppler, get_code_phase
    export get_correlator_sample_shifts
    export get_early_late_sample_spacing
    export get_early_late_index_shift
    export get_secondary_code_or_bit_found
    export get_correlator_carrier_phase, get_correlator_carrier_frequency
    export get_state
    export get_system
    export get_cn0
    export track
    export get_bits
    export get_num_bits

    export TrackingState
    export NumAnts
    export NumAccumulators
    export MomentsCN0Estimator
    export AbstractCN0Estimator
    export EarlyPromptLateCorrelator
    #export VeryEarlyPromptLateCorrelator
    export SecondaryCodeOrBitDetector
    export GainControlledSignal

    struct NumAnts{x}
    end

    NumAnts(x) = NumAnts{x}()

    struct NumAccumulators{x}
    end

    NumAccumulators(x) = NumAccumulators{x}()

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
    include("boc.jl")
    include("downconvert_and_correlate.jl")
end
