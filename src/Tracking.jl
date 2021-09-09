module Tracking
    using
        DocStringExtensions,
        GNSSSignals,
        StaticArrays,
        TrackingLoopFilters,
        StructArrays,
        LoopVectorization
    using Unitful: upreferred, Hz, dBHz, ms
    import Base.zero, Base.length, Base.resize!

    export
        get_early,
        get_prompt,
        get_late,
        get_correlator,
        get_accumulators,
        get_num_accumulators,
        get_early_index,
        get_prompt_index,
        get_late_index,
        get_accumulator,
        get_carrier_doppler,
        get_carrier_phase,
        get_code_doppler,
        get_code_phase,
        get_correlator_sample_shifts,
        get_early_late_sample_spacing,
        get_early_late_index_shift,
        get_secondary_code_or_bit_found,
        get_correlator_carrier_phase,
        get_correlator_carrier_frequency,
        get_state,
        get_system,
        get_cn0,
        track,
        TrackingState,
        NumAnts,
        NumAccumulators,
        MomentsCN0Estimator,
        AbstractCN0Estimator,
        get_bits,
        get_num_bits,
        EarlyPromptLateCorrelator,
        #VeryEarlyPromptLateCorrelator,
        SecondaryCodeOrBitDetector,
        GainControlledSignal

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
