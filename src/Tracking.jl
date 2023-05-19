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

    export get_early, get_prompt, get_late,
        get_prn,
        get_code_phase,
        get_code_doppler,
        get_carrier_phase,
        get_carrier_doppler,
        get_integrated_samples,
        get_correlator,
        get_last_fully_integrated_correlator,
        get_last_fully_integrated_filtered_prompt,
        get_sample_of_last_fully_integrated_correlator,
        get_secondary_code_or_bit_detector,
        get_prompts_buffer,
        get_bit_buffer,
        get_bits,
        get_accumulators,
        get_early_late_sample_spacing,
        track,
        NumAnts,
        NumAccumulators,
        MomentsCN0Estimator,
        AbstractCN0Estimator,
        EarlyPromptLateCorrelator,
        VeryEarlyPromptLateCorrelator,
        SecondaryCodeOrBitDetector,
        GainControlledSignal,
        AbstractPostCorrFilter,
        SatState,
        SystemSatsState,
        CPUDownconvertAndCorrelator,
        GPUDownconvertAndCorrelator,
        ConventionalPLLAndDLL,
        ConventionalPLLsAndDLLs,
        DefaultPostCorrFilter,
        TrackState,
        add_sats!,
        remove_sats!,
        get_sat_states,
        get_sat_state,
        get_system

    struct NumAnts{x}
    end

    NumAnts(x) = NumAnts{x}()

    struct NumAccumulators{x}
    end

    NumAccumulators(x) = NumAccumulators{x}()

    TupleLike{T <: Tuple} = Union{T, NamedTuple{<:Any, T}}

    struct DopplersAndFilteredPrompt
        carrier_doppler::typeof(1.0Hz)
        code_doppler::typeof(1.0Hz)
        filtered_prompt::ComplexF64
    end

    """
    $(SIGNATURES)

    Get the number of samples in the signal.
    """
    @inline function get_num_samples(signal)
        length(signal)
    end

    @inline function get_num_samples(signal::AbstractMatrix)
        size(signal, 1)
    end

    include("code_replica.jl")
    include("carrier_replica.jl")
    include("downconvert.jl")
    include("cn0_estimation.jl")
    include("bit_buffer.jl")
    include("correlator.jl")
    include("discriminators.jl")
    include("post_corr_filter.jl")
    include("secondary_code_or_bit_detector.jl")
    include("gpsl1.jl")
    include("gpsl5.jl")
    include("galileo_e1b.jl")
    include("sat_state.jl")
    include("sample_parameters.jl")
    include("update_sat_state.jl")
    include("downconvert_and_correlate.jl")
    include("gpu_downconvert_and_correlate.jl")
    include("conventional_pll_and_dll.jl")
    include("tracking_state.jl")
    include("track.jl")
end
