module Tracking

using CUDA
using DocStringExtensions
using GNSSSignals
using LoopVectorization
using StaticArrays
using StructArrays
using TrackingLoopFilters
using Dictionaries
using Accessors

using Acquisition: AcquisitionResults
using Unitful: upreferred, Hz, dBHz, ms
import Base.zero, Base.length, Base.resize!, LinearAlgebra.dot

export get_early,
    get_prompt,
    get_late,
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
    get_bit_buffer,
    get_bits,
    get_accumulators,
    get_early_late_sample_spacing,
    get_system_sats_state,
    track,
    NumAnts,
    NumAccumulators,
    MomentsCN0Estimator,
    AbstractCN0Estimator,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    SecondaryCodeOrBitDetector,
    AbstractPostCorrFilter,
    SatState,
    SystemSatsState,
    CPUSatDownconvertAndCorrelator,
    CPUSystemDownconvertAndCorrelator,
    GPUSatDownconvertAndCorrelator,
    GPUSystemDownconvertAndCorrelator,
    SystemConventionalPLLAndDLL,
    TrackingConventionalPLLAndDLL,
    NoSatPostProcess,
    NoSystemPostProcess,
    NoTrackingPostProcess,
    ConventionalPLLAndDLL,
    DefaultPostCorrFilter,
    TrackState,
    merge_sats,
    filter_out_sats,
    get_sat_states,
    get_sat_state,
    get_system,
    estimate_cn0,
    get_default_correlator,
    convert_code_to_texture_memory

const Maybe{T} = Union{T,Nothing}

const StructVecOrMat{T} = Union{StructVector{T},StructArray{T,2}}

struct NumAnts{x} end

NumAnts(x) = NumAnts{x}()

struct NumAccumulators{x} end

NumAccumulators(x) = NumAccumulators{x}()

TupleLike{T<:Tuple} = Union{T,NamedTuple{<:Any,T}}

abstract type AbstractSatDopplerEstimator end
abstract type AbstractSystemDopplerEstimator end
abstract type AbstractTrackingDopplerEstimator end

"""
$(SIGNATURES)

Abstract downconverter and correlator type. Structs for
downconversion and correlation must have this abstract type as a
parent.
"""
abstract type AbstractSatDownconvertAndCorrelator end
abstract type AbstractSystemDownconvertAndCorrelator end

abstract type AbstractSatPostProcess end
abstract type AbstractSystemPostProcess end
abstract type AbstractTrackingPostProcess end

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

struct TrackState{
    S<:MultipleSystemSatsState{
        N,
        <:AbstractGNSS,
        <:SatState{
            <:AbstractCorrelator,
            <:AbstractSatDopplerEstimator,
            <:AbstractSatDownconvertAndCorrelator,
            <:AbstractSatPostProcess,
        },
        <:AbstractSystemDopplerEstimator,
        <:AbstractSystemDownconvertAndCorrelator,
        <:AbstractSystemPostProcess,
    } where {N},
    DE<:AbstractTrackingDopplerEstimator,
    P<:AbstractTrackingPostProcess,
}
    multiple_system_sats_state::S
    doppler_estimator::DE
    post_process::P
    num_samples::Int
end

include("sample_parameters.jl")
include("doppler_estimator.jl")
include("downconvert_and_correlate.jl")
include("downconvert_and_correlate_cpu.jl")
include("downconvert_and_correlate_gpu.jl")
include("conventional_pll_and_dll.jl")
include("tracking_state.jl")
include("post_process.jl")
include("track.jl")
end
