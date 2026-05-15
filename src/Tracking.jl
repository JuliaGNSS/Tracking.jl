module Tracking

using DocStringExtensions
using FastSinCos
using GNSSSignals
using SIMD
using StaticArrays
using TrackingLoopFilters
using Dictionaries
using Accessors
using ConstructionBase
using Polyester

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
    get_signal_start_sample,
    get_correlator,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_bit_buffer,
    get_bits,
    get_num_bits,
    get_accumulators,
    get_early_late_sample_spacing,
    get_num_ants,
    has_bit_or_secondary_code_been_found,
    track,
    track!,
    NumAnts,
    NumAccumulators,
    MomentsCN0Estimator,
    AbstractCN0Estimator,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    AbstractPostCorrFilter,
    TrackedSignal,
    TrackedSat,
    get_signal,
    get_signals,
    get_doppler_estimator_state,
    max_code_length,
    AbstractDopplerEstimator,
    init_estimator_state,
    update_estimator_on_handoff,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    DefaultPostCorrFilter,
    TrackState,
    add_satellite!,
    add_satellite,
    merge_sats,
    filter_out_sats,
    get_sat_states,
    get_sat_state,
    get_system,
    estimate_cn0,
    get_default_correlator,
    AbstractCorrelator,
    AbstractDownconvertAndCorrelator,
    SatelliteDicts,
    get_num_accumulators,
    get_correlator_sample_shifts,
    calc_signal_samples_to_integrate,
    update,
    get_code_frequency,
    get_code_length,
    get_codes,
    get_modulation,
    get_secondary_code,
    update_accumulator

const Maybe{T} = Union{T,Nothing}

"""
$(SIGNATURES)

Type parameter wrapper for specifying the number of antennas in the system.
Use `NumAnts(n)` to create an instance.
"""
struct NumAnts{x} end

NumAnts(x) = NumAnts{x}()

"""
$(SIGNATURES)

Type parameter wrapper for specifying the number of correlator accumulators.
Use `NumAccumulators(n)` to create an instance.
"""
struct NumAccumulators{x} end

NumAccumulators(x) = NumAccumulators{x}()

TupleLike{T<:Tuple} = Union{T,NamedTuple{<:Any,T}}

"""
$(SIGNATURES)

Abstract supertype for doppler estimators. Concrete subtypes carry estimator
configuration (and any cross-satellite or cross-system shared state). The
per-satellite state used by the estimator lives in each [`TrackedSat`](@ref)
wrapper — see [`init_estimator_state`](@ref) for the extension point.
"""
abstract type AbstractDopplerEstimator end

"""
$(SIGNATURES)

Abstract downconverter and correlator type. Structs for
downconversion and correlation must have this abstract type as a
parent.
"""
abstract type AbstractDownconvertAndCorrelator end

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
include("correlators/correlator.jl")
include("correlators/early_prompt_late.jl")
include("correlators/very_early_prompt_late.jl")
include("discriminators.jl")
include("post_corr_filter.jl")
include("gpsl1ca.jl")
include("gpsl5i.jl")
include("galileo_e1b.jl")
include("sat_state.jl")


"""
$(SIGNATURES)

Main tracking state container holding satellite states for multiple GNSS systems
and the Doppler estimator (e.g., PLL/DLL). This is the primary struct used for
tracking operations.

`signal_groups` is a NamedTuple that mirrors the keys of `satellites` and
holds the per-capability signal-instance tuple. `add_satellite!` reads it
to build new sats with the right signal shape (signal instances can't be
reconstructed from the `TrackedSat` value type alone — concrete signal
subtypes like `GPSL1CA{Matrix{Int16}}` carry the code-matrix instance as
a field, so re-instantiation requires the original instance).
"""
struct TrackState{S<:SatelliteDicts,SG<:TupleLike,DE<:AbstractDopplerEstimator}
    satellites::S
    signal_groups::SG
    doppler_estimator::DE
end

include("sample_parameters.jl")
include("downconvert_and_correlate.jl")
include("downconvert_and_correlate_fused.jl")
include("downconvert_and_correlate_cpu.jl")
include("conventional_pll_and_dll.jl")
include("tracking_state.jl")
include("track.jl")

end
