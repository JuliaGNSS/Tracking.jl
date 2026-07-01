module Tracking

using BitIntegers
using DocStringExtensions
using FastSinCos
using GNSSSignals
using SIMD
using SpecialFunctions: erfinv
using StaticArrays
using TrackingLoopFilters
using Dictionaries
using Accessors
using Polyester

# 1800-bit exact-width unsigned for the GPS L1C-P overlay-code search.
# Defined once at module load. Benchmarked at ~71 μs for the full
# 1800-phase Hamming-distance sweep, ~1.5× faster than a padded
# UInt1856 variant because no mask is needed on shift/XOR (see the
# sync-detection-redesign plan in docs/plans for the comparison).
BitIntegers.@define_integers 1800

using Unitful: upreferred, uconvert, Hz, dBHz, ms, s
import Base.zero, Base.length, Base.resize!

export get_early,
    get_prompt,
    get_late,
    get_prn,
    get_code_phase,
    get_code_doppler,
    get_carrier_phase,
    get_carrier_doppler,
    get_integrated_samples,
    get_preferred_num_code_blocks_to_integrate,
    set_preferred_num_code_blocks_to_integrate!,
    reset_loop_filters!,
    get_signal_start_sample,
    get_correlator,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_filtered_prompt,
    get_filtered_prompts,
    get_bit_buffer,
    get_bits,
    get_num_bits,
    get_soft_bits,
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
    current_code_wrap,
    AbstractDopplerEstimator,
    init_estimator_state,
    update_estimator_on_handoff,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    OneBitDownconvertAndCorrelator,
    OneBitThreadedDownconvertAndCorrelator,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    DefaultPostCorrFilter,
    TrackState,
    add_satellite!,
    add_satellite,
    remove_satellite!,
    remove_satellite,
    merge_sats,
    get_sat_states,
    get_sat_state,
    estimate_cn0,
    get_default_correlator,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    AbstractCorrelator,
    AbstractDownconvertAndCorrelator,
    SatelliteDicts,
    SignalGroup,
    SignalGroups,
    BandMeasurement,
    BandMeasurements,
    band_key,
    band_keys,
    get_samples,
    get_sampling_frequency,
    get_intermediate_frequency,
    get_num_accumulators,
    get_correlator_sample_shifts,
    calc_signal_samples_to_integrate,
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

include("band_measurement.jl")
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
include("gpsl1c_d.jl")
include("gpsl1c_p.jl")
include("gpsl5i.jl")
include("galileo_e1b.jl")
include("sat_state.jl")

"""
$(SIGNATURES)

Main tracking state container holding satellite states for multiple GNSS
systems and the Doppler estimator (e.g., PLL/DLL). This is the primary
struct used for tracking operations.

`groups` is a NamedTuple of [`SignalGroup`](@ref)s. Each group bundles its
per-group `satellites` dictionary, signal-instance tuple, band, and
antenna count.
"""
struct TrackState{G<:SignalGroups,DE<:AbstractDopplerEstimator}
    groups::G
    doppler_estimator::DE
end

include("sample_parameters.jl")
include("downconvert_and_correlate.jl")
include("downconvert_and_correlate_fused.jl")
include("downconvert_and_correlate_cpu.jl")
include("downconvert_and_correlate_onebit.jl")
include("conventional_pll_and_dll.jl")
include("tracking_state.jl")
include("track.jl")

end
