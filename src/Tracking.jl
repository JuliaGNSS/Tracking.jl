module Tracking

using BitIntegers
using DocStringExtensions
using FastSinCos
using GNSSSignals
using SIMD
using SinCosLUT
using SpecialFunctions: erfinv
using StaticArrays
using TrackingLoopFilters
using Dictionaries
using Accessors
using Polyester
using Random: AbstractRNG, Xoshiro

# 1800-bit exact-width unsigned for the 1800-chip overlay-code searches of
# GPS L1C-P and BeiDou B1C-P.
# Defined once at module load. Benchmarked at ~71 μs for the full
# 1800-phase Hamming-distance sweep, ~1.5× faster than a padded
# UInt1856 variant because no mask is needed on shift/XOR (see the
# sync-detection-redesign plan in docs/plans for the comparison).
BitIntegers.@define_integers 1800

using Unitful: upreferred, uconvert, ustrip, dimension, NoUnits, Hz, dBHz, ms, s
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
    get_last_fully_integrated_num_code_blocks,
    get_last_fully_integrated_integration_time,
    get_filtered_prompts,
    get_correlator_outputs,
    append_correlator_output!,
    CorrelatorOutput,
    get_bit_buffer,
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
    NWPRCN0Estimator,
    NoCN0Estimator,
    NoiseRefCN0Estimator,
    AbstractCN0Estimator,
    CN0UpdateContext,
    requires_noise_density,
    AbstractNoiseEstimator,
    CorrelatorNoiseEstimator,
    NoiseObservation,
    noise_observation,
    noise_observation_from_correlator,
    noise_observation_from_samples,
    append_noise_observation!,
    update_noise!,
    get_noise_density,
    EarlyPromptLateCorrelator,
    VeryEarlyPromptLateCorrelator,
    AbstractPostCorrFilter,
    get_weights,
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
    Int16DownconvertAndCorrelator,
    Int16ThreadedDownconvertAndCorrelator,
    OneBitDownconvertAndCorrelator,
    OneBitThreadedDownconvertAndCorrelator,
    TwoBitDownconvertAndCorrelator,
    TwoBitThreadedDownconvertAndCorrelator,
    ConventionalPLLAndDLL,
    ConventionalAssistedPLLAndDLL,
    VectorPLLAndDLL,
    SatVectorPLLAndDLL,
    enable_vt!,
    disable_vt!,
    reset_code_discr_acc!,
    reset_carrier_discr_acc!,
    mean_code_discr,
    mean_carrier_discr,
    set_code_freq_updates!,
    set_carrier_freq_updates!,
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
    get_band_id,
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

The per-sat correlation loop, per-group body, and public
`downconvert_and_correlate(!)` entry points are defined once on this abstract
type (see `downconvert_and_correlate_cpu.jl`); a subtype customises behaviour by
overriding the dispatch hooks it needs — `_despread_one_signal!` (the one
correlation primitive, which both the per-satellite path and the noise reference
go through), `_correlate_signals` / `_scratch_buffers` (multi-signal kernel +
scratch), `_threading` (serial vs. Polyester `@batch`, default serial), and
`_check_sample_type` (per-backend sample-type check, default no-op). A subtype
that overrides none inherits the single-threaded CPU plumbing rather than getting
a `MethodError`; one that overrides only `_despread_one_signal!` gets a working
single-signal path, satellites and noise measurement alike.
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
# `cn0_estimators/` after `bit_buffer.jl`: the CN0 estimators' update context
# carries the navigation-bit state (`BitBuffer`) and reads the signal's
# blocks-per-bit trait from there. Within the folder the shared file comes
# first, since every concrete estimator subtypes `AbstractCN0Estimator`.
include("bit_buffer.jl")
include("cn0_estimators/cn0_estimator.jl")
include("cn0_estimators/moments.jl")
include("cn0_estimators/no_cn0.jl")
include("cn0_estimators/nwpr.jl")
include("cn0_estimators/noise_ref.jl")
include("correlators/correlator.jl")
include("correlators/early_prompt_late.jl")
include("correlators/very_early_prompt_late.jl")
# `noise_estimators/` after `correlators/`: the software source despreads
# through a real correlator, and `TrackState`'s `NoiseEstimators` field type
# needs `AbstractNoiseEstimator` to exist before the struct is defined below.
include("noise_estimators/noise_estimator.jl")
include("noise_estimators/correlator.jl")
include("discriminators.jl")
include("post_corr_filter.jl")
include("gps/l1ca.jl")
include("gps/l1c_d.jl")
include("gps/l1c_p.jl")
include("gps/l2c.jl")
include("gps/l5.jl")
include("galileo/e1b.jl")
include("galileo/e1c.jl")
include("galileo/e5a.jl")
include("galileo/e5b.jl")
include("galileo/e6.jl")
include("beidou/b1i.jl")
include("beidou/b3i.jl")
include("beidou/b2a.jl")
include("beidou/b2b.jl")
include("beidou/b1c.jl")
include("sat_state.jl")

"""
$(SIGNATURES)

Main tracking state container holding satellite states for multiple GNSS
systems and the Doppler estimator (e.g., PLL/DLL). This is the primary
struct used for tracking operations.

`groups` is a NamedTuple of [`SignalGroup`](@ref)s. Each group bundles its
per-group `satellites` dictionary, signal-instance tuple, band, and
antenna count.

`noise_estimators` is a NamedTuple of [`AbstractNoiseEstimator`](@ref)s keyed by
**signal** id (`GNSSSignals.get_signal_id` — `:GPSL1CA`, `:GalileoE1B`, …), the
same NamedTuple idiom [`BandMeasurements`](@ref) uses for bands, so a lookup
folds to a compile-time constant. Keyed by signal and not by band because the
floor a record divides by is the *post-correlation* one, which depends on the
despreading modulation (see [`AbstractNoiseEstimator`](@ref)). A signal gets an
entry only where its C/N₀ estimator reads a noise density (see
[`requires_noise_density`](@ref)); signals with no such estimator get none, and
then the noise measurement costs exactly nothing. Each estimator averages **in
place**, so `TrackState` itself is never rebuilt for a noise update.

`noise_descriptor` is per-call **scratch**, not state: one reusable heap cell in
which the correlate step parks the chunk's noise descriptors so that a threaded
backend's parallel loop can reach them through a single pointer instead of a
by-value copy (see `_park_noise_items!`). It is shared — not copied — by every
`TrackState` derived from this one, exactly as the per-satellite scratch vectors
are, and nothing outside one `downconvert_and_correlate!` call reads it.
"""
struct TrackState{G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    groups::G
    doppler_estimator::DE
    noise_estimators::NE
    noise_descriptor::Base.RefValue{Any}
end

# Three-argument construction: the descriptor cell is scratch, so a freshly built
# `TrackState` starts with an empty one and fills it on its first correlate call.
# Every *derived* state (`TrackState(track_state; …)` and the in-place mutators)
# threads the existing cell through instead, which is what keeps the box inside it
# alive across `track!`'s copies and the steady state allocation-free.
TrackState(
    groups::SignalGroups,
    doppler_estimator::AbstractDopplerEstimator,
    noise_estimators::NoiseEstimators,
) = TrackState(groups, doppler_estimator, noise_estimators, Base.RefValue{Any}(nothing))

include("sample_parameters.jl")
include("downconvert_and_correlate.jl")
include("downconvert_and_correlate_fused.jl")
include("downconvert_and_correlate_cpu.jl")
include("downconvert_and_correlate_int16.jl")
include("downconvert_and_correlate_onebit.jl")
include("downconvert_and_correlate_twobit.jl")
include("conventional_pll_and_dll.jl")
include("vector_pll_and_dll.jl")
include("tracking_state.jl")
include("track.jl")

end
