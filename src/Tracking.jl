module Tracking

using DocStringExtensions
using FastSinCos
using GNSSSignals
using SIMD
using SinCosLUT
using StaticArrays
using TrackingLoopFilters
using Dictionaries
using Accessors
using Polyester
using Random: AbstractRNG, Xoshiro

# The device-independent loop core — correlator records, discriminators, the
# bit buffer and its per-signal sync detectors, the C/N₀ estimators and the
# noise window, the loop-filter rules and the Doppler estimators' per-record
# `step_loop` — lives in TrackingLoops.jl, so that a hardware correlator's loop
# process runs the very same code without this package's sample-domain half.
# Its API is its own: users load it next to Tracking to configure a correlator,
# an estimator or a filter, and this package does not re-export it — except the
# functions a user reads this package's results with, below the export list.
# Only what the code here calls, extends or re-exports is imported, by name, so
# that a binding TrackingLoops renames or drops fails at precompile rather than
# on the first call that reaches it.
import TrackingLoops
import TrackingLoops:
    AbstractCN0Estimator,
    AbstractCorrelator,
    AbstractDopplerEstimator,
    AbstractNoiseEstimator,
    AbstractPostCorrFilter,
    append_noise_observation!,
    BitBuffer,
    buffer,
    _calc_num_code_blocks_that_form_a_bit,
    calc_num_code_blocks_to_integrate,
    ConventionalAssistedPLLAndDLL,
    CorrelatorNoiseEstimator,
    CorrelatorOutput,
    default_cn0_estimator,
    default_num_code_blocks_to_integrate,
    DefaultPostCorrFilter,
    EarlyPromptLateCorrelator,
    estimate_cn0,
    FixedNCOWord,
    fold_record,
    get_accumulators,
    get_code_block_buffer_type,
    get_correlator_sample_shifts,
    get_default_correlator,
    get_early,
    get_late,
    get_num_ants,
    get_prompt,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    init_estimator_state,
    LoopRecord,
    Maybe,
    MomentsCN0Estimator,
    _next_noise_prn,
    NO_LANDING_SAMPLE,
    _noise_density_and_ready,
    noise_density_type,
    _noise_observation,
    _noise_window_filling,
    NoiseEstimators,
    NoiseObservation,
    NoiseUpdateContext,
    _num_ants,
    _num_ants_of_density_type,
    _num_ants_val,
    NumAnts,
    _pool_taps,
    requires_noise_density,
    reset,
    reset_estimator_state,
    SatConventionalPLLAndDLL,
    step_loop,
    update,
    update_accumulator,
    update_noise!

using Unitful: upreferred, uconvert, ustrip, dimension, NoUnits, Hz, dBHz, ms, s
import Base.zero, Base.length, Base.resize!

export get_prn,
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
    get_bit_buffer,
    get_num_bits,
    track,
    track!,
    TrackedSignal,
    TrackedSat,
    get_signal,
    get_signals,
    get_doppler_estimator_state,
    max_code_length,
    current_code_wrap,
    update_estimator_on_handoff,
    CPUDownconvertAndCorrelator,
    CPUThreadedDownconvertAndCorrelator,
    Int16DownconvertAndCorrelator,
    Int16ThreadedDownconvertAndCorrelator,
    OneBitDownconvertAndCorrelator,
    OneBitThreadedDownconvertAndCorrelator,
    TwoBitDownconvertAndCorrelator,
    TwoBitThreadedDownconvertAndCorrelator,
    TrackState,
    add_satellite!,
    add_satellite,
    remove_satellite!,
    remove_satellite,
    merge_sats,
    get_sat_states,
    get_sat_state,
    AbstractDownconvertAndCorrelator,
    SatelliteDicts,
    SignalGroup,
    SignalGroups,
    BandMeasurement,
    BandMeasurements,
    band_keys,
    get_samples,
    get_sampling_frequency,
    get_intermediate_frequency,
    calc_signal_samples_to_integrate

# TrackingLoops functions a user of this package reads its results with: the
# ones this package implements for its states, satellites and signals, and the
# accessors of the correlators it hands back. Re-exported so reading a C/N₀, a
# prompt or the decoded bits needs no second package.
export estimate_cn0,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    get_num_ants,
    append_noise_observation!,
    get_prompt,
    get_early,
    get_late,
    get_accumulators

TupleLike{T<:Tuple} = Union{T,NamedTuple{<:Any,T}}

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
# The software fill path of the noise reference: despreading an untracked PRN
# through this package's own kernels (the window itself is TrackingLoops').
include("noise_estimators/correlator.jl")
include("sat_state.jl")

"""
$(SIGNATURES)

Main tracking state container holding satellite states for multiple GNSS
systems and the Doppler estimator (e.g., PLL/DLL). This is the primary
struct used for tracking operations.

`groups` is a NamedTuple of [`SignalGroup`](@ref)s. Each group bundles its
per-group `satellites` dictionary, signal-instance tuple, band, and
antenna count.

`noise_estimators` is a NamedTuple of [`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)s keyed by
**signal** id (`GNSSSignals.get_signal_id` — `:GPSL1CA`, `:GalileoE1B`, …), the
same NamedTuple idiom [`BandMeasurements`](@ref) uses for bands, so a lookup
folds to a compile-time constant. Keyed by signal and not by band because the
floor a record divides by is the *post-correlation* one, which depends on the
despreading modulation (see [`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)). A signal gets an
entry only where its C/N₀ estimator reads a noise density (see
[`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:AbstractCN0Estimator}})); signals with no such estimator get none, and
then the noise measurement costs exactly nothing. Each estimator averages **in
place**, so `TrackState` itself is never rebuilt for a noise update.

`noise_descriptor` is per-call **scratch**, not state: one reusable heap cell in
which the correlate step parks the chunk's noise descriptors so that a threaded
backend's parallel loop can reach them through a single pointer instead of a
by-value copy (see `_park_noise_items!`). It is shared — not copied — by every
`TrackState` derived from this one, exactly as the per-satellite scratch vectors
are, and nothing outside one `downconvert_and_correlate!` call reads it.

`sample_offsets` counts, per band (keyed by `get_band_id`), the samples every
[`track!`](@ref) call before this one consumed: the origin of the current
measurement on a sample grid that runs on across calls. The correlate phase
indexes its records from the current measurement's start, and the estimate phase
adds the band's offset when it hands a record to the Doppler estimator, so an
estimator that keeps per-satellite state of its own (TrackingLoops'
`VectorPLLAndDLL`) sees every satellite of the band on one time grid. The
in-place mutators share it with the state they return; [`track`](@ref) and the
non-mutating constructors ([`merge_sats`](@ref), [`add_satellite`](@ref),
[`remove_satellite`](@ref)) give their result a copy, so tracking the result
leaves the input's grid where it was.

A fresh `TrackState` starts every band at `0`. To carry an estimator that keeps
state of its own over into a rebuilt `TrackState` (to track another band with
the same `VectorPLLAndDLL`, say), pass the old state's offsets as the
constructor's `sample_offsets` keyword. Otherwise the estimator receives records
that lie in its past. A band the old state did not track starts on the same
time grid at `elapsed_time * sampling_frequency` (in samples of the new band).
"""
struct TrackState{G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    groups::G
    doppler_estimator::DE
    noise_estimators::NE
    noise_descriptor::Base.RefValue{Any}
    sample_offsets::Dict{Symbol,Int}
end

# Three-argument construction: the descriptor cell is scratch, so a freshly built
# `TrackState` starts with an empty one and fills it on its first correlate call.
# Every *derived* state (`TrackState(track_state; …)` and the in-place mutators)
# threads the existing cell through instead, which is what keeps the box inside it
# alive across `track!`'s copies and the steady state allocation-free.
TrackState(
    groups::SignalGroups,
    doppler_estimator::AbstractDopplerEstimator,
    noise_estimators::NoiseEstimators;
    sample_offsets = nothing,
) = TrackState(
    groups,
    doppler_estimator,
    noise_estimators,
    Base.RefValue{Any}(nothing),
    _initial_sample_offsets(groups, sample_offsets),
)

# Every band of `groups` at `0`, overridden by the `sample_offsets` the caller
# carries over (a `NamedTuple` or `AbstractDict` keyed by `get_band_id`). The
# result is always a fresh `Dict`, so it never aliases the caller's.
function _initial_sample_offsets(groups::SignalGroups, sample_offsets)
    offsets = Dict{Symbol,Int}(get_band_id(group.band) => 0 for group in groups)
    isnothing(sample_offsets) && return offsets
    for (band, offset) in pairs(sample_offsets)
        haskey(offsets, band) || throw(
            ArgumentError(
                string(
                    "`sample_offsets` names band `",
                    band,
                    "`, which no group of this TrackState tracks (bands: ",
                    join(keys(offsets), ", "),
                    ").",
                ),
            ),
        )
        offset >= 0 || throw(
            ArgumentError("`sample_offsets` must be non-negative, got $band => $offset."),
        )
        offsets[band] = offset
    end
    offsets
end

include("sample_parameters.jl")
include("downconvert_and_correlate.jl")
include("downconvert_and_correlate_fused.jl")
include("downconvert_and_correlate_cpu.jl")
include("downconvert_and_correlate_int16.jl")
include("downconvert_and_correlate_onebit.jl")
include("downconvert_and_correlate_twobit.jl")
include("conventional_pll_and_dll.jl")
include("tracking_state.jl")
include("track.jl")
include("precompile.jl")

end
