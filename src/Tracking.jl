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

# The device-independent loop core (correlator records, discriminators, bit
# buffer and sync, C/N₀ and noise estimation, loop filters, the Doppler
# estimators' `step_loop`) lives in TrackingLoops.jl, so a hardware correlator's
# loop runs the same code. Users load it next to Tracking; it is not re-exported,
# except the result readers below the export list. Imported by name so a binding
# TrackingLoops renames or drops fails at precompile.
import TrackingLoops
import TrackingLoops:
    AbstractCN0Estimator,
    AbstractCorrelator,
    AbstractDopplerEstimator,
    AbstractNoiseEstimator,
    AbstractPostCorrFilter,
    append_noise_observation!,
    BitBuffer,
    _calc_num_code_blocks_that_form_a_bit,
    calc_num_code_blocks_to_integrate,
    carrier_phase_correction,
    ConventionalAssistedPLLAndDLL,
    CorrelatorNoiseEstimator,
    CorrelatorOutput,
    default_cn0_estimator,
    default_num_code_blocks_to_integrate,
    DefaultPostCorrFilter,
    EarlyPromptLateCorrelator,
    estimate_cn0,
    FixedNCOWord,
    apply_record,
    get_accumulators,
    get_correlator_sample_shifts,
    get_default_correlator,
    get_early,
    get_late,
    get_num_ants,
    get_prompt,
    get_soft_bits,
    has_bit_or_secondary_code_been_found,
    init_estimator_state,
    is_dataless_with_secondary_code,
    LoopRecord,
    SignalLoopState,
    Maybe,
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
    step_loop,
    update,
    update_accumulator,
    update_noise!

# `dBHz` serves the manual: its custom C/N₀ estimator doctest reads `Tracking.dBHz`.
using Unitful: upreferred, uconvert, ustrip, dimension, NoUnits, Hz, dBHz, s
# Extended with the methods that hand tracked satellites over for a PVT solve; see
# `position_velocity_time.jl`.
import PositionVelocityTime
using GNSSDecoder: GNSSDecoderState
import Base.zero, Base.length, Base.resize!

export get_prn,
    get_code_phase,
    get_code_doppler,
    get_carrier_phase,
    get_carrier_phase_correction,
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

# TrackingLoops functions that read this package's results, re-exported so
# reading them needs no second package.
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

Supertype of the downconvert-and-correlate backends.

The `downconvert_and_correlate(!)` loops are defined once on this type (see
`downconvert_and_correlate_cpu.jl`); a subtype implements the hooks:
`_despread_one_signal!` (the correlation primitive, for satellites and the noise
reference alike; alone it gives a working single-signal backend),
`_correlate_signals` / `_scratch_buffers` (multi-signal kernel and scratch),
`_threading` (default serial) and `_check_sample_type` (default no-op).
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
# Software fill path of the noise reference: despreading an untracked PRN.
include("noise_estimators/correlator.jl")
include("sat_state.jl")

"""
$(SIGNATURES)

The tracking state: the satellites of every signal group and the Doppler
estimator.

`groups` is a NamedTuple of [`SignalGroup`](@ref)s.

`noise_estimators` is a NamedTuple of
[`AbstractNoiseEstimator`](@extref TrackingLoops.AbstractNoiseEstimator)s keyed by
**signal** id (`GNSSSignals.get_signal_id`), since the post-correlation floor
depends on the despreading modulation. Only signals whose C/N₀ estimator reads a
noise density (see
[`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:AbstractCN0Estimator}}))
get one, so the others pay nothing. They average **in place**.

`noise_descriptor` is per-call **scratch**: a heap cell where the correlate step
parks the chunk's noise descriptors, so a threaded backend reaches them through
one pointer (see `_park_noise_items!`). Shared by every derived `TrackState`.

`sample_offsets` counts, per band (keyed by `get_band_id`), the samples earlier
[`track!`](@ref) calls consumed: the origin of the current measurement on a
sample grid that runs on across calls. The estimate phase adds it to a record's
index, so an estimator with per-satellite state of its own (TrackingLoops'
`VectorPLLAndDLL`) sees every satellite on one time grid. The in-place mutators
share it; [`track`](@ref), [`merge_sats`](@ref), [`add_satellite`](@ref) and
[`remove_satellite`](@ref) give their result a copy.

A fresh `TrackState` starts every band at `0`. To carry such an estimator over
into a rebuilt `TrackState` (e.g. to add a band), pass the old state's offsets
as the `sample_offsets` keyword, or the estimator receives records from its
past. A band the old state did not track starts at
`elapsed_time * sampling_frequency` samples of the new band.
"""
struct TrackState{G<:SignalGroups,DE<:AbstractDopplerEstimator,NE<:NoiseEstimators}
    groups::G
    doppler_estimator::DE
    noise_estimators::NE
    noise_descriptor::Base.RefValue{Any}
    sample_offsets::Dict{Symbol,Int}
end

# A fresh state gets an empty descriptor cell; derived states pass theirs on,
# which keeps the steady state allocation-free.
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

# Every band at `0`, overridden by the caller's `sample_offsets`, in a fresh `Dict`.
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
include("position_velocity_time.jl")
include("precompile.jl")

end
