module Tracking

using DocStringExtensions
using FastSinCos
using GNSSSignals
using SIMD
using SinCosLUT
using StaticArrays
using Dictionaries
using Accessors
using Polyester

# The device-independent loop core (post-correlation processing, discriminators,
# bit synchronization, C/N₀ and noise estimation, loop filters, the tracking
# loops' `step_loop!`) lives in TrackingLoops.jl, so a hardware correlator's
# loop runs the same code. Users load it next to Tracking; it is not re-exported,
# except the result readers below the export list. Imported by name so a binding
# TrackingLoops renames or drops fails at precompile.
import TrackingLoops
# The carrier filter of the default tracking loop, `ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)`.
using TrackingLoopFilters: ThirdOrderAssistedBilinearLF
import TrackingLoops:
    AbstractCN0Estimator,
    AbstractCorrelatorBank,
    AbstractNoiseDensityEstimator,
    AbstractBeamformer,
    AbstractTrackingLoop,
    append_despread_noise!,
    append_noise_observation!,
    BitSynchronizer,
    calc_carrier_phase_correction,
    calc_correlator_sample_offsets,
    calc_num_code_periods_to_integrate,
    ConstantNCOCommand,
    CorrelatorNoiseEstimator,
    CorrelatorBankOutput,
    default_cn0_estimator,
    default_num_code_periods_to_integrate,
    DefaultBeamformer,
    drain_soft_bits!,
    EarlyPromptLateCorrelatorBank,
    estimate_cn0,
    get_correlator_outputs,
    get_default_correlator_bank,
    get_early,
    get_late,
    get_num_antennas,
    get_num_code_periods_into_bit,
    get_prompt,
    get_secondary_code_phase,
    get_soft_bits,
    init_loop_state,
    is_bit_or_secondary_code_synced,
    is_dataless_with_secondary_code,
    map_correlator_outputs,
    noise_correlator_bank,
    noise_density_and_ready,
    noise_despreads,
    noise_window_filling,
    NoiseEstimators,
    NoiseObservation,
    NoiseUpdateContext,
    NumAntennas,
    PostCorrelationState,
    PredetectionIntegration,
    process_correlator_bank_output!,
    requires_noise_density,
    reseed_loop_state,
    ScalarTrackingLoop,
    step_loop!,
    update_noise!,
    with_correlator_outputs

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
    get_carrier_doppler,
    get_num_integrated_samples,
    get_preferred_num_code_periods_to_integrate,
    set_preferred_num_code_periods_to_integrate!,
    reset_loop_filters!,
    get_next_sample_index,
    get_correlator_bank,
    get_last_fully_integrated_correlator,
    get_last_fully_integrated_beamformed_prompt,
    get_last_fully_integrated_num_code_blocks,
    get_last_fully_integrated_integration_time,
    get_beamformed_prompts,
    get_correlator_bank_outputs,
    append_correlator_bank_output!,
    get_bit_synchronizer,
    get_num_bits,
    track,
    track!,
    TrackedSignal,
    TrackedSat,
    get_signal,
    get_signals,
    get_loop_state,
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
    get_tracked_sats,
    get_tracked_sat,
    SignalGroup,
    BandMeasurement,
    get_band_ids,
    get_signal_observations

# TrackingLoops functions that read this package's results, re-exported so
# reading them needs no second package.
export estimate_cn0,
    calc_carrier_phase_correction,
    get_soft_bits,
    is_bit_or_secondary_code_synced,
    get_num_antennas,
    append_noise_observation!,
    get_prompt,
    get_early,
    get_late,
    get_correlator_outputs

# The default tracking loop: a PLL assisted by an FLL until phase lock, and a
# carrier-aided DLL.
_default_tracking_loop() = ScalarTrackingLoop(ThirdOrderAssistedBilinearLF)

TupleLike{T<:Tuple} = Union{T,NamedTuple{<:Any,T}}
const Maybe{T} = Union{T,Nothing}

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

Get the number of samples per antenna in `samples`.
"""
@inline function get_num_samples(samples)
    length(samples)
end

@inline function get_num_samples(samples::AbstractMatrix)
    size(samples, 1)
end

include("band_measurement.jl")
include("code_replica.jl")
include("carrier_replica.jl")
include("downconvert.jl")
# Software fill path of the noise reference: despreading an untracked PRN.
include("noise_estimators/correlator.jl")
include("tracked_sat.jl")

"""
$(SIGNATURES)

The tracking state: the satellites of every signal group and the tracking loop.

`groups` is a NamedTuple of [`SignalGroup`](@ref)s.

`noise_estimators` is a NamedTuple of
[`AbstractNoiseDensityEstimator`](@extref TrackingLoops.AbstractNoiseDensityEstimator)s keyed by
**signal** id (`GNSSSignals.get_signal_id`), since the post-correlation floor
depends on the despreading modulation. Only signals whose C/N₀ estimator reads a
noise density (see
[`requires_noise_density`](@extref TrackingLoops.requires_noise_density-Tuple{Type{<:TrackingLoops.AbstractCN0Estimator}}))
get one, so the others pay nothing. They average **in place**.

`noise_descriptor` is per-call **scratch**: a heap cell where the correlate step
parks the chunk's noise descriptors, so a threaded backend reaches them through
one pointer (see `_park_noise_items!`). Shared by every derived `TrackState`.

`sample_index_offsets` counts, per band (keyed by `get_band_id`), the samples earlier
[`track!`](@ref) calls consumed: the origin of the current measurement on a
sample grid that runs on across calls. [`step_tracking_loops!`](@ref) adds it to
an integration's sample index, so a tracking loop with shared state of its own (TrackingLoops'
`VectorTrackingLoop`) sees every satellite on one time grid. The in-place mutators
share it; [`track`](@ref), [`merge_sats`](@ref), [`add_satellite`](@ref) and
[`remove_satellite`](@ref) give their result a copy.

A fresh `TrackState` starts every band at `0`. To carry such a tracking loop over
into a rebuilt `TrackState` (e.g. to add a band), pass the old state's offsets
as the `sample_index_offsets` keyword, or the tracking loop receives integrations
from its past. A band the old state did not track starts at
`elapsed_time * sampling_frequency` samples of the new band.
"""
struct TrackState{G<:SignalGroups,TL<:AbstractTrackingLoop,NE<:NoiseEstimators}
    groups::G
    tracking_loop::TL
    noise_estimators::NE
    noise_descriptor::Base.RefValue{Any}
    sample_index_offsets::Dict{Symbol,Int}
end

# A fresh state gets an empty descriptor cell; derived states pass theirs on,
# which keeps the steady state allocation-free.
TrackState(
    groups::SignalGroups,
    tracking_loop::AbstractTrackingLoop,
    noise_estimators::NoiseEstimators;
    sample_index_offsets = nothing,
) = TrackState(
    groups,
    tracking_loop,
    noise_estimators,
    Base.RefValue{Any}(nothing),
    _initial_sample_index_offsets(groups, sample_index_offsets),
)

# Every band at `0`, overridden by the caller's `sample_index_offsets`, in a fresh `Dict`.
function _initial_sample_index_offsets(groups::SignalGroups, sample_index_offsets)
    offsets = Dict{Symbol,Int}(get_band_id(group.band) => 0 for group in groups)
    isnothing(sample_index_offsets) && return offsets
    for (band, offset) in pairs(sample_index_offsets)
        haskey(offsets, band) || throw(
            ArgumentError(
                string(
                    "`sample_index_offsets` names band `",
                    band,
                    "`, which no group of this TrackState tracks (bands: ",
                    join(keys(offsets), ", "),
                    ").",
                ),
            ),
        )
        offset >= 0 || throw(
            ArgumentError(
                "`sample_index_offsets` must be non-negative, got $band => $offset.",
            ),
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
include("tracking_loops.jl")
include("tracking_state.jl")
include("track.jl")
include("position_velocity_time.jl")
include("precompile.jl")

end
