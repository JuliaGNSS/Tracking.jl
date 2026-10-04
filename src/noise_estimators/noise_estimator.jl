"""
$(SIGNATURES)

Abstract supertype for per-signal noise estimators — the source of the noise
**density** `N₀` that [`NoiseRefCN0Estimator`](@ref) divides each record's prompt
power by.

One instance is held per **signal** (not per satellite or RF band) in
[`TrackState`](@ref)'s `noise_estimators`; see the Noise Estimator docs page for
why it is per signal and a density.

# The interface

Three required methods, all with a default on this abstract type:

  - [`update_noise!`](@ref) — the **software** fill path, called only inside
    `downconvert_and_correlate!`.
  - [`append_noise_observation!`](@ref) — the **hardware** fill path, parallel to
    [`append_correlator_output!`](@ref).
  - [`get_noise_density`](@ref) — the signal's current density, or `nothing`
    while nothing has been measured yet.

A fourth, `Tracking.noise_window_looks`, is optional and only matters to a source
that reports a **covariance** (see `noise_window_looks`).

The two fill paths live on disjoint call graphs, so one concrete type,
[`CorrelatorNoiseEstimator`](@ref), serves both.
"""
abstract type AbstractNoiseEstimator end

"""
Type alias for a NamedTuple of [`AbstractNoiseEstimator`](@ref)s keyed by signal
id — the shape [`TrackState`](@ref) holds them in. Not a `Dictionary`, whose
value type would turn abstract (and allocate per chunk) once two signals hold
different estimator types.
"""
const NoiseEstimators = NamedTuple{<:Any,<:Tuple{Vararg{AbstractNoiseEstimator}}}

"""
$(SIGNATURES)

One noise measurement, the producer-side value parallel to
[`CorrelatorOutput`](@ref): the already-normalised **density** plus the
bookkeeping a sliding window needs. Build one with [`noise_observation`](@ref),
[`noise_observation_from_correlator`](@ref) or
[`noise_observation_from_samples`](@ref) rather than by hand.

Fields:

  - `noise_density` — `N₀`, of dimension `1/Hz`. For an antenna array it is the
    `M×M` spatial covariance `R̂` instead, of the same dimension elementwise: the
    diagonal is each antenna's own `N₀` and the off-diagonals their noise
    correlation. A satellite reduces it to its own scalar floor through its
    beamforming weights, `wᴴR̂w` (see [`AbstractPostCorrFilter`](@ref)).
  - `num_sub_integrations` — `M`, the number of **independent looks**, and so
    the statistical weight (relative variance `1/M`). Not the sample count; see
    docs/src/noise_estimator.md.
  - `duration` — the wall-clock span of samples covered, which bounds the
    sliding window. Not `M · N / f_s` in general: the software source's looks
    are simultaneous taps.
  - `prn` — the PRN whose code measured it; also carries the software source's
    PRN-rotation position (see `_next_noise_prn`).
"""
struct NoiseObservation{D,T}
    noise_density::D
    num_sub_integrations::Int
    duration::T
    prn::Int16
end

# Canonical density type. Every builder converts to it, so a band's window is
# concretely typed however the caller spelled their sampling frequency.
const NoiseDensity = typeof(1.0 / 1.0Hz)

# The same, per element of a multi-antenna window's spatial covariance `R̂`. The
# diagonal carries each antenna's own `N₀`; the off-diagonals carry the antennas'
# noise correlation, which is exactly what a beamformer's weights see.
const NoiseCovarianceElement = typeof((1.0 + 0.0im) / 1.0Hz)

# What a window measures, as a function of how many antennas it despreads. A bare
# `SMatrix`, not `Hermitian`, so it stays isbits and `NoiseWindowTotals` stays
# inline in its `Ref` (allocation-free).
_density_type_for_num_ants(::NumAnts{1}) = NoiseDensity
_density_type_for_num_ants(::NumAnts{M}) where {M} =
    SMatrix{M,M,NoiseCovarianceElement,M * M}

# The inverse, on the density *type*, so `TrackState` can check any
# `AbstractNoiseEstimator`'s antenna count through `noise_density_type` alone.
_num_ants_of_density_type(::Type{NoiseDensity}) = NumAnts(1)
_num_ants_of_density_type(::Type{<:SMatrix{M,M}}) where {M} = NumAnts(M)

# Retype an observation onto a window's own `{D,T}`, for one assembled by hand
# (the builders already emit the canonical pair). The identity method is explicit
# because the retyping method below is more specific than Base's
# `convert(::Type{T}, ::T)`; with it, the steady state converts nothing by dispatch.
Base.convert(::Type{NoiseObservation{D,T}}, o::NoiseObservation{D,T}) where {D,T} = o
Base.convert(::Type{NoiseObservation{D,T}}, o::NoiseObservation) where {D,T} =
    NoiseObservation{D,T}(
        convert(D, o.noise_density),
        o.num_sub_integrations,
        convert(T, o.duration),
        o.prn,
    )

"""
$(SIGNATURES)

The concrete type `get_noise_density` returns for `estimator` once its window
holds anything, so the fold's density argument stays a plain scalar rather than
a `Union{Nothing,…}`.

Defaults to `typeof(1.0/1.0Hz)`, which every shipped builder produces.
"""
noise_density_type(::AbstractNoiseEstimator) = NoiseDensity

"""
$(SIGNATURES)

Build a [`NoiseObservation`](@ref) from **one dump given as a raw complex
accumulation** `Σ x·c` over `integrated_samples` samples — the most faithful
form, since Tracking does the squaring, scaling and averaging.

`code_amplitude` is the RMS amplitude of the sampled replica (see
`Tracking.normalize`); leave it at `1` for a ±1 code, pass the code's RMS for a
multi-level one (CBOC). `prn` records which code measured it.

For an antenna array pass an `SVector` of `M` per-antenna accumulations; the
observation then carries the spatial covariance `b·bᴴ`. Report every element:
each satellite's floor is `wᴴR̂w` under its own weights.
"""
noise_observation(
    accumulation::Complex,
    integrated_samples,
    sampling_frequency;
    code_amplitude = 1,
    prn::Integer = 0,
) = _noise_observation(
    abs2(accumulation),
    1,
    integrated_samples,
    sampling_frequency,
    code_amplitude,
    prn,
    integrated_samples / sampling_frequency,
)

noise_observation(
    accumulation::StaticVector{M,<:Complex},
    integrated_samples,
    sampling_frequency;
    code_amplitude = 1,
    prn::Integer = 0,
) where {M} = _noise_observation(
    accumulation * accumulation',
    1,
    integrated_samples,
    sampling_frequency,
    code_amplitude,
    prn,
    integrated_samples / sampling_frequency,
)

"""
$(SIGNATURES)

Build a [`NoiseObservation`](@ref) from `M = num_sub_integrations` dumps
**pre-summed on chip** as `Σ_m |B_m|²`, covering `total_samples = M · N` samples
in all.

Pre-summing must be **incoherent** (`Σ|B_m|²`, never `|Σ B_m|²`). `M > 1` only
cuts host traffic; one observation per dump is what [`noise_observation`](@ref)
covers.

`duration` defaults to `total_samples / sampling_frequency`, right for
consecutive sub-integrations. Pass it explicitly when they are **simultaneous**
(e.g. correlator taps over the same `N` samples: `N / sampling_frequency`).

For an antenna array, pass `Σ_m B_m·B_mᴴ` as an `SMatrix` (sum the outer
products, never the outer product of the sum).
"""
noise_observation_from_correlator(
    accumulated_power,
    num_sub_integrations,
    total_samples,
    sampling_frequency;
    code_amplitude = 1,
    prn::Integer = 0,
    duration = total_samples / sampling_frequency,
) = _noise_observation(
    accumulated_power,
    num_sub_integrations,
    total_samples,
    sampling_frequency,
    code_amplitude,
    prn,
    duration,
)

noise_observation_from_correlator(
    accumulated_power::StaticMatrix{M,M},
    num_sub_integrations,
    total_samples,
    sampling_frequency;
    code_amplitude = 1,
    prn::Integer = 0,
    duration = total_samples / sampling_frequency,
) where {M} = _noise_observation(
    accumulated_power,
    num_sub_integrations,
    total_samples,
    sampling_frequency,
    code_amplitude,
    prn,
    duration,
)

"""
$(SIGNATURES)

Build a [`NoiseObservation`](@ref) from a front-end / AGC **power monitor**:
`Σ|x|²` over `num_samples` raw samples.

This is [`noise_observation_from_correlator`](@ref) at a one-sample
sub-integration, so `M == num_samples` and no `code_amplitude` applies. Unlike a
despread it weights the spectrum flat, so it misses coloured interference.

For an antenna array, pass `Σ x·xᴴ` as an `SMatrix`.
"""
noise_observation_from_samples(
    accumulated_power,
    num_samples,
    sampling_frequency;
    prn::Integer = 0,
) = _noise_observation(
    accumulated_power,
    num_samples,
    num_samples,
    sampling_frequency,
    1,
    prn,
    num_samples / sampling_frequency,
)

noise_observation_from_samples(
    accumulated_power::StaticMatrix{M,M},
    num_samples,
    sampling_frequency;
    prn::Integer = 0,
) where {M} = _noise_observation(
    accumulated_power,
    num_samples,
    num_samples,
    sampling_frequency,
    1,
    prn,
    num_samples / sampling_frequency,
)

# Shared core of the three builders: `N₀ = power / (total_samples · A_c² · f_s)`.
#
# Both fields are `convert`ed to the canonical types (not just units), so a
# window's element type does not depend on how the caller spelled their sampling
# frequency; a `Float32` observation would otherwise miss the window's method and
# be dropped silently.
#
# The canonical type is picked by dispatching on the value, not by passing a
# `Type` argument: on Julia 1.10 that defeated optimisation of `update_noise!`'s
# loop and boxed temporaries.
@inline _canonical_density(density::Number) = convert(NoiseDensity, density)
@inline _canonical_density(density::StaticMatrix{M,M}) where {M} =
    convert(SMatrix{M,M,NoiseCovarianceElement,M * M}, density)

@inline function _noise_observation(
    accumulated_power,
    num_sub_integrations,
    total_samples,
    sampling_frequency,
    code_amplitude,
    prn,
    duration,
)
    density = _canonical_density(
        accumulated_power / (total_samples * code_amplitude^2 * sampling_frequency),
    )
    NoiseObservation(
        density,
        Int(num_sub_integrations),
        convert(typeof(1.0s), duration),
        Int16(prn),
    )
end

"""
$(SIGNATURES)

Measure one signal's noise on its band's samples and append the resulting
observations to `estimator`'s window, returning `estimator`.

`measurement` is the band's [`BandMeasurement`](@ref) (the samples are per band;
only the despreading code is per signal). `first_sample:last_sample` is the slice
this call may consume. `context` is a [`NoiseUpdateContext`](@ref).

This is the **software** fill path, called only inside
`downconvert_and_correlate!`. The default is a no-op, for sources fed from
outside. Implementations mutate the window in place and return the same struct,
since `TrackState` is immutable.
"""
update_noise!(
    estimator::AbstractNoiseEstimator,
    measurement,
    first_sample,
    last_sample,
    context,
) = estimator

"""
$(SIGNATURES)

Append one externally built [`NoiseObservation`](@ref) to `estimator`'s sliding
window and return `estimator`. This is the **hardware** fill path, parallel to
[`append_correlator_output!`](@ref). The window is mutated in place, so this
works through the immutable `TrackState`.

The [`TrackState`](@ref) form selects the signal:

```julia
append_noise_observation!(track_state, obs)             # single-signal TrackState
append_noise_observation!(track_state, obs, :GPSL1CA)   # explicit signal id
append_noise_observation!(track_state, obs, GPSL1CA)    # or the signal type
```
"""
append_noise_observation!(estimator::AbstractNoiseEstimator, ::NoiseObservation) = estimator

"""
$(SIGNATURES)

The signal's current noise density `N₀`, or `nothing` while the window holds
nothing yet.

A **read, not a drain**: the window keeps sliding across chunks and `track!`
calls. The returned quantity has dimension `1/Hz` (see
[`AbstractNoiseEstimator`](@ref)).
"""
get_noise_density(::AbstractNoiseEstimator) = nothing

"""
$(SIGNATURES)

Per-call side information handed to [`update_noise!`](@ref): what a software
noise source may need and a [`BandMeasurement`](@ref) does not carry.

Fields:

  - `signal` — **the signal being measured** (the estimator's key). Its code
    period sets the sub-integration length and its code family is what the
    reference despreads with, so the measurement carries that signal's spectral
    weighting.
  - `chunk_index` — the index of the chunk being measured, on the same grid
    `downconvert_and_correlate!` uses.
  - `downconvert_and_correlator` — the backend running this call. A software
    source **must** despread on it: the one- and two-bit accumulators are
    popcount counts, so another kernel would put `|P|²` and `N̂₀` on different
    scales and drop the quantisation loss from the measurement.

One struct, so adding a field later is not a signature change. It deliberately
carries no satellite state; see [`CorrelatorNoiseEstimator`](@ref).
"""
struct NoiseUpdateContext{S<:AbstractGNSSSignal,DC<:AbstractDownconvertAndCorrelator}
    signal::S
    chunk_index::Int
    downconvert_and_correlator::DC
end

# The signal's density as a plain scalar plus a "ready" flag, splitting
# `get_noise_density`'s `Union{Nothing,D}` once per signal per chunk so everything
# below stays monomorphic. Implemented by `_noise_density_and_ready` below.
#
# "Ready" also requires a finite, positive floor. A zero floor is reachable (a
# front-end dropout gives all-zero samples) and would give a `NaN` C/N₀, which
# clears every lock threshold (see `NoCN0Estimator`). Not ready means the fold
# skips the update and `estimate_cn0` reports `-Inf dB-Hz`.
#
# Dispatched, because neither `isfinite` nor `>` is defined on an `SMatrix`.
@inline _finite_density(d::Number) = isfinite(d)
@inline _finite_density(R::StaticMatrix) = all(isfinite, R)

# "Positive" for a covariance means a positive diagonal (each antenna measured
# some power), not positive definiteness: `R̂` is never inverted, and strongly
# correlated antennas legitimately give a near-singular one.
@inline _positive_density(d::Number) = d > zero(d)
@inline _positive_density(R::StaticMatrix) =
    all(i -> real(R[i, i]) > zero(real(R[i, i])), axes(R, 1))

"""
$(SIGNATURES)

How many independent looks the estimator's current density is averaged from, or
`nothing` if the source does not report it.

Only the rank gate reads this: an `M×M` covariance averaged from fewer than `M`
looks is rank-deficient by construction, so the fold withholds it until there are
enough. A source returning `nothing` is never gated.
"""
noise_window_looks(::AbstractNoiseEstimator) = nothing

# Is the estimate built from enough looks to span its own dimensions? Scalars
# never gate. The software reference pools three rank-1 outer products per
# observation, so at `M > 3` the first `⌈M/3⌉` observations are held back.
# Without this, weights near the unmeasured subspace read a floor far too low
# (down to 2 % at `M = 4`, ≈16 dB optimistic), with the diagonal still positive.
@inline _sufficient_looks(::Number, looks) = true
@inline _sufficient_looks(::StaticMatrix, ::Nothing) = true
@inline _sufficient_looks(::StaticMatrix{M,M}, looks::Integer) where {M} = looks >= M

@inline function _noise_density_and_ready(estimator::AbstractNoiseEstimator)
    density = get_noise_density(estimator)
    D = noise_density_type(estimator)
    isnothing(density) && return (zero(D), false)
    d = density::D
    _finite_density(d) &&
    _positive_density(d) &&
    _sufficient_looks(d, noise_window_looks(estimator)) || return (zero(D), false)
    (d, true)
end

# Not ready *only* because the look count is short (window still filling), so a
# zero floor still reaches the fold's warning. Used only to keep that warning
# quiet during a normal multi-antenna startup.
@inline function _noise_window_filling(estimator::AbstractNoiseEstimator)
    density = get_noise_density(estimator)
    isnothing(density) && return false
    d = density::noise_density_type(estimator)
    _finite_density(d) &&
        _positive_density(d) &&
        !_sufficient_looks(d, noise_window_looks(estimator))
end
