"""
$(SIGNATURES)

Abstract supertype for per-band noise estimators — the source of the noise
**density** `N₀` that [`NoiseRefCN0Estimator`](@ref) divides each record's prompt
power by.

One instance is held per RF band in [`TrackState`](@ref)'s `noise_estimators`
NamedTuple (keyed by `GNSSSignals.get_band_id`), never per satellite and never
per signal: a band has exactly one noise floor, and averaging it once per band
is what makes the reference's own variance negligible against the per-record
prompt statistics.

# Why a density and not a power

For a correlator normalised the way [`normalize`](@ref) normalises — by
`integrated_samples * code_amplitude` — white input noise of per-sample variance
`σ²` gives `E|P|² = σ²/N = N₀/T`. So `N₀ = σ²/f_s` is **independent of the
integration time**, which is what lets one per-band figure serve records of any
length. It is stored as a Unitful quantity of dimension `1/Hz`, so the
consumer's `⟨|P|²⟩/N₀ − 1/T` is dimension-checked rather than trusted.

# The interface

Three methods, all with a default on this abstract type:

  - [`update_noise!`](@ref) — measure one band's slice of samples and append the
    resulting observations. This is the **software** fill path; it has exactly
    one call site, inside `downconvert_and_correlate!`.
  - [`append_noise_observation!`](@ref) — append an observation built elsewhere.
    This is the **hardware** fill path (FPGA/ASIC correlator or a front-end
    power monitor), parallel to [`append_correlator_output!`](@ref).
  - [`get_noise_density`](@ref) — the band's current density, or `nothing` while
    nothing has been measured yet. A read, not a drain: the window keeps
    sliding.

The two fill paths live on **disjoint call graphs** — a hardware producer never
calls `downconvert_and_correlate!` — so one concrete type serves both and no
type parameter distinguishes them. See [`CorrelatorNoiseEstimator`](@ref), the
only shipped implementation.

A source that must not measure from samples (a pure power-monitor integration,
say) subtypes this and gives `update_noise!` a no-op method; that is why the
interface is exported.
"""
abstract type AbstractNoiseEstimator end

"""
Type alias for a NamedTuple of [`AbstractNoiseEstimator`](@ref)s keyed by band id
— the shape [`TrackState`](@ref) holds them in. A NamedTuple and not a
`Dictionary`, because a dictionary would need an abstract value type as soon as
two bands hold different estimator types, which is type-unstable and would
allocate on every chunk.
"""
const NoiseEstimators = NamedTuple{<:Any,<:Tuple{Vararg{AbstractNoiseEstimator}}}

"""
$(SIGNATURES)

One noise measurement, the producer-side value parallel to
[`CorrelatorOutput`](@ref). It carries the already-normalised **density** plus
the two pieces of bookkeeping a sliding window needs, so the internal contract
is one scalar and every unit conversion happens at the boundary where the caller
knows their hardware.

Build one with [`noise_observation`](@ref),
[`noise_observation_from_correlator`](@ref) or
[`noise_observation_from_samples`](@ref) rather than by hand — they are the
three shapes a producer actually has, and all three reduce to the same `N₀` on
the same white input.

Fields:

  - `noise_density` — `N₀`, of dimension `1/Hz`.
  - `num_sub_integrations` — `M`, the number of **independent looks** that went
    into it, and therefore its statistical weight. Not the sample count: for
    `Σ_m |B_m|²` over `M` sub-integrations of `N` samples each, the density needs
    only `M·N`, but its relative variance is `1/M`. A producer handing back one
    1 ms coherent despread and one handing back 64 16-chip despreads report the
    same total samples and wildly different precision, so the window is averaged
    weighted by `M`.
  - `duration` — the wall-clock span of samples the observation covers, which is
    what bounds the sliding window. It is *not* `M · N / f_s` in general: the
    software source's `M` looks are the correlator's taps, which are
    simultaneous rather than consecutive.
  - `prn` — the PRN whose code was used to measure it. Legitimate metadata (the
    observation really was measured with that code), and it is what carries the
    software source's PRN-rotation position without any scalar state.
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

"""
$(SIGNATURES)

The concrete type `get_noise_density` returns for `estimator` once its window
holds anything. Used to keep the fold's density argument a plain scalar — never
a `Union{Nothing,…}` — so `update` stays devirtualised and allocation-free even
while the window is still empty.

Defaults to `typeof(1.0/1.0Hz)`, which every shipped builder produces.
"""
noise_density_type(::AbstractNoiseEstimator) = NoiseDensity

"""
$(SIGNATURES)

Build a [`NoiseObservation`](@ref) from **one dump given as a raw complex
accumulation** `Σ x·c` over `integrated_samples` samples — the most faithful
form a producer can report, since Tracking does the squaring, the scaling and
all the averaging itself.

`code_amplitude` is the RMS amplitude of the sampled replica (see
[`normalize`](@ref)); leave it at `1` for a ±1 code, pass the code's RMS for a
multi-level one (CBOC). `prn` records which code measured it.
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

"""
$(SIGNATURES)

Build a [`NoiseObservation`](@ref) from `M = num_sub_integrations` dumps
**pre-summed on chip** as `Σ_m |B_m|²`, covering `total_samples = M · N` samples
in all.

Pre-summing must be **incoherent** (`Σ|B_m|²`, never `|Σ B_m|²`), which is why
this builder takes a power rather than a complex value. `M > 1` exists only so a
high-rate producer can cut host traffic; the natural granularity is one
observation per dump with `M = 1`, which [`noise_observation`](@ref) covers.

`duration` defaults to `total_samples / sampling_frequency`, which is right when
the `M` sub-integrations are consecutive in time. Pass it explicitly when they
are **simultaneous** — the software source's `M` looks are the reference
correlator's taps, which all integrate the same `N` samples, so their span is
`N / sampling_frequency` and not `M` times that.
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

"""
$(SIGNATURES)

Build a [`NoiseObservation`](@ref) from a front-end / AGC **power monitor**:
`Σ|x|²` over `num_samples` raw samples.

This is [`noise_observation_from_correlator`](@ref) at a sub-integration of one
sample — the replica degenerates to a single ±1 value, `|x·c|² = |x|²`, and the
window mean is exactly `σ̂²` — so `M == num_samples` and no `code_amplitude`
applies. It is the minimum-variance, flat-weighting end of the same continuum
the despread source sits at the spectral-fidelity end of; a tilted front end is
the one thing that tells them apart.
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

# Shared core of the three builders: `N₀ = power / (total_samples · A_c² · f_s)`.
# Both the density and the duration are converted to their canonical units, so a
# band's window stays concretely typed no matter how the caller spelled the
# sampling frequency.
@inline function _noise_observation(
    accumulated_power,
    num_sub_integrations,
    total_samples,
    sampling_frequency,
    code_amplitude,
    prn,
    duration,
)
    density = uconvert(
        Hz^-1,
        accumulated_power / (total_samples * code_amplitude^2 * sampling_frequency),
    )
    NoiseObservation(density, Int(num_sub_integrations), uconvert(s, duration), Int16(prn))
end

"""
$(SIGNATURES)

Measure the noise on one band's samples and append the resulting observations to
`estimator`'s window, returning `estimator`.

`measurement` is the band's [`BandMeasurement`](@ref); `first_sample` and
`last_sample` bound the slice of it this call may consume (the current chunk).
`context` is a [`NoiseUpdateContext`](@ref), carrying what a source may need but
a `BandMeasurement` does not have — the band's reference signal, which PRNs are
currently tracked on it, and the chunk index.

This is the **software** fill path and it has exactly one call site, inside
`downconvert_and_correlate!`. A source whose observations arrive from outside
(see [`append_noise_observation!`](@ref)) never reaches it, so the default here
is a no-op: it mutates the window in place and leaves the struct alone, which is
what lets `TrackState` hold it immutably.
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
[`append_correlator_output!`](@ref) — see there for how the two differ.

The window is mutated in place and the struct is not rebuilt, so this is
allocation-free in steady state and works through the immutable `TrackState`.

The [`TrackState`](@ref) form selects the band:

```julia
append_noise_observation!(track_state, obs)         # single-band TrackState
append_noise_observation!(track_state, obs, :L1)    # explicit band id
```
"""
append_noise_observation!(estimator::AbstractNoiseEstimator, ::NoiseObservation) = estimator

"""
$(SIGNATURES)

The band's current noise density `N₀`, or `nothing` while the window holds
nothing yet.

A **read, not a drain**: the window keeps sliding across chunks and across
`track!` calls, so every record of every satellite on the band divides by the
same figure. The returned quantity has dimension `1/Hz` — see
[`AbstractNoiseEstimator`](@ref) for why a density rather than a power.
"""
get_noise_density(::AbstractNoiseEstimator) = nothing

"""
$(SIGNATURES)

Per-call side information handed to [`update_noise!`](@ref): what a software
noise source may need and a [`BandMeasurement`](@ref) does not carry.

Fields:

  - `signal` — the band's **reference signal**: the first signal of the first
    group on that band, in `groups` order. Its primary code period sets the
    sub-integration length and its code family is the one the reference
    despreads with. Deterministic, needs no configuration, and for white noise
    the choice is immaterial (`N̂₀ = σ²/f_s` is code-independent once `A_c` is
    divided out) — it only selects *which* spectral weighting a tilted front end
    is measured with.
  - `tracked_prn_sets` — a tuple of the key sets of the band's groups'
    satellite dictionaries, so a source can pick a PRN nobody is tracking
    without allocating a merged set.
  - `chunk_index` — the index of the chunk being measured, on the same grid
    `downconvert_and_correlate!` uses.
  - `downconvert_and_correlator` — the backend running this call. A software
    source **must** despread on it rather than on a kernel of its own choosing:
    the one- and two-bit accumulators are popcount counts rather than sample
    sums, so a float-kernel reference would compare `|P|²` and `N̂₀` on
    incompatible scales, and the quantisation loss would drop out of the
    measurement instead of being carried by it.

Kept as one struct so that adding a field later is not a signature change.
"""
struct NoiseUpdateContext{
    S<:AbstractGNSSSignal,
    K<:Tuple,
    DC<:AbstractDownconvertAndCorrelator,
}
    signal::S
    tracked_prn_sets::K
    chunk_index::Int
    downconvert_and_correlator::DC
end

# Is `prn` tracked by any of the band's groups? Tuple recursion so the walk over
# a heterogeneous tuple of key sets stays concrete and allocation-free.
@inline _is_prn_tracked(::Tuple{}, prn) = false
@inline _is_prn_tracked(sets::Tuple, prn) =
    (prn in first(sets)) || _is_prn_tracked(Base.tail(sets), prn)

@inline _is_prn_tracked(context::NoiseUpdateContext, prn) =
    _is_prn_tracked(context.tracked_prn_sets, prn)

# The band's density as a plain scalar plus a "is it meaningful yet" flag —
# the union-free form the fold threads down to `_apply_correlator_output`.
# `get_noise_density`'s `Union{Nothing,D}` is split here, once per band per
# chunk, so that everything below it stays monomorphic.
@inline function _noise_density_and_ready(estimator::AbstractNoiseEstimator)
    density = get_noise_density(estimator)
    D = noise_density_type(estimator)
    isnothing(density) && return (zero(D), false)
    (density::D, true)
end
