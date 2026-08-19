"""
$(SIGNATURES)

Abstract supertype for per-signal noise estimators — the source of the noise
**density** `N₀` that [`NoiseRefCN0Estimator`](@ref) divides each record's prompt
power by.

One instance is held per **signal** in [`TrackState`](@ref)'s `noise_estimators`
NamedTuple (keyed by `GNSSSignals.get_signal_id`), never per satellite: every
satellite of a signal shares one floor, and averaging it once per signal is what
makes the reference's own variance negligible against the per-record prompt
statistics.

# Why per signal and not per RF band

What a record actually divides by is the **post-correlation** density. For an
interferer of PSD `S_I(f)` on top of the thermal floor that is

```
N₀,eff = N₀,thermal + ∫ S_I(f) · |G(f)|² df
```

— the spectral separation coefficient, weighted by the *despreading modulation's*
own spectrum `|G(f)|²`. Two signals sharing one band, one antenna and one front
end therefore see different floors the moment the interference is not white, and
the difference is not small: BPSK(1) has its peak at DC and a null at
±1.023 MHz, BOC(1,1) is the reverse, so a CW tone at band centre is rejected by
Galileo E1B and lands squarely in GPS L1 C/A, while a tone at ±1.023 MHz does the
opposite. Front-end tilt, filter roll-off at the band edge and adjacent-band
leakage all colour the floor the same way, more quietly.

Because [`CorrelatorNoiseEstimator`](@ref) *despreads* rather than metering
power, keying by signal makes that integral **measured rather than modelled**:
the reference runs the consumer's own code, so its spectral weighting is the
consumer's by construction. This is the same argument that makes the reference
backend-free — it traverses the identical path as the prompt — extended from the
quantiser to the interference environment.

It also matches the hardware. A noise channel is a tracking channel with a wrong
PRN, and a tracking channel is configured with a code; per-signal is what an FPGA
would build anyway.

# Why a density and not a power

For a correlator normalised the way `Tracking.normalize` normalises — by
`integrated_samples * code_amplitude` — white input noise of per-sample variance
`σ²` gives `E|P|² = σ²/N = N₀/T`. So `N₀ = σ²/f_s` is **independent of the
integration time**, which is what lets one per-signal figure serve records of any
length. It is stored as a Unitful quantity of dimension `1/Hz`, so the
consumer's `⟨|P|²⟩/N₀ − 1/T` is dimension-checked rather than trusted.

# The interface

Three methods, all with a default on this abstract type:

  - [`update_noise!`](@ref) — measure one signal's slice of its band's samples
    and append the resulting observations. This is the **software** fill path; it
    has exactly one call site, inside `downconvert_and_correlate!`.
  - [`append_noise_observation!`](@ref) — append an observation built elsewhere.
    This is the **hardware** fill path (FPGA/ASIC correlator or a front-end
    power monitor), parallel to [`append_correlator_output!`](@ref).
  - [`get_noise_density`](@ref) — the signal's current density, or `nothing`
    while nothing has been measured yet. A read, not a drain: the window keeps
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
Type alias for a NamedTuple of [`AbstractNoiseEstimator`](@ref)s keyed by signal
id — the shape [`TrackState`](@ref) holds them in. A NamedTuple and not a
`Dictionary`, because a dictionary would need an abstract value type as soon as
two signals hold different estimator types, which is type-unstable and would
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

  - `noise_density` — `N₀`, of dimension `1/Hz`. For an antenna array it is the
    `M×M` spatial covariance `R̂` instead, of the same dimension elementwise: the
    diagonal is each antenna's own `N₀` and the off-diagonals their noise
    correlation. A satellite reduces it to its own scalar floor through its
    beamforming weights, `wᴴR̂w` (see [`AbstractPostCorrFilter`](@ref)).
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

# The same, per element of a multi-antenna window's spatial covariance `R̂`. The
# diagonal carries each antenna's own `N₀`; the off-diagonals carry the antennas'
# noise correlation, which is exactly what a beamformer's weights see.
const NoiseCovarianceElement = typeof((1.0 + 0.0im) / 1.0Hz)

# What a window measures, as a function of how many antennas it despreads.
#
# A bare `SMatrix`, not a `Hermitian` wrapper: `zero`, `+` and `/` are all the
# window's totals need, StaticArrays supplies them, and the result stays isbits —
# which is what keeps `NoiseWindowTotals` inline in its `Ref` and the whole
# update path allocation-free (see `NoiseWindowTotals`).
_density_type_for_num_ants(::NumAnts{1}) = NoiseDensity
_density_type_for_num_ants(::NumAnts{M}) where {M} =
    SMatrix{M,M,NoiseCovarianceElement,M * M}

# The inverse, read off a density *type*. Deliberately phrased on the type rather
# than on the estimator, so it answers for **any** `AbstractNoiseEstimator` —
# including a user's own hardware source — through `noise_density_type` alone.
# That is what lets `TrackState` check an explicitly-passed estimator's antenna
# count against its signal group's without knowing the estimator's concrete type.
_num_ants_of_density_type(::Type{NoiseDensity}) = NumAnts(1)
_num_ants_of_density_type(::Type{<:SMatrix{M,M}}) where {M} = NumAnts(M)

# Retype an observation onto a window's own `{D,T}`. The builders already emit the
# canonical pair, so this is for an observation assembled by hand: its field types
# are then whatever the caller's arithmetic produced, and a window that will not
# take it has no way to say so except by dropping it.
#
# The identity method is written out rather than left to Base's
# `convert(::Type{T}, ::T)`, because Base's does *not* win here — binding `{D,T}`
# in the first argument makes the retyping method specific enough to be chosen
# even when the second argument already matches. It reduces to a no-op either way
# (`NoiseObservation` is isbits and `convert(D, ::D)` is identity, so the append
# path measures zero bytes), but only with this method is "the steady state
# converts nothing" true of the dispatch rather than of the optimiser.
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
`Tracking.normalize`); leave it at `1` for a ±1 code, pass the code's RMS for a
multi-level one (CBOC). `prn` records which code measured it.

Pass an `SVector` of `M` per-antenna accumulations for an antenna array, and the
observation carries the spatial covariance `b·bᴴ` instead of the scalar `|b|²`.
Report every element the front end has: a satellite's floor is `wᴴR̂w` under its
own beamforming weights, and an array collapsed to one number cannot answer that.
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

Pre-summing must be **incoherent** (`Σ|B_m|²`, never `|Σ B_m|²`), which is why
this builder takes a power rather than a complex value. `M > 1` exists only so a
high-rate producer can cut host traffic; the natural granularity is one
observation per dump with `M = 1`, which [`noise_observation`](@ref) covers.

`duration` defaults to `total_samples / sampling_frequency`, which is right when
the `M` sub-integrations are consecutive in time. Pass it explicitly when they
are **simultaneous** — the software source's `M` looks are the reference
correlator's taps, which all integrate the same `N` samples, so their span is
`N / sampling_frequency` and not `M` times that.

For an antenna array, pass the pre-summed spatial covariance `Σ_m B_m·B_mᴴ` as an
`SMatrix`. The incoherence requirement is the same one term for term: sum the
outer products, never form the outer product of the sum.
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

This is [`noise_observation_from_correlator`](@ref) at a sub-integration of one
sample — the replica degenerates to a single ±1 value, `|x·c|² = |x|²`, and the
window mean is exactly `σ̂²` — so `M == num_samples` and no `code_amplitude`
applies. It is the minimum-variance, flat-weighting end of the same continuum
the despread source sits at the spectral-fidelity end of; a tilted front end is
the one thing that tells them apart.

For an antenna array, pass `Σ x·xᴴ` as an `SMatrix` — the array's raw spatial
covariance, which is what a beamformer's weights are reduced against.
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
# Both fields are `convert`ed to the canonical types — `NoiseDensity` and
# `typeof(1.0s)` — rather than merely `uconvert`ed to canonical *units*, so that a
# window's element type is fixed by this function and not by how the caller
# spelled their sampling frequency. `uconvert` alone preserves the input's number
# type, so a producer reporting `4f6Hz` built a `NoiseObservation{Quantity{Float32,
# …},…}`, which then matched no `append_noise_observation!` method for the
# `Float64` window and fell through to the abstract no-op — the observation was
# silently dropped and the window stayed empty forever. Canonicalising here is
# what makes the documented hardware fill path work for any spelling; the
# converting method on `CorrelatorNoiseEstimator` covers a hand-built observation
# that never came through a builder.
#
# The canonical type is chosen by dispatching on the measured value itself — a
# scalar power lands on `NoiseDensity`, an `M×M` covariance on the matching
# `SMatrix` — rather than by threading the type in as an argument. That is not a
# style preference: a `Type` passed by value is opaque to Julia 1.10's optimiser,
# and inside `update_noise!`'s per-sub-integration loop it was enough to stop the
# whole frame being optimised, boxing isbits temporaries several calls deeper.
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

`measurement` is the band's [`BandMeasurement`](@ref) — the samples are a band
property, one front end feeding every signal on it; only the despreading code,
and therefore the measured floor, is per signal. `first_sample` and `last_sample`
bound the slice of it this call may consume (the current chunk). `context` is a
[`NoiseUpdateContext`](@ref), carrying what a source may need but a
`BandMeasurement` does not have — the signal being measured, which PRNs are
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

A **read, not a drain**: the window keeps sliding across chunks and across
`track!` calls, so every record of every satellite tracking that signal divides
by the same figure. The returned quantity has dimension `1/Hz` — see
[`AbstractNoiseEstimator`](@ref) for why a density rather than a power, and why
per signal rather than per band.
"""
get_noise_density(::AbstractNoiseEstimator) = nothing

"""
$(SIGNATURES)

Per-call side information handed to [`update_noise!`](@ref): what a software
noise source may need and a [`BandMeasurement`](@ref) does not carry.

Fields:

  - `signal` — **the signal being measured**, i.e. the one this estimator's key
    names. Its primary code period sets the sub-integration length and its code
    family is the one the reference despreads with, so the measurement carries
    that signal's own spectral weighting — which is the point of keying the
    estimator by signal (see [`AbstractNoiseEstimator`](@ref)). Nothing here is
    chosen: it is the consumer's signal, not a per-band stand-in for it.
  - `chunk_index` — the index of the chunk being measured, on the same grid
    `downconvert_and_correlate!` uses.
  - `downconvert_and_correlator` — the backend running this call. A software
    source **must** despread on it rather than on a kernel of its own choosing:
    the one- and two-bit accumulators are popcount counts rather than sample
    sums, so a float-kernel reference would compare `|P|²` and `N̂₀` on
    incompatible scales, and the quantisation loss would drop out of the
    measurement instead of being carried by it.

Kept as one struct so that adding a field later is not a signature change.

Note what is **not** here: which PRNs are currently tracked. The reference used
to need it, to skip tracked PRNs and so keep a fixed code phase 0 off their
correlation peaks. [`CorrelatorNoiseEstimator`](@ref) now randomises the phase
instead, so it rotates over the whole family and reads nothing from satellite
state — which is what makes "open-loop" literal rather than approximate.
"""
struct NoiseUpdateContext{S<:AbstractGNSSSignal,DC<:AbstractDownconvertAndCorrelator}
    signal::S
    chunk_index::Int
    downconvert_and_correlator::DC
end

# The signal's density as a plain scalar plus a "is it meaningful yet" flag —
# the union-free form the fold threads down to `_apply_correlator_output`.
# `get_noise_density`'s `Union{Nothing,D}` is split here, once per signal per
# chunk, so that everything below it stays monomorphic.
#
# "Meaningful" is more than "present". A measured floor of zero is not a floor to
# divide by, and it is reachable without any misuse: a front-end dropout or a
# buffer underrun delivers all-zero samples, every tap despreads to zero and the
# window fills with `N̂₀ = 0`. `|P|²/0` is then `Inf`, or — since that same buffer
# gives a zero prompt — `NaN`, and a `NaN dB-Hz` compares `>=` **true** against
# every lock threshold, so a dead signal would read as locked. A non-finite
# density (from a producer's own arithmetic via `append_noise_observation!`) is
# rejected for the same reason. Not-ready is the honest answer and needs no new
# machinery: the fold already skips the update and warns, and `estimate_cn0`
# already reports `-Inf dB-Hz` — the house convention for a missing estimate.
#
# Both tests are dispatched, because a multi-antenna window's density is an
# `SMatrix` and neither `isfinite` nor `>` is defined on one.
@inline _finite_density(d::Number) = isfinite(d)
@inline _finite_density(R::StaticMatrix) = all(isfinite, R)

# "Positive" for a covariance means **each antenna measured some power**, i.e. the
# diagonal is positive — deliberately not positive definiteness. A near-singular
# `R̂` is legitimate, not a fault: strongly correlated antennas (in the limit, one
# signal fed to every element) produce a rank-1 covariance that is still the
# right answer for every `wᴴR̂w` anyone will ask of it. `R̂` is never inverted, so
# rank never matters; only a genuinely dead element does, and that shows on the
# diagonal.
@inline _positive_density(d::Number) = d > zero(d)
@inline _positive_density(R::StaticMatrix) =
    all(i -> real(R[i, i]) > zero(real(R[i, i])), axes(R, 1))

@inline function _noise_density_and_ready(estimator::AbstractNoiseEstimator)
    density = get_noise_density(estimator)
    D = noise_density_type(estimator)
    isnothing(density) && return (zero(D), false)
    d = density::D
    _finite_density(d) && _positive_density(d) || return (zero(D), false)
    (d, true)
end
