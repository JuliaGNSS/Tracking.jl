"""
$(SIGNATURES)

Abstract post correlation filter for the prompt signal — the linear combiner that
reduces a correlator's per-antenna taps to the single complex value the
discriminators and the C/N₀ estimator see.

An implementation must provide two methods:

  - [`get_weights`](@ref)`(filter, ::NumAnts{M})` — the combining weights `w` this
    filter currently applies. A plain `ComplexF64` at `M == 1`, an
    `SVector{M,ComplexF64}` above.
  - [`update`](@ref)`(filter, prompt)` — return the filter to use for the next
    record, which is where an adaptive beamformer recomputes its weights.

The contract is deliberately **linear in the antennas**: the combined tap is
`wᴴ·b`. That is what lets the C/N₀ path stay honest without any per-satellite
noise state. The noise reference measures one spatial covariance `R̂` per signal
(see [`AbstractNoiseEstimator`](@ref)), and each satellite reduces it to its own
scalar floor through its own weights, `N₀ = wᴴR̂w` — exact for any fixed `w`,
because `E[|wᴴn|²] = wᴴRw`. A filter that combined its antennas non-linearly
would have no such closed form, and its C/N₀ would silently be a ratio of two
different noise scales.
"""
abstract type AbstractPostCorrFilter end

"""
$(SIGNATURES)

This is the default post correlation filter. For a single antenna channel it
returns the prompt value as is. For multi antenna systems it selects the **last**
antenna's value — the noise reference's covariance reduces to that antenna's own
floor under these weights, so the two sides of the C/N₀ ratio agree.
"""
struct DefaultPostCorrFilter <: AbstractPostCorrFilter end

"""
$(SIGNATURES)

Update the post-correlation filter with the latest prompt and return the
new filter (immutable update). This is the extension point for custom
[`AbstractPostCorrFilter`](@ref)s — e.g. a beamformer adapting its weights.
"""
update(filter::DefaultPostCorrFilter, prompt) = filter

"""
$(SIGNATURES)

The linear combining weights `w` that `filter` currently applies to a correlator
tap, as a `ComplexF64` for `NumAnts(1)` and an `SVector{M,ComplexF64}` for
`NumAnts(M)`. The combined tap is `wᴴ·b`, and the matching post-combination noise
floor is `wᴴR̂w` — see [`AbstractPostCorrFilter`](@ref).

Required for every [`AbstractPostCorrFilter`](@ref); there is no fallback,
because a filter whose weights are unknown cannot be given a correct C/N₀.

# Examples

A beamformer that averages the antenna elements:

```julia
struct MyBeamformer <: AbstractPostCorrFilter end
Tracking.update(f::MyBeamformer, prompt) = f
Tracking.get_weights(::MyBeamformer, ::NumAnts{1}) = 1.0 + 0.0im
Tracking.get_weights(::MyBeamformer, ::NumAnts{M}) where {M} =
    SVector{M,ComplexF64}(ntuple(_ -> 1 / M + 0.0im, M))
```
"""
function get_weights end

# The scalar case is a plain `ComplexF64` rather than an `SVector{1}`, because a
# single-antenna correlator's accumulators are plain `ComplexF64` too (see
# `type_for_num_ants`). `one(ComplexF64)` is exact, so `_combine_antennas` and
# `_reduce_noise_density` are both bit-identical no-ops at `M == 1`.
get_weights(::DefaultPostCorrFilter, ::NumAnts{1}) = one(ComplexF64)
get_weights(::DefaultPostCorrFilter, ::NumAnts{M}) where {M} =
    SVector{M,ComplexF64}(ntuple(i -> i == M ? one(ComplexF64) : zero(ComplexF64), M))

# Combine one correlator tap through the weights: `wᴴ·b`. The conjugation is on
# `w`, matching `_reduce_noise_density`'s `wᴴRw` — the two must agree or the
# C/N₀ ratio's numerator and denominator would use different combiners.
@inline _combine_antennas(w::Number, tap::Number) = conj(w) * tap
@inline _combine_antennas(w::StaticVector, tap::StaticVector) = w'tap

# The same, for every tap of a correlator at once: an `M`-antenna correlator
# reduced to the single-channel one the discriminators read.
#
# Its own function, rather than the closure written inline at the call site, and
# that is load-bearing rather than tidy: `_apply_correlator_output` is large
# enough that a closure *capturing* the weights is built on the heap there, which
# costs an allocation per record on the one path that must have none (see the
# `track!` allocation guards in `test/track_in_place.jl`). Built here, in a
# function small enough to inline, the capture stays in registers.
@inline _combine_correlator(correlator, weights) = update_accumulator(
    correlator,
    map(tap -> _combine_antennas(weights, tap), get_accumulators(correlator)),
)

# Reduce a measured noise floor to the scalar one that *this* combiner sees.
#
# `nothing` means no noise estimator is configured for the signal at all, a
# static property of the setup: it must pass straight through so the C/N₀
# context's type parameter stays `Nothing` and `NoiseRefCN0Estimator` throws the
# wiring-mistake error at the first record (see `_signal_noise_density`).
#
# The scalar case is `|w|²·N₀`, which for `DefaultPostCorrFilter`'s exact
# `1.0+0.0im` is bit-identical to the unreduced density. The matrix case is
# `wᴴRw`, the exact post-combination floor for any fixed `w` — no `‖w‖²` factor
# on top, because the covariance already carries the antennas' relative scale.
@inline _reduce_noise_density(::Nothing, w) = nothing
@inline _reduce_noise_density(density::Number, w::Number) = abs2(w) * density
@inline _reduce_noise_density(R::StaticMatrix, w::StaticVector) = real(w' * R * w)
