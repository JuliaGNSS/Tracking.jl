"""
$(SIGNATURES)

A CN0 estimator that estimates nothing: it keeps no state, does no work per
record, and reports `-Inf dB-Hz`. Use it to say **"this signal's C/N₀ is not
measured"** rather than publish a number that cannot be trusted — for a
non-driver signal whose C/N₀ nobody reads, or as [`NWPRCN0Estimator`](@ref)'s
`fallback` in place of the [`MomentsCN0Estimator`](@ref)'s ~27.6 dB-Hz noise
floor where no coherent window exists.

`-Inf` and not `NaN` because, with Unitful's `Level` comparison,
`NaN dB-Hz >= threshold` is `true` for every threshold, so a `NaN` would clear
every lock detector. `-Inf` compares `false` against any finite threshold. This
is the package-wide convention for a missing estimate.
"""
struct NoCN0Estimator <: AbstractCN0Estimator end

"""
$(SIGNATURES)

Ignore the prompt and return the estimator unchanged — [`NoCN0Estimator`](@ref)
holds no state.
"""
update(estimator::NoCN0Estimator, prompt) = estimator

"""
$(SIGNATURES)

Always `-Inf dB-Hz`; see [`NoCN0Estimator`](@ref) for why that value and not
`NaN`.
"""
estimate_cn0(::NoCN0Estimator, integration_time) = dBHz(0.0 / integration_time)
