"""
$(SIGNATURES)

The band's noise reference, measured by **despreading an untracked PRN** — the
only [`AbstractNoiseEstimator`](@ref) Tracking ships, and the one to configure on
a hardware path too (you simply fill it with
[`append_noise_observation!`](@ref) instead of letting [`update_noise!`](@ref)
fill it).

# Why a despread and not a power meter

The reference traverses the **identical** quantise → downconvert → despread →
accumulate path as the prompt, reusing the same kernel and the same replica
generator. So the measured `N₀` already contains every imperfection of that
chain — the 1-bit / 2-bit quantisation loss, the quantiser's operating point
under load, the input scaling an AGC step moves, the code amplitude of a CBOC
replica — with no per-backend model and no closed-form correction. A `Σ|x|²`
power meter would need one line per backend and would still miss the
second-order coupling where a strong signal shifts that operating point.

# It is open-loop

The reference has **no feedback of any kind**: it runs at the band's nominal IF
with zero Doppler and an off-peak code phase from a rotating PRN, so there is no
discriminator, no loop filter and no NCO update ever written to it. It is
therefore correct with **zero satellites tracked**, and on an FPGA a noise
channel is a strict subset of a tracking channel rather than an addition to one.

# Fields / configuration

  - `window_duration` — how far back the sliding window reaches (1 s by
    default). The one genuine tunable, because it trades: longer means lower
    variance, but a longer smear across AGC changes. At the default and one
    sub-integration per code period the window holds `K_n ≈ 1000` looks per tap,
    which costs ≤0.08 dB against a variance-free reference at every C/N₀ — below
    which no further tuning is warranted.
  - `tap_code_shift` — the reference correlator's tap spacing in chips (1.5 by
    default). Taps are only worth having if they are statistically independent,
    and at the *tracking* default of ±0.5 chip they are not: for jointly circular
    Gaussian accumulations `corr(|Bᵢ|²,|Bⱼ|²) = |ρᵢⱼ|²`, and the sampled-code
    correlation at half a chip is ≈0.5, so three taps are worth 2.25 independent
    looks. At ≥1 chip they are worth 2.98. 1.5 chips sits in the autocorrelation
    null and clear of both the 1- and 2-chip sidelobe values.
  - `buffered` — the sliding window itself, a **length-managed FIFO** written in
    place. The `Vector`'s own length is the position, so there is no ring index
    to write back and the struct is never rebuilt — which is what lets per-band
    state live in an immutable [`TrackState`](@ref).
  - `code_replica` — scratch for the software despread, grown once and reused.
    Held here rather than borrowed from the backend's `ScratchBuffers` so the
    reference can never alias the per-signal replica path.

The sub-integration length is deliberately **not** a field: it is the primary
code period of the band's reference signal, derived rather than configured.
Coherent integration buys nothing for noise-power estimation — one dump carries
100 % relative error however long it is, because `Var(|B|²) = (E|B|²)²` — so all
the information lives in the *number* of looks, and shortening the dump is the
only way to buy more of them. Three things say not to: SIMD (64 kernel calls per
1 ms chunk instead of one, and a shorter dump would need a new kernel variant,
forfeiting the bit-identical-arithmetic property the whole approach rests on),
DC balance (a full-period Gold code is balanced, `E[(Σc)²] ≈ 1`; any sub-period
window behaves like iid signs, ≈16 at 16 chips, so an ADC offset biases a short
despread and not a full-period one), and cadence parity with a hardware
correlator, which naturally dumps on the code epoch. The variance is bought back
with `window_duration` instead, which is nearly free — a 1 s window is ~1000
`Float64` densities per band.
"""
struct CorrelatorNoiseEstimator{D,T} <: AbstractNoiseEstimator
    window_duration::typeof(1.0s)
    tap_code_shift::Float64
    buffered::Vector{NoiseObservation{D,T}}
    code_replica::Vector{Int8}
end

"""
$(SIGNATURES)

Construct a [`CorrelatorNoiseEstimator`](@ref) averaging over the last
`window_duration` of observations, with the reference correlator's taps spaced
`tap_code_shift` chips apart. See the type's docstring for what the two
parameters buy and why nothing else is configurable.

The window is `sizehint!`-ed to four times the number of code-period
observations it expects to hold. The headroom is what makes the FIFO's
`push!`/`popfirst!` pair measure **exactly** zero bytes: at 1× or 2× Julia
periodically shifts the front offset back and reallocates.
"""
function CorrelatorNoiseEstimator(;
    window_duration = 1.0s,
    tap_code_shift = 1.5,
    num_ants::NumAnts = NumAnts(1),
)
    window_duration > zero(window_duration) ||
        throw(ArgumentError("window_duration must be positive, got $window_duration"))
    tap_code_shift > 0 ||
        throw(ArgumentError("tap_code_shift must be positive, got $tap_code_shift"))
    buffered = NoiseObservation{NoiseDensity,typeof(1.0s)}[]
    # Four times the count a 1 ms sub-integration would put in the window; see
    # the docstring for why the headroom rather than an exact fit.
    sizehint!(buffered, 4 * max(1, round(Int, window_duration / 1.0ms)) + 1)
    CorrelatorNoiseEstimator(
        uconvert(s, float(window_duration)),
        Float64(tap_code_shift),
        buffered,
        Int8[],
    )
end

"""
$(SIGNATURES)

Append `observation` to the band's sliding window, dropping entries off the
front while the remainder still spans `window_duration`. Returns `estimator`.

The window is bounded in **time**, not in observation count, which is what lets
one producer report 16-chip accumulations and another a single pre-averaged
0.2 s figure under the same configuration and with no scalar state — the span
is computable from the FIFO itself.
"""
function append_noise_observation!(
    estimator::CorrelatorNoiseEstimator{D,T},
    observation::NoiseObservation{D,T},
) where {D,T}
    buffered = estimator.buffered
    push!(buffered, observation)
    _trim_noise_window!(buffered, estimator.window_duration)
    estimator
end

# Drop entries off the front while what is left still spans `window_duration`,
# keeping the window minimal but never shorter than the configured span (and
# never empty, so a single observation longer than the window still counts).
# The span is summed once and decremented as entries go, so this is O(K) per
# append with no allocation.
@inline function _trim_noise_window!(buffered::Vector{<:NoiseObservation}, window_duration)
    total = zero(window_duration)
    @inbounds for i in eachindex(buffered)
        total += buffered[i].duration
    end
    @inbounds while length(buffered) > 1 && total - buffered[1].duration >= window_duration
        total -= buffered[1].duration
        popfirst!(buffered)
    end
    nothing
end

"""
$(SIGNATURES)

The window's `M`-weighted mean density, or `nothing` while it is empty.

Weighted by `num_sub_integrations` because that is the number of independent
looks each entry represents: the density itself needs only the sample count, but
its relative variance is `1/M`, so a one-dump entry and a 64-dump entry combine
correctly only when weighted this way.
"""
function get_noise_density(estimator::CorrelatorNoiseEstimator{D}) where {D}
    buffered = estimator.buffered
    isempty(buffered) && return nothing
    weighted = zero(D)
    total_looks = 0
    @inbounds for i in eachindex(buffered)
        observation = buffered[i]
        weighted += observation.num_sub_integrations * observation.noise_density
        total_looks += observation.num_sub_integrations
    end
    total_looks == 0 && return nothing
    weighted / total_looks
end

noise_density_type(::CorrelatorNoiseEstimator{D}) where {D} = D

"""
$(SIGNATURES)

Number of observations currently in the band's window. Diagnostic only — the
window is bounded in time, so this varies with the producer's dump cadence.
"""
Base.length(estimator::CorrelatorNoiseEstimator) = length(estimator.buffered)
