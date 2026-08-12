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

"""
$(SIGNATURES)

Measure the band's noise over samples `first_sample:last_sample` of
`measurement` and append the resulting observations, returning `estimator`.

The slice is split into `num_sub` **equal** sub-integrations,

```julia
num_sub = max(1, round(Int, slice_duration / code_period))
```

with `code_period` the primary code period of the band's reference signal — see
[`CorrelatorNoiseEstimator`](@ref) for why the sub-integration length is derived
rather than configured. Equal slices rather than fixed-length ones so no
remainder is wasted and every window entry is statistically identical; the
`max(1, …)` matters, because a band whose reference code period exceeds the
chunk (a 10 ms code with 1 ms chunks) would otherwise yield no observations at
all, forever. Each observation is pushed **individually**: pre-averaging would
put one entry per chunk in the window, re-welding its time span to the Doppler
update rate.

Each sub-integration starts its replica at code phase 0 and runs for its slice.
**Nothing here follows a code-period boundary**, and that is deliberate: the
reference despreads with a wrong PRN at a wrong phase, so there is no signal to
accumulate coherently and `N̂₀ = |B|²/(N·A_c²·f_s)` is unbiased for any `N` and
any starting phase. Following boundaries would be actively harmful — a chunk
rarely holds a whole number of code periods, so a boundary-aligned reference
would have to carry a partial accumulator across calls, which is exactly the
scalar state that has no per-band anchor. Successive observations are
independent because their *sample* ranges are disjoint; the replica repeating is
irrelevant.

All of the correlator's taps are used and pooled into one power sum, at
`tap_code_shift` chips apart so they are genuinely independent looks (see the
type's docstring). The despread runs on the caller's backend kernel — the same
one the prompt goes through — which is what makes the measurement model-free.

The reference is **open-loop**: the band's nominal IF as the carrier frequency,
zero Doppler and zero carrier phase. ±5 kHz of Doppler is 0.5 % of a 1 MHz code
band, so no satellite's state is needed and the reference is correct with
nothing tracked at all.
"""
function update_noise!(
    estimator::CorrelatorNoiseEstimator,
    measurement::BandMeasurement,
    first_sample::Integer,
    last_sample::Integer,
    context::NoiseUpdateContext,
)
    num_samples = Int(last_sample) - Int(first_sample) + 1
    num_samples > 0 || return estimator
    signal_type = context.signal
    sampling_frequency = measurement.sampling_frequency
    # Nominal chip rate: the reference runs at zero Doppler.
    code_frequency = get_code_frequency(signal_type)
    code_period = get_code_length(signal_type) / code_frequency
    slice_duration = num_samples / sampling_frequency
    num_sub = max(1, round(Int, uconvert(NoUnits, slice_duration / code_period)))
    num_sub = min(num_sub, num_samples)          # never fewer than one sample per slice

    # One antenna only, and it must be the one `DefaultPostCorrFilter` selects
    # (`last`) — the ratio would otherwise compare two different channels. A
    # post-corr filter with a noise gain other than unity invalidates the ratio
    # outright; that is documented as a requirement rather than corrected.
    samples = _noise_reference_samples(measurement.samples)
    correlator = EarlyPromptLateCorrelator(;
        preferred_early_late_to_prompt_code_shift = estimator.tap_code_shift,
    )
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    num_taps = length(sample_shifts)
    code_amplitude = get_code_amplitude(signal_type)

    prn, code_phase = _next_noise_prn(estimator, signal_type, context)
    # `gen_code_replica!` writes *at* `start_sample` while the kernel reads the
    # replica offset by `start_sample - 1`, so the buffer spans the whole
    # measurement rather than one slice.
    replica_size =
        get_num_samples(measurement) + maximum(sample_shifts) - minimum(sample_shifts)
    length(estimator.code_replica) < replica_size &&
        resize!(estimator.code_replica, replica_size)

    for sub = 1:num_sub
        # Equal slices: the `sub`-th covers samples `⌊(sub−1)N/num_sub⌋+1 … ⌊subN/num_sub⌋`.
        slice_start = Int(first_sample) + div((sub - 1) * num_samples, num_sub)
        slice_stop = Int(first_sample) + div(sub * num_samples, num_sub) - 1
        slice_samples = slice_stop - slice_start + 1
        slice_samples > 0 || continue
        accumulators = _correlate_noise_reference!(
            context.downconvert_and_correlator,
            estimator.code_replica,
            samples,
            zero(correlator),
            signal_type,
            prn,
            sample_shifts,
            code_phase,
            code_frequency,
            measurement.intermediate_frequency,
            sampling_frequency,
            slice_start,
            slice_samples,
        )
        # Pool the taps. At ≥1 chip spacing they are independent looks at one
        # scalar, so nothing about their relative values means anything — the
        # opposite of the prompt path, where the taps' *differences* are the code
        # and carrier discriminants.
        power = 0.0
        for k = 1:num_taps
            power += abs2(accumulators[k])
        end
        append_noise_observation!(
            estimator,
            noise_observation_from_correlator(
                power,
                num_taps,
                num_taps * slice_samples,
                sampling_frequency;
                code_amplitude,
                prn,
                # The taps are simultaneous, not consecutive: they all integrate
                # the same `slice_samples`, so that — and not `num_taps` times it
                # — is what the time-bounded window must count.
                duration = slice_samples / sampling_frequency,
            ),
        )
    end
    estimator
end

# The antenna the reference despreads. `DefaultPostCorrFilter` is `last(x)`, so
# an antenna array's noise must be measured on its last column too.
@inline _noise_reference_samples(samples::AbstractVector) = samples
@inline _noise_reference_samples(samples::AbstractMatrix) =
    view(samples, :, size(samples, 2))

# Which PRN to borrow a code from, and at which phase. Advancing a PRN needs a
# position, and a position is scalar state with no per-band anchor — so it is
# carried in the window instead: `NoiseObservation` records the PRN it was
# measured with, and the next one steps on from `last(buffered).prn`. An empty
# window starts at the lowest untracked PRN.
#
# Rotation is worth this much and no more. The leakage term is already a sum over
# the whole sky, so it is averaged over that many cross-correlation draws; what a
# *fixed* PRN would cost is a slowly drifting bias (≈±0.28 dB on the ≈0.9 dB
# whole-sky leakage, moving as the geometry does). Rotation replaces that drift
# with a stable mean, which — since the leakage is deliberately uncorrected — is
# the whole objective. It does not touch the self-leakage term at all: every
# wrong PRN collects the tracked satellite's own power to the same expected
# degree.
#
# Preferring an untracked PRN is what makes code phase 0 safe. If every PRN of
# the family is tracked, reuse one and offset by half a code period, which is far
# outside the correlation peak.
@inline function _next_noise_prn(
    estimator::CorrelatorNoiseEstimator,
    signal_type::AbstractGNSSSignal,
    context::NoiseUpdateContext,
)
    num_prns = size(get_codes(signal_type), 2)
    previous = isempty(estimator.buffered) ? 0 : Int(last(estimator.buffered).prn)
    for step = 1:num_prns
        prn = mod(previous - 1 + step, num_prns) + 1
        _is_prn_tracked(context, prn) || return prn, 0.0
    end
    prn = mod(previous, num_prns) + 1
    prn, get_code_length(signal_type) / 2
end
