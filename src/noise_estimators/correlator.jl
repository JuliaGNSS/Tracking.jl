"""
$(SIGNATURES)

The signal's noise reference, measured by **despreading an untracked PRN** — the
only [`AbstractNoiseEstimator`](@ref) Tracking ships, and the one to configure on
a hardware path too (you simply fill it with
[`append_noise_observation!`](@ref) instead of letting [`update_noise!`](@ref)
fill it).

# Why a despread and not a power meter

The reference traverses the **identical** quantise → downconvert → despread →
accumulate path as the prompt, reusing the same kernel, the same replica
generator and — because the estimator is keyed by signal — the same *code*. So
the measured `N₀` already contains every imperfection of that chain — the 1-bit /
2-bit quantisation loss, the quantiser's operating point under load, the input
scaling an AGC step moves, the code amplitude of a CBOC replica — with no
per-backend model and no closed-form correction. A `Σ|x|²` power meter would need
one line per backend and would still miss the second-order coupling where a
strong signal shifts that operating point.

It is also what makes the spectral weighting right rather than assumed. Sharing
the consumer's code means the despread evaluates that signal's own
`∫ S_I(f)·|G(f)|² df` — the coloured-interference case a power meter, or a
reference despread with some *other* signal's code, both get wrong. See
[`AbstractNoiseEstimator`](@ref).

# It is open-loop

The reference has **no feedback of any kind**: it runs at the band's nominal IF,
at a randomly dithered offset from zero Doppler, with a random code phase and a
rotating PRN. There is no discriminator, no loop filter and no NCO update ever
written to it — and, since the randomisation removed the need to know which PRNs
are in use, it reads **nothing** from satellite state, not even the tracked key
set. It is therefore correct with **zero satellites tracked**, and on an FPGA a
noise channel is a strict subset of a tracking channel rather than an addition to
one.

# Why the code phase and the Doppler are randomised

A reference at a *fixed* code phase and *exactly* zero Doppler is a stationary
target, and the failure that matters is not an attacker — it is that a hit never
goes away. The replica is re-anchored to code phase 0 on a grid running at the
**nominal** chip rate, so the relative phase against any incoming signal at zero
Doppler is frozen: whatever it lands on, it stays. A signal that is present but
untracked — a spoofed PRN, or simply a visible satellite the receiver has not
acquired — has a `3 × 1.5 / 1023 ≈ 0.44 %` chance per PRN of sitting inside one
of the three taps, and if it does it contributes `T·(C/N₀)/3 ≈ 10.5·N₀` at
45 dB-Hz to every observation on that PRN, **indefinitely** (≈+1.5 dB on `N̂₀`
after the rotation's dilution). Worse, the geometry is perverse: relative phase
drifts only at the signal's own code Doppler (`f_d/1540`), so *low* Doppler is
both what makes a hit possible and what makes it permanent.

Drawing a fresh code phase and a fresh carrier offset per sub-integration turns
that standing bias into an independent per-observation trial. For a full sky at
45 dB-Hz the residual is `≈0.07 dB` — a handful of 1 %-sized outliers scattered
through a 1000-entry window, rather than a permanent shift of the floor — and,
because a hit now needs *both* draws to land, the two randomisations multiply.

It also makes the untracked-PRN restriction unnecessary, which is why the
reference now rotates over **all** PRNs of the family: a random phase lands
within ±1 chip of a tracked satellite's peak with probability `≈0.6 %`, worth
`10.5/1000 ≈ 0.045 dB` for the one observation it touches — an order of magnitude
under the whole-sky leakage the estimator already carries uncorrected. A constant
32-PRN pool also stops the rotation from shortening (and the dilution from
worsening) exactly as the receiver acquires more satellites.

Neither draw biases the measurement. `N̂₀ = |B|²/(N·A_c²·f_s)` is unbiased for
any starting phase, and for the carrier the reference measures the same noise
power wherever it sits; `carrier_dither` of ±5 kHz smears the spectral weighting
`∫ S_I(f)·|G(f)|² df` by 0.25 % of a 2 MHz main lobe, which no coloured
interferer resolves. On an FPGA an arbitrary code phase is *easier* than phase 0
— a free-running code generator gives you one, where phase 0 needs a reset.

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
  - `carrier_dither` — half-width of the uniform offset added to the band's
    nominal IF, per sub-integration (5 kHz by default, i.e. the terrestrial GNSS
    Doppler spread). Zero pins the reference at the IF exactly, which restores
    half of the stationary target described above — useful to isolate the
    code-phase draw in a test, and not otherwise.
  - `rng` — the source of the code-phase and carrier draws, **seeded and
    therefore reproducible** by default (`Xoshiro(0)`, one stream per estimator,
    advanced in place like `buffered` so the struct is never rebuilt). What the
    randomisation has to defeat is a *stationary* reference, not a reader: an
    attacker cannot observe the chunk grid the phase would be measured against,
    so scattering the draws is the whole requirement and unpredictability buys
    nothing on top. Reproducibility, on the other hand, buys a lot in a library
    whose other guarantees are phrased as bit-identical arithmetic. Pass
    `rng = Random.default_rng()` for a task-local, non-reproducible stream (also
    the one to use if you ever share an estimator across threads).
  - `buffered` — the sliding window itself, a **length-managed FIFO** written in
    place. The `Vector`'s own length is the position, so there is no ring index
    to write back and the struct is never rebuilt — which is what lets per-signal
    state live in an immutable [`TrackState`](@ref).
  - `code_replica` — scratch for the software despread, grown once and reused.
    Held here rather than borrowed from the backend's `ScratchBuffers` so the
    reference can never alias the per-signal replica path.

The sub-integration length is deliberately **not** a field: it is the primary
code period of the signal it measures, derived rather than configured.
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
`Float64` densities per signal.
"""
struct CorrelatorNoiseEstimator{D,T,R<:AbstractRNG} <: AbstractNoiseEstimator
    window_duration::typeof(1.0s)
    tap_code_shift::Float64
    carrier_dither::typeof(1.0Hz)
    buffered::Vector{NoiseObservation{D,T}}
    code_replica::Vector{Int8}
    rng::R
end

"""
$(SIGNATURES)

Construct a [`CorrelatorNoiseEstimator`](@ref) averaging over the last
`window_duration` of observations, with the reference correlator's taps spaced
`tap_code_shift` chips apart and its carrier dithered by up to `carrier_dither`
either side of the band's nominal IF. See the type's docstring for what the
parameters buy and why nothing else is configurable.

`rng` is drawn from for the per-sub-integration code phase and carrier offset. It
defaults to a **seeded** `Xoshiro(0)`, so a run reproduces exactly; pass
`Random.default_rng()` for a task-local, unpredictable stream. See the type's
docstring for why scattering the draws — rather than making them unguessable — is
the whole requirement.

The window is `sizehint!`-ed to four times the number of code-period
observations it expects to hold. The headroom is what makes the FIFO's
`push!`/`popfirst!` pair measure **exactly** zero bytes: at 1× or 2× Julia
periodically shifts the front offset back and reallocates.
"""
function CorrelatorNoiseEstimator(;
    window_duration = 1.0s,
    tap_code_shift = 1.5,
    carrier_dither = 5000.0Hz,
    num_ants::NumAnts = NumAnts(1),
    rng::AbstractRNG = Xoshiro(0),
)
    window_duration > zero(window_duration) ||
        throw(ArgumentError("window_duration must be positive, got $window_duration"))
    tap_code_shift > 0 ||
        throw(ArgumentError("tap_code_shift must be positive, got $tap_code_shift"))
    carrier_dither >= zero(carrier_dither) ||
        throw(ArgumentError("carrier_dither must not be negative, got $carrier_dither"))
    buffered = NoiseObservation{NoiseDensity,typeof(1.0s)}[]
    # Four times the count a 1 ms sub-integration would put in the window; see
    # the docstring for why the headroom rather than an exact fit.
    sizehint!(buffered, 4 * max(1, round(Int, window_duration / 1.0ms)) + 1)
    CorrelatorNoiseEstimator(
        uconvert(s, float(window_duration)),
        Float64(tap_code_shift),
        uconvert(Hz, float(carrier_dither)),
        buffered,
        Int8[],
        rng,
    )
end

"""
$(SIGNATURES)

Append `observation` to the signal's sliding window, dropping entries off the
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

Number of observations currently in the signal's window. Diagnostic only — the
window is bounded in time, so this varies with the producer's dump cadence.
"""
Base.length(estimator::CorrelatorNoiseEstimator) = length(estimator.buffered)

"""
$(SIGNATURES)

Measure this signal's noise over samples `first_sample:last_sample` of
`measurement` and append the resulting observations, returning `estimator`.

The slice is split into `num_sub` **equal** sub-integrations,

```julia
num_sub = max(1, round(Int, slice_duration / code_period))
```

with `code_period` the primary code period of the signal being measured — see
[`CorrelatorNoiseEstimator`](@ref) for why the sub-integration length is derived
rather than configured. Equal slices rather than fixed-length ones so no
remainder is wasted and every window entry is statistically identical; the
`max(1, …)` matters, because a signal whose code period exceeds the
chunk (a 10 ms code with 1 ms chunks) would otherwise yield no observations at
all, forever. Each observation is pushed **individually**: pre-averaging would
put one entry per chunk in the window, re-welding its time span to the Doppler
update rate.

Each sub-integration draws a **fresh random code phase** and a **fresh random
carrier offset** within ±`carrier_dither` of the band's nominal IF, then runs for
its slice. **Nothing here follows a code-period boundary**, and that is
deliberate: the reference despreads with a wrong PRN at a wrong phase, so there
is no signal to accumulate coherently and `N̂₀ = |B|²/(N·A_c²·f_s)` is unbiased
for any `N` and any starting phase. Following boundaries would be actively
harmful — a chunk rarely holds a whole number of code periods, so a
boundary-aligned reference would have to carry a partial accumulator across
calls, which is exactly the scalar state that has no per-signal anchor.
Successive observations are independent because their *sample* ranges are
disjoint; the replica repeating is irrelevant.

The draws are per sub-integration rather than per chunk so that every window
entry is an independent trial — see [`CorrelatorNoiseEstimator`](@ref) for why a
*stationary* phase and Doppler turn a chance alignment with a present-but-
untracked signal into a permanent bias, and the randomisation turns it back into
an occasional single-observation outlier.

All of the correlator's taps are used and pooled into one power sum, at
`tap_code_shift` chips apart so they are genuinely independent looks (see the
type's docstring). The despread runs on the caller's backend kernel — the same
one the prompt goes through — which is what makes the measurement model-free.

The reference is **open-loop**: no discriminator, no loop filter, no NCO update,
and nothing read from satellite state at all. The PRN rotates over the whole
family, tracked or not, because a random phase makes avoiding the tracked ones
unnecessary.
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

    prn = _next_noise_prn(estimator, signal_type)
    # `gen_code_replica!` writes *at* `start_sample` while the kernel reads the
    # replica offset by `start_sample - 1`, so the buffer spans the whole
    # measurement rather than one slice.
    replica_size =
        get_num_samples(measurement) + maximum(sample_shifts) - minimum(sample_shifts)
    length(estimator.code_replica) < replica_size &&
        resize!(estimator.code_replica, replica_size)

    rng = estimator.rng
    code_length = get_code_length(signal_type)
    carrier_dither = estimator.carrier_dither
    intermediate_frequency = measurement.intermediate_frequency

    for sub = 1:num_sub
        # Equal slices: the `sub`-th covers samples `⌊(sub−1)N/num_sub⌋+1 … ⌊subN/num_sub⌋`.
        slice_start = Int(first_sample) + div((sub - 1) * num_samples, num_sub)
        slice_stop = Int(first_sample) + div(sub * num_samples, num_sub) - 1
        slice_samples = slice_stop - slice_start + 1
        slice_samples > 0 || continue
        # One draw each, per sub-integration. Both are free — the kernels already
        # take a code phase and a carrier frequency — and both are what keep a
        # chance alignment with a present-but-untracked signal from freezing into
        # a standing bias. See the type's docstring.
        code_phase = rand(rng) * code_length
        carrier_frequency = intermediate_frequency + (2 * rand(rng) - 1) * carrier_dither
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
            carrier_frequency,
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
        # Straight to the builders' shared core rather than through
        # `noise_observation_from_correlator`'s keyword form: the keywords cost a
        # per-call allocation on Julia 1.10 that 1.11+ elides, and this is the one
        # place it runs per chunk. The arguments are the same ones that builder
        # would forward — note the duration, which is `slice_samples` and not
        # `num_taps` times it, because the taps are simultaneous rather than
        # consecutive: they all integrate the same samples.
        append_noise_observation!(
            estimator,
            _noise_observation(
                power,
                num_taps,
                num_taps * slice_samples,
                sampling_frequency,
                code_amplitude,
                prn,
                slice_samples / sampling_frequency,
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

# Which PRN to borrow a code from. Advancing a PRN needs a position, and a
# position is scalar state with no per-signal anchor — so it is carried in the
# window instead: `NoiseObservation` records the PRN it was measured with, and
# the next one steps on from `last(buffered).prn`. An empty window starts at 1.
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
# The rotation runs over the **whole family**, tracked PRNs included. Skipping
# the tracked ones was what made a fixed code phase 0 safe; a random phase makes
# it unnecessary, and keeping the skip would have cost more than it bought — it
# is the only thing the reference read from satellite state, and it shrank the
# pool (worsening the dilution of any one bad draw) exactly as the receiver
# acquired more satellites. Landing within ±1 chip of a tracked peak now has
# probability ≈0.6 % and is worth ≈0.045 dB for the single observation it
# touches, against the ≈0.9 dB of whole-sky leakage already carried uncorrected.
@inline function _next_noise_prn(
    estimator::CorrelatorNoiseEstimator,
    signal_type::AbstractGNSSSignal,
)
    num_prns = size(get_codes(signal_type), 2)
    previous = isempty(estimator.buffered) ? 0 : Int(last(estimator.buffered).prn)
    mod(previous, num_prns) + 1
end
