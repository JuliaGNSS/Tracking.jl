# Running sums over a `CorrelatorNoiseEstimator`'s window, so that neither
# appending to it nor reading it has to walk it. `span` answers "does the window
# still cover `window_duration`?" for the trim, and `weighted_density / looks` is
# `get_noise_density` outright.
#
# `stale` counts appends since the last exact recomputation. Incremental
# add/subtract on `Float64` drifts, and `append_noise_observation!` is public — a
# producer may legitimately mix a 0.2 s pre-averaged entry with 1 ms ones, so the
# magnitudes are not guaranteed comparable and the cancellation is not purely
# hypothetical. Rebuilding once per window's worth of appends bounds the drift to
# one window of operations while staying O(1) amortised.
#
# An ordinary immutable value, rebuilt whole on every change and stored in the
# estimator's single `Ref`. Immutable rather than a mutable struct because
# nothing in `src` is one — that is what lets per-signal state live in an
# immutable `TrackState` — and one `Ref` of a four-field value rather than four
# `Ref`s of a field each: it is isbits, so the box holds it inline and an update
# is one write to one cache line instead of four chases to four heap objects.
#
# The estimator itself is still never rebuilt. Only this cache is, and it is a
# cache of `buffered`'s contents — derived state with an anchor, recomputed from
# it exactly (see `stale`) — not the free-floating scalar the design rules out.
struct NoiseWindowTotals{D,T}
    span::T
    weighted_density::D
    looks::Int
    stale::Int
end

NoiseWindowTotals{D,T}() where {D,T} = NoiseWindowTotals{D,T}(zero(T), zero(D), 0, 0)

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

  - `rng` — the source of the code-phase and carrier draws, **seeded** by default
    (`Xoshiro(0)`, one stream per estimator, advanced in place like `buffered` so
    the struct is never rebuilt). What the randomisation has to defeat is a
    *stationary* reference, not a reader: an attacker cannot observe the chunk
    grid the phase would be measured against, so scattering the draws is the whole
    requirement and unpredictability buys nothing on top. A seeded default keeps a
    run repeatable, which is worth having in a library whose other guarantees are
    phrased as bit-identical arithmetic. Pass `rng = Random.default_rng()` for a
    task-local stream (also the one to use if you ever share an estimator across
    threads).

    Repeatable **on one Julia version**, and no further: `Xoshiro`'s stream is not
    part of Julia's compatibility guarantee and does change across releases. Do
    not build a tolerance around a particular draw — size the window so the
    assertion holds for *any* draw. Every `N̂₀` here is a `1/√(3K)` estimate over
    `K` observations, and that is the number to design against.

  - `buffered` — the sliding window itself, a **length-managed FIFO** written in
    place. The `Vector`'s own length is the position, so there is no ring index
    to write back and the struct is never rebuilt — which is what lets per-signal
    state live in an immutable [`TrackState`](@ref).

  - `totals` — the window's running sums (span, `M`-weighted density, looks),
    maintained as entries are pushed and dropped. They are what keeps both
    `append_noise_observation!` and `get_noise_density` **O(1)**; recomputing
    them per call made each an O(K) scan, and since the whole design buys its
    variance by making `K` large (≈2500 entries at a 1 s window and a 0.4 ms
    chunk), that put a per-chunk cost directly proportional to the accuracy
    asked for. Written in place through `Ref` cells, exactly as `buffered` is,
    so the "never rebuilt" property still holds — and a *cache* of `buffered`
    rather than state of its own, recomputed from it exactly often enough that
    the incremental arithmetic cannot drift.

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
    totals::Base.RefValue{NoiseWindowTotals{D,T}}
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
defaults to a **seeded** `Xoshiro(0)`, so a run repeats on a given Julia version;
pass `Random.default_rng()` for a task-local stream. See the type's docstring for
why scattering the draws — rather than making them unguessable — is the whole
requirement, and for why the seed is not something to calibrate against.

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
        Ref(NoiseWindowTotals{NoiseDensity,typeof(1.0s)}()),
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
is a property of the FIFO rather than of a configured count.

**O(1)** per call, amortised, whatever the window holds — see `totals` in the
type's docstring for why that matters here rather than being a micro-optimisation.

Any [`NoiseObservation`](@ref) is accepted and retyped onto the window's own field
types, which costs nothing for one the builders produced (they already emit the
canonical pair) and is what keeps a hand-assembled or `Float32` one from matching
no method here and falling through to the abstract no-op — where it would be
**dropped silently** and leave the window empty forever.
"""
function append_noise_observation!(
    estimator::CorrelatorNoiseEstimator{D,T},
    observation::NoiseObservation,
) where {D,T}
    observation = convert(NoiseObservation{D,T}, observation)
    buffered = estimator.buffered
    totals = estimator.totals
    push!(buffered, observation)
    _add_observation!(totals, observation)
    _trim_noise_window!(buffered, totals, estimator.window_duration)
    _refresh_totals_if_stale!(buffered, totals)
    estimator
end

@inline function _add_observation!(totals::Base.RefValue{<:NoiseWindowTotals}, observation)
    t = totals[]
    totals[] = typeof(t)(
        t.span + observation.duration,
        t.weighted_density + observation.num_sub_integrations * observation.noise_density,
        t.looks + observation.num_sub_integrations,
        t.stale,
    )
    nothing
end

@inline function _drop_observation!(totals::Base.RefValue{<:NoiseWindowTotals}, observation)
    t = totals[]
    totals[] = typeof(t)(
        t.span - observation.duration,
        t.weighted_density - observation.num_sub_integrations * observation.noise_density,
        t.looks - observation.num_sub_integrations,
        t.stale,
    )
    nothing
end

# Drop entries off the front while what is left still spans `window_duration`,
# keeping the window minimal but never shorter than the configured span (and
# never empty, so a single observation longer than the window still counts).
# The span comes from `totals` rather than being re-summed, so an append that
# drops one entry does O(1) work rather than walking the whole window.
@inline function _trim_noise_window!(
    buffered::Vector{<:NoiseObservation},
    totals::Base.RefValue{<:NoiseWindowTotals},
    window_duration,
)
    @inbounds while length(buffered) > 1 &&
                    totals[].span - buffered[1].duration >= window_duration
        _drop_observation!(totals, buffered[1])
        popfirst!(buffered)
    end
    nothing
end

# Recompute the totals exactly, once per window's worth of appends. See
# `NoiseWindowTotals` for why the incremental sums cannot simply be trusted
# forever; the O(K) walk every K appends is O(1) amortised.
@inline function _refresh_totals_if_stale!(
    buffered::Vector{<:NoiseObservation},
    totals::Base.RefValue{<:NoiseWindowTotals},
)
    t = totals[]
    totals[] = typeof(t)(t.span, t.weighted_density, t.looks, t.stale + 1)
    t.stale + 1 < length(buffered) && return nothing
    _refresh_totals!(buffered, totals)
end

function _refresh_totals!(
    buffered::Vector{<:NoiseObservation},
    totals::Base.RefValue{NoiseWindowTotals{D,T}},
) where {D,T}
    span = zero(T)
    weighted_density = zero(D)
    looks = 0
    @inbounds for i in eachindex(buffered)
        observation = buffered[i]
        span += observation.duration
        weighted_density += observation.num_sub_integrations * observation.noise_density
        looks += observation.num_sub_integrations
    end
    totals[] = NoiseWindowTotals{D,T}(span, weighted_density, looks, 0)
    nothing
end

"""
$(SIGNATURES)

The window's `M`-weighted mean density, or `nothing` while it is empty.

Weighted by `num_sub_integrations` because that is the number of independent
looks each entry represents: the density itself needs only the sample count, but
its relative variance is `1/M`, so a one-dump entry and a 64-dump entry combine
correctly only when weighted this way.

Read straight off the window's running totals, so this is **O(1)**. It is called
once per signal per chunk from the Doppler-estimator fold, which is why it may
not walk the window: doing so charged every chunk for the window length, i.e.
for the very `K` the estimator's accuracy is bought with.
"""
function get_noise_density(estimator::CorrelatorNoiseEstimator)
    totals = estimator.totals[]
    totals.looks == 0 && return nothing
    totals.weighted_density / totals.looks
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
    # measurement rather than one slice. Only the software backends read it at
    # all; `_despread_one_signal!` sizes it from this on those that do and ignores
    # it on the ones that pack the code plane in-kernel.
    replica_size =
        get_num_samples(measurement) + maximum(sample_shifts) - minimum(sample_shifts)

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
        # The same despread primitive the per-satellite path runs, which is what
        # keeps the measurement model-free: `carrier_phase = 0.0` because the
        # reference is open-loop and has no NCO to continue from, and
        # `use_band_cache = false` because this runs before any group has packed
        # the bit-wise backends' shared sign planes.
        accumulators = get_accumulators(
            _despread_one_signal!(
                context.downconvert_and_correlator,
                zero(correlator),
                samples,
                signal_type,
                prn,
                sample_shifts,
                code_phase,
                0.0,
                code_frequency,
                carrier_frequency,
                sampling_frequency,
                slice_start,
                slice_samples,
                replica_size,
                false,
            ),
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
