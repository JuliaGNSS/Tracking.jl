# Running sums over a `CorrelatorNoiseEstimator`'s window, so neither appending
# nor reading walks it: `span` drives the trim, `weighted_density / looks` is
# `get_noise_density`.
#
# `stale` counts appends since the last exact recomputation. Incremental Float64
# add/subtract drifts (a producer may mix 0.2 s entries with 1 ms ones), so the
# sums are rebuilt once per window's worth of appends: bounded drift, O(1)
# amortised.
#
# Immutable and isbits, held in the estimator's single `Ref` (one inline write
# per update); a cache derived from `buffered`, so the estimator itself is never
# rebuilt.
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

It is open-loop, model-free and randomised per sub-integration; see "The
software source" in the docs for why.

Neither draw biases the measurement: `N̂₀ = |B|²/(N·A_c²·f_s)` is unbiased for
any starting phase, and ±5 kHz of `carrier_dither` smears the spectral weighting
by 0.25 % of a 2 MHz main lobe, which no coloured interferer resolves.

# Fields / configuration

  - `window_duration` — how far back the sliding window reaches (1 s by
    default): lower variance against a longer smear across AGC changes. The
    default holds `K_n ≈ 1000` looks per tap, ≤0.08 dB from a variance-free
    reference at every C/N₀.

  - `tap_code_shift` — tap spacing in chips (1.5 by default). Taps are only
    independent looks when their codes decorrelate: since
    `corr(|Bᵢ|²,|Bⱼ|²) = |ρᵢⱼ|²`, three taps at ±0.5 chip are worth 2.25 looks,
    at ≥1 chip 2.98. 1.5 chips sits in the autocorrelation null, clear of both
    the 1- and 2-chip sidelobe values.

  - `carrier_dither` — half-width of the uniform offset added to the band's
    nominal IF per sub-integration (5 kHz by default, the terrestrial Doppler
    spread). Zero disables the carrier draw; useful only in tests.

  - `rng` — source of the code-phase and carrier draws, seeded `Xoshiro(0)` by
    default so a run repeats. The draws only need to be scattered, not
    unpredictable. Pass `Random.default_rng()` for a task-local stream (also if
    sharing an estimator across threads). The stream changes across Julia
    versions, so do not tune a tolerance to a particular draw; design against the
    `1/√(3K)` relative error over `K` observations.

  - `buffered` — the sliding window, a FIFO written in place, so the struct is
    never rebuilt and can live in an immutable [`TrackState`](@ref).

  - `totals` — running sums over `buffered` (see `NoiseWindowTotals`), keeping
    `append_noise_observation!` and `get_noise_density` **O(1)** even though
    the window is deliberately large (≈2500 entries at 1 s and 0.4 ms chunks).

The sub-integration length is deliberately **not** a field: it is the signal's
primary code period.
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
`tap_code_shift` chips apart, its carrier dithered by up to `carrier_dither`
either side of the band's nominal IF, and `rng` for the per-sub-integration
draws. See the type's docstring for each.

The window is `sizehint!`-ed to four times the code-period observations it
expects; with less headroom the FIFO's `push!`/`popfirst!` periodically
reallocates.

`num_ants` must match the antenna count of the signal group this estimator is
keyed to. Above one antenna the window holds the `M×M` spatial covariance `R̂`,
which each satellite reduces to its own floor `wᴴR̂w` (see
[`AbstractPostCorrFilter`](@ref)). `TrackState` provisions the right count when
`noise_estimators` is left at `nothing`, and rejects a mismatch otherwise.
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
    D = _density_type_for_num_ants(num_ants)
    buffered = NoiseObservation{D,typeof(1.0s)}[]
    # Four times the 1 ms sub-integration count; see the docstring.
    sizehint!(buffered, 4 * max(1, round(Int, window_duration / 1.0ms)) + 1)
    CorrelatorNoiseEstimator(
        uconvert(s, float(window_duration)),
        Float64(tap_code_shift),
        uconvert(Hz, float(carrier_dither)),
        buffered,
        Ref(NoiseWindowTotals{D,typeof(1.0s)}()),
        rng,
    )
end

# The antenna count this estimator despreads, derived from its density type so
# the window and the reference correlator cannot disagree.
@inline _num_ants(::CorrelatorNoiseEstimator{D}) where {D} = _num_ants_of_density_type(D)

"""
$(SIGNATURES)

Append `observation` to the signal's sliding window, dropping entries off the
front while the remainder still spans `window_duration`. Returns `estimator`.

The window is bounded in **time**, not in observation count, so producers with
different dump lengths share one configuration. **O(1)** per call, amortised.

Any [`NoiseObservation`](@ref) is retyped onto the window's field types (free for
builder output), so a hand-assembled or `Float32` one is not silently dropped by
the abstract no-op.
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

# Drop entries off the front while what is left still spans `window_duration`
# (never emptying it, so a single observation longer than the window counts).
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

# Recompute the totals exactly once per window's worth of appends; see
# `NoiseWindowTotals`.
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

Weighted by `num_sub_integrations`, the independent looks per entry (see
[`NoiseObservation`](@ref)). **O(1)**: read off the running totals, since it is
called once per signal per chunk.
"""
function get_noise_density(estimator::CorrelatorNoiseEstimator)
    totals = estimator.totals[]
    totals.looks == 0 && return nothing
    totals.weighted_density / totals.looks
end

noise_density_type(::CorrelatorNoiseEstimator{D}) where {D} = D

# `looks` is the running sum of `num_sub_integrations`, i.e. the look count.
noise_window_looks(estimator::CorrelatorNoiseEstimator) = estimator.totals[].looks

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

with `code_period` the primary code period of the signal being measured. Equal
slices so every window entry is statistically identical; the `max(1, …)` keeps a
code period longer than the chunk from yielding no observations ever. Each
observation is pushed **individually**, so the window's span is not tied to the
chunk rate.

Each sub-integration draws a fresh random code phase and carrier offset (see
[`CorrelatorNoiseEstimator`](@ref)). **Nothing follows a code-period boundary**:
with a wrong PRN there is nothing to accumulate coherently, `N̂₀` is unbiased for
any `N` and starting phase, and aligning to boundaries would need a partial
accumulator carried across calls. Observations are independent because their
sample ranges are disjoint.

All taps are pooled into one power (see `_pool_taps`), on the caller's backend
kernel. With `num_ants > 1` every antenna column is despread and the taps pool
into the spatial covariance `R̂ = Σ_k b_k·b_kᴴ`; each satellite reduces it to its
own floor `wᴴR̂w` at C/N₀ time.
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

    # Every antenna column, despread through a correlator of the estimator's own
    # width, so each satellite can ask for the floor its own beamformer sees.
    samples = measurement.samples
    num_ants = _num_ants(estimator)
    correlator = EarlyPromptLateCorrelator(;
        num_ants,
        preferred_early_late_to_prompt_code_shift = estimator.tap_code_shift,
    )
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    num_taps = length(sample_shifts)
    code_amplitude = get_code_amplitude(signal_type)

    prn = _next_noise_prn(estimator, signal_type)
    # `gen_code_replica!` writes *at* `start_sample` while the kernel reads the
    # replica offset by `start_sample - 1`, so the buffer spans the whole
    # measurement. Backends that pack the code plane in-kernel ignore it.
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
        # Fresh draws per sub-integration; see the type's docstring.
        code_phase = rand(rng) * code_length
        carrier_frequency = intermediate_frequency + (2 * rand(rng) - 1) * carrier_dither
        # The per-satellite despread primitive. `carrier_phase = 0.0`: open-loop,
        # no NCO to continue. `use_band_cache = false`: the group's packed sign
        # planes may belong to another band or buffer (see `_despread_one_signal!`).
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
        power = _pool_taps(accumulators, num_ants)
        # The builders' core directly: the keyword form allocates per call on
        # Julia 1.10. The duration is one slice, not `num_taps` of them, because
        # the taps integrate the same samples.
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

# Pool the taps, which at ≥1 chip spacing are independent looks at the same
# noise: `Σ_k |b_k|²` for one antenna, `Σ_k b_k·b_kᴴ` for an array (diagonal =
# per-antenna power, off-diagonals = noise correlation).
#
# Iterate rather than index over `1:num_taps`: a runtime index into an `SVector`
# defeats the unroll and, at `M > 1`, makes Julia 1.10 spill to the heap.
@inline function _pool_taps(accumulators, ::NumAnts{1})
    power = 0.0
    for tap in accumulators
        power += abs2(tap)
    end
    power
end

@inline function _pool_taps(accumulators, ::NumAnts{M}) where {M}
    covariance = zero(SMatrix{M,M,ComplexF64,M * M})
    for tap in accumulators
        covariance += tap * tap'
    end
    covariance
end

# Which PRN to borrow a code from. The rotation position is carried in the window
# (`last(buffered).prn`) rather than as scalar state; an empty window starts at 1.
#
# Rotation turns a fixed PRN's slowly drifting cross-correlation bias (≈±0.28 dB
# on ≈0.9 dB of whole-sky leakage) into a stable mean; it does not affect the
# self-leakage. It covers the whole family, tracked PRNs included: with a random
# code phase, landing within ±1 chip of a tracked peak costs only ≈0.045 dB for
# that one observation (≈0.6 % of draws).
@inline function _next_noise_prn(
    estimator::CorrelatorNoiseEstimator,
    signal_type::AbstractGNSSSignal,
)
    num_prns = size(get_codes(signal_type), 2)
    previous = isempty(estimator.buffered) ? 0 : Int(last(estimator.buffered).prn)
    mod(previous, num_prns) + 1
end
