"""
$(SIGNATURES)

C/N₀ against a **measured noise reference**: per record,

```
Ĉ/N₀ = ⟨|P|²⟩ / N̂₀ − 1/T
```

with `N̂₀` the signal's noise density from its [`AbstractNoiseEstimator`](@ref) and
`T` that record's own integration time. The ring averages the per-record terms
and [`estimate_cn0`](@ref) converts the mean once.

Unlike [`NWPRCN0Estimator`](@ref) this is **non-coherent**: no bit sync, no
window, no `M`, no fallback, no saturation ceiling, and immune to residual
carrier-phase error. It therefore works uniformly on **every** signal, including
GPS L1C-D, Galileo E1B and secondary-coded signals before sync (issue #217).

# What it needs

A noise density for its **own signal** (see [`AbstractNoiseEstimator`](@ref) for
why per signal). On the sample-driven path [`TrackState`](@ref) provisions a
[`CorrelatorNoiseEstimator`](@ref) automatically (see
[`requires_noise_density`](@ref)); on a correlator-ingest path fill the same type
with [`append_noise_observation!`](@ref). With **no** source configured, `update`
throws: a silent substitution would hide a backend that never populates the
reference. While a configured window is still **empty**, the fold skips the
update (warning once per signal) and [`estimate_cn0`](@ref) reports `-Inf dB-Hz`.

# Deliberate properties

  - **No `fallback`.** One record plus a density is a valid, if noisy, estimate
    on every signal.
  - **`estimate_cn0`'s `integration_time` is ignored**; `T` is applied per record
    in `update`.
  - **Only the average is floored, never the per-record term.** Negative terms
    at low C/N₀ keep the mean unbiased; clamping them would reintroduce a floor
    like [`MomentsCN0Estimator`](@ref)'s.

# The one bias it carries

The reference despreads with a *wrong* PRN, so it also collects the tracked
satellite's own power, `ε_self/N₀ = C/f_chip`, which NWPR cancels in its ratio:
C/N₀ reads low by ≈0.13 dB at 45 dB-Hz and ≈0.40 at 50 on GPS L1 C/A. It is not
corrected, since that would couple every satellite's C/N₀ to every other's. See
"The one bias it carries" in docs/src/cn0_estimator.md.

# Fields / configuration

  - `num_records` — how many records the estimate averages over, i.e. the memory
    of the estimator (100 by default, ~100 ms at GPS L1 C/A).
  - `buffered_cn0`, `current_index`, `filled_length` — the ring of per-record
    C/N₀ terms in **linear Hz**, written in place, plus its position and fill.
"""
struct NoiseRefCN0Estimator <: AbstractCN0Estimator
    num_records::Int
    buffered_cn0::Vector{Float64}
    current_index::Int
    filled_length::Int
end

"""
$(SIGNATURES)

Construct a fresh [`NoiseRefCN0Estimator`](@ref) averaging over the last
`num_records` records.
"""
function NoiseRefCN0Estimator(; num_records::Int = 100)
    num_records >= 1 ||
        throw(ArgumentError("num_records must be at least 1, got $num_records"))
    NoiseRefCN0Estimator(num_records, zeros(Float64, num_records), 0, 0)
end

length(estimator::NoiseRefCN0Estimator) = estimator.filled_length
get_current_index(estimator::NoiseRefCN0Estimator) = estimator.current_index

"""
$(SIGNATURES)

`true`: a signal carrying this estimator is provisioned a
[`CorrelatorNoiseEstimator`](@ref); see [`requires_noise_density`](@ref).
"""
requires_noise_density(::Type{NoiseRefCN0Estimator}) = true

"""
$(SIGNATURES)

Fold one record's `prompt` into the ring: `|P|²/N̂₀ − 1/T` in linear Hz, with the
density and the record's own integration time taken from `context`. The term is
**not** clamped (see [`NoiseRefCN0Estimator`](@ref)).
"""
function update(estimator::NoiseRefCN0Estimator, prompt, context::CN0UpdateContext)
    cn0 = ustrip(
        Hz,
        uconvert(Hz, abs2(prompt) / context.noise_density - 1 / context.integration_time),
    )
    buffer_length = Base.length(estimator.buffered_cn0)
    next_index = mod(estimator.current_index, buffer_length) + 1
    estimator.buffered_cn0[next_index] = cn0
    NoiseRefCN0Estimator(
        estimator.num_records,
        estimator.buffered_cn0,
        next_index,
        min(estimator.filled_length + 1, buffer_length),
    )
end

# A missing density or integration time is `Nothing` in the context's type, so
# both are caught by dispatch on the first record. The third method resolves the
# ambiguity of `{…,Nothing,Nothing}` and reports the missing source, the root
# cause.
@noinline update(
    ::NoiseRefCN0Estimator,
    prompt,
    ::CN0UpdateContext{<:AbstractGNSSSignal,<:Unsigned,Nothing},
) = _throw_noise_ref_needs_density()
@noinline update(
    ::NoiseRefCN0Estimator,
    prompt,
    ::CN0UpdateContext{<:AbstractGNSSSignal,<:Unsigned,Nothing,Nothing},
) = _throw_noise_ref_needs_density()
@noinline update(
    ::NoiseRefCN0Estimator,
    prompt,
    ::CN0UpdateContext{<:AbstractGNSSSignal,<:Unsigned,<:Any,Nothing},
) = _throw_noise_ref_needs_integration_time()

@noinline _throw_noise_ref_needs_density() = throw(
    ArgumentError(
        "NoiseRefCN0Estimator needs a noise density for its signal, but no " *
        "AbstractNoiseEstimator is configured for it. Pass " *
        "`cn0_estimator = NWPRCN0Estimator()` if you are feeding externally " *
        "supplied correlator outputs without a noise observation, or give " *
        "`TrackState` a `noise_estimators` entry for this signal and fill it " *
        "with `append_noise_observation!`.",
    ),
)

@noinline _throw_noise_ref_needs_integration_time() = throw(
    ArgumentError(
        "NoiseRefCN0Estimator has this signal's noise density but no integration " *
        "time: it applies `T` per record (`|P|²/N̂₀ - 1/T`), so it cannot fold one " *
        "without it. Build the context with both — " *
        "`CN0UpdateContext(signal, bit_buffer, num_code_blocks; noise_density, " *
        "integration_time)` — which is what the tracking loop does; " *
        "`integration_time` is that record's own " *
        "`integrated_samples / sampling_frequency`.",
    ),
)

"""
$(SIGNATURES)

A bare prompt stream carries no noise density and no integration time, so this
estimator has no two-argument form — use the three-argument
[`update(::AbstractCN0Estimator, ::Any, ::CN0UpdateContext)`](@ref), which is
what the tracking loop calls, or [`MomentsCN0Estimator`](@ref) /
[`NWPRCN0Estimator`](@ref) for folding a captured prompt stream by hand.
"""
@noinline update(::NoiseRefCN0Estimator, prompt) = throw(
    ArgumentError(
        "NoiseRefCN0Estimator cannot be updated from a bare prompt: it needs the " *
        "signal's noise density and the record's integration time, both of which " *
        "live in the CN0UpdateContext. Call the three-argument `update`.",
    ),
)

"""
$(SIGNATURES)

Mean of the buffered per-record terms, converted once with `dBHz`.

`integration_time` is **ignored**: `T` was applied per record in `update`, so a
record lengthened by [`set_preferred_num_code_blocks_to_integrate!`](@ref) is
handled correctly beside shorter ones. The argument stays for interface
uniformity.

An empty ring, a mean that has not cleared zero, and a **non-finite** mean all
report `-Inf dB-Hz`, never `NaN` (see [`NoCN0Estimator`](@ref) for why). The
explicit `isfinite` matters because a `NaN` passes `mean_cn0 <= 0`;
`_noise_density_and_ready` keeps a zero floor out of the ring, this catches a
`NaN` prompt.
"""
function estimate_cn0(estimator::NoiseRefCN0Estimator, integration_time)
    filled = length(estimator)
    filled == 0 && return dBHz(0.0Hz)
    total = 0.0
    @inbounds for i = 1:filled
        total += estimator.buffered_cn0[i]
    end
    mean_cn0 = total / filled
    isfinite(mean_cn0) && mean_cn0 > 0 || return dBHz(0.0Hz)
    dBHz(mean_cn0 * Hz)
end
