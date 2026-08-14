"""
$(SIGNATURES)

C/N₀ against a **measured noise reference**: per record,

```
Ĉ/N₀ = ⟨|P|²⟩ / N̂₀ − 1/T
```

with `N̂₀` the signal's noise density from its [`AbstractNoiseEstimator`](@ref) and
`T` that record's own integration time. The ring averages the per-record terms
and [`estimate_cn0`](@ref) converts the mean once.

Unlike [`NWPRCN0Estimator`](@ref) this is **non-coherent**: it needs no bit sync,
no window length `M`, no coherent window and no fallback, it is immune to
residual carrier-phase error, and it has no saturation ceiling. It is therefore
the one estimator that works uniformly on **every** signal — including GPS
L1C-D and Galileo E1B (`blocks_per_bit == 1`) and any secondary-coded signal
before sync, where NWPR admits no coherent window at all and defers permanently
to a fallback with a different bias (issue #217).

# What it needs

A noise density for its **own signal** — not for its RF band, because what this
divides by is the post-correlation floor and that depends on the despreading
modulation (see [`AbstractNoiseEstimator`](@ref)). On the sample-driven path it
is automatic: [`TrackState`](@ref) provisions a
[`CorrelatorNoiseEstimator`](@ref) for every signal whose estimator asks for one
(see [`requires_noise_density`](@ref)), and `track!` measures before the fold
reads. On a correlator-ingest path you configure the same type and fill it with
[`append_noise_observation!`](@ref) per signal instead. With **no** source
configured for the signal, `update` throws — a wiring mistake, and a silent
substitution would hide a backend that fails to populate the reference. While a
configured source's window is merely still **empty**, the update is skipped
(the fold warns once per signal) and [`estimate_cn0`](@ref) reports `-Inf dB-Hz`
until it fills.

# Three deliberate properties

  - **No `fallback` field.** The three reasons NWPR needs one — warm-up, no
    admissible window pre-sync, and no window ever at one symbol per code block
    — all vanish. One record plus a density is a valid, if noisy, estimate on
    every signal.
  - **`estimate_cn0`'s `integration_time` argument is ignored**, because `T` is
    applied per record at update time, where it is actually known. The signature
    stays for interface uniformity; this is not an oversight.
  - **Only the average is floored, never the per-record term.** Individual terms
    go negative at low C/N₀ and that is exactly what makes the mean unbiased;
    clamping per record would reintroduce a noise floor of the kind
    [`MomentsCN0Estimator`](@ref) has.

# The one bias it carries

The reference despreads with a *wrong* PRN, so besides the thermal floor and the
other satellites' interference — both of which NWPR sees identically — it also
collects **the tracked satellite's own power**, `ε_self/N₀ = C/f_chip`. NWPR does
not: it despreads with the correct code, where the satellite's power appears in
both `NBP` and `WBP` and cancels in the ratio. On GPS L1 C/A that is 0.013 dB at
35 dB-Hz, 0.042 at 40, 0.13 at 45 and 0.40 at 50 (≈10× smaller on the 10.23 Mcps
signals), always reading **low**. It is carried rather than corrected, because
the correction `N̂₀ ← N̂₀ − Σᵢ Ĉᵢ/f_chip` would make every satellite's C/N₀ a
function of every other satellite's estimate — the cross-satellite feedback loop
this design exists without. Below 40 dB-Hz it sits under NWPR's own +0.05 dB
bias, and that is the range where lock and loss decisions are made.

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

This estimator reads its signal's noise density, so a signal carrying it is
provisioned with a [`CorrelatorNoiseEstimator`](@ref). Declared on the type, so
the provisioning decision folds out of a group's slot type — see
[`requires_noise_density`](@ref).
"""
requires_noise_density(::Type{NoiseRefCN0Estimator}) = true

"""
$(SIGNATURES)

Fold one record's `prompt` into the ring: `|P|²/N̂₀ − 1/T` in linear Hz, with the
density and the record's own integration time taken from `context`.

The term is **not** clamped — at low C/N₀ individual terms are negative, and
that is what keeps the mean unbiased.
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

# No noise source configured for this signal: a static property of the setup, so it
# is caught at the first record and loudly. A silent substitution would hide a
# backend that never populates the reference — the whole point of the estimator.
@noinline function update(
    ::NoiseRefCN0Estimator,
    prompt,
    ::CN0UpdateContext{<:AbstractGNSSSignal,<:Unsigned,Nothing},
)
    throw(
        ArgumentError(
            "NoiseRefCN0Estimator needs a noise density for its signal, but no " *
            "AbstractNoiseEstimator is configured for it. Pass " *
            "`cn0_estimator = NWPRCN0Estimator()` if you are feeding externally " *
            "supplied correlator outputs without a noise observation, or give " *
            "`TrackState` a `noise_estimators` entry for this signal and fill it " *
            "with `append_noise_observation!`.",
        ),
    )
end

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

`integration_time` is **ignored**: `T` was applied per record in `update`, where
each record's own value was known — a record lengthened by
[`set_preferred_num_code_blocks_to_integrate!`](@ref) is therefore handled
correctly even when it sits in the ring beside shorter ones. The argument stays
for interface uniformity with the other estimators.

An empty ring, and a mean that has not cleared zero, both report `-Inf dB-Hz` —
the house convention that a missing estimate is `-Inf` and never `NaN` (see
[`NoCN0Estimator`](@ref) for why).
"""
function estimate_cn0(estimator::NoiseRefCN0Estimator, integration_time)
    filled = length(estimator)
    filled == 0 && return dBHz(0.0Hz)
    total = 0.0
    @inbounds for i = 1:filled
        total += estimator.buffered_cn0[i]
    end
    mean_cn0 = total / filled
    mean_cn0 <= 0 && return dBHz(0.0Hz)
    dBHz(mean_cn0 * Hz)
end
